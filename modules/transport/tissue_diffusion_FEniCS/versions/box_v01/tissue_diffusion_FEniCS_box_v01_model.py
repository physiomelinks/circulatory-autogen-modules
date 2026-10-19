"""
Diffusion of one solute in a 3D box of tissue, solved with FEniCSx (dolfinx 0.8-0.9), coupled to
0D modules (capillaries, a varicosity, ...) through libcuflynx.coupling.

    dC/dt = sigma_diff lap(C) + sum_k chi_k J_k / V_k + M C / (C + C50)      in the box
    no flux through the box's faces

Each connected 0D module k owns an exchange region: one cell of a grid_nx x grid_ny x grid_nz
block grid of the box, placed along x at the middle of y and z (neighbour k of N at
x-cell floor((k + 1/2) grid_nx / N)). It sends J_k (J_s, mol/s, e.g. a capillary's O2 flux or a
varicosity's NE release) into its region as a uniform source, and receives the region's mean
concentration C_k (C_t, mM). These are the variables of the module's capillary_to_flux_port,
the same port the finite-volume CellML module tissue_diffusion_volume has, so a model can swap
one for the other.

Two discretisations (parameter fv_scheme):

* fv_scheme = 0: Q1 finite elements on the grid refined ``refine`` times per direction;
* fv_scheme = 1: one DG0 value per grid cell with two-point face fluxes
  sigma_diff A_f (C_a - C_b) / |x_a - x_b| -- the same discrete equations as a grid of
  tissue_diffusion_volume cells joined by tissue_diffusion_face faces, which makes the two
  directly comparable.

Time stepping: theta method (theta = 1/2 is Crank-Nicolson) for diffusion, with the source and
the reaction from the start of the step. The matrix is assembled and factorised once per step
size. Under MPI the mesh is distributed and the region means are summed over the ranks, so every
rank returns the same values.

Field snapshots (optional, for plots): with ``field_times`` (a list of times) among the
parameters, write(t) records the mean concentration in every mesh element at the first output
time at or after each, as an array (nx, ny, nz) of the mesh (the grid times ``refine``), in
``self.fields``; close() saves them to <output_dir>/<row>_fields.npz (``times``, ``fields``,
``L`` the box size in m, ``grid`` and ``refine``).
"""
import os

import numpy as np

try:
    import ufl
    from mpi4py import MPI
    from petsc4py import PETSc
    from dolfinx import fem
    from dolfinx import mesh as dmesh
    import dolfinx.fem.petsc as fem_petsc
except ImportError as e:  # pragma: no cover - reported to the user
    raise ImportError(f'tissue_diffusion_FEniCS needs FEniCSx (dolfinx 0.8 or 0.9): '
                      f'mamba install -c conda-forge fenics-dolfinx=0.9 ({e})') from e


def _functionspace(mesh, element):
    make = getattr(fem, 'functionspace', None) or getattr(fem, 'FunctionSpace')
    return make(mesh, element)


def _petsc_vec(function):
    return function.x.petsc_vec if hasattr(function.x, 'petsc_vec') else function.vector


class DiffusionBox:
    """The class the module's api block names; see the module docstring."""

    def __init__(self, params, neighbours, comm=None, info=None):
        p = dict(params)
        self.info = dict(info or {})
        self.names = list(neighbours['C_t'])
        if list(neighbours.get('J_s', self.names)) != self.names:
            raise ValueError('C_t and J_s must connect to the same modules (one capillary_to_flux_port)')
        n_regions = len(self.names)
        # an mpi4py communicator if we have one (under mpiexec), else this process alone
        self.comm = comm if isinstance(comm, MPI.Comm) else MPI.COMM_SELF

        L = np.array([p['Lx'], p['Ly'], p['Lz']], dtype=float)
        grid = np.array([p['grid_nx'], p['grid_ny'], p['grid_nz']], dtype=int)
        self.fv = p.get('fv_scheme', 0.0) > 0.5
        refine = 1 if self.fv else max(1, int(round(p.get('refine', 2))))
        self.theta = float(p.get('theta', 0.5))
        self.M, self.C50 = float(p['M']), float(p['C50'])

        # solved in lengths scaled by the box's largest side: in metres the matrix entries are
        # ~1e-15, below PETSc's LU zero-pivot tolerance
        self.Lref = float(L.max())
        L = L / self.Lref
        self.mesh = dmesh.create_box(self.comm, [np.zeros(3), L], list(grid * refine), dmesh.CellType.hexahedron)
        self.V = _functionspace(self.mesh, ('DG', 0) if self.fv else ('Lagrange', 1))
        self.Q = _functionspace(self.mesh, ('DG', 0))

        # exchange regions: one grid cell each
        h = L / grid
        self.regions, self.chi, self.volumes = [], [], []
        for k in range(n_regions):
            idx = np.array([int((k + 0.5) * grid[0] / n_regions), grid[1] // 2, grid[2] // 2])
            lo, hi = idx * h, (idx + 1) * h
            chi = fem.Function(self.Q)
            chi.interpolate(lambda x, lo=lo, hi=hi: np.all(
                [(x[i] > lo[i] - 1e-12 * L[i]) & (x[i] < hi[i] + 1e-12 * L[i]) for i in range(3)], axis=0
            ).astype(float))
            self.regions.append((lo, hi))
            self.chi.append(chi)
            self.volumes.append(float(np.prod(hi - lo)))   # scaled; the physical volume is * Lref**3

        self.u_n = fem.Function(self.V)
        self.u_n.x.array[:] = float(p['C_init'])
        self.uh = fem.Function(self.V)
        self.src = fem.Function(self.Q)      # sum_k chi_k J_k / V_k, mM/s
        self.D = fem.Constant(self.mesh, PETSc.ScalarType(p['sigma_diff'] / self.Lref ** 2))
        self.dt = fem.Constant(self.mesh, PETSc.ScalarType(1e-3))
        u, v = ufl.TrialFunction(self.V), ufl.TestFunction(self.V)
        reaction = self.M * self.u_n / (self.u_n + self.C50)
        if self.fv:
            # two-point flux across each interior face, as tissue_diffusion_face
            X = _functionspace(self.mesh, ('DG', 0, (3,)))
            centres = fem.Function(X)
            centres.interpolate(lambda x: x)
            d = centres('+') - centres('-')
            dist = ufl.sqrt(ufl.dot(d, d))

            def diffusion(w, q):
                return self.D * (w('+') - w('-')) * (q('+') - q('-')) / dist * ufl.dS
        else:
            def diffusion(w, q):
                return self.D * ufl.inner(ufl.grad(w), ufl.grad(q)) * ufl.dx
        a = u * v * ufl.dx + self.theta * self.dt * diffusion(u, v)
        L_form = (self.u_n * v * ufl.dx - (1.0 - self.theta) * self.dt * diffusion(self.u_n, v)
                  + self.dt * (self.src + reaction) * v * ufl.dx)
        self.a_form, self.L_form = fem.form(a), fem.form(L_form)
        self.mean_forms = [fem.form(self.u_n * chi * ufl.dx) for chi in self.chi]
        self.total_form = fem.form(self.u_n * ufl.dx)
        self.b = fem_petsc.create_vector(self.L_form)
        self._solvers = {}
        self.history = []   # (t, region means, total amount) at the output times

        # field snapshots: element means, laid out (nx, ny, nz) on the mesh
        self.field_times = sorted(float(t) for t in (p.get('field_times') or []))
        self.fields = []     # (t, array (nx, ny, nz))
        self.grid, self.refine, self.L = grid, refine, L * self.Lref
        if self.field_times:
            n_local = self.Q.dofmap.index_map.size_local
            centres = self.Q.tabulate_dof_coordinates()[:n_local]
            n_mesh = grid * refine
            self._element_index = np.clip(np.floor(centres / (L / n_mesh)).astype(int), 0, n_mesh - 1)
            self._element_mean = fem.Function(self.Q)

    # --- helpers ----------------------------------------------------------------------------------
    def _solver(self, dt):
        key = round(dt, 15)
        if key not in self._solvers:
            self.dt.value = dt
            A = fem_petsc.assemble_matrix(self.a_form)
            A.assemble()
            ksp = PETSc.KSP().create(self.mesh.comm)
            ksp.setOperators(A)
            ksp.setType(PETSc.KSP.Type.PREONLY)
            ksp.getPC().setType(PETSc.PC.Type.LU)
            self._solvers[key] = (A, ksp)
        return self._solvers[key][1]

    def region_means(self):
        local = np.array([fem.assemble_scalar(f) for f in self.mean_forms])
        return self.comm.allreduce(local, op=MPI.SUM) / np.array(self.volumes)

    def total_amount(self):
        """sum C V over the box (mol), for checking conservation."""
        return self.comm.allreduce(fem.assemble_scalar(self.total_form), op=MPI.SUM) * self.Lref ** 3

    # --- the coupling contract ---------------------------------------------------------------------
    def initial_outputs(self):
        return {'C_t': self.region_means()}

    def step(self, t, dt, inputs):
        J = np.asarray(inputs['J_s'], dtype=float)
        self.src.x.array[:] = 0.0
        for chi, Jk, Vk in zip(self.chi, J, self.volumes):
            self.src.x.array[:] += chi.x.array * (Jk / (Vk * self.Lref ** 3))   # mol/s / m3 = mM/s
        ksp = self._solver(dt)
        self.dt.value = dt
        with self.b.localForm() as b_local:
            b_local.set(0.0)
        fem_petsc.assemble_vector(self.b, self.L_form)
        self.b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
        ksp.solve(self.b, _petsc_vec(self.uh))
        self.uh.x.scatter_forward()
        self.u_n.x.array[:] = self.uh.x.array
        return {'C_t': self.region_means()}

    def element_means(self):
        """The mean concentration in each mesh element, as an array (nx, ny, nz) of the mesh (on
        rank 0; None on the others). For Q1 the element mean of a trilinear field is its value
        at the element's centre."""
        self._element_mean.interpolate(self.u_n)
        n_local = self.Q.dofmap.index_map.size_local
        parts = self.comm.gather((self._element_index, self._element_mean.x.array[:n_local].copy()), root=0)
        if self.comm.rank != 0:
            return None
        out = np.full(tuple(self.grid * self.refine), np.nan)
        for idx, values in parts:
            out[idx[:, 0], idx[:, 1], idx[:, 2]] = values
        return out

    def write(self, t):
        self.history.append((t, self.region_means(), self.total_amount()))
        while self.field_times and t >= self.field_times[0] - 1e-12:
            self.field_times.pop(0)
            field = self.element_means()
            if field is not None:
                self.fields.append((t, field))

    def snapshot(self):
        self._saved = self.u_n.x.array.copy()

    def restore(self):
        self.u_n.x.array[:] = self._saved

    def close(self):
        out_dir = self.info.get('output_dir')
        if self.fields and out_dir and self.comm.rank == 0:
            os.makedirs(out_dir, exist_ok=True)
            np.savez(os.path.join(out_dir, f"{self.info.get('row', 'tissue')}_fields.npz"),
                     times=np.array([t for t, _ in self.fields]), fields=np.array([f for _, f in self.fields]),
                     L=self.L, grid=self.grid, refine=self.refine)
        for A, ksp in self._solvers.values():
            ksp.destroy()
            A.destroy()
        self._solvers = {}
