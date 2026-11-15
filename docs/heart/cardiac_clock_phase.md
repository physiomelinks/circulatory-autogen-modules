# The cardiac clocks' phase, computed from the beat phase

The heart versions in this library (the flat hearts `vp`, `vp_Ca`, `vp_LVprop`, `vp_devel`, `vp_new_valve`,
`vp_wCont`, `vp_wCont_ASD`, `vp_wCont_nonstiff` and the chambers `Liang2009_v01`, `Liang2009NoVisco_v01`,
`Liang2009NoVisco_v02` that the heart supermodules `Argus2026_v01`, `controlled_Argus2026_v01` and
`nonstiff_Argus2026_v01` are built from) drive each chamber's activation from a phase $\psi\in[0,1)$. This
page derives the form they use, $\psi = P(s)$ computed from the beat phase $s$, from the form they used
before, a phase state $\chi$ integrated with a piecewise rate. The two give the same trajectories; the
second has no state that waits on a switch, which is what made about one run in ten start an atrium
35-76 ms early or late.

`tests/test_cardiac_clock_phase.py` checks every statement here that the code depends on: the closed
form against the integrated law, the output $\chi$ against the integrated state, the two forms of each
version against each other in a simulation, and the onset's insensitivity to the solver.

Notation: in the CellML, the atrium's clock is `chi_a`, `chi_afloor`, `chi_afloor_final`, `t_astart`,
`T_ac`, `T_ar`; the ventricle's `chi_v`, `chi_vfloor`, `chi_vfloor_final`, `t_vstart`, `T_vc`, `T_vr`;
a chamber's `chi`, `chi_floor`, `chi_floor_final`, `t_start`, `T_c`, `T_r`. Below, one clock is written
$\chi$, $\psi$ (`chi_floor`), $t_s$, $T_c$, $T_r$ for any of them.

## 1. The original law

The beat phase $s$ is a state with

$$\frac{ds}{dt} = \frac{1}{T_w},\qquad s(0) = 0,\qquad m = s - \lfloor s\rfloor \quad(\texttt{mt}),$$

where $T_w$ is the instantaneous period: $T_w = T$ in the constant-period hearts (`vp`, `vp_Ca`, ...,
`cardiac_clock/Liang2009_v01`), and $T_w = \max(T + \Delta T, T_{\min})$ (`T_wCont`, `T_eff`) in the
controlled ones (`vp_wCont`, `vp_wCont_ASD`, `vp_wCont_nonstiff`, `cardiac_clock/controlled`), where a
controller sets $\Delta T$ at any time.

Each chamber's phase $\chi$ was a state, $\chi(0) = 0$, with the wrapped phase
$\psi = \chi - \lfloor\chi\rfloor$ and the rate

$$
\frac{d\chi}{dt} = \begin{cases}
\dfrac{1}{4T_{c,w}} & t_s' \le m \le t_s' + \varepsilon_1 \ \text{ and }\ \psi \le \tfrac14 & \text{(start window)}\\[6pt]
\dfrac{1}{4T_{c,w}} & t_s' + \varepsilon_1 > 1 \ \text{ and }\ m \le t_s' + \varepsilon_1 - 1 \ \text{ and }\ \psi \le \tfrac14 & \text{(start window, wrapped)}\\[6pt]
\dfrac{1}{4T_{c,w}} & \varepsilon_2 < \psi \le \tfrac14 & \text{(contraction)}\\[6pt]
\dfrac{1}{4T_{r,w}} & \tfrac14 \le \psi < \tfrac12 \quad(\chi_v:\ \tfrac14 < \psi < \tfrac12) & \text{(relaxation)}\\[6pt]
\dfrac{1}{2\,(T_w - T_{c,w} - T_{r,w})} & \psi \ge \tfrac12 \quad(\chi_v:\ \psi > \tfrac12) & \text{(rest)}\\[6pt]
0 & \text{otherwise} & \text{(wait)}
\end{cases}
$$

with $t_s' = t_s/T$ (`t_astart_norm`, `t_vstart_norm`, `t_start_norm`), $\varepsilon_1 = 0.07$,
$\varepsilon_2 = 0.02$, and the durations scaled with the period:

$$T_{c,w} = \frac{T_c\,T_w}{T},\qquad T_{r,w} = \frac{T_r\,T_w}{T}$$

(`T_ac_wCont`, `T_c_eff`, ...; in the constant-period hearts $T_w = T$ and $T_{c,w} = T_c$). The
pieces are taken in order: the first whose condition holds gives the rate.

Boundary conventions. The atrial clock and the chambers use $\ge \tfrac14$ for relaxation and
$\ge \tfrac12$ for rest; the ventricular clock `chi_v` of the flat hearts uses strict $>$ at both. At
$\psi = \tfrac14$ both forms are in the contraction piece ($\psi \le \tfrac14$ comes first). At
$\psi = \tfrac12$ the strict form falls to "otherwise": a rate of 0 at that single point (the review of
2026-10, `reviews/heart_review.yaml`, noted it). A solution passes through such a point without
stopping unless it lands on it exactly, so it changes no trajectory; section 3 shows the closed form is
continuous there, so the closed form takes the same value whichever piece holds at the point.

The activation follows from $\psi$ as before, in both forms:

$$\psi_f = \begin{cases} 2\psi & \psi \le \tfrac12\\ 0 & \text{otherwise}\end{cases}\quad(\texttt{chi\_floor\_final}),\qquad
e = \tfrac12\bigl(1 - \cos 2\pi\psi_f\bigr).$$

What the window and the wait do: $\chi$ rests at an integer (the wait, rate 0, for $\psi \le \varepsilon_2$
outside the window) until $m$ reaches $t_s'$, then contracts. The dead band $\psi\le\varepsilon_2$ lets
$\chi$ leave the integer only inside the window.

## 2. The rate in the beat phase does not depend on the period

Every piece of the rate is a constant divided by $T_w$:

$$\frac{1}{4T_{c,w}} = \frac{T}{4T_c}\cdot\frac{1}{T_w},\qquad
\frac{1}{4T_{r,w}} = \frac{T}{4T_r}\cdot\frac{1}{T_w},\qquad
\frac{1}{2(T_w - T_{c,w} - T_{r,w})} = \frac{T}{2(T - T_c - T_r)}\cdot\frac{1}{T_w}.$$

Dividing by $ds/dt = 1/T_w$, wherever $ds/dt > 0$,

$$\frac{d\chi}{ds} = g(\psi, m),\qquad
g = \begin{cases}
\dfrac{T}{4T_c} & \text{start window or contraction}\\[4pt]
\dfrac{T}{4T_r} & \text{relaxation}\\[4pt]
\dfrac{T}{2(T - T_c - T_r)} & \text{rest}\\[4pt]
0 & \text{wait}
\end{cases}$$

with the same conditions as above, and $m = s - \lfloor s\rfloor$ a function of $s$. $T_w$ has gone:
$\chi$ as a function of $s$ is the solution of an equation that contains no $T_w$, so it is the same
function of $s$ whatever $T_w(t)$ does — constant, controlled, or changed in steps at any moment of the
beat. Only the map from $t$ to $s$ depends on the period.

This needs the durations to scale with the period ($T_{c,w}/T_w = T_c/T$). It fails for absolute
durations (section 5).

## 3. The closed form

Integrate $g$ in $s$ from $\chi(0) = 0$, $s(0) = 0$. Write

$$p = s - t_s',\qquad \phi = p - \lfloor p \rfloor \in [0, 1),\qquad a = \frac{T_c}{T},\qquad b = \frac{T_c + T_r}{T}.$$

**Before the first onset** ($s < t_s'$, possible when $t_s' > 0$): $\psi = 0 \le \varepsilon_2$ and
$m = s < t_s'$, so the window is closed and the rate is 0 (the wait): $\chi = 0$.

**First beat**, from $s = t_s'$: the window opens ($m = t_s'$), so $\chi$ contracts at
$g = T/(4T_c)$. It leaves the dead band after $\Delta s = 4\varepsilon_2 a$, inside the window when
$4\varepsilon_2 T_c/T < \varepsilon_1$ (0.0136 for the atrium's $T_c/T = 0.17$, 0.027 for the ventricle's
0.34, against 0.07; the condition holds for $T_c < 0.875\,T$). Then:

| piece | $\psi$ from → to | rate $g$ | $\Delta s$ taken | $\phi$ at its end |
|---|---|---|---|---|
| contraction | $0 \to \tfrac14$ | $T/(4T_c)$ | $\tfrac14 \cdot 4T_c/T = a$ | $a$ |
| relaxation | $\tfrac14 \to \tfrac12$ | $T/(4T_r)$ | $\tfrac14 \cdot 4T_r/T = b - a$ | $b$ |
| rest | $\tfrac12 \to 1$ | $T/(2(T - T_c - T_r))$ | $\tfrac12 \cdot 2(T - T_c - T_r)/T = 1 - b$ | $1$ |

So $\psi$ is linear in $\phi$ on each piece:

$$
P(s) = \begin{cases}
0 & s < t_s'\\[4pt]
\dfrac{\phi}{4a} = \dfrac{T\phi}{4T_c} & T\phi \le T_c\\[8pt]
\dfrac14 + \dfrac{\phi - a}{4(b - a)} = \dfrac14 + \dfrac{T\phi - T_c}{4T_r} & T_c < T\phi < T_c + T_r\\[8pt]
\dfrac12 + \dfrac{\phi - b}{2(1 - b)} = \dfrac12 + \dfrac{T\phi - T_c - T_r}{2(T - T_c - T_r)} & \text{otherwise}
\end{cases}
$$

$P$ is continuous: at $T\phi = T_c$ the first two pieces give $\tfrac14$, at $T\phi = T_c + T_r$ the last
two give $\tfrac12$. So whichever piece a comparison puts the boundary point in, the value is the same;
the closed form is the same for the strict ($\chi_v$) and non-strict comparisons.

**One cycle of $\chi$ is one beat.** The three pieces take $a + (b - a) + (1 - b) = 1$ in $s$: $\chi$
reaches its next integer exactly at $s = t_s' + 1$, when $\phi$ wraps to 0, and that is exactly when
$m = t_s'$ again, the moment the window opens. So from the first onset on, $\chi$ never waits: arriving
at the integer, it finds the window open and starts the next contraction at once. The start window and
the wait act only from the initial state, before the first onset. (In a numerical solution $\chi$
drifts against $s$ and arrives slightly early or late; the window and the wait absorb that drift. The
closed form has no drift to absorb.)

**The activation.** Since $e$ is computed from $\psi$, $e = \tfrac12(1 - \cos 2\pi\psi_f(P(s)))$ is a
function of $s$ alone. It is $C^1$ in $s$: $e$ and $de/ds$ are 0 on both sides of the onset and at the end
of relaxation, and $de/ds = 0$ at the peak ($\psi = \tfrac14$); its second derivative jumps at the four
phase boundaries. The CellML computes $\phi$ from $s$ with one `floor`, which $s$ crosses at rate
$1/T_w$, never waiting on it.

In the CellML (each chamber, and each of the two clocks of a flat heart):

```
phi       = (s - t_start_norm) - floor(s - t_start_norm)
chi_floor = 0                                                      if s < t_start_norm
            0.25 phi T / T_c                                       if phi T <= T_c
            0.25 + 0.25 (phi T - T_c) / T_r                        if phi T <  T_c + T_r
            0.5 + 0.5 (phi T - T_c - T_r) / (T - T_c - T_r)        otherwise
chi       = 0                                                      if s < t_start_norm
            floor(s - t_start_norm) + chi_floor                    otherwise
```

The chambers read $s$ from the cardiac clock (its `cardiac_clock_port` carries `mt`, `T_eff` and `s`).
$T_{c,w}$, $T_{r,w}$ (`T_ac_wCont`, `T_c_eff`, ...) and `mt` are still computed, as outputs.
$\varepsilon_1$ and $\varepsilon_2$ are gone: they only acted in the wait.

## 4. The output $\chi$

$\chi$ stays an output, so the verification output maps and plots that name it (`heart/chi_a`,
`heart/ra/chi`, ...) keep working:

$$\chi = \begin{cases} 0 & s < t_s'\\ \lfloor s - t_s'\rfloor + P(s) & \text{otherwise.}\end{cases}$$

Its integer part counts the onsets: it goes up by one at each onset $s = k + t_s'$, as the integrated
state's did when $\psi$ wrapped from 1 to 0. So for the initial state every version had,
$\chi(0) = 0$ and $s(0) = 0$, $\chi$ equals the old state at every time: 0 up to the first onset, then
$\lfloor p\rfloor + P$. (For the ventricle, $t_s' = 0$ and this is $\lfloor s\rfloor + P(s)$.) `chi_floor`
is the wrapped phase $P(s)$, as before.

When it differs. The old state could start anywhere; the new output cannot:
- an initial $\chi$ with another integer part (a model warmed up for $n$ beats and restarted with
  $\chi = n + \psi_0$ but $s$ reset, or a ventricle stored some cycles behind its atrium): the output
  follows $\lfloor s - t_s' \rfloor$, not the stored integer. Only the integer part differs; `chi_floor`
  and $e$ are unaffected;
- an initial phase that is not $P(s_0)$ (a stored state with $\chi$ lagging or leading $s$, or $s(0)$
  set without $\chi$): the old law ran that lag out over the first beat, through the wait (a lagging
  $\chi$) or the rest piece (a leading one); the new form has no lag, so the first beat differs until the
  old state would have re-synchronised. From the next onset on they agree;
- $s(0) > t_s'$ with $\chi(0) = 0$: the old law waited for the next window, the new form starts at
  $P(s_0)$.

No version, parameterisation or system model in this library sets an initial $\chi$ or $s$ other than 0.

## 5. When it holds

The derivation uses only:
1. durations proportional to the period, $T_{c,w} = T_c T_w/T$ and $T_{r,w} = T_r T_w/T$, which makes
   $d\chi/ds$ independent of $T_w$ (section 2);
2. $ds/dt = 1/T_w > 0$;
3. $T_c + T_r < T$, and the first contraction leaving the dead band inside the window,
   $4\varepsilon_2 T_c/T < \varepsilon_1$;
4. the initial state $\chi(0) = 0$, $s(0) = 0$ (section 4).

Every version listed at the top meets all four, so the change is in place: the same versions, the same
parameters, the same trajectories.

It fails with absolute durations under a varying period: contraction and relaxation lasting $T_c$, $T_r$
seconds whatever $T_w$ is (the 2022 design, circulatory_autogen `ec32e8f3`, "allows dynamic changes in
period and contractile and relaxation period"). Then $d\chi/ds = T_w/(4T_c)$ depends on $T_w$, $\chi$ is not
a function of $s$, and the wait does real work: it holds $\chi$ at the integer until the window. In the
discrete-period scenarios below, that form ("chiabs") keeps its own design but its phase leaves $P(s)$
by up to 0.5 (half a cycle) as soon as the period changes. Such a version would need the state (or event
handling); it cannot use this closed form.

## 6. Why, and the evidence

**The fault.** The integrated law switches on $\lfloor\chi\rfloor$ of a state that waits at an integer.
The atrium waits at exactly $\chi_a = 0$ for the first 0.836 s ($t_s' = 0.8$ of a 1.045 s beat), then the
rate jumps from 0 to 1.37 s⁻¹. A Newton iterate a hair below the integer ($\chi_a = -9.4\times10^{-31}$)
reads $\psi \approx 1$, the rest piece, so the finite-difference Jacobian taken there spans the jump:
$J_{\chi\chi} = -9.0\times10^{17}$. CVODE's Newton correction for $\chi_a$ is then about zero and it accepts
the predictor's extrapolation, $\chi_a = 0.093$ one step later against an exact $5.6\times10^{-4}$. Round-off
decides which runs do this, so it happens at any tolerance, in libOpenCOR, Myokit and a C replica alike:
the atrium of the next beat starts 34-63 ms early (healed at the next onset, since $\chi_a$ waits for the
window) or 60 ms late (carried into the next beat), and with it the mitral pulse. Source: the 2026-11-15
investigation of `control_parasymp` (heart `vp_wCont`), C replica at tol 1e-10, RHS log
`valve-flip/runs/rhslog_1e-10.txt` (line 695); report `valve-flip/report.md`.

**The sensitivity sweep** (2 s of `control_parasymp`, libCellML C code, CVODE 7.7.0 BDF with a dense
difference-quotient Jacobian; 25 tolerances $10^{-10}$-$10^{-7}$ with rtol = atol, each as is, with the
initial states $\times(1 + 10^{-12})$ and with a maximum step of 0.999 ms; and 20 s runs at 9 tolerances
as is and perturbed, with $T_w$ moving 1.00-1.47 s under parasympathetic control):

| form | bad runs, 2 s (atrial onset off by 35-76 ms) | solver failures | largest onset shift | 20 s runs with the clock > 1 ms off |
|---|---|---|---|---|
| integrated $\chi$ | 7 / 75 | 1 | 76.3 ms | 3 / 18 (77, 63, 45 ms) |
| closed form $P(s)$ | 0 / 75 | 0 | 0.0065 ms | 0 / 18 |

The closed form's median run differs from the reference (the current model by CasADi CVODES at
$10^{-12}$/$10^{-14}$) by $4.9\times10^{-7}$ in $e_a$ and $6.0\times10^{-7}$ (relative) in left-atrial
volume, against $4.0\times10^{-6}$ and $5.2\times10^{-7}$ for the integrated form's own median. Source:
`clock-fix/runs/summary_base_A_C_AC_B_D_E8_E16.json`, `clock-fix/runs/long_table.txt`; scripts
`clock-fix/sens.sh`, `long.sh`, `analyse.py`, `analyse_long.py`; option (b) in `clock-fix/h/patch.py`.

**Discrete period changes** (`control_parasymp` with $T_w$ overridden: a step 1.0 → 0.6 s or
0.6 → 1.2 s at eight points of the beat — during the wait, just before the first onset, just before a
later onset, inside the start window, in atrial contraction, relaxation and rest, and in ventricular
systole — a 0.3 s ramp, a period held per beat from a sequence or from the controller, and the
controller as written; 20 scenarios). Criteria fixed before the runs: every onset after the first
within 0.5 ms of $s = k + t_s'$, one onset per beat, $\psi$ within $10^{-3}$ of $P(s)$, no jump in $e_a$.

| form | scenarios passing | largest onset error | largest $\lvert\psi - P(s)\rvert$ |
|---|---|---|---|
| integrated $\chi$ (tol $3\times10^{-10}$ and $10^{-9}$) | 20 / 20 | 0.0033 ms | $1.1\times10^{-6}$ |
| closed form $P(s)$ | 20 / 20 | 0.0035 ms | 0 |
| explicit modes with root-finding and reset (reference) | 20 / 20 | 0.0035 ms | $1.3\times10^{-12}$ |
| absolute durations ("chiabs") | 11 / 20 against its own design | — | 0.50 |

So under any period change both forms put each phase at the same fraction of the instantaneous beat.
Source: `clock-fix/test_discrete_period.py`, results `clock-fix/runs/T/results.json`.

**The library versions.** `tests/test_cardiac_clock_phase.py` simulates four version families in both
forms (the converter in the test rebuilds the original law from a version's CellML, or the closed form from
the original; its round trip reproduces every version file) with CVODE at rtol = atol = $10^{-12}$ for 5.5 s,
with the period stepped three times (in atrial contraction, atrial relaxation and ventricular systole) for the
controlled ones. Largest difference over the states and activations, relative to each one's maximum:

| family (version) | largest difference | in $e$ (activations) |
|---|---|---|
| flat, constant period (`heart/vp`) | $1.0\times10^{-7}$ (`zeta_miv_pre`) | $\le 8.6\times10^{-9}$ |
| flat, controlled, period steps (`heart/vp_wCont`) | $1.6\times10^{-7}$ (`zeta_miv_pre`) | $\le 1.7\times10^{-8}$ |
| chamber in its clock harness (`chamber/Liang2009NoVisco_v02`) | $2.9\times10^{-10}$ ($e$) | $2.9\times10^{-10}$ |
| supermodule, controlled, period steps (`heart/controlled_Argus2026_v01`) | $5.5\times10^{-8}$ (`puv/zeta_pre`) | $\le 8.1\times10^{-9}$ |

At $10^{-10}$ the integrated form's own error after a period step is larger ($e_v$ $9\times10^{-6}$,
`zeta_miv_pre` $10^{-4}$ between $10^{-10}$ and $10^{-12}$), while the closed form moves by $3\times10^{-12}$:
the closed form is the more accurate of the two at a given tolerance. Just before the first atrial onset, a
change of $-10^{-12}$ in the integrated `chi_a` (or a chamber's `chi`) moves its rate by $7.6\times10^{11}$
times the change; in the closed form no rate moves by more than $10^3$ times it. In the library's V&V suite
the change turns the `stability_test` of `heart/vp`, `vp_Ca` and `vp_wCont` from failing (CVODE at rtol
$10^{-6}$ against the reference: 0.013) to passing. Myokit on a `-march=native` SUNDIALS 7.7.0 build did not
show the onset shift in 500 runs of `heart/vp_wCont`'s test network (tolerances $10^{-10}$-$10^{-6}$, as is
and with the initial states perturbed by $\pm10^{-12}$), which is why the tolerance sweep in the test is an
expected failure that may pass for the integrated form.

All `clock-fix/` and `valve-flip/` paths are the investigation's scratch directories
(`~/.cache/claude-tmp/` on the machine it ran on); the options page that compared five fixes is
"Atrial clock fix options" (2026-11-15).
