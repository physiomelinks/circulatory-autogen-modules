"""References for the extracellular matrix: its cost to the cell (tools/build_adaptive_cell.py, energy_metabolism) and
the space between cells (tools/build_extracellular_space.py). Research notes: AdaptOrDeity
docs/research/food-uptake-and-extracellular-space.md, "Cost of making the matrix (EPS)"."""

MATRIX_BIB = {
    'sutherland2001biofilm': dict(
        type='article', author='Sutherland, Ian W.', title='Biofilm exopolysaccharides: a strong and sticky framework',
        journal='Microbiology', year='2001', volume='147', number='1', pages='3--9', doi='10.1099/00221287-147-1-3',
        note='Synthase-dependent and Wzx/Wzy exopolysaccharide pathways from nucleotide sugars; no ATP per gram given'),
    'welman2003exopolysaccharides': dict(
        type='article', author='Welman, Alan D. and Maddox, Ian S.',
        title='Exopolysaccharides from lactic acid bacteria: perspectives and challenges', journal='Trends in Biotechnology',
        year='2003', volume='21', number='6', pages='269--274', doi='10.1016/S0167-7799(03)00107-0',
        note='Lactic acid bacteria make tens to hundreds of mg/L exopolysaccharide, a few percent of the sugar fermented'),
    'zhang2003biodegradability': dict(
        type='article', author='Zhang, Xiaoqi and Bishop, Paul L.',
        title='Biodegradability of biofilm extracellular polymeric substances', journal='Chemosphere', year='2003',
        volume='50', number='1', pages='63--69', doi='10.1016/S0045-6535(02)00319-3',
        note='EPS is biodegradable and serves as a substrate when nutrients run short'),
    'flemming2010biofilm': dict(
        type='article', author='Flemming, Hans-Curt and Wingender, Jost', title='The biofilm matrix',
        journal='Nature Reviews Microbiology', year='2010', volume='8', number='9', pages='623--633',
        doi='10.1038/nrmicro2415', note='The matrix is often more than 90% of a biofilm\'s dry mass'),
    'stewart2003diffusion': dict(
        type='article', author='Stewart, Philip S.', title='Diffusion in biofilms', journal='Journal of Bacteriology',
        year='2003', volume='185', number='5', pages='1485--1491', doi='10.1128/JB.185.5.1485-1491.2003',
        note='Table 1: diffusion coefficients in water at 25 C (O2 2.0e-9, glucose 6.7e-10 m2/s); in biofilms De/Daq 0.6 for '
             'light gases and 0.25 for most organic solutes'),
    'mackie1955diffusion': dict(
        type='article', author='Mackie, J. S. and Meares, P.',
        title='The diffusion of electrolytes in a cation-exchange resin membrane I. Theoretical',
        journal='Proceedings of the Royal Society of London A', year='1955', volume='232', number='1191',
        pages='498--509', doi='10.1098/rspa.1955.0234',
        note='Obstruction by polymer of volume fraction phi: D/D0 = ((1 - phi)/(1 + phi))^2'),
    'amsden1998solute': dict(
        type='article', author='Amsden, Brian', title='Solute diffusion within hydrogels. Mechanisms and models',
        journal='Macromolecules', year='1998', volume='31', number='23', pages='8382--8395', doi='10.1021/ma980765f',
        note='Review of obstruction, free-volume and hydrodynamic models of diffusion in hydrogels'),
}
