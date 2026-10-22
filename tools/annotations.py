"""
The biological identity of the bond-graph modules' species and reactions, written onto their CellDL documents'
elements by tools/celldl.py (BG-RDF annotations: bqbiol:is a ChEBI term and bqbiol:isPartOf a GO compartment on
each species, bqbiol:isVersionOf a Rhea reaction, or bqbiol:hasPart one for each reaction a lumped one combines).
cam_testing's bondgraph_reactions_test compares each annotated reaction with Rhea's.

ChEBI terms are the charge states Rhea uses (major microspecies at pH 7.3). Rhea terms are its undirected (master)
reactions; `python -m cam_testing.rhea find CHEBI:... ` lists candidates for a reaction's participants.
"""

# a species' label in the diagram (tools/celldl.py species_label) -> its ChEBI term
SPECIES_CHEBI = {
    'glc': 4167,        # D-glucose
    'ATP': 30616, 'ADP': 456216, 'AMP': 456215,
    'Pi': 43474,        # phosphate (hydrogenphosphate)
    'PPi': 33019,       # diphosphate(3-)
    'H': 15378, 'H2O': 15377, 'CO2': 16526, 'O2': 15379,
    'G6P': 61548,       # D-glucose 6-phosphate
    'G1P': 58601,       # alpha-D-glucose 1-phosphate
    'F6P': 57634,       # beta-D-fructose 6-phosphate
    'F16BP': 32966,     # beta-D-fructose 1,6-bisphosphate
    'DHAP': 57642,      # dihydroxyacetone phosphate
    'GAP': 59776,       # D-glyceraldehyde 3-phosphate
    'BPG13': 57604,     # (2R)-3-phospho-glyceroyl phosphate
    'PG3': 58272,       # (2R)-3-phosphoglycerate
    'PG2': 58289,       # (2R)-2-phosphoglycerate
    'PEP': 58702,       # phosphoenolpyruvate
    'Pyr': 15361,       # pyruvate
    'lac': 16651,       # (S)-lactate
    'NAD': 57540, 'NADH': 57945,
    'Gly3P': 57597,     # sn-glycerol 3-phosphate
    'glycerol': 17754,
    'ADPglc': 57498,    # ADP-alpha-D-glucose
    'Gly': 15444,       # glycogen, as Rhea's [(1->4)-alpha-D-glucosyl](n): its unit's ChEBI term
    'glyc': 15444,      # the lumped cytosol's glycogen
    'mx': 18246,        # the lumped cytosol's matrix, as Rhea's (1->4)-beta-D-glucan: its unit's ChEBI term
    'AcCoA': 57288, 'CoA': 57287,
    'Cit': 16947,       # citrate
    'Icit': 15562,      # D-threo-isocitrate
    'OG': 16810,        # 2-oxoglutarate
    'SucCoA': 57292,    # succinyl-CoA
    'Suc': 30031,       # succinate
    'Fum': 29806,       # fumarate
    'Mal': 15589,       # (S)-malate
    'OAA': 16452,       # oxaloacetate
    'Q': 16389,         # a ubiquinone
    'QH2': 17976,       # a ubiquinol
    'cytcox': 29034,    # Fe(III)-[cytochrome c]
    'cytcred': 29033,   # Fe(II)-[cytochrome c]
}

# a compartment (bgf:hasLocation) -> its GO cellular component
COMPARTMENT_GO = {
    'c': '0005829',     # cytosol
    'm': '0005759',     # mitochondrial matrix
    'p': '0005758',     # mitochondrial intermembrane space
}

BIOCHEMICAL, TRANSPORT = '0000176', '0000655'
# glycolysis's steps, PGI to LDH, which the lumped cytosol's glycolysis combines
_GLYCOLYSIS = ['11816', '16109', '14729', '18585', '10300', '14801', '15901', '10164', '18157', '23444']

# module -> reaction name -> {'rhea': id} or {'parts': [ids]} (a lumped reaction), optionally 'sbo'
REACTIONS = {
    'glycolysis': {
        'HK': {'rhea': '17825'}, 'PGI': {'rhea': '11816'}, 'PFK': {'rhea': '16109'}, 'ALD': {'rhea': '14729'},
        'TPI': {'rhea': '18585'}, 'GAPDH': {'rhea': '10300'}, 'PGK': {'rhea': '14801'}, 'PGM': {'rhea': '15901'},
        'ENO': {'rhea': '10164'}, 'PK': {'rhea': '18157'}, 'LDH': {'rhea': '23444'},
    },
    'gluconeogenesis': {
        'FBPASE': {'rhea': '11064'}, 'G3PDH': {'rhea': '11092'}, 'G6PASE': {'rhea': '16689'}, 'GK': {'rhea': '21644'},
    },
    'glycogen': {
        'PGM': {'rhea': '23536'},       # phosphoglucomutase: G6P = G1P
        'GLGC': {'rhea': '12120'}, 'GLGA': {'rhea': '18189'}, 'GLGP': {'rhea': '41732'}, 'PPA': {'rhea': '24576'},
    },
    'mitochondrion': {
        'PDH': {'rhea': '28042'}, 'CS': {'rhea': '16845'}, 'ACO': {'rhea': '10336'}, 'IDH': {'rhea': '23632'},
        'OGDH': {'rhea': '27786'}, 'SCS': {'rhea': '17661'}, 'CII': {'rhea': '13713'}, 'FUM': {'rhea': '12460'},
        'MDH': {'rhea': '21432'},
        'CI': {'rhea': '29091', 'sbo': TRANSPORT}, 'CIII': {'rhea': '11484', 'sbo': TRANSPORT},
        'CIV': {'rhea': '11436', 'sbo': TRANSPORT}, 'ATPS': {'rhea': '57720', 'sbo': TRANSPORT},
        'ANT': {'rhea': '34999', 'sbo': TRANSPORT}, 'PIC': {'rhea': '29939', 'sbo': TRANSPORT},
        'PYRC': {'rhea': '64720', 'sbo': TRANSPORT}, 'LEAK': {'rhea': '34979', 'sbo': TRANSPORT},
    },
    'energy_metabolism_adaptive': {
        'HK': {'rhea': '17825'},
        'HYD': {'rhea': '13065'},       # ATP hydrolysis
        'LG': {'parts': _GLYCOLYSIS},   # G6P to 2 lactate
        'MOB': {'parts': ['41732', '23536']},                  # glycogen phosphorylase, phosphoglucomutase
        'STORE': {'parts': ['23536', '12120', '18189', '24576']},   # to glycogen through ADP-glucose
        # matrix: to a glucan through UDP-glucose (UTP remade from ATP), and one ATP for the export
        'MX': {'parts': ['23536', '19889', '24576', '19929', '25098', '13065']},
    },
}
IDENTIFIERS = 'https://identifiers.org/'


def species_terms(label, compartment):
    '''The BG-RDF lines annotating a species element.'''
    lines = []
    if label in SPECIES_CHEBI:
        lines.append(f'bqbiol:is <{IDENTIFIERS}CHEBI:{SPECIES_CHEBI[label]}>')
    if compartment in COMPARTMENT_GO:
        lines.append(f'bqbiol:isPartOf <{IDENTIFIERS}GO:{COMPARTMENT_GO[compartment]}>')
    return lines


def reaction_terms(module, name):
    '''The BG-RDF lines annotating a reaction element of ``module``.'''
    terms = REACTIONS.get(module, {}).get(name)
    if not terms:
        return []
    lines = [f'bqbiol:isVersionOf <{IDENTIFIERS}rhea:{terms["rhea"]}>'] if 'rhea' in terms else []
    lines += [f'bqbiol:hasPart <{IDENTIFIERS}rhea:{p}>' for p in terms.get('parts', [])]
    if terms.get('sbo'):
        lines.append(f'bqbiol:isVersionOf <{IDENTIFIERS}SBO:{terms["sbo"]}>')
    return lines
