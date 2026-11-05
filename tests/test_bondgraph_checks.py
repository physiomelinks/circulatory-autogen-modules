"""The bond-graph checks (cam_testing/bondgraph.py): a version's CellDL document opens in the editor, and its
annotated reactions have Rhea's participants and stoichiometry. Documents are built here; Rhea's reactions come
from the package's cache (cam_testing/data/rhea)."""
import os

import pytest

from cam_testing import bondgraph, checks, rhea

pytest.importorskip('rdflib')

CHEBI = {'Glc': 4167, 'ATP': 30616, 'G6P': 61548, 'ADP': 456216, 'AMP': 456215, 'H': 15378, 'H2O': 15377,
         'Pi': 43474, 'Na': 29101, 'K': 29103, 'ADPglc': 57498, 'Gly': 15444}
PREFIXES = '''@prefix : <#> .
@prefix bgf: <https://bg-rdf.org/ontologies/bondgraph-framework#> .
@prefix celldl: <http://celldl.org/ontologies/celldl#> .
@prefix bqbiol: <http://biomodels.net/biology-qualifiers/> .
'''


class Version:
    '''What the checks use of a library version: its files.'''
    def __init__(self, directory, stem='mod_v01'):
        self.dir, self.stem = str(directory), stem
        self.config, self.spec, self.is_supermodule = {}, {}, False

    def path(self, suffix):
        return os.path.join(self.dir, f'{self.stem}_{suffix}')


def species(name, location='c', cls='ChemicalPotential', chebi=True):
    return (f'{name}_{location}', cls, name, location, chebi)


def reaction(name, rhea_ids=(), parts=()):
    return (name, 'ChemicalReactionNode', name, 'c', (tuple(rhea_ids), tuple(parts)))


def write(directory, elements, bonds, stem='mod_v01', drawn_extra=(), undrawn=(), junctions=()):
    '''A CellDL document of ``elements`` (from species / reaction) and ``bonds`` ((source, target, count)).'''
    turtle = [PREFIXES, f'<> a celldl:Document, bgf:BondgraphModel ; bgf:hasBondElement '
              + ', '.join(f':{e[0]}' for e in elements) + ' .']
    if junctions:
        turtle.append('<> bgf:hasJunctionStructure ' + ', '.join(f':{j}' for j in junctions) + ' .')
    for j in junctions:
        turtle.append(f':{j} a celldl:Component, bgf:ZeroNode .')
    for eid, cls, name, location, notes in elements:
        lines = [f':{eid} a celldl:Component, bgf:{cls} ; bgf:hasSpecies "{name}" ; bgf:hasLocation "{location}"']
        if cls == 'ChemicalReactionNode':
            rhea_ids, parts = notes
            lines += [f'bqbiol:isVersionOf <https://identifiers.org/rhea:{r}>' for r in rhea_ids]
            lines += [f'bqbiol:hasPart <https://identifiers.org/rhea:{r}>' for r in parts]
        elif notes:
            lines += [f'bqbiol:is <https://identifiers.org/CHEBI:{CHEBI[name]}>',
                      'bqbiol:isPartOf <https://identifiers.org/GO:0005829>']
        turtle.append(' ;\n  '.join(lines) + ' .')
    for source, target, count in bonds:
        turtle.append(f':{source}.{target} a celldl:Connection ; celldl:hasSource :{source} ; celldl:hasTarget :{target} ;'
                      f' bgf:hasSource :{source} ; bgf:hasTarget :{target}' + (f' ; bgf:bondCount {count}' if count != 1 else '') + ' .')
    drawn = [e[0] for e in elements if e[0] not in undrawn] + list(junctions) + list(drawn_extra)
    svg = ('<svg xmlns="http://www.w3.org/2000/svg"><metadata id="celldl-rdf-metadata" data-content-type="text/turtle">'
           '<![CDATA[' + '\n'.join(turtle) + ']]></metadata><g id="celldl-diagram-layer">'
           + ''.join(f'<g id="{d}" class="celldl-Component"/>' for d in drawn)
           + ''.join(f'<path id="{s}.{t}" class="celldl-Connection"/>' for s, t, _ in bonds) + '</g></svg>')
    version = Version(directory, stem)
    with open(version.path('bondgraph.svg'), 'w') as f:
        f.write(svg)
    return version


def hexokinase(directory, rhea_ids=('17825',), substrates=(('Glc', 1), ('ATP', 1)),
               products=(('G6P', 1), ('ADP', 1), ('H', 1)), chebi=True, **kw):
    elements = [reaction('Re_HK', rhea_ids)] + [species(n, chebi=chebi) for n, _ in substrates + products]
    bonds = [(f'{n}_c', 'Re_HK', c) for n, c in substrates] + [('Re_HK', f'{n}_c', c) for n, c in products]
    return write(directory, elements, bonds, **kw)


def only(result):
    assert len(result.details) == 1, result.details
    return result.details[0]


# ---- not a bond graph ---------------------------------------------------------------------------

def test_a_version_without_a_bond_graph_is_not_applicable(tmp_path):
    version = Version(tmp_path)
    for check in (bondgraph.celldl_check, bondgraph.reactions_check):
        result = check(version)
        assert (result.status, result.message) == (checks.NOT_APPLICABLE, bondgraph.NOT_A_BONDGRAPH)


def test_the_report_shows_both_columns_as_not_applicable_without_a_bond_graph(tmp_path):
    from cam_testing import report
    assert set(bondgraph.BONDGRAPH_TESTS) <= set(report.REPORT_TESTS)
    for test in bondgraph.BONDGRAPH_TESTS:
        assert report.not_applicable_reason(Version(tmp_path), test) == bondgraph.NOT_A_BONDGRAPH
        assert report.TEST_SHORT[test] and report.TEST_ABOUT[test]
    hexokinase(tmp_path)
    assert report.not_applicable_reason(Version(tmp_path), 'bondgraph_celldl_test') is None


# ---- opens in CellDL ----------------------------------------------------------------------------

def test_a_document_the_editor_opens_passes(tmp_path, monkeypatch):
    monkeypatch.delenv('BG_TOOLS', raising=False)
    result = bondgraph.celldl_check(hexokinase(tmp_path))
    assert result.status == checks.PASSED, result.details


def test_what_keeps_the_editor_from_opening_a_document_fails(tmp_path, monkeypatch):
    monkeypatch.delenv('BG_TOOLS', raising=False)
    version = hexokinase(tmp_path, undrawn=('ATP_c',), drawn_extra=('stray',))
    problems = bondgraph.celldl_check(version).details
    assert any('ATP_c is not drawn' in p for p in problems)
    assert any('stray is drawn but not in the RDF' in p for p in problems)
    # the .ttl beside it must be the same graph
    os.makedirs(tmp_path / 'b')
    version = hexokinase(tmp_path / 'b')
    with open(version.path('bondgraph.ttl'), 'w') as f:
        f.write(PREFIXES + '<> a bgf:BondgraphModel .\n')
    result = bondgraph.celldl_check(version)
    assert result.status == checks.FAILED and 'bondgraph.ttl is not the document' in result.details[0]
    # not CellDL at all
    with open(version.path('bondgraph.svg'), 'w') as f:
        f.write('<svg xmlns="http://www.w3.org/2000/svg"/>')
    assert 'not a CellDL document' in bondgraph.celldl_check(version).message


# ---- reactions against Rhea ----------------------------------------------------------------------

def test_an_annotated_reaction_with_rheas_stoichiometry_matches(tmp_path):
    result = bondgraph.reactions_check(hexokinase(tmp_path))
    assert result.status == checks.PASSED, result.details
    assert 'matches (as written)' in only(result)
    # written backwards, every coefficient doubled: still Rhea's reaction
    os.makedirs(tmp_path / 'b')
    result = bondgraph.reactions_check(hexokinase(tmp_path / 'b', substrates=(('G6P', 2), ('ADP', 2), ('H', 2)),
                                                  products=(('Glc', 2), ('ATP', 2))))
    assert result.status == checks.PASSED and 'reversed, scaled by 2' in only(result)


def test_a_difference_only_in_h_plus_or_water_is_a_warning(tmp_path):
    result = bondgraph.reactions_check(hexokinase(tmp_path, products=(('G6P', 1), ('ADP', 1))))
    assert result.status == checks.PASSED
    assert result.metrics['warn'] == 1 and 'H+ (c): 0 where Rhea has 1' in only(result)


def test_a_wrong_or_missing_participant_fails(tmp_path):
    result = bondgraph.reactions_check(hexokinase(tmp_path, products=(('G6P', 1), ('AMP', 1), ('H', 1))))
    assert result.status == checks.FAILED
    detail = only(result)
    assert 'missing ADP (CHEBI:456216, c)' in detail and 'extra AMP (CHEBI:456215, c)' in detail
    os.makedirs(tmp_path / 'b')
    result = bondgraph.reactions_check(hexokinase(tmp_path / 'b', substrates=(('Glc', 2), ('ATP', 1))))
    assert result.status == checks.FAILED and 'D-glucose (CHEBI:4167, c): 2 where Rhea has 1' in only(result)


def test_unannotated_reactions_leave_it_pending_and_unannotated_species_fail(tmp_path):
    result = bondgraph.reactions_check(hexokinase(tmp_path, rhea_ids=()))
    assert result.status == checks.PENDING and only(result) == 'Re_HK: not annotated'
    os.makedirs(tmp_path / 'b')
    result = bondgraph.reactions_check(hexokinase(tmp_path / 'b', chebi=False))
    assert result.status == checks.FAILED and 'no ChEBI term on' in only(result)


def test_a_reaction_not_in_the_cache_fails_with_the_command_that_fetches_it(tmp_path):
    result = bondgraph.reactions_check(hexokinase(tmp_path), load=lambda rid: None)
    assert result.status == checks.FAILED and 'python -m cam_testing.rhea fetch 17825' in only(result)


def test_a_transport_reaction_maps_rheas_in_and_out_to_its_compartments(tmp_path):
    # Na+/K+-ATPase: 3 Na+ out of the cell and 2 K+ in, for an ATP hydrolysed inside
    elements = [reaction('Re_NKA', ['18353']), species('Na', 'c'), species('Na', 'e'), species('K', 'c'),
                species('K', 'e'), species('ATP'), species('H2O'), species('ADP'), species('Pi'), species('H')]
    bonds = [('Na_c', 'Re_NKA', 3), ('K_e', 'Re_NKA', 2), ('ATP_c', 'Re_NKA', 1), ('H2O_c', 'Re_NKA', 1),
             ('Re_NKA', 'Na_e', 3), ('Re_NKA', 'K_c', 2), ('Re_NKA', 'ADP_c', 1), ('Re_NKA', 'Pi_c', 1), ('Re_NKA', 'H_c', 1)]
    result = bondgraph.reactions_check(write(tmp_path, elements, bonds))
    assert result.status == checks.PASSED, result.details
    # K+ moved the wrong way: no reading of in and out fits
    os.makedirs(tmp_path / 'b')
    bonds = [b if 'K' not in b[0] + b[1] else (('K_c', 'Re_NKA', 2) if b[1] == 'Re_NKA' else ('Re_NKA', 'K_e', 2)) for b in bonds]
    assert bondgraph.reactions_check(write(tmp_path / 'b', elements, bonds)).status == checks.FAILED


def test_a_lumped_reaction_is_compared_by_its_net_equation(tmp_path):
    # hexokinase then glucose-6-phosphatase: a futile cycle, ATP + H2O = ADP + Pi + H+
    elements = [reaction('Re_cycle', parts=['17825', '16689'])] + [species(n) for n in ('ATP', 'H2O', 'ADP', 'Pi', 'H')]
    bonds = [('ATP_c', 'Re_cycle', 1), ('H2O_c', 'Re_cycle', 1), ('Re_cycle', 'ADP_c', 1), ('Re_cycle', 'Pi_c', 1),
             ('Re_cycle', 'H_c', 1)]
    result = bondgraph.reactions_check(write(tmp_path, elements, bonds))
    assert result.status == checks.PASSED, result.details
    assert '17825+16689' in only(result)


def test_a_lumped_reaction_may_take_a_part_more_than_once(tmp_path):
    # two hexokinase steps (2 Glc + 2 ATP = 2 G6P + 2 ADP + 2 H+) and one glucose-6-phosphatase back:
    # Glc + 2 ATP + H2O = G6P + 2 ADP + Pi + 2 H+
    elements = [reaction('Re_lump', parts=['17825', '16689'])] + [species(n) for n in ('Glc', 'ATP', 'H2O', 'G6P', 'ADP', 'Pi', 'H')]
    bonds = [('Glc_c', 'Re_lump', 1), ('ATP_c', 'Re_lump', 2), ('H2O_c', 'Re_lump', 1), ('Re_lump', 'G6P_c', 1),
             ('Re_lump', 'ADP_c', 2), ('Re_lump', 'Pi_c', 1), ('Re_lump', 'H_c', 2)]
    result = bondgraph.reactions_check(write(tmp_path, elements, bonds))
    assert result.status == checks.PASSED, result.details
    # one ADP too many: no multiples of the parts give it
    os.makedirs(tmp_path / 'b')
    bonds[4] = ('Re_lump', 'ADP_c', 3)
    assert bondgraph.reactions_check(write(tmp_path / 'b', elements, bonds)).status == checks.FAILED


def test_an_electrical_element_with_a_species_name_is_not_a_participant(tmp_path):
    # a membrane potential, written with a species name, bonded both ways (it nets to nothing anyway)
    elements = [reaction('Re_HK', ['17825'])] + [species(n) for n in ('Glc', 'ATP', 'G6P', 'ADP', 'H')] \
        + [('psi_m', 'ElectricalCapacitiveNode', 'psi', 'm', False)]
    bonds = [('Glc_c', 'Re_HK', 1), ('ATP_c', 'Re_HK', 1), ('Re_HK', 'G6P_c', 1), ('Re_HK', 'ADP_c', 1),
             ('Re_HK', 'H_c', 1), ('psi_m', 'Re_HK', 1), ('Re_HK', 'psi_m', 1)]
    result = bondgraph.reactions_check(write(tmp_path, elements, bonds))
    assert result.status == checks.PASSED, result.details


def test_a_polymer_counts_the_units_it_gains(tmp_path):
    # glycogen synthase: [glucosyl](n) + ADP-glucose = [glucosyl](n+1) + ADP + H+, a store of glucosyl units
    elements = [reaction('Re_GS', ['18189']), species('ADPglc'), species('Gly', cls='ChemicalStorageNode'),
                species('ADP'), species('H')]
    bonds = [('ADPglc_c', 'Re_GS', 1), ('Re_GS', 'Gly_c', 1), ('Re_GS', 'ADP_c', 1), ('Re_GS', 'H_c', 1)]
    result = bondgraph.reactions_check(write(tmp_path, elements, bonds))
    assert result.status == checks.PASSED, result.details


def test_participants_are_found_through_junctions(tmp_path):
    # ATP reaches the reaction through a 0-junction shared with nothing else
    elements = [reaction('Re_HK', ['17825'])] + [species(n) for n in ('Glc', 'ATP', 'G6P', 'ADP', 'H')]
    bonds = [('Glc_c', 'Re_HK', 1), ('ATP_c', 'J0', 1), ('J0', 'Re_HK', 1), ('Re_HK', 'G6P_c', 1),
             ('Re_HK', 'ADP_c', 1), ('Re_HK', 'H_c', 1)]
    result = bondgraph.reactions_check(write(tmp_path, elements, bonds, junctions=('J0',)))
    assert result.status == checks.PASSED, result.details


def test_the_cache_holds_the_tests_reactions():
    for rid in ('17825', '16689', '18353', '18189'):
        entry = rhea.load(rid)
        assert entry and entry['id'] == rid and entry['participants']
    assert rhea.load('18353')['is_transport'] and not rhea.load('17825')['is_transport']


def test_rheas_equation_and_chebi_list_give_its_participants():
    # the REST API's form: the distinct compounds' ChEBI terms in the order the equation names them
    parts = rhea.parse_equation('2 K(+)(out) + 3 Na(+)(in) + ATP + H2O = 2 K(+)(in) + 3 Na(+)(out) + ADP + phosphate + H(+)',
                                ['CHEBI:29103', 'CHEBI:29101', 'CHEBI:30616', 'CHEBI:15377', 'CHEBI:456216', 'CHEBI:43474',
                                 'CHEBI:15378'])
    entry = {'is_transport': True, 'participants': parts}
    assert bondgraph.rhea_vector(entry) == bondgraph.rhea_vector(rhea.load('18353'))
    parts = rhea.parse_equation('[(1->4)-alpha-D-glucosyl](n) + ADP-alpha-D-glucose = [(1->4)-alpha-D-glucosyl](n+1) + ADP + H(+)',
                                ['CHEBI:15444', 'CHEBI:57498', 'CHEBI:456216', 'CHEBI:15378'])
    assert bondgraph.rhea_vector({'participants': parts}) == bondgraph.rhea_vector(rhea.load('18189'))


def test_a_species_written_more_specifically_than_rheas_matches_by_chebi_is_a(tmp_path, monkeypatch):
    # the model's G6P is a (made-up) more specific term, which ChEBI says is_a Rhea's D-glucose 6-phosphate
    monkeypatch.setitem(CHEBI, 'G6Pa', 99999)
    monkeypatch.setattr(rhea, '_ancestors', {'CHEBI:99999': ['CHEBI:61548']})
    result = bondgraph.reactions_check(hexokinase(tmp_path, products=(('G6Pa', 1), ('ADP', 1), ('H', 1))))
    assert result.status == checks.PASSED, result.details
    assert 'by ChEBI is_a' in only(result)
    monkeypatch.setattr(rhea, '_ancestors', {})
    assert bondgraph.reactions_check(hexokinase(tmp_path, products=(('G6Pa', 1), ('ADP', 1), ('H', 1)))).status == checks.FAILED


def test_a_transport_reactions_unlocated_participants_match_wherever_the_model_has_them(tmp_path):
    # Na+/K+-ATPase with its ATP, water and products in a third compartment x: Rhea gives them no location
    elements = [reaction('Re_NKA', ['18353']), species('Na', 'c'), species('Na', 'e'), species('K', 'c'),
                species('K', 'e')] + [species(n, 'x') for n in ('ATP', 'H2O', 'ADP', 'Pi', 'H')]
    bonds = [('Na_c', 'Re_NKA', 3), ('K_e', 'Re_NKA', 2), ('ATP_x', 'Re_NKA', 1), ('H2O_x', 'Re_NKA', 1),
             ('Re_NKA', 'Na_e', 3), ('Re_NKA', 'K_c', 2), ('Re_NKA', 'ADP_x', 1), ('Re_NKA', 'Pi_x', 1), ('Re_NKA', 'H_x', 1)]
    result = bondgraph.reactions_check(write(tmp_path, elements, bonds))
    assert result.status == checks.PASSED, result.details


LEAK = {'id': '34979', 'equation': 'H(+)(in) = H(+)(out)', 'is_transport': True,
        'participants': [{'chebi': 'CHEBI:15378', 'side': 'L', 'coefficient': '1', 'location': 'in', 'name': 'H(+)'},
                         {'chebi': 'CHEBI:15378', 'side': 'R', 'coefficient': '1', 'location': 'out', 'name': 'H(+)'}]}


def test_a_reaction_of_protons_alone_compares_them_strictly(tmp_path):
    # a proton leak, H+(in) = H+(out): with nothing else in it, H+ isn't let off with a warning
    elements = [reaction('Re_LEAK', ['34979']), species('H', 'c'), species('H', 'm')]
    result = bondgraph.reactions_check(write(tmp_path, elements, [('H_c', 'Re_LEAK', 1), ('Re_LEAK', 'H_m', 1)]),
                                       load=lambda rid: LEAK)
    assert result.status == checks.PASSED and 'matches' in only(result)
    os.makedirs(tmp_path / 'b')
    result = bondgraph.reactions_check(write(tmp_path / 'b', elements, [('H_c', 'Re_LEAK', 2), ('Re_LEAK', 'H_m', 1)]),
                                       load=lambda rid: LEAK)
    assert result.status == checks.FAILED
