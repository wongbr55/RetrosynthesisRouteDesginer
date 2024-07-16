"""
File containing code for template matching a core to substrate and fragment and rule retrieval
"""

from rdkit import Chem
from rdkit.Chem.rdchem import Mol, Atom, Bond, BondType
from rdkit.Chem.rdChemReactions import ChemicalReaction
from typing import Set, Tuple, List

import networkx as nx
import matplotlib.pyplot as plt
from utils import node_match, edge_match, highlight_reaction_core
import core_extraction as ce

from classes.partial_molecule import ReactionCore, Fragment

##################################################################################################
# Template matching/Fragment retrieval
##################################################################################################

def get_reactants_for_substrate(substrate: str, reactant_core: ReactionCore, product_core: ReactionCore) -> Set[Mol]:
    """
    Matches the core template to a fixed substrate and returns proper reactants
    """
    substrate_mol = Chem.MolFromSmiles(substrate)
    counter = 1
    for atom in substrate_mol.GetAtoms():
        atom.SetAtomMapNum(counter)
        counter += 1
    
    highlight_reaction_core(substrate_mol, set(), set(), "substrate.png")

    substrate_graph = nx.Graph()
    for atom in substrate_mol.GetAtoms():
        substrate_graph.add_node(atom.GetAtomMapNum(), label=atom.GetAtomicNum())
    for bond in substrate_mol.GetBonds():
        substrate_graph.add_edge(bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum(), weight=bond.GetBondTypeAsDouble())
    
    core_graph = nx.Graph()
    for atom in product_core.atoms:
        core_graph.add_node(atom.GetAtomMapNum(), label=atom.GetAtomicNum())
    for bond in product_core.bonds:
        core_graph.add_edge(bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum(), weight=bond.GetBondTypeAsDouble())
        
    matcher = nx.algorithms.isomorphism.GraphMatcher(substrate_graph, core_graph, node_match=node_match, edge_match=edge_match)
    subgraphs = [subgraph for subgraph in matcher.subgraph_isomorphisms_iter()]
    
    # assume only one match, deal with later
    # now we have the reaction core in the substrate molecule
    # we can figure out where to break bonds to generate fragments in core area
    core_fragment_edges = set()
    for bond in product_core.bonds:
        atom1, atom2 = bond.GetBeginAtom(), bond.GetEndAtom()
        map1, map2 = atom1.GetAtomMapNum(), atom2.GetAtomMapNum()
        if not reactant_core.check_bond(bond) and reactant_core.check_atom(atom1) and reactant_core.check_atom(atom2): # and reactant_core.check_atom(atom1) and reactant_core.check_atom(atom2):
             core_fragment_edges.add((bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()))
    core_to_substrate = {subgraphs[0][key]: key for key in subgraphs[0]}
    substrate_fragment_edges = {(core_to_substrate[edge[0]], core_to_substrate[edge[1]]) for edge in core_fragment_edges}
    # print(substrate_fragment_edges)
    fragment = Fragment({atom for atom in substrate_mol.GetAtoms()}, {bond for bond in substrate_mol.GetBonds()})
    fragments = fragment.fragment_with_multiple_edges(substrate_fragment_edges)
    # use the reaction core to make changes to different fragments
        
    next_index_to_add = max({atom.GetAtomMapNum() for atom in substrate_mol.GetAtoms()})
    reactants = set()
    for fragment in fragments:
        next_index_to_add = fragment.transform(reactant_core, core_to_substrate, core_fragment_edges, next_index_to_add)
        reactants.add(fragment.get_mol())
    
    return reactants


if __name__ == "__main__":
    # core = ce.get_reaction_core(["O=C1C(C(OCC)=O)CCC1", "C=CC(C[R1])=C"], ["O=C1C(C(OCC)=O)(CCC(C[R1])=C)CCC1"])
    
    # wiley2_table_3
    # core = ce.get_reaction_core(["O=C(NC)C1=CC=CC=C1/C=C([R1])/[R2]"], ["O=C(N1C)C2=CC=CC=C2C1C([R1])[R2]"])
    # product_core = ce.extend_reaction_core(core[1], core[2], core[0])
    # # print(product_core)
    # highlight_reaction_core(core[0].get_mol(), set(), set(), "reactant_core.png")
    # highlight_reaction_core(product_core.get_mol(), set(), set(), "product_core.png")
    # mols = get_reactants_for_substrate("O=C(N1C)C2=CC=CC=C2C1CC3=CC=C(C(C)=O)C=C3", core[0], product_core)
    
    # wiley26_scheme_2
    # core = ce.get_reaction_core(["[R]O"], ["[R]OC1CCCO1", "[H][H]"])
    # product_core = ce.extend_reaction_core(core[1], core[2], core[0])
    # mols = get_reactants_for_substrate("COC1=CC=C(CCOC2CCCO2)C=C1", core[0], product_core)
    
    # wiley38_table_2
    # core = ce.get_reaction_core(["IC1CCOCC1", "OC(C1=CC([R])=CC=C1)=O"], ["O=C(C1CCOCC1)C2=CC=CC([R])=C2"])
    # product_core = ce.extend_reaction_core(core[1], core[2], core[0])
    
    # highlight_reaction_core(core[0].get_mol(), set(), set(), "reactant_core.png")
    # highlight_reaction_core(product_core.get_mol(), set(), set(), "product_core.png")
    # mols = get_reactants_for_substrate("O=C(C1CCOCC1)C2=CC=C(C(C)(C)C)C=C2", core[0], product_core)
    
    # rsc4_scheme_2
    # core = ce.get_reaction_core(["[R2]C(C1=CC([R1])=CC=C1)=C", "[R4]C(C1=CC([R3])=CC=C1)=C", "[R]O"], ["[R2]C(O[R])(C1=CC([R1])=CC=C1)CCC(O[R])(C2=CC=CC([R3])=C2)[R4]"])
    # product_core = ce.extend_reaction_core(core[1], core[2], core[0])
    # # print(product_core)
    # highlight_reaction_core(core[0].get_mol(), set(), set(), "reactant_core.png")
    # highlight_reaction_core(product_core.get_mol(), set(), set(), "product_core.png")
    # mols = get_reactants_for_substrate("COC(C1=CC=CC=C1)CCC(OC)C2=CC=CC=C2", core[0], product_core)
    
    # wiley27_scheme_2
    core = ce.get_reaction_core(["[R1]N(C1=C(/C2=N\[H])C=CC([R2])=C1)C(N2C3=CC=C([R3])C=C3)=O"], ["O=C(N([R1])C(C=C([R2])C=C1)=C1C2=N3)N2C4=C3C=C([R3])C=C4", "[H][H]"])
    product_core = ce.extend_reaction_core(core[1], core[2], core[0])
    # # print(product_core)
    highlight_reaction_core(core[0].get_mol(), set(), set(), "reactant_core.png")
    highlight_reaction_core(product_core.get_mol(), set(), set(), "product_core.png")
    mols = get_reactants_for_substrate("O=C(N(CC1=CC=CC=C1)C(C=CC=C2)=C2C3=N4)N3C5=C4C=C(CCO)C=C5", core[0], product_core)   

    # product_core = ce.extend_reaction_core(core[1], core[2], core[0])
    # highlight_reaction_core(core[0].get_mol(), set(), set(), "extended_reactant_core.png")
    # highlight_reaction_core(product_core.get_mol(), set(), set(), "extended_product_core.png")
    # mols = get_reactants_for_substrate("O=C(C1CCOCC1)C2=CC=C(C)C=C2", core[0], product_core)
    # mols = get_reactants_for_substrate("C1=CC=CC=C1", core[0], product_core)
    index = 1
    for mol in mols:
        highlight_reaction_core(mol, set(), set(), "new_reactant" + str(index) + ".png")
        print(Chem.MolToSmiles(mol))
        index += 1
    