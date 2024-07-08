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

    substrate_graph = nx.Graph()
    for atom in substrate_mol.GetAtoms():
        substrate_graph.add_node(atom.GetAtomMapNum(), label=atom.GetAtomicNum())
    for bond in substrate_mol.GetBonds():
        if bond.GetBondType() == BondType.SINGLE:
            bond_type = 1
        elif bond.GetBondType() == BondType.DOUBLE:
            bond_type = 2
        else:
            bond_type = 3
        substrate_graph.add_edge(bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum(), weight=bond_type)
    
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
        if not check_if_bond_in_core(bond, reactant_core):
             core_fragment_edges.add((bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()))
    core_to_substrate = {subgraphs[0][key]: key for key in subgraphs[0]}
    substrate_fragment_edges = {(core_to_substrate[edge[0]], core_to_substrate[edge[1]]) for edge in core_fragment_edges}
    fragment = Fragment({atom for atom in substrate_mol.GetAtoms()}, {bond for bond in substrate_mol.GetBonds()})
    fragments = fragment.fragment_with_multiple_edges(substrate_fragment_edges)

    # use the reaction core to make changes to different fragments
    for fragment in fragments:
        fragment.transform(reactant_core, core_to_substrate, core_fragment_edges)
    
    return {fragment.get_mol() for fragment in fragments}
    

def check_if_bond_in_core(bond: Bond, other_core: ReactionCore) -> bool:
    """
    Checks whether or not a bond exists in a set of molecules
    """
    map_num1, map_num2 = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
    for possible_bond in other_core.bonds:
        possible1, possible2 = possible_bond.GetBeginAtom().GetAtomMapNum(), possible_bond.GetEndAtom().GetAtomMapNum()
        if (map_num1 == possible1 and map_num2 == possible2) or (map_num1 == possible2 and map_num2 == possible1):
            return True
    
    return False


if __name__ == "__main__":
    core = ce.get_reaction_core(["O=C1C(C(OCC)=O)CCC1", "C=CC(C[R1])=C"], ["O=C1C(C(OCC)=O)(CCC(C[R1])=C)CCC1"])
    # core = ce.get_reaction_core(["[R]C1=CC=CC=C1", "C"], ["C1=CC=CC=C1", "[R]C"])
    product_core = ce.extend_reaction_core(core[1], core[2], core[0])
    mols = get_reactants_for_substrate("O=C1C(C(OCC)=O)(CCC(C[C])=C)CCC1", core[0], product_core)
    # mols = get_reactants_for_substrate("C1=CC=CC=C1", core[0], product_core)
    index = 1
    for mol in mols:
        highlight_reaction_core(mol, set(), set(), "new_reactant" + str(index) + ".png")
        index += 1
    