"""
File containing code for template matching a core to substrate and fragment and rule retrieval
"""

from rdkit import Chem
from rdkit.Chem.rdchem import Mol, Atom, Bond, BondType
from rdkit.Chem.rdChemReactions import ChemicalReaction
from typing import Set, Tuple, List

import networkx as nx

from classes.partial_molecule import ReactionCore, Fragment
from classes.rules import Rule, ModifyBond, AddBond
from utils import find_atom


##################################################################################################
# Template matching/Fragment retrieval
##################################################################################################

def get_reactants_for_substrate(substrate: Mol, reactant_core: ReactionCore, product_core: ReactionCore):
    """
    Matches the core template to a fixed substrate
    """
    
    # use nx to do subgraph matching
    substrate_graph = nx.Graph()
    for bond in substrate.GetBonds():
        if bond.GetBondType() == BondType.SINGLE:
            bond_type = 1
        elif bond.GetBondType() == BondType.DOUBLE:
            bond_type = 2
        else:
            bond_type = 3
        substrate_graph.add_edge(bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtomMapNum().GetAtomMapNum(), bond_type=bond_type)
    
    core_graph = nx.Graph()
    for bond in product_core.bonds:
        if bond.GetBondType() == BondType.SINGLE:
            bond_type = 1
        elif bond.GetBondType() == BondType.DOUBLE:
            bond_type = 2
        else:
            bond_type = 3
        core_graph.add_edge(bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtomMapNum().GetAtomMapNum(), bond_type=bond_type)
    
    matcher = nx.algorithms.isomorphism.GraphMatcher(core_graph, substrate_graph)
    subgraphs = [subgraph for subgraph in matcher.subgraph_isomorphisms_iter()]
    
    # assume only one match, deal with later
    
    # now we have the reaction core in the substrate molecule
    # we can figure out where to break bonds to generate fragments in core area
    core_fragment_edges = set()
    for bond in reactant_core.bonds:
        if not check_if_bond_in_core(bond, product_core):
            core_fragment_edges.add((bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()))
            
    substrate_fragment_edges = {(subgraphs[0][edge[0]], subgraphs[0][edge[1]]) for edge in core_fragment_edges}
    fragment = Fragment({atom for atom in substrate.GetAtoms()}, {bond for bond in substrate.GetBonds()})
    fragments = fragment.fragment_with_multiple_edges(substrate_fragment_edges)
    
    # use the reaction core to make changes to different fragments
    for fragment in fragments:
        fragment.transform(reactant_core, subgraphs[0])
    
    return fragments
    

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
               