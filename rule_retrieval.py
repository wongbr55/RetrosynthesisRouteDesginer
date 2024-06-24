"""
File containing code for getting reactants given a reaction core
"""

from rdkit import Chem
from rdkit.Chem.rdchem import Mol, Atom, Bond, BondType
from rdkit.Chem.rdChemReactions import ChemicalReaction
from typing import Set, Tuple

from partial_molecule import ReactionCore, Fragment

def get_fragments(substrate: Mol, core: ReactionCore, reactant_template: Set[Mol], product_template: Set[Mol]):
    """
    Given a reaction core of a given molecule "substrate", finds the corrosponding reactants that made that
    Substrate
    """
    fragment_map_nums = set()
    # find the atoms where bonds have broken
    for mol in reactant_template:
        for atom in mol.GetAtoms():
            for bond in atom.GetBonds():
                if not check_if_bond_in_mols(bond, product_template):
                    fragment_map_nums.add((bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()))     

    # with complete set of "edges" for the fragments, we break up the molecule
    return get_fragments(substrate, fragment_map_nums)
    

def check_if_bond_in_mols(bond: Bond, molecules: Set[Mol]) -> bool:
    """
    Checks whether or not a bond exists in a set of molecules
    """
    map_num1, map_num2 = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
    for mol in molecules:
        for possible_bond in mol.GetBonds():
            possible1, possible2 = possible_bond.GetBeginAtom().GetAtomMapNum(), possible_bond.GetEndAtom().GetAtomMapNum()
            print("Comparing with: " + str(possible1, possible2))
            if (map_num1 == possible1 and map_num2 == possible2) or (map_num1 == possible2 and map_num2 == possible1):
                return True
    
    return False
    

def get_fragments(substrate: Mol, fragment_edges: Set[Tuple[int]]):
    """
    Gets all of the fragments from the substrate
    """
    
    fragments = set()
    curr_fragment = Fragment(substrate.GetAtoms(), substrate.GetBonds())
    fragments.add(curr_fragment)
    
    fragments = curr_fragment.fragment_with_multiple_edges(fragment_edges)
    
    # for edge in fragment_edges:
    #     # find fragment that has the edge (since fragments has all of the starting Atom and Bond objects)
    #     # this will always work
    #     for fragment in fragments:
    #         # once we have found, generate new fragments
    #         if fragment.check_bond_map_num(edge[0], edge[1]):
    #             new_fragment = fragment.fragment(edge)
    #             if new_fragment is not None:
    #                 fragments.add(new_fragment)
    #             break
    
    return fragments     

