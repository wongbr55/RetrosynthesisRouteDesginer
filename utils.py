"""
File that contains common utility functions
"""

from rdkit import Chem
from rdkit.Chem.rdchem import Mol, Atom, Bond, BondType
from rdkit.Chem import Draw
from rdkit.Chem.rdChemReactions import ChemicalReaction
from typing import List, Set, Tuple, Iterable

import networkx as nx
from networkx import Graph

REMOVE_BOND_FLAG = "REMOVE"
ADD_BOND_FLAG = "ADD"
MODIFY_BOND_FLAG = "MODIFY"

R_GROUP_PLACEHOLDER = 30


VALENCY_TABLE = {7: 3}

def generate_graph(mol: Mol) -> Graph:
    """
    Generates Mol as a Graph object
    Does not include bond 
    """
    graph = Graph()
    for atom in mol.GetAtoms():
        graph.add_node(atom.GetAtomMapNum(), label=atom.GetAtomicNum())
    for bond in mol.GetBonds():
        graph.add_edge(bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum(), bond_type=bond.GetBondTypeAsDouble())
    
    return graph

def generate_graph_with_indicies_as_labels(mol: Mol) -> Graph:
    """
    Generates graph and uses indicies of atoms as labels instead of map numbers
    """
    graph = Graph()
    for atom in mol.GetAtoms():
        graph.add_node(atom.GetIdx(), label=atom.GetAtomicNum())
    for bond in mol.GetBonds():
        graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond_type=bond.GetBondTypeAsDouble())
    
    return graph

def create_atom(atom: Atom) -> Atom:
    """
    Creates a new Atom object from atom with the same properties
    This is done to clear all bonds and neighbours from the Atom object
    """
    new_atom = Atom(atom.GetAtomicNum())
    new_atom.SetAtomMapNum(atom.GetAtomMapNum())
    new_atom.SetFormalCharge(atom.GetFormalCharge())
    new_atom.SetHybridization(atom.GetHybridization())
    
    return new_atom

def find_atom(atom_map_num: int, mols: Iterable[Mol]) -> Atom:
    """
    Finds Atom object that has the same atom map number as atom_map_num
    """
    for mol in mols:
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() == atom_map_num:
                return atom
    return None

def find_atom_index(atom_map_num: int, mols: Iterable[Mol]) -> Atom:
    """
    Finds Atom object that has the same atom map number as atom_map_num
    """
    for i in range(0, len(mols)):
        for atom in mols[i].GetAtoms():
            if atom.GetAtomMapNum() == atom_map_num:
                return i
    return -1

def find_atom_index_in_mol(atom_map_num: int, mols: Iterable[Mol]) -> Atom:
    """
    Finds Atom object that has the same atom map number as atom_map_num
    """
    for i in range(0, len(mols)):
        for atom in mols[i].GetAtoms():
            if atom.GetAtomMapNum() == atom_map_num:
                return atom.GetIdx()
    return -1

def find_bond_set(index1: int, index2: int, bonds: Iterable[Bond]) -> Bond:
    for bond in bonds:
        possible1, possible2 = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
        if (index1 == possible1 and index2 == possible2) or (index1 == possible2 and index2 == possible1):
            return bond
    return None

def find_bond(index1: int, index2: int, mols: Iterable[Mol]) -> Bond:
    """
    Finds Bond object that has same bond map numbers as index1 and index2
    """
    
    for mol in mols:
        for bond in mol.GetBonds():
            possible1, possible2 = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
            if (index1 == possible1 and index2 == possible2) or (index1 == possible2 and index2 == possible1):
                return bond
    return None


def highlight_reaction_core(mol, changing_atoms, changing_bonds, filename):
    """
    Draws highlted reaction core on reactant and product
    Only used for debugging
    :param mol:
    :param changing_atoms:
    :param changing_bonds:
    :return:
    """
    atom_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomMapNum() in changing_atoms]
    bond_indices = [bond.GetIdx() for bond in mol.GetBonds() if
                    (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) in changing_bonds or (
                        bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()) in changing_bonds]
    return Draw.MolToFile(mol, filename, highlightAtoms=atom_indices, highlightBonds=bond_indices)


def print_atom_and_bond(reaction_core: Tuple[Set[Atom], Set[Bond]]) -> None:
    """
    Prints the map numbers of Atom and Bond objects in a reaction core
    Used for debugging purposes
    :param reaction_core:
    :return:
    """
    print("## ATOMS ##")
    for atom in reaction_core[0]:
        print("Atom :" + str(atom.GetAtomMapNum()))
    print("## BONDS ##")
    for bond in reaction_core[1]:
        print("From " + str(bond.GetBeginAtom().GetAtomMapNum()) + " to " + str(bond.GetEndAtom().GetAtomMapNum()))
        

# below are used to networkx graph matching
def node_match(n1, n2):
    return n1['label'] == n2['label']


def edge_match(e1, e2):
    return e1['bond_type'] == e2['bond_type']
