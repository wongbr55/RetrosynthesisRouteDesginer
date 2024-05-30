"""
Rule extraction. Used to create database/graph of possible product types to reactants
Implemented Herustics from Route Designer
"""
import rdkit
from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import rdChemReactions


def get_reaction_core_for_single_reactant(reaction_molecule: Mol, products: list[Mol]):
    """
    Gets the reaction core for a single molecule
    Returns a tuple containing the changed atoms and bonds that represent the reaction core
    :param reaction_molecule:
    :param products:
    :return: tuple[set]
    """
    reactant_atom_set = set()
    product_atom_set = set()

    for atom in reaction_molecule.GetAtoms():
        if atom.GetAtomMapNum() > 0:
            reactant_atom_set.add(atom)
    for product in products:
        for atom in product.GetAtoms():
            if atom.GetAtomMapNum() > 0:
                product_atom_set.add(atom)

    changed_atoms = reactant_atom_set.symmetric_difference(product_atom_set)

    reactant_bond_set = set()
    product_bond_set = set()

    for bond in reaction_molecule.GetBonds():
        if bond.GetBeginAtom() in changed_atoms or bond.GetEndAtom() in changed_atoms:
            reactant_bond_set.add((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))

    for product in products:
        for bond in product.GetBonds():
            if bond.GetBeginAtom() in changed_atoms or bond.GetEndAtom() in changed_atoms:
                product_bond_set.add((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))

    changed_bonds = reactant_bond_set.symmetric_difference(product_bond_set)

    return changed_atoms, changed_bonds


def get_reaction_core(reactant_molecule_list: list[Mol], rxn_smarts: str):
    """
    Identifies the reaction core of the reactants
    :param reactant_molecule_list:
    :return:
    """

    reaction = rdChemReactions.ReactionFromSmarts(rxn_smarts)
    product_mols = reaction.RunReactants(reactant_molecule_list)

    # get all the reaction cores
    reaction_cores = []
    for reactant in reactant_molecule_list:
        new_core = get_reaction_core_for_single_reactant(reactant, product_mols)
        reaction_cores.append(new_core)

    # find intersection of all the cores
    changed_atoms = set()
    changed_bonds = set()
    for core in reaction_cores:
        changed_atoms = changed_atoms.symmetric_difference(core[0])
        changed_bonds = changed_bonds.symmetric_difference(core[1])

    return changed_atoms, changed_bonds
