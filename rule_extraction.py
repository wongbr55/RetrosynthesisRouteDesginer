"""
Rule extraction. Used to create database/graph of possible product types to reactants
Implemented Herustics from Route Designer
"""
from rdkit import Chem
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import rdChemReactions
from rdkit.Chem.rdChemReactions import ChemicalReaction


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


def get_reaction_core_helper(reactant_molecule_list: list[Mol], reaction: ChemicalReaction):
    """
    Identifies the reaction core of the reactants
    Assume RXN property for reactant and product are set
    :param reactant_molecule_list:
    :return:
    """

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


def get_reaction_core_for_single_reaction(reactant_smiles: list[str], product_smiles: list[str]):
    """
    Runs the full reaction core extraction
    :param reactant_smiles:
    :param product_smiles:
    :return:
    """
    # setup Mol objects
    reactant_mol = []
    product_mol = []
    for smile in reactant_smiles:
        new_mol = Chem.MolFromSmiles(smile)
        new_mol.SetProp("RXN", "reactant")
        reactant_mol.append(new_mol)
    for smile in product_smiles:
        new_mol = Chem.MolFromSmiles(smile)
        new_mol.SetProp("RXN", "product")
        product_mol.append(new_mol)

    # setup ChemicalReaction
    reaction = Chem.rdChemReactions.ReactionFromMolecule(reactant_mol[0])
    for i in range(1, len(reactant_mol)):
        reaction.AddReactantTemplate(reactant_mol[i])
    for product in product_mol:
        reaction.AddProductTemplate(product)

    # get core
    # we create a new Mol object using the changed Atoms and Bonds
    changed_props = get_reaction_core_helper(reactant_mol, reaction)
    atom_conversion = {}
    index = 0
    for atom in changed_props[0]:

        atom_conversion[atom] = index
        index += 1
