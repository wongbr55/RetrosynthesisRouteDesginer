"""
Core extraction. Code provided extracts core from reactions
Implemented Herustics from Route Designer
"""
from rdkit import Chem
from rdkit.Chem.rdchem import Mol, Atom, Bond
from rdkit.Chem.rdchem import BondType
from rdkit.Chem import rdChemReactions, Draw
from rdkit.Chem.rdChemReactions import ChemicalReaction

##################################################################################################
# CORE EXTRACTION
##################################################################################################


def get_reaction_core_for_single_reaction(reactant_smiles: list[str], product_smiles: list[str]) -> \
        tuple[tuple[set[Atom], set[Bond]], list[Mol], list[Mol]]:
    """
    Runs the full reaction core extraction
    Returns a tuple:
    first index is a tuple with new reaction core atoms and bonds
    second index is list of reactant in Mol objects
    third is list of products in Mol objects
    :param reactant_smiles:
    :param product_smiles:
    :return:
    """
    # setup Mol objects
    reactant_mol = []
    product_mol = []
    reaction = ChemicalReaction()
    for smile in reactant_smiles:
        new_mol = Chem.MolFromSmiles(smile)
        # new_mol.SetProp("RXN", "reactant")
        reactant_mol.append(new_mol)
    for smile in product_smiles:
        new_mol = Chem.MolFromSmiles(smile)
        # new_mol.SetProp("RXN", "product")
        product_mol.append(new_mol)

    reactant_counter = 0
    for reactant in reactant_mol:
        for num in range(0, reactant.GetNumAtoms()):
            reactant.GetAtomWithIdx(num).SetAtomMapNum(reactant_counter)
            reactant_counter += 1
        reaction.AddReactantTemplate(reactant)

    product_counter = 0
    for product in product_mol:
        for num in range(0, product.GetNumAtoms()):
            product.GetAtomWithIdx(num).SetAtomMapNum(product_counter)
            product_counter += 1
        reaction.AddProductTemplate(product)

    # get core
    changed_props = get_reaction_core_helper(reactant_mol, reaction)
    return changed_props, reactant_mol, product_mol


def get_reaction_core_helper(reactant_molecule_list: list[Mol], reaction: ChemicalReaction) -> tuple[set[Atom], set[Bond]]:
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
        new_core = get_reaction_core_for_single_reactant(reactant, product_mols[0])
        reaction_cores.append(new_core)
    # find intersection of all the cores
    changed_atoms = set()
    changed_bonds = set()
    for core in reaction_cores:
        changed_atoms = changed_atoms.symmetric_difference(core[0])
        changed_bonds = changed_bonds.symmetric_difference(core[1])

    return changed_atoms, changed_bonds


def get_reaction_core_for_single_reactant(reaction_molecule: Mol, products: list[Mol]) -> tuple[set[Atom], set[Bond]]:
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


def highlight_reaction_core(mol, changing_atoms, changing_bonds):
    """
    Draws highlted reaction core on reactant and product
    Only used for double-checking work
    :param mol:
    :param changing_atoms:
    :param changing_bonds:
    :return:
    """
    atom_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomMapNum() in changing_atoms]
    bond_indices = [bond.GetIdx() for bond in mol.GetBonds() if
                    (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) in changing_bonds or (
                        bond.GetEndAtomIdx(), bond.GetBeginAtomIdx()) in changing_bonds]
    return Draw.MolToImage(mol, highlightAtoms=atom_indices, highlightBonds=bond_indices)


##################################################################################################
# CORE EXTENSION
##################################################################################################


def extend_reaction_core(reactant_mol: list[Mol], product_mol: list[Mol], reaction_core: tuple[set[Atom], set[Bond]]) -> None:
    """
    Takes the current core in terms of atoms and bonds and extends the reaction core
    Heurstics found in section 2.2 of Route Desiginer paper
    Returns the new reaction core
    :param reactant_mol:
    :param product_mol:
    :param reaction_core:
    :return:
    """

    # extend core
    for atom in reaction_core[0]:
        # we first employ some special heursitics for primary bonds
        extend_primary_bonds(atom, reaction_core)
        # then we apply the general heurstics
        add_atoms_to_core(atom, reaction_core)
    # get leaving groups and add them to the reaction core
    get_leaving_group(reactant_mol, product_mol, reaction_core)


def get_leaving_group(reactant_mol: list[Mol], product_mol: list[Mol], reaction_core: tuple[set[Atom], set[Bond]]) -> None:
    """
    Gets leaving group Atoms and Bonds and mutates reaction_core to add them
    :param reactant_mol:
    :param product_mol:
    :param reaction_core:
    :return:
    """

    # check for leaving groups and add them to reaction core
    reactant_atom_set = set()
    product_atom_set = set()
    reactant_bond_set = set()
    product_bond_set = set()

    for reactant in reactant_mol:
        for atom in reactant.GetAtoms():
            if atom.GetAtomMapNum() > 0:
                reactant_atom_set.add(atom)
        for bond in reactant.GetBonds():
            if bond.GetBeginAtom() in reactant_atom_set or bond.GetEndAtom() in reactant_atom_set:
                reactant_bond_set.add(bond)
    for product in product_mol:
        for atom in product.GetAtoms():
            if atom.GetAtomMapNum() > 0:
                product_atom_set.add(atom)
        for bond in product.GetBonds():
            if bond.GetBeginAtom() in product_atom_set or bond.GetEndAtom() in product_atom_set:
                product_bond_set.add(bond)

    leaving_group_atoms = reactant_atom_set.symmetric_difference(product_atom_set)
    leaving_group_bonds = reactant_bond_set.symmetric_difference(product_bond_set)
    for leaving_atom in leaving_group_atoms:
        reaction_core[0].add(leaving_atom)
    for leaving_bond in leaving_group_bonds:
        if (leaving_bond.GetBeginAtom() not in reaction_core[0] and leaving_bond.GetEndAtom() not in reaction_core[0])\
                and (leaving_bond.GetBeginAtom() in reactant_atom_set and
                     leaving_bond.GetEndAtom() in reactant_atom_set):
            reaction_core[1].add(leaving_bond)


def extend_primary_bonds(atom: Atom, reaction_core: tuple[set[Atom], set[Bond]]) -> None:
    """
    Adds new atoms to reaction core
    Assume that atom is a member of the reaction core
    We are checking primary bonds as they employ different heurstics than regular bonds
    Mutates reaction_core with new atoms and bonds
    Preconditions
        - atom in reaction_core[0]
    :param atom:
    :param reaction_core:
    :return:
    """
    have_searched_aromatic = False
    for bond in atom.GetBonds():
        if atom == bond.GetEndAtom():
            new_atom = bond.GetBeginAtom()
        else:
            new_atom = bond.GetEndAtom()

        # check for secoondary double or triple bonds
        for secondary_bond in new_atom.GetBonds():
            if secondary_bond.GetBondType() == BondType.DOUBLE or secondary_bond.GetBondType() == BondType.TRIPLE:
                if new_atom == secondary_bond.GetBeginAtom():
                    newest_atom = secondary_bond.GetEndAtom()
                else:
                    newest_atom = secondary_bond.GetBeginAtom()
                reaction_core[0].add(new_atom)
                reaction_core[0].add(newest_atom)
                reaction_core[1].add(bond)
                reaction_core[1].add(secondary_bond)

        # check for aromatic bonds
        if bond.GetIsAromatic() and not have_searched_aromatic:
            get_aromatic_ring(atom, reaction_core)
            have_searched_aromatic = True


def get_aromatic_ring(atom: Atom, reaction_core: tuple[set[Atom], set[Bond]]) -> None:
    """
    Gets all the atoms and bonds in an aromatic ring
    Mutates the
    :param atom:
    :param reaction_core:
    :return:
    """

    ownning_mol = atom.GetOwningMol()
    aromatic_rings = Chem.GetSymmSSSR(ownning_mol)

    # find the aromatic ring containing atom
    aromatic_ring_containing_atom = None
    for ring in aromatic_rings:
        if atom.GetIdx() in ring:
            aromatic_ring_containing_atom = ring
            break

    # retrieve bonds connected to atom in the ring
    for bond_idx in aromatic_ring_containing_atom:
        bond = ownning_mol.GetBondWithIdx(bond_idx)
        reaction_core[1].add(bond)
        # we can add both atom objects as reaction_core[0] is a set, will not hold duplicates
        reaction_core[0].add(bond.GetBeginAtomIdx)
        reaction_core[0].add(bond.GetEndAtomIdx)


def add_atoms_to_core(atom: Atom, reaction_core: tuple[set[Atom], set[Bond]]) -> None:
    """
    Adds new atoms to the reaction core
    Essentially employs a BFS search on molecule starting at atom
    Mutates reaction core for extension
    :param atom:
    :param reaction_core:
    :return:
    """

    queue = [atom]
    seen_so_far = set()
    while len(queue) != 0:
        curr_atom = queue.pop(0)
        for bond in curr_atom.GetBonds():
            if curr_atom == bond.GetEndAtom():
                new_atom = bond.GetBeginAtom()
            else:
                new_atom = bond.GetEndAtom()
            if not check_external_nonaromatic_bond(bond) and new_atom not in reaction_core[0] and \
                    bond not in reaction_core[1] and new_atom not in seen_so_far:
                reaction_core[0].add(new_atom)
                reaction_core[1].add(bond)
                queue.append(new_atom)
        seen_so_far.add(curr_atom)


def check_external_nonaromatic_bond(bond: Bond) -> bool:
    """
    Checks if current bond is an external nonaromatic carbon-carbon bond
    Returns True of False upon completion
    :param bond:
    :return: boolean
    """
    atom1 = bond.GetBeginAtom()
    atom2 = bond.GetEndAtom()

    # check to see if atoms are carbon
    if atom1.GetAtomicNum() != 6 or atom2.GetAtomicNum() != 6:
        return False
    # check if the bond is aromatic
    if bond.GetIsAromatic():
        return False
    # check if the bond is external (not part of a ring)
    if atom1.IsInRing() or atom2.IsInRing():
        return False
    # otherwise we are good
    return True
