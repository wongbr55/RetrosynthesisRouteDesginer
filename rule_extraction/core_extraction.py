"""
Core extraction. Code provided extracts core from reactions
Implemented Herustics from Route Designer
"""
from rdkit import Chem
from rdkit.Chem.rdchem import Mol, Atom, Bond, BondType
from rdkit.Chem import Draw
from rdkit.Chem.rdChemReactions import ChemicalReaction
from typing import List, Set, Tuple


##################################################################################################
# CORE EXTRACTION
##################################################################################################


def get_reaction_core(reaction_smarts: str) -> Tuple[Tuple[Set[Atom], Set[tuple]], List[Mol], List[Mol]]:
    """
    Runs the full reaction core
    Returns a tuple:
    first index is a tuple with new reaction core atoms and bonds
    second index is list of reactant in Mol objects
    third is list of products in Mol objects
    :return:
    """

    reaction = Chem.rdChemReactions.ReactionFromSmarts(reaction_smarts)

    # setup Mol objects
    reactant_mol = []
    product_mol = []
    for mol in reaction.GetReactants():
        reactant_mol.append(mol)
    for mol in reaction.GetProducts():
        product_mol.append(mol)

    # get core
    changed_props = get_reaction_core_helper(reactant_mol, product_mol)
    return changed_props, reactant_mol, product_mol


def get_reaction_core_helper(reactant_molecule_list: List[Mol], product_mols: List[Mol]) -> Tuple[Set[Atom], Set[Bond]]:
    """
    Identifies the reaction core of the reactants
    Assume RXN property for reactant and product are set
    :param reactant_molecule_list:
    :return:
    """

    # get all the reaction cores
    reaction_cores = []
    for reactant in reactant_molecule_list:
        new_core = get_reaction_core_for_single_reactant(reactant, product_mols)
        reaction_cores.append(new_core)
    # find intersection of all the cores
    changed_atoms = set()
    changed_bonds = set()

    for core in reaction_cores:
        atom_no_duplicates = set()
        for atom in core[0]:
            if not any(atom.GetAtomMapNum() == other_atom.GetAtomMapNum() for other_atom in changed_atoms):
                atom_no_duplicates.add(atom)
        changed_atoms = changed_atoms.union(atom_no_duplicates)

        bonds_no_duplicates = set()
        for bond in core[1]:
            if not any(bond.Match(other_bond) for other_bond in changed_bonds):
                bonds_no_duplicates.add(bond)
        changed_bonds = changed_bonds.union(bonds_no_duplicates)

    return changed_atoms, changed_bonds


def get_reaction_core_for_single_reactant(reaction_molecule: Mol, products: List[Mol]) -> Tuple[Set[Atom], Set[Bond]]:
    """
    Gets the reaction core for a single molecule
    Returns a tuple containing the changed atoms and bonds that represent the reaction core
    :param reaction_molecule:
    :param products:
    :return: tuple[set]
    """

    changed_atoms = set()
    changed_bonds = set()
    changed_atom_map_num = set()

    for atom in reaction_molecule.GetAtoms():
        index = atom.GetAtomMapNum()

        for product in products:
            for prod_atom in product.GetAtoms():
                # once we found the matching atom, we compare attributes
                if prod_atom.GetAtomMapNum() == index and not compare_props(atom, prod_atom):
                    changed_atoms.add(atom)
                    changed_atom_map_num.add(atom.GetAtomMapNum())

    for bond in reaction_molecule.GetBonds():
        if bond.GetBeginAtom().GetAtomMapNum() in changed_atom_map_num and bond.GetEndAtom().GetAtomMapNum() in changed_atom_map_num:
            # changed_bonds.add((bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()))
            changed_bonds.add(bond)

    return changed_atoms, changed_bonds


def find_matching_product_atom(reactant_atom: Atom, products: List[Mol]):
    """
    Finds the matching product atom for a given reaction atom reactant_atom
    Returns the matching product atom
    :param reactant_atom:
    :param products:
    :return:
    """

    for product in products:
        for atom in product.GetAtoms():
            if reactant_atom.GetAtomMapNum() == atom.GetAtomMapNum():
                return atom
    return None


def compare_props(reactant_atom: Atom, product_atom: Atom):
    """
    Compares a fixed atom when it is in the reactant and when it is in the product
    Returns boolean (whether atoms are the same or different)
    :param reactant_atom:
    :param product_atom:
    :return:
    """

    if reactant_atom.GetDegree() != product_atom.GetDegree():
        return False
    if reactant_atom.GetFormalCharge() != product_atom.GetFormalCharge():
        return False
    if reactant_atom.GetNumRadicalElectrons() != product_atom.GetNumRadicalElectrons():
        return False

    reactant_atom_neighbors = {atom.GetAtomMapNum() for atom in reactant_atom.GetNeighbors()}
    product_atom_neighbors = {atom.GetAtomMapNum() for atom in product_atom.GetNeighbors()}
    for reactant_num in reactant_atom_neighbors:
        if reactant_num not in product_atom_neighbors:
            return False

    for bond in reactant_atom.GetBonds():
        if not compare_bond(reactant_atom, product_atom, bond):
            return False

    return True


def compare_bond(reactant_atom: Atom, product_atom: Atom, reactant_bond: Bond):
    """

    :param reactant_atom:
    :param product_atom:
    :param reactant_bond:
    :return:
    """
    end_atom = reactant_bond.GetOtherAtom(reactant_atom)
    for product_bond in product_atom.GetBonds():
        if product_bond.GetOtherAtom(product_atom).GetAtomMapNum() == end_atom.GetAtomMapNum():
            if product_bond.GetBondType() != reactant_bond.GetBondType():
                return False

            return True
    return False


def highlight_reaction_core(mol, changing_atoms, changing_bonds):
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
    return Draw.MolToImage(mol, highlightAtoms=atom_indices, highlightBonds=bond_indices)


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

##################################################################################################
# CORE EXTENSION
##################################################################################################


def extend_reaction_core(reactant_mol: List[Mol], product_mol: List[Mol], reaction_core: Tuple[Set[Atom], Set[Bond]]) -> None:
    """
    Takes the current core in terms of atoms and bonds and extends the reaction core
    Heurstics found in section 2.2 of Route Designer paper
    Returns the new reaction core
    :param reactant_mol:
    :param product_mol:
    :param reaction_core:
    :return:
    """

    # create copy of reaction_core[0] to iterate through
    copy_core = {atom for atom in reaction_core[0]}
    atom_map_nums = {atom.GetAtomMapNum() for atom in reaction_core[0]}
    # extend core
    for atom in copy_core:
        # we first employ some special heursitics for primary bonds
        extend_primary_bonds(atom, reaction_core, atom_map_nums)
        # then we apply the general heurstics
        add_atoms_to_core(atom, reaction_core, atom_map_nums)
    # get leaving groups and add them to the reaction core
    get_leaving_group(reactant_mol, product_mol, reaction_core, atom_map_nums)


def extend_primary_bonds(atom: Atom, reaction_core: Tuple[Set[Atom], Set[Bond]], atom_map_nums: Set[int]) -> None:
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
    for bond in atom.GetBonds():
        new_atom = bond.GetOtherAtom(atom)
        # check for secondary double or triple bonds
        for secondary_bond in new_atom.GetBonds():
            newest_atom = secondary_bond.GetOtherAtom(new_atom)
            if secondary_bond.GetBondType() == BondType.DOUBLE or secondary_bond.GetBondType() == BondType.TRIPLE and \
                    new_atom.GetAtomMapNum() not in atom_map_nums and newest_atom.GetAtomMapNum() not in atom_map_nums:
                reaction_core[0].add(new_atom)
                reaction_core[0].add(newest_atom)
                reaction_core[1].add(bond)
                reaction_core[1].add(secondary_bond)
                atom_map_nums.add(new_atom.GetAtomMapNum())
                atom_map_nums.add(new_atom.GetAtomMapNum())

        # check for aromatic bonds
        if bond.GetIsAromatic():
            get_aromatic_ring(atom, reaction_core, atom_map_nums)


def get_aromatic_ring(atom: Atom, reaction_core: Tuple[Set[Atom], Set[Bond]], atom_map_nums: Set[int]) -> None:
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

        if bond.GetBeginAtom().GetAtomMapNum() not in atom_map_nums and bond.GetEndAtom().GetAtomMapNum() not in atom_map_nums:
            reaction_core[1].add(bond)
            reaction_core[0].add(bond.GetBeginAtom())
            reaction_core[0].add(bond.GetEndAtom())
            atom_map_nums.add(bond.GetBeginAtom().GetAtomMapNum())
            atom_map_nums.add(bond.GetEndAtom().GetAtomMapNum())


def add_atoms_to_core(atom: Atom, reaction_core: Tuple[Set[Atom], Set[Bond]], atom_map_nums: Set[int]) -> None:
    """
    Adds new atoms to the reaction core
    Essentially employs a BFS search on molecule starting at atom
    Mutates reaction core for extension
    :param atom:
    :param reaction_core:
    :return:
    """

    queue = [atom]
    atom_seen_so_far = set()
    while len(queue) != 0:
        curr_atom = queue.pop()
        for bond in curr_atom.GetBonds():
            new_atom = bond.GetOtherAtom(curr_atom)
            if not check_external_nonaromatic_bond(bond) and new_atom not in reaction_core[0] and \
                    bond not in reaction_core[1] and new_atom.GetAtomMapNum() not in atom_seen_so_far and \
                    new_atom.GetAtomMapNum() not in atom_map_nums:
                reaction_core[0].add(new_atom)
                reaction_core[1].add(bond)
                atom_map_nums.add(new_atom.GetAtomMapNum())
                queue.append(new_atom)
            elif check_external_nonaromatic_bond(bond):
                atom_seen_so_far.add(curr_atom.GetAtomMapNum())
        atom_seen_so_far.add(curr_atom.GetAtomMapNum())


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


def get_leaving_group(reactant_mol: List[Mol], product_mol: List[Mol], reaction_core: Tuple[Set[Atom], Set[Bond]],
                      atom_map_nums: Set[int]) -> None:
    """
    Gets leaving group Atoms and Bonds and mutates reaction_core to add them
    :param reactant_mol:
    :param product_mol:
    :param reaction_core:
    :return:
    """

    # we check for bonds being broken in reaction core
    # the only place leaving groups can occur is from bonds in reaction core as their properties change
    # such as neighbors, degree, etc.
    for core_atom in reaction_core[0]:
        # find leaving group
        product_atom = find_matching_product_atom(core_atom, product_mol)
        reactant_neighbors = {atom.GetAtomMapNum() for atom in core_atom.GetNeighbors()}
        product_neighbors = {atom.GetAtomMapNum() for atom in product_atom.GetNeighbors()}
        leaving_map_num = {map_num for map_num in reactant_neighbors if map_num not in product_neighbors and
                           map_num not in atom_map_nums}
        # add leaving group to reaction core
        # add atoms
        leaving_atoms = {find_atom(num, reactant_mol) for num in leaving_map_num}
        for atom in leaving_atoms:
            reaction_core[0].add(atom)
        atom_map_nums = atom_map_nums.union(leaving_map_num)
        # add bonds
        for atom in leaving_atoms:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomMapNum() in leaving_map_num and neighbor.GetAtomMapNum() not in atom_map_nums:
                    reaction_core[0].add(neighbor)
                    atom_map_nums.add(neighbor.GetAtomMapNum())


def find_atom(atom_map_num: int, mols: List[Mol]) -> Atom:
    """
    Finds Atom object that has the same atom map number as atom_map_num
    :param atom_map_num:
    :param mols:
    :return:
    """
    for mol in mols:
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum == atom_map_num:
                return atom
    return None
