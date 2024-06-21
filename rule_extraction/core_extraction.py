"""
Core extraction. Code provided extracts core from reactions
Implemented Herustics from Route Designer
"""
from rdkit import Chem
from rdkit.Chem.rdchem import Mol, Atom, Bond, BondType
from rdkit.Chem import Draw
from rdkit.Chem.rdChemReactions import ChemicalReaction
from typing import List, Set, Tuple
from rxnmapper import RXNMapper

from utils import ReactionCore

##################################################################################################
# CORE EXTRACTION
##################################################################################################


def get_reaction_core(reactants_smiles: List[str], products_smiles: List[str]) ->Tuple[ReactionCore, List[Mol], List[Mol]]:
    """
    Runs the full reaction core
    Returns a tuple:
    first index is a ReactionCore object
    second index is list of reactant in Mol objects
    third is list of products in Mol objects
    """
    reactants_smarts = '.'.join([smiles for smiles in reactants_smiles])
    products_smarts = '.'.join([smiles for smiles in products_smiles])
    reaction_smarts_without_atom_map = reactants_smarts + ">>" + products_smarts
    
    # add atom map nums to reaction smarts
    rxn_mapper = RXNMapper()
    reaction_smarts = rxn_mapper.get_attention_guided_atom_maps([reaction_smarts_without_atom_map])[0]["mapped_rxn"]  
    return get_reaction_core_with_smarts(reaction_smarts)


def get_reaction_core_with_smarts(reaction_smarts: str) -> Tuple[ReactionCore, List[Mol], List[Mol]]:
    """
    Gets the reaction core given the reaction smarts
    Essentially just creates ChemicalReaction object and passes to get_reaction_core_helper
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
    core = get_reaction_core_helper(reactant_mol, product_mol)
    return core, reactant_mol, product_mol


def get_reaction_core_helper(reactant_molecule_list: List[Mol], product_mols: List[Mol]) -> ReactionCore:
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
    
    core = ReactionCore()
    core.add_atoms(changed_atoms)
    core.add_bonds(changed_bonds)

    return core


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
    return Draw.MolToFile(mol, "testing.png", highlightAtoms=atom_indices, highlightBonds=bond_indices)


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


def extend_reaction_core(reactant_mol: List[Mol], product_mol: List[Mol], reaction_core: ReactionCore) -> str:
    """
    Takes the current core in terms of atoms and bonds and extends the reaction core
    Heurstics found in section 2.2 of Route Designer paper
    Returns the reaction SMARTS for the reaction core
    """

    # create copy of reaction_core[0] to iterate through
    copy_core = {atom for atom in reaction_core.atoms}
    atom_map_nums = {atom_map_num for atom_map_num in reaction_core.atom_map_nums}
    # extend core
    for atom in copy_core:
        # we first employ some special heursitics for primary bonds
        extend_primary_bonds(atom, reaction_core, atom_map_nums)
        # then we apply the general heurstics
        add_atoms_to_core(atom, reaction_core, atom_map_nums)
    # get leaving groups and add them to the reaction core
    get_leaving_group(reactant_mol, product_mol, reaction_core, atom_map_nums)
    
    # make identical reaction core for products and get SMARTS
    product_core = ReactionCore()
    duplicate_reaction_core_to_product_core(reaction_core, product_core, product_mol)
    return reaction_core.get_smarts() + ">>" + product_core.get_smarts()


def extend_primary_bonds(atom: Atom, reaction_core: ReactionCore, atom_map_nums: Set[int]) -> None:
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
            if (secondary_bond.GetBondType() == BondType.DOUBLE or secondary_bond.GetBondType() == BondType.TRIPLE) and \
                    new_atom.GetAtomMapNum() not in atom_map_nums and newest_atom.GetAtomMapNum() not in atom_map_nums:
                reaction_core.add_atom(new_atom)
                reaction_core.add_atom(newest_atom)
                # print("Added " + str((atom.GetAtomMapNum(), new_atom.GetAtomMapNum())) + " through extension")
                # print("Added " + str((new_atom.GetAtomMapNum(), newest_atom.GetAtomMapNum())) + " through extension")
                reaction_core.add_bond(bond)
                reaction_core.add_bond(secondary_bond)
                atom_map_nums.add(new_atom.GetAtomMapNum())
                atom_map_nums.add(new_atom.GetAtomMapNum())

        # check for aromatic bonds
        if bond.GetIsAromatic():
            get_aromatic_ring(atom, reaction_core, atom_map_nums)


def get_aromatic_ring(atom: Atom, reaction_core: ReactionCore, atom_map_nums: Set[int]) -> None:
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
            reaction_core.add_bond(bond)
            reaction_core.add_atom(bond.GetBeginAtom())
            reaction_core.add_atom(bond.GetEndAtom())
            atom_map_nums.add(bond.GetBeginAtom().GetAtomMapNum())
            atom_map_nums.add(bond.GetEndAtom().GetAtomMapNum())


def add_atoms_to_core(atom: Atom, reaction_core: ReactionCore, atom_map_nums: Set[int]) -> None:
    """
    Adds new atoms to the reaction core
    Essentially employs a DFS search on molecule starting at atom
    Mutates reaction core for extension
    :param atom:
    :param reaction_core:
    :return:
    """

    queue = [atom]
    atom_seen_so_far = set()
    while len(queue) != 0:
        curr_atom = queue.pop()
        print("Checking atom " + str(curr_atom.GetAtomMapNum()))
        for bond in curr_atom.GetBonds():
            new_atom = bond.GetOtherAtom(curr_atom)
            print("Checking bond with " + str(new_atom.GetAtomMapNum()))
            # if never seen before add to core
            if not check_external_nonaromatic_bond(bond, new_atom) and new_atom.GetAtomMapNum() not in atom_seen_so_far and \
                    new_atom.GetAtomMapNum() not in atom_map_nums:
                reaction_core.add_atom(new_atom)
                reaction_core.add_bond(bond)
                atom_map_nums.add(new_atom.GetAtomMapNum())
                print("Appending " + str(new_atom.GetAtomMapNum()))
                queue.append(new_atom)
            # if added previously by primary core extension still add to queue
            elif not check_external_nonaromatic_bond(bond, new_atom) and new_atom.GetAtomMapNum() in atom_map_nums \
                        and new_atom.GetAtomMapNum() not in atom_seen_so_far:
                queue.append(new_atom)
                atom_seen_so_far.add(new_atom)
            # otherwise do not add
            elif check_external_nonaromatic_bond(bond, new_atom):
                atom_seen_so_far.add(new_atom.GetAtomMapNum())
                
        atom_seen_so_far.add(curr_atom.GetAtomMapNum())


def check_external_nonaromatic_bond(bond: Bond, new_atom: Atom) -> bool:
    """
    Checks if current bond is an external nonaromatic carbon-carbon bond
    Returns True of False upon completion
    :param bond:
    :return: boolean
    """
    
    # check if the atom is aromatic
    if new_atom.GetIsAromatic():
        # print("Bond between " + str(atom1.GetAtomMapNum()) + " and " + str(atom2.GetAtomMapNum()) + " fails on aromatic")
        print("Atom " + str(new_atom.GetAtomMapNum()) + " fails on aromatic")
        return True
    # check if it is in a ring
    if new_atom.IsInRing():
        print("Atom " + str(new_atom.GetAtomMapNum()) + " fails on ring")
        return True
    # TODO check if it is external
    if new_atom.GetDegree() == 1:
        return True
    
    # print("Bond between " + str(atom1.GetAtomMapNum()) + " and " + str(atom2.GetAtomMapNum()) + " passes")
    return False


def get_leaving_group(reactant_mol: List[Mol], product_mol: List[Mol], reaction_core: ReactionCore,
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
    for core_atom in reaction_core.atoms:
        # find leaving group
        product_atom = find_atom(core_atom.GetAtomMapNum(), product_mol)
        reactant_neighbors = {atom.GetAtomMapNum() for atom in core_atom.GetNeighbors()}
        product_neighbors = {atom.GetAtomMapNum() for atom in product_atom.GetNeighbors()}
        leaving_map_num = {map_num for map_num in reactant_neighbors if map_num not in product_neighbors and
                           map_num not in atom_map_nums}
        # add leaving group to reaction core
        # add atoms
        leaving_atoms = {find_atom(num, reactant_mol) for num in leaving_map_num}
        for atom in leaving_atoms:
            reaction_core.add_atom(atom)
        atom_map_nums = atom_map_nums.union(leaving_map_num)
        # add bonds
        for atom in leaving_atoms:
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetAtomMapNum() in leaving_map_num and neighbor.GetAtomMapNum() not in atom_map_nums:
                    reaction_core.add_atom(neighbor)
                    reaction_core.add_bond(bond)
                    atom_map_nums.add(neighbor.GetAtomMapNum())


def duplicate_reaction_core_to_product_core(reactant_core: ReactionCore, product_core: ReactionCore, products: List[Mol]) -> None:
    """
    Takes the reactant core and duplicates it to the product core with the proper bonds and such in the product
    """
    for atom in reactant_core.atoms:
        product_atom = find_atom(atom.GetAtomMapNum(), products)
        product_core.add_atom(product_atom)
    
    for bond in reactant_core.bonds:
        product_bond = find_bond(bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum(), products)
        if product_core.check_atom(product_bond.GetBeginAtom()) and product_core.check_atom(product_bond.GetEndAtom()):
            product_core.add_bond(product_bond)


def find_atom(atom_map_num: int, mols: List[Mol]) -> Atom:
    """
    Finds Atom object that has the same atom map number as atom_map_num
    :param atom_map_num:
    :param mols:
    :return:
    """
    for mol in mols:
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() == atom_map_num:
                return atom
    return None


def find_bond(index1: int, index2: int, mols: List[Mol]) -> Bond:
    """
    Finds Bond object that has same bond map numbers as index1 and index2
    """
    
    for mol in mols:
        for bond in mol.GetBonds():
            possible1, possible2 = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
            if (index1 == possible1 and index2 == possible2) or (index1 == possible2 and index2 == possible1):
                return bond
    return None
