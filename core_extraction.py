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

import re

from classes.partial_molecule import ReactionCore
from utils import find_atom, find_bond, highlight_reaction_core

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
    
    replaced_r_group_reactant = [re.sub(r'R\d+', 'Zn', smile) for smile in reactants_smiles]
    replaced_r_group_product = [re.sub(r'R\d+', 'Zn', smile) for smile in products_smiles]
    if any("R" in smiles for smiles in replaced_r_group_reactant) or any ("R" in smiles for smiles in replaced_r_group_product):
        replaced_r_group_reactant = [re.sub(r'R+', 'Zn', smile) for smile in replaced_r_group_reactant]
        replaced_r_group_product = [re.sub(r'R+', 'Zn', smile) for smile in replaced_r_group_product]
    
    reactants_smarts = '.'.join([mol for mol in replaced_r_group_reactant])
    products_smarts = '.'.join([mol for mol in replaced_r_group_product])
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

    # if we have atoms that are left unmapped we need to give them atom map numbers
    next_index_to_add = 0
    for mol in reactant_mol:
        next_index_to_add = max(next_index_to_add, max({atom.GetAtomMapNum() for atom in mol.GetAtoms()}))
    for mol in product_mol:
        next_index_to_add = max(next_index_to_add, max({atom.GetAtomMapNum() for atom in mol.GetAtoms()}))
    # print("Next index to add: " + str(next_index_to_add))
    counter = 1
    for mol in reactant_mol:
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() == 0:
                atom.SetAtomMapNum(next_index_to_add + 1)
                next_index_to_add += 1
                # print("Next index to add: " + str(next_index_to_add))
        highlight_reaction_core(mol, set(), set(), "reactant_template" + str(counter) + ".png")
        counter += 1
    counter = 1
    for mol in product_mol:
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() == 0:
                atom.SetAtomMapNum(next_index_to_add + 1)
                next_index_to_add += 1
        highlight_reaction_core(mol, set(), set(), "product_template" + str(counter) + ".png")
        counter += 1
    
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


##################################################################################################
# CORE EXTENSION
##################################################################################################


def extend_reaction_core(reactant_mol: List[Mol], product_mol: List[Mol], reaction_core: ReactionCore) -> ReactionCore:
    """
    Takes the current core in terms of atoms and bonds and extends the reaction core
    Heurstics found in section 2.2 of Route Designer paper
    Returns the reaction core that corrosponds to the products
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
    
    # in an aptempt to make the reaction core less generic we try to expand by an extra 1 or 2 bonds
    add_extra_bonds(reaction_core, reactant_mol)
    
    # make identical reaction core for products and get SMARTS
    product_core = ReactionCore()
    duplicate_reaction_core_to_product_core(reaction_core, product_core, product_mol)
    # print(reaction_core)
    return product_core

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
        # if new_atom.IsInRing():
        #     get_aromatic_ring(atom, reaction_core, atom_map_nums)


def get_aromatic_ring(atom: Atom, reaction_core: ReactionCore, atom_map_nums: Set[int]) -> None:
    """
    Gets all the atoms and bonds in an aromatic ring
    Mutates the
    :param atom:
    :param reaction_core:
    :return:
    """

    ownning_mol = atom.GetOwningMol()
    aromatic_rings = ownning_mol.GetRingInfo().AtomRings()
    # print(aromatic_rings)

    # find the aromatic ring containing atom
    aromatic_ring_containing_atom = []
    for ring in aromatic_rings:
        if atom.GetIdx() in ring:
            aromatic_ring_containing_atom.append(ring)

    # retrieve bonds connected to atom in the ring
    for ring in aromatic_ring_containing_atom:
        for i, atom_idx in enumerate(ring):
            next_atom_idx = ring[(i + 1) % len(ring)]
            bond = ownning_mol.GetBondBetweenAtoms(atom_idx, next_atom_idx)
            if bond is not None:
                if not reaction_core.check_atom(bond.GetBeginAtom()):
                    atom_map_nums.add(bond.GetBeginAtom().GetAtomMapNum())
                    reaction_core.add_atom(bond.GetBeginAtom())
                # print((bond.GetBeginAtom(), bond.GetEndAtom()))
                if not reaction_core.check_atom(bond.GetEndAtom()):
                    atom_map_nums.add(bond.GetEndAtom().GetAtomMapNum())
                    reaction_core.add_atom(bond.GetEndAtom())
                # print((bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()))
                reaction_core.add_bond(bond)
            
    # for ring in aromatic_ring_containing_atom:
    #     for bond_idx in ring:
    #         bond = ownning_mol.GetBondWithIdx(bond_idx)

    #         if bond.GetBeginAtom().GetAtomMapNum() not in atom_map_nums and bond.GetEndAtom().GetAtomMapNum() not in atom_map_nums:
    #             reaction_core.add_bond(bond)
    #             reaction_core.add_atom(bond.GetBeginAtom())
    #             reaction_core.add_atom(bond.GetEndAtom())
    #             atom_map_nums.add(bond.GetBeginAtom().GetAtomMapNum())
    #             atom_map_nums.add(bond.GetEndAtom().GetAtomMapNum())


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
        for bond in curr_atom.GetBonds():
            new_atom = bond.GetOtherAtom(curr_atom)
            # if never seen before add to core
            if not check_external_nonaromatic_bond(new_atom) and new_atom.GetAtomMapNum() not in atom_seen_so_far and \
                    new_atom.GetAtomMapNum() not in atom_map_nums:
                reaction_core.add_atom(new_atom)
                reaction_core.add_bond(bond)
                atom_map_nums.add(new_atom.GetAtomMapNum())
                queue.append(new_atom)
            # if added previously by primary core extension still add to queue
            elif not check_external_nonaromatic_bond(new_atom) and new_atom.GetAtomMapNum() in atom_map_nums \
                        and new_atom.GetAtomMapNum() not in atom_seen_so_far:
                queue.append(new_atom)
                atom_seen_so_far.add(new_atom)
            # otherwise do not add
            elif check_external_nonaromatic_bond(new_atom):
                atom_seen_so_far.add(new_atom.GetAtomMapNum())
                
        atom_seen_so_far.add(curr_atom.GetAtomMapNum())
    
def add_extra_bonds(reaction_core: ReactionCore, reactant_mol: Mol) -> None:
    """
    Extends reaction core by 1 bond for each atom in reactant core
    Then makes sure to connect any unconnected atoms that should be connected in the reaction core
    """
    # extend by 1 bond
    atoms_to_add = set()
    bonds_to_add = set()
    for atom in reaction_core.atoms:
        for bond in atom.GetBonds():
            other_atom = bond.GetOtherAtom(atom)
            # check against 30 to make sure it is not zinc (R group place holder)
            if not reaction_core.check_bond(bond) and other_atom.GetAtomicNum() != 30:
                bonds_to_add.add(bond)
                if not reaction_core.check_atom(other_atom):
                    atoms_to_add.add(other_atom)
            # for other_bond in other_atom.GetBonds():
            #     other_other_atom = other_bond.GetOtherAtom(other_atom)
            #     if not reaction_core.check_bond(other_bond) and other_other_atom.GetAtomicNum() != 30:
            #         bonds_to_add.add(bond)
            #         if not reaction_core.check_atom(other_other_atom):
            #             atoms_to_add.add(other_other_atom)
    
    reaction_core.add_atoms(atoms_to_add)
    reaction_core.add_bonds(bonds_to_add)
    
    atoms_to_add = set()
    bonds_to_add = set()
    # adds bonds for atoms that are neighbors in reaction core
    for atom_map in reaction_core.atom_map_nums:
        atom = reaction_core.get_atom_from_atom_map_num(atom_map)
        if atom is not None:
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                new_num = neighbor.GetAtomMapNum()
                if reaction_core.check_atom(neighbor) and \
                    not reaction_core.check_bond_map_num(atom_map, neighbor.GetAtomMapNum()):
                        atoms_to_add.add(neighbor)
                        bonds_to_add.add(bond)
                        
    reaction_core.add_atoms(atoms_to_add)
    reaction_core.add_bonds(bonds_to_add)
    # print(reaction_core)

def check_external_nonaromatic_bond(new_atom: Atom) -> bool:
    """
    Checks if current bond is an external nonaromatic carbon-carbon bond
    Returns True of False upon completion
    """
    
    # check if the atom is aromatic
    if new_atom.GetIsAromatic():
        return True
    # check if it is in a ring
    if new_atom.IsInRing():
        return True
    # check if it is external
    # we check 6 for carbon or 30 for zinc (placeholder for R group)
    if new_atom.GetDegree() == 1 and (new_atom.GetAtomicNum() == 6 or new_atom.GetAtomicNum() == 30):
        return True

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
    to_add_atoms = set()
    to_add_bonds = set()
    for core_atom in reaction_core.atoms:
        # find leaving group
        product_atom = find_atom(core_atom.GetAtomMapNum(), product_mol)
        if product_atom is not None:
            reactant_neighbors = {atom.GetAtomMapNum() for atom in core_atom.GetNeighbors()}
            product_neighbors = {atom.GetAtomMapNum() for atom in product_atom.GetNeighbors()}
            leaving_map_num = {map_num for map_num in reactant_neighbors if map_num not in product_neighbors and
                            map_num not in atom_map_nums}
            # add leaving group to reaction core
            # add atoms
            leaving_atoms = {find_atom(num, reactant_mol) for num in leaving_map_num}
            for atom in leaving_atoms:
                to_add_atoms.add(atom)
            atom_map_nums = atom_map_nums.union(leaving_map_num)
            # add bonds
            for atom in leaving_atoms:
                for bond in atom.GetBonds():
                    neighbor = bond.GetOtherAtom(atom)
                    if neighbor.GetAtomMapNum() in leaving_map_num and neighbor.GetAtomMapNum() not in atom_map_nums:
                        to_add_atoms.add(neighbor)
                        to_add_bonds.add(bond)
                        atom_map_nums.add(neighbor.GetAtomMapNum())
    
    # # check to see if some are just simply absent
    # new_atoms_to_add = set()
    # already_added = {atom.GetAtomMapNum() for atom in to_add_atoms}
    # for mol in reactant_mol:
    #     for atom in mol.GetAtoms():
    #         map_num = atom.GetAtomMapNum()
    #         if not find_atom(map_num, product_mol) and map_num not in already_added:
    #             new_atoms_to_add.add(atom)
    # # check again for bonds
    # new_bonds_to_add = set()
    # for mol in reactant_mol:
    #     for bond in mol.GetBonds():
    #         endpoint1, endpoint2 = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
    #         if find_atom(endpoint1, new_atoms_to_add) and find_atom(endpoint2, new_atoms_to_add):
    #             new_bonds_to_add.add(bond)
    
    # to_add_atoms = to_add_atoms.union(new_atoms_to_add)
    # to_add_bonds = to_add_bonds.union(new_bonds_to_add)
    
    for atom, bond in zip(to_add_atoms, to_add_bonds):
        reaction_core.add_atom(atom)
        reaction_core.add_bond(bond)


def duplicate_reaction_core_to_product_core(reactant_core: ReactionCore, product_core: ReactionCore, products: List[Mol]) -> None:
    """
    Takes the reactant core and duplicates it to the product core with the proper bonds and such in the product
    """
    for atom in reactant_core.atoms:
        product_atom = find_atom(atom.GetAtomMapNum(), products)
        if product_atom is not None:
            product_core.add_atom(product_atom)
    
    for product in products:
        for bond in product.GetBonds():
            if product_core.check_atom(bond.GetBeginAtom()) and product_core.check_atom(bond.GetEndAtom()):
                product_core.add_bond(bond)
