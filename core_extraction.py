"""
Core extraction. Code provided extracts core from reactions
Implemented Herustics from Route Designer
"""
from rdkit import Chem
from rdkit.Chem.rdchem import Mol, Atom, Bond, BondType, EditableMol
from rdkit.Chem import Draw
from rdkit.Chem.rdChemReactions import ChemicalReaction

from typing import List, Set, Tuple

from rxnmapper import RXNMapper
import networkx as nx
from networkx import Graph

import re

from classes.partial_molecule import ReactionCore
from utils import find_atom, find_bond, highlight_reaction_core, create_atom, find_atom_index, REMOVE_BOND_FLAG, ADD_BOND_FLAG, MODIFY_BOND_FLAG

##################################################################################################
# CORE EXTRACTION
##################################################################################################


def get_reaction_core(reactants_smiles: List[str], products_smiles: List[str]) ->Tuple[List[Mol], List[Mol], List[Mol]]:
    """
    Runs the full reaction core
    Returns a tuple:
    first index is a reactant side of reaction core object
    second index is list of reactant in Mol objects
    third is list of products in Mol objects
    """
    
    replaced_r_group_reactant = [re.sub(r'R\d+', 'Zn', smile) for smile in reactants_smiles]
    replaced_r_group_product = [re.sub(r'R\d+', 'Zn', smile) for smile in products_smiles]
    if any("R" in smiles for smiles in replaced_r_group_reactant) or any ("R" in smiles for smiles in replaced_r_group_product):
        replaced_r_group_reactant = [re.sub(r'R+', 'Zn', smile) for smile in replaced_r_group_reactant]
        replaced_r_group_product = [re.sub(r'R+', 'Zn', smile) for smile in replaced_r_group_product]
    
    replaced_x_group_reactant = [re.sub(r'X+', 'Cu', smile) for smile in replaced_r_group_reactant]
    replaced_x_group_product = [re.sub(r'X+', 'Cu', smile) for smile in replaced_r_group_product]
    
    reactants_smarts = '.'.join([mol for mol in replaced_x_group_reactant])
    products_smarts = '.'.join([mol for mol in replaced_x_group_product])
    reaction_smarts_without_atom_map = reactants_smarts + ">>" + products_smarts
    
    # add atom map nums to reaction smarts
    rxn_mapper = RXNMapper()
    results = rxn_mapper.get_attention_guided_atom_maps([reaction_smarts_without_atom_map])
    reaction_smarts = ""
    highest_confidence = 0
    for result in results:
        if result["confidence"] > highest_confidence:
            reaction_smarts = result["mapped_rxn"]
            highest_confidence = result["confidence"]
                
    reaction_smarts = reaction_smarts.replace("Cu", "*")
    return get_reaction_core_with_smarts(reaction_smarts)


def get_reaction_core_with_smarts(reaction_smarts: str) -> Tuple[List[Mol], List[Mol], List[Mol]]:
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
        # highlight_reaction_core(mol, set(), set(), "reactant_template" + str(counter) + ".png")
        counter += 1
    counter = 1
    for mol in product_mol:
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() == 0:
                atom.SetAtomMapNum(next_index_to_add + 1)
                next_index_to_add += 1
        # highlight_reaction_core(mol, set(), set(), "product_template" + str(counter) + ".png")
        counter += 1
    
    # get core
    core = get_reaction_core_helper(reactant_mol, product_mol)
    return core, reactant_mol, product_mol


def get_reaction_core_helper(reactant_molecule_list: List[Mol], product_mols: List[Mol]) -> List[Mol]:
    """
    Identifies the reaction core of the reactants
    Assume RXN property for reactant and product are set
    :param reactant_molecule_list:
    :return:
    """

    # get all the reaction cores
    reaction_cores = []
    for reactant in reactant_molecule_list:
        new_reactant_mol = get_reaction_core_for_single_reactant(reactant, product_mols)
        reaction_cores.append(new_reactant_mol)

    return [mol.GetMol() for mol in reaction_cores]


def get_reaction_core_for_single_reactant(reaction_molecule: Mol, products: List[Mol]) -> EditableMol:
    """
    Gets the reaction core for a single molecule
    Returns a tuple containing the changed atoms and bonds that represent the reaction core
    :param reaction_molecule:
    :param products:
    :return: tuple[set]
    """

    # changed_atoms = set()
    # changed_bonds = set()
    changed_atom_map_num = set()
    new_reactant_mol = EditableMol(Mol())
    reaction_to_core = {}
    for atom in reaction_molecule.GetAtoms():
        index = atom.GetAtomMapNum()
        for product in products:
            for prod_atom in product.GetAtoms():
                # once we found the matching atom, we compare attributes
                if prod_atom.GetAtomMapNum() == index and not compare_props(atom, prod_atom):
                    new_atom = create_atom(atom)
                    reaction_to_core[prod_atom.GetAtomMapNum()] = new_reactant_mol.AddAtom(new_atom)
                    changed_atom_map_num.add(atom.GetAtomMapNum())

    for bond in reaction_molecule.GetBonds():
        if bond.GetBeginAtom().GetAtomMapNum() in changed_atom_map_num and bond.GetEndAtom().GetAtomMapNum() in changed_atom_map_num:
            # changed_bonds.add((bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()))
            index1, index2 = reaction_to_core[bond.GetBeginAtom().GetAtomMapNum()], reaction_to_core[bond.GetEndAtom().GetAtomMapNum()]
            new_reactant_mol.AddBond(index1, index2, bond.GetBondType())

    return new_reactant_mol


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


def extend_reaction_core(reactant_mol: List[Mol], product_mol: List[Mol], reaction_core: List[Mol]) -> Tuple[List[Mol], List[Mol]]:
    """
    Takes the current core in terms of atoms and bonds and extends the reaction core
    Heurstics found in section 2.2 of Route Designer paper
    Returns the reaction core that corrosponds to the products
    """

    # create copy of reaction_core[0] to iterate through
    copy_core = set()
    for mol in reaction_core:
        copy_core = copy_core.union({atom for atom in mol.GetAtoms()})
    
    atom_map_nums = {atom.GetAtomMapNum() for atom in copy_core}
    # extend core
    for atom in copy_core:
        atom_from_reactants = find_atom(atom.GetAtomMapNum(), reactant_mol)
        # find the reaction core Mol instance that is changing
        core_mol = find_atom(atom.GetAtomMapNum(), reaction_core).GetOwningMol()
        core_mol_index = find_atom_index(atom.GetAtomMapNum(), reaction_core)
        core_mol = EditableMol(core_mol)
        atom_index = find_atom(atom.GetAtomMapNum(), reaction_core).GetIdx()
        
        # we first employ some special heursitics for primary bonds
        extend_primary_bonds(atom_from_reactants, core_mol, atom_map_nums, atom_index)
        # then we apply the general heurstics
        add_atoms_to_core(atom, core_mol, atom_map_nums, atom_index)
        reaction_core[core_mol_index] = core_mol.GetMol()
    # get leaving groups and add them to the reaction core
    get_leaving_group(reactant_mol, product_mol, reaction_core, atom_map_nums)
    
    # if the number of reaction core's is not equal to the number of reactants, we connect them
    smiles = ""
    for mol in reaction_core:
        smiles += "." + Chem.MolToSmiles(mol)
    smiles = smiles[1:]
    if smiles.count(".") != len(reactant_mol) - 1:
        reaction_core = [Chem.MolFromSmiles(smile, sanitize=False) for smile in smiles.split(".")]
        combine_cores(reaction_core, reactant_mol)
    # make identical reaction core for products and get SMARTS
    product_core = duplicate_reaction_core_to_product_core(reaction_core, product_mol)
    return reaction_core, product_core


def extend_primary_bonds(atom: Atom, core_mol: EditableMol, atom_map_nums: Set[int], atom_index: int) -> None:
    """
    Adds new atoms to reaction core, or more specifically core_mol
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
            if (secondary_bond.GetBondType() == BondType.DOUBLE or secondary_bond.GetBondType() == BondType.TRIPLE or secondary_bond.GetBondType() == BondType.AROMATIC) and \
                    new_atom.GetAtomMapNum() not in atom_map_nums and newest_atom.GetAtomMapNum() not in atom_map_nums:
                atom1, atom2 = create_atom(new_atom), create_atom(newest_atom)
                index1, index2 = core_mol.AddAtom(atom1), core_mol.AddAtom(atom2)
                core_mol.AddBond(atom_index, index1, bond.GetBondType())
                core_mol.AddBond(index1, index2, secondary_bond.GetBondType())

                atom_map_nums.add(new_atom.GetAtomMapNum())
                atom_map_nums.add(newest_atom.GetAtomMapNum())

        # check for aromatic bonds
        if (new_atom.IsInRing() or new_atom.GetIsAromatic()) and new_atom.GetAtomMapNum() not in atom_map_nums:
            atom1 = create_atom(new_atom)
            index1 = core_mol.AddAtom(atom1)
            core_mol.AddBond(atom_index, index1, bond.GetBondType())
            atom_map_nums.add(new_atom.GetAtomMapNum())
            

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


def add_atoms_to_core(atom: Atom, core_mol: EditableMol, atom_map_nums: Set[int], atom_index: int) -> None:
    """
    Adds new atoms to the reaction core
    Essentially employs a DFS search on molecule starting at atom
    Mutates reaction core for extension
    :param atom:
    :param reaction_core:
    :return:
    """

    queue = [atom]
    queue_indicies = [atom_index]
    atom_seen_so_far = set()
    while len(queue) != 0:
        curr_atom = queue.pop()
        curr_index = queue_indicies.pop()
        for bond in curr_atom.GetBonds():
            new_atom = bond.GetOtherAtom(curr_atom)
            # if never seen before add to core
            if not check_external_nonaromatic_bond(new_atom) and new_atom.GetAtomMapNum() not in atom_seen_so_far and \
                    new_atom.GetAtomMapNum() not in atom_map_nums:
                # reaction_core.add_atom(new_atom)
                atom_to_add = create_atom(new_atom)
                index_to_add = core_mol.AddAtom(atom_to_add)
                core_mol.AddBond(curr_index, index_to_add, bond.GetBondType())
                # reaction_core.add_bond(bond)
                atom_map_nums.add(new_atom.GetAtomMapNum())
                queue.append(new_atom)
                queue_indicies.append(curr_index)
            # if added previously by primary core extension still add to queue
            elif not check_external_nonaromatic_bond(new_atom) and new_atom.GetAtomMapNum() in atom_map_nums \
                        and new_atom.GetAtomMapNum() not in atom_seen_so_far:
                queue.append(new_atom)
                queue_indicies.append(find_atom(new_atom.GetAtomMaoNum(), [core_mol.GetMol()]).GetIdx())
                atom_seen_so_far.add(new_atom.GetAtomMapNum())
            # otherwise do not add
            elif check_external_nonaromatic_bond(new_atom):
                atom_seen_so_far.add(new_atom.GetAtomMapNum())
                
        atom_seen_so_far.add(curr_atom.GetAtomMapNum())
    

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


def get_leaving_group(reactant_mol: List[Mol], product_mol: List[Mol], reaction_core: List[Mol],
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
    for mol in reactant_mol:
        for core_atom in mol.GetAtoms():
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

        matching_core = reaction_core[reactant_mol.index(mol)]
        atom_map_to_core = {}
        for atom in to_add_atoms:
            atom_map_to_core = matching_core.AddAtom(create_atom(atom))
        for bond in to_add_bonds:
            index1, index2 = atom_map_to_core[bond.GetBeginAtom().GetAtomMapNum()], atom_map_to_core[bond.GetEndAtom().GetAtomMapNum()]
            matching_core.AddBond(index1, index2, bond.GetBondType())
        to_add_atoms, to_add_bonds = set(), set()

def duplicate_reaction_core_to_product_core(reactant_core: List[Mol], products: List[Mol]) -> List[Mol]:
    """
    Takes the reactant core and duplicates it to the product core with the proper bonds and such in the product
    """
    product_core = EditableMol(Mol())
    all_atoms = set()    
    for mol in reactant_core:
        all_atoms = all_atoms.union({atom for atom in mol.GetAtoms()})
    
    atom_map_to_index = {}
    for atom in all_atoms:
        product_atom = find_atom(atom.GetAtomMapNum(), products)
        atom_map_to_index[product_atom.GetAtomMapNum()] = product_core.AddAtom(create_atom(product_atom))
    
    for product in products:
        for bond in product.GetBonds():
            if bond.GetBeginAtom().GetAtomMapNum() in atom_map_to_index and bond.GetEndAtom().GetAtomMapNum() in atom_map_to_index:
                product_core.AddBond(atom_map_to_index[bond.GetBeginAtom().GetAtomMapNum()], atom_map_to_index[bond.GetEndAtom().GetAtomMapNum()], bond.GetBondType())
    
    return [product_core.GetMol()]


def combine_cores(reaction_cores: List[Mol], reactants: List[Mol]):
    """
    Checks to see if we can combine core1 and core2 with any of the reactants
    """
    for core1 in reaction_cores:
        for core2 in reaction_cores:
            result = combine_2_fixed_cores(core1, core2, reactants)
            if result is not None:
                reaction_cores.remove(core1)
                reaction_cores.remove(core2)
                reaction_cores.append(result)
                return


def combine_2_fixed_cores(core1: Mol, core2: Mol, reactants: List[Mol]) -> Mol:
    if core1 != core2:
        for reactant in reactants:
            result = connect_two_subgraphs(reactant, core1, core2)
            if result is not None:
                return result
    return None

def connect_two_subgraphs(overall_mol: Mol, mol1: Mol, mol2: Mol) -> Mol:
    """
    Connects mol1 and mol2 if they are apart of the same molecule
    Returns mol1 and mol2 combined with their shortest path
    Returns None if they are not part of the same graph
    """
    graph = generate_graph(overall_mol)
    graph1, graph2 = generate_graph(mol1), generate_graph(mol2)
    min_distance = float("inf")
    min_path = [[]]
    
    for u in graph1.nodes():
        for v in graph2.nodes():
            try:
                distance = nx.shortest_path_length(graph, source=u, target=v)
                if distance < min_distance:
                    min_distance = distance
                    min_path[0] = nx.shortest_path(graph, source=u, target=v)
            except nx.NetworkXNoPath:
                continue
    # if no path is found, return none
    if min_path == []:
        return None
    # first combine mol1 and mol2
    final_mol = EditableMol(Mol())
    map_nums_to_index = {}
    for atom in mol1.GetAtoms():
        map_nums_to_index[atom.GetAtomMapNum()] = final_mol.AddAtom(create_atom(atom))
    for bond in mol1.GetBonds():
        index1, index2 = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
        final_mol.AddBond(map_nums_to_index[index1], map_nums_to_index[index2], bond.GetBondType())
    for atom in mol2.GetAtoms():
        map_nums_to_index[atom.GetAtomMapNum()] = final_mol.AddAtom(create_atom(atom))
    for bond in mol2.GetBonds():
        index1, index2 = bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum()
        final_mol.AddBond(map_nums_to_index[index1], map_nums_to_index[index2], bond.GetBondType())
    
    # now add new path
    for atom_map in min_path[0]:
        if find_atom(atom_map, [final_mol.GetMol()]) is None:
            mathcing_atom = find_atom(atom_map, [overall_mol])
            map_nums_to_index[atom_map] = final_mol.AddAtom(create_atom(mathcing_atom))
            
    
    for i in range(0, len(min_path[0]) - 1):
        index1, index2 = min_path[0][i], min_path[0][i + 1]
        bond = find_bond(index1, index2, [overall_mol])
        final_mol.AddBond(map_nums_to_index[index1], map_nums_to_index[index2], bond.GetBondType())
    
    return final_mol.GetMol()
    
    
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

if __name__ == "__main__":
    
    # wiley2_table_3
    core = get_reaction_core(["O=C(NC)C1=CC=CC=C1/C=C([R1])/[R2]"], ["O=C(N1C)C2=CC=CC=C2C1C([R1])[R2]"])
    reactant_core, product_core = extend_reaction_core(core[1], core[2], core[0])
    highlight_reaction_core(core[1][0], set(), set(), "og_reactant.png")
    highlight_reaction_core(core[2][0], set(), set(), "og_product.png")
    highlight_reaction_core(reactant_core[0], set(), set(), "reactant_core.png")
    highlight_reaction_core(product_core[0], set(), set(), "product_core.png")