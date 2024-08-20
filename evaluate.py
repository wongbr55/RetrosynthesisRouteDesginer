"""
File that contains code to evalute algorithim
Contains functionaility for both exact string mathcing using canonical SMILES and Tanimoto similiarity
"""

from rdkit import Chem, DataStructs
from typing import Set, Tuple, List, Dict
from core_extraction import get_reaction_core, extend_reaction_core
from template_rule_retrieval import get_reactants_for_substrate
from reverse_reaction import get_reactants

import json
import csv


def evaluate_reverse_reaction(ground_truth_csv: str, reaction_template_csv: str) -> Dict:
    """
    Evaluates reverser_reaction.get_reactants
    """
    # get the reaction templates
    reaction_templates = {}
    filenames = []
    with open(reaction_template_csv, "r") as csvfile:
        templates = csv.reader(csvfile)
        next(templates)
        for row in templates:
            image_num = row[0]
            filenames.append(row[1])
            products = []
            reactants = []
            if row[2] != "NULL": products.append(row[2])
            if row[3] != "NULL": products.append(row[3])
            if row[4] != "NULL": reactants.append(row[4])
            if row[5] != "NULL": reactants.append(row[5])
            if row[6] != "NULL": reactants.append(row[6])
            reaction_templates[image_num] = (reactants, products)
    # get the evalutaions for each substrate
    evaluations = {}
    with open(ground_truth_csv, "r") as csvfile:
        ground_truths = csv.reader(csvfile)
        next(ground_truths)
        for row in ground_truths:
            template = []
            image_num = row[0]
            template.append(reaction_templates[image_num][0])
            template.append(reaction_templates[image_num][1])
            
            substrate = row[3]
            ground_truth_reactants = []
            if row[5] != "NULL": ground_truth_reactants.extend(row[5].split("."))
            if row[6] != "NULL": ground_truth_reactants.extend(row[6].split("."))
            if row[7] != "NULL": ground_truth_reactants.extend(row[7].split("."))
            try:
                generated_reactants = get_reactants(substrate, template[0], template[1])
                # if image_num == "16":
                    
                evaluation = evaluate(ground_truth_reactants, generated_reactants)
                evaluations[row[1] + " (" + row[2] + ")"] = evaluation
            except:
                evaluations[row[1] + " (" + row[2] + ")"] = "ERROR"

    return evaluations, filenames

def evaluate_csv(ground_truth_csv: str, reaction_template_csv: str) -> Dict[str, Tuple[List[List[int]], List[float]]]:
    """
    Evaluates a csv file using exact string matching and tanimoto accuracy heurstics from evaluate function
    Assume csv format is of form: Figure, Product 1 SMILES, Product 2 SMILES, Reactant 1 SMILES, Reactant 2 SMILES, Reactant 3 SMILES
    Returns a dict of evaluation metrics for each substrate
    """
    
    # reaction cores maps a file to its reactant and product cores
    # to be used when evalauting files in the ground truth
    reaction_cores = {}
    with open(reaction_template_csv, "r") as csvfile:
        reaction_templates = csv.reader(csvfile)
        next(reaction_templates)
        for row in reaction_templates:
            filename = row[0]
            products = []
            reactants = []
            if row[1] != "": products.append(row[1])
            if row[2] != "": products.append(row[2])
            if row[3] != "": reactants.append(row[3])
            if row[4] != "": reactants.append(row[4])
            if row[5] != "": reactants.append(row[5])
            try:
                core = get_reaction_core(reactants, products)
                product_core = extend_reaction_core(core[1], core[2], core[0])
                reaction_cores[filename] = (core[0], product_core)
            except:
                reaction_cores[filename] = "ERROR UPON CORE EXTRACTION"
    
    evaluations = {}
    with open(ground_truth_csv, "r") as csvfile:
        ground_truths = csv.reader(csvfile)
        next(ground_truths)
        for row in ground_truths:
            # get the reaction cores:
            cores = []
            filename = row[0]
            for file in reaction_cores:
                if file in filename and reaction_cores[file] is not str:
                    cores.append(reaction_cores[file][0])
                    cores.append(reaction_cores[file][1])
                    break
                elif file in filename and reaction_cores[file] is str:
                    break
            if cores == []:
                evaluations[filename] = "UNABLE TO RETREIEVE CORE"
            else:
            # get generated reactants and evaluate againsy ground truth
                substrate = row[1]
                ground_truth_reactants = []
                if row[3] != "": ground_truth_reactants.append(row[3])
                if row[4] != "": ground_truth_reactants.append(row[4])
                if row[5] != "": ground_truth_reactants.append(row[5])
                try:
                    generated_reactants = get_reactants_for_substrate(substrate, cores[0], cores[1])
                    generated_reactant_smiles = [Chem.MolToSmiles(mol) for mol in generated_reactants]
                    # print(generated_reactant_smiles)
                    evaluation = evaluate(ground_truth_reactants, generated_reactant_smiles)
                    evaluations[filename] = evaluation
                except:
                    evaluations[filename] = "ERROR UPON FRAGMENTATION/REPAIRMENT"
    
    with open("evaluations.json", "w") as file:
        json.dump(evaluations, file)
    return evaluations

def evaluate(ground_truth_reactants: List[str], generated_reactants: List[str]) -> Tuple[List[List[int]], List[float]]:
    """
    Uses exact canonical string match and tanimoto similarity to evalute algorithim for a single reaction
    Returns a tuple of similarity matricies for exact string matching and tanimoto similarity (in that order)
    """
    # get canonical strings for both ground truth and generated
    canonical_ground_truth = []
    canonical_generated = []
    for smiles in ground_truth_reactants:
        canonical_ground_truth.append(canonicalize_smiles(smiles))
    for smiles in generated_reactants:
        canonical_generated.append(canonicalize_smiles(smiles))
    # print(canonical_generated)
    
    # get smilarity matricies
    exact_match_matrix = []
    tanimoto_matrix = []
    for generated in canonical_generated:
        curr_exact_match = []
        curr_tanimoto = 0
        for ground_truth in canonical_ground_truth:
            if generated == ground_truth:
                curr_exact_match.append(1)
            else:
                curr_exact_match.append(0)
            similarity = tanimoto_similarity(generated, ground_truth)
            if similarity == 1:
                curr_tanimoto += 1
        exact_match_matrix.append(curr_exact_match)
        tanimoto_matrix.append(curr_tanimoto / len(canonical_ground_truth))
    
    return {"canonical match": exact_match_matrix, "tanimoto similarity": tanimoto_matrix, "generated": canonical_generated, "ground truth": canonical_ground_truth}
            

def canonicalize_smiles(smiles: str) -> str:
    """
    Returns the canonicalized form of smiles
    Returns "" upon error
    """
    # print("SMILES going in: " + smiles)
    try:
        canon_smiles = Chem.CanonSmiles(smiles)
    except:
        canon_smiles = ""
    return canon_smiles
    

def tanimoto_similarity(smiles1: str, smiles2: str) -> float:
    """
    Returns the tanimoto similarity of smiles1 and smiles2
    Returns -1 upon error
    """
    try:
        mol1, mol2 = Chem.MolFromSmiles(smiles1), Chem.MolFromSmiles(smiles2)
        fp1, fp2 = Chem.RDKFingerprint(mol1), Chem.RDKFingerprint(mol2)
        return DataStructs.FingerprintSimilarity(fp1, fp2)
    except:
        return -1
    
if __name__ == "__main__":
    results, file_names = evaluate_reverse_reaction("Ground_truth_evaluation set_substrate_scope.csv", "Ground_truth_evaluation_set_reaction_diagram.csv")
    canonical_correct = 0
    canonical_total = 0
    tanimoto_correct = 0
    tanimoto_total = 0
    num_errors = 0
    
    
    file_accuracy = {key: {"canonical correct": 0,
                           "canonical total": 0, 
                           "tanimoto correct": 0, 
                           "tanimoto total": 0,
                           "num errors": 0,
                           "total": 0} for key in file_names}
    
    for key in results:
        
        file_for_curr = ""
        for file in file_accuracy:
            if file in key:
                file_for_curr = file
                break

        if results[key] != "ERROR":
            for match in results[key]["canonical match"]:
                if 1 in match:
                    file_accuracy[file_for_curr]["canonical correct"] += 1
                    canonical_correct += 1
                file_accuracy[file_for_curr]["canonical total"] += 1
                canonical_total += 1
            for match in results[key]["tanimoto similarity"]:
                if match >= 1 / len(results[key]["tanimoto similarity"]):
                    file_accuracy[file_for_curr]["tanimoto correct"] += 1
                    tanimoto_correct += 1
                file_accuracy[file_for_curr]["tanimoto total"] += 1
                tanimoto_total += 1
        else:
            file_accuracy[file_for_curr]["num errors"] += 1
            num_errors += 1
        file_accuracy[file_for_curr]["total"] += 1
            
    
    for key in file_accuracy:
        if file_accuracy[key]["canonical total"] > 0:
            canonical_accuracy = str(file_accuracy[key]["canonical correct"] / file_accuracy[key]["canonical total"])
        else:
            canonical_accuracy = "0"
        if file_accuracy[key]["tanimoto total"] > 0:
            tanimoto_accuracy = str(file_accuracy[key]["tanimoto correct"] / file_accuracy[key]["tanimoto total"])
        else:
            tanimoto_accuracy = "0"
        percent_errors = str(file_accuracy[key]["num errors"] / file_accuracy[key]["total"])
        print(key + ": canonical accuracy: " + canonical_accuracy + \
            ", tanimoto accuracy: " + tanimoto_accuracy + ", percent errors: " + percent_errors)
    
    print("##### RESULTS #####")
    print("Canonical accuracy: " + str(canonical_correct / canonical_total))
    print("Tanimoto accuracy: " + str(tanimoto_correct / tanimoto_total))
    print("Total num with errors: " + str(num_errors))
        
    with open("evaluations_for_reverse_reaction.json", "w") as file:
        json.dump(results, file)