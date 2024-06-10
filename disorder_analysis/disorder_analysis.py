import argparse
import os
import subprocess
import uuid

import requests
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from Bio import ExPASy, SwissProt
from pymol import cmd


def get_seq_from_pdb(pdb_file):
    print("Reading PDB file...")

    parser = PDBParser()
    structure = parser.get_structure('X', pdb_file)
    sequences = {}
    for model in structure:
        for chain in model:
            seq = ''
            for residue in chain:
                if residue.id[0] == ' ':
                    seq += seq1(residue.resname)
            sequences[chain.id] = seq
    return sequences


def get_sequence_from_uniprot(uniprot_id):
    print(f"Downloading sequence from Uniprot for {uniprot_id}...")
    handle = ExPASy.get_sprot_raw(uniprot_id)
    record = SwissProt.read(handle)
    return str(record.sequence)


def download_pdb(uniprot_id):
    print("Downloading PDB structure...")
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    response = requests.get(url)
    output_file = "disorder_analysis/pdb/" + uniprot_id + ".pdb"
    if response.status_code == 200:
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        with open(output_file, 'wb') as f:
            f.write(response.content)
    else:
        raise Exception(f"Failed to fetch PDB file for UniProt ID {uniprot_id}: {response.status_code}")

    return output_file


def run_iupred(sequences):
    print("Running IUPred...")
    file = "tmp" + str(uuid.uuid4()) + ".fasta"
    result = {}
    for key in list(sequences.keys()):
        with open(file, "w") as tmp:
            tmp.write(f'>{key}\n')
            tmp.write(f'{sequences[key]}\n')

        iupred_command = ["python", "disorder_analysis/iupred2a/iupred2a.py", file, "short"]
        out = subprocess.run(iupred_command, capture_output=True, encoding="utf-8").stdout
        out = out.split("POS\tRES\tIUPRED2\n")[1].split("\n")
        disorder = ""
        for line in out:
            if line != "":
                float_val = float(line.split("\t")[2])
                if float_val < 0.5:
                    disorder += "-"
                else:
                    disorder += "D"
        result[key] = disorder
    os.remove(file)
    return result


def color_pdb_by_disorder(pdb_file, iupred_result, output_file):
    print("Colouring structure...")
    print("\n")
    cmd.load(pdb_file)
    for chain_id, dis in iupred_result.items():
        for i, d in enumerate(dis):
            if d == "D":
                resi = i + 1
                selection = f"chain {chain_id} and resi {resi}"
                cmd.color("red", selection)
    cmd.save(output_file)


def main(input):
    print("Disorder Analysis and Visualization")
    print("\n")

    result_dir_path = "disorder_analysis/result"
    if not os.path.isdir(result_dir_path):
        os.makedirs(result_dir_path)

    if input.endswith('.pdb'):
        if not os.path.exists(input):
            print(f"PDB file {input} does not exist.")
            exit(1)
        sequences = get_seq_from_pdb(input)
        pdb_file = input
        if '/' in input:
            output = f"{result_dir_path}/{input.split('/')[-1].split('.')[0]}.pse"
        else:
            output = f"{result_dir_path}/{input.split('.')[0]}.pse"
    else:
        sequences = {'A': get_sequence_from_uniprot(input)}
        pdb_file = download_pdb(input)
        output = f"{result_dir_path}/{input}.pse"

    # Run IUPred
    iupred_result = run_iupred(sequences)

    # Color PDB
    color_pdb_by_disorder(pdb_file, iupred_result, output)

    print("Success!")
    print(f"Result file: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Disorder prediction and visualization for provided Uniprot ID/PDB file.')
    parser.add_argument('input', help='Input (UniProt ID or path to PDB file)')

    args = parser.parse_args()
    main(args.input)
