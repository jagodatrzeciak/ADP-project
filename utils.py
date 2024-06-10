import requests
import os
import subprocess

from Bio import SeqIO


def rmsd_calculation(pdb_files):
    """
    Open the Pymol windows with the pdb files in silent mode,
    make superposition and calculate RMSD if more than one file.
    :param pdb_files: List of PDB file paths
    :return: RMSD value if more than one file, None otherwise
    """
    script_content = ""
    for i, pdb in enumerate(pdb_files):
        script_content += f"load {pdb}, model{i + 1}\n"

    if len(pdb_files) > 1:
        script_content += "super " + ", ".join([f"model{i + 1}" for i in range(len(pdb_files))]) + "\n"
        script_content += "print('RMSD:', cmd.rms_cur('model1', 'model2'))\n"

    script_path = "temp_pymol_script.pml"
    with open(script_path, 'w') as script_file:
        script_file.write(script_content)

    result = subprocess.run(["pymol", "-cq", script_path], capture_output=True, text=True)
    os.remove(script_path)
    rmsd = None
    if "RMSD: " in result.stdout:
        rmsd = float(result.stdout.split("RMSD: ")[1].strip())
        print(f"RMSD: {rmsd}")
    return rmsd


def check_fasta_correctness(fasta_sequence):
    """
    Check if the FASTA sequence is correctly formatted. ('>' at the beginning, only standars amino acids)
    :param fasta_sequence:
    :return: 1 if the FASTA sequence is correctly formatted, 0 if it contains non-standard amino acids, -1 if there is uncorrect formatting
    """
    fasta_sequence = fasta_sequence.strip().split('\n')
    if len(fasta_sequence) < 2 or fasta_sequence[0][0] != '>':
        print(fasta_sequence[0], len(fasta_sequence))
        return -1
    seq = "".join(fasta_sequence[1:]).upper()
    for aa in seq:
        if aa not in 'ACDEFGHIKLMNPQRSTVWY':
            return 0
    return 1


def fasta_from_pdb(pdb_file_path):
    """
    Extract the FASTA sequence from the PDB file.
    :param pdb_file_path:
    :return: path tp fasta file
    """
    fasta_sequence = ''
    try:
        with open(pdb_file_path, 'r') as pdb_file:
            record = SeqIO.parse(pdb_file, 'pdb-atom')
            record = next(record)
            fasta_sequence = f">{pdb_file_path.split('_')[0]}\n{str(record.seq)}"
        with open(f"results/{pdb_file_path.split('.')[0]}.fasta", "w") as fasta_file:
            fasta_file.write(fasta_sequence)
        return f"results/{pdb_file_path.split('.')[0]}.fasta"
    except:
        return None


def write_fasta_file_to_temp(fasta_sequence):
    """
    Write the FASTA sequence to a temporary file.
    :param fasta_sequence:
    :return: The path to the temporary file
    """
    header = fasta_sequence.strip().split('\n')[0][1:]
    with open(f"results/{header}.fasta", "w") as fasta_file:
        fasta_file.write(fasta_sequence)
    return f"results/{header}.fasta"


def pdb_correctness(pdb_id):
    """
    Process the PDB ID input.
    :param pdb_id:
    :return: path_to_fasta_file if the PDB ID is correct, None otherwise
    """
    pdb_id = pdb_id.strip()
    try:
        fasta = requests.get(f'https://www.rcsb.org/fasta/entry/{pdb_id}')
        with open(f"results/{pdb_id}.fasta", "w") as fasta_file:
            header = fasta.text.split('\n')[0].strip().split("|")[0]
            to_write = f"{header}\n" + "".join(fasta.text.split('\n')[1:])
            fasta_file.write(to_write)
        return f"results/{pdb_id}.fasta"
    except:
        return None


def uniprot_correctness(uniprot_id):
    """
    Process the UniProt ID input.
    :param uniprot_id:
    :return: path_to_fasta_file if the UniProt ID is correct, None otherwise
    """
    try:
        uniprot_id = uniprot_id.strip()
        fasta = requests.get(f'https://www.uniprot.org/uniprot/{uniprot_id}.fasta')
        with open(f"results/{uniprot_id}.fasta", "w") as fasta_file:
            header = fasta.text.split('\n')[0].strip().split("|")[1]
            to_write = f">{header}\n" + "".join(fasta.text.split('\n')[1:])
            fasta_file.write(to_write)
        return f"results/{uniprot_id}.fasta"
    except:
        return None


def get_sequence_length(path_to_fasta):
    """
    Get the length of the sequence from the FASTA file.
    :param path_to_fasta:
    :return: length of the sequence
    """
    for record in SeqIO.parse(path_to_fasta, "fasta"):
        return len(record.seq)
    return 0
