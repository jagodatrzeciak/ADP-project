import requests
# copilot please change wget to requests


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


def write_fasta_file_to_temp(fasta_sequence):
    """
    Write the FASTA sequence to a temporary file.
    :param fasta_sequence:
    :return: The path to the temporary file
    """
    header = fasta_sequence.strip().split('\n')[0][1:]
    with open(f"{header}.fasta", "w") as fasta_file:
        fasta_file.write(fasta_sequence)
    return fasta_file


def pdb_correctness(pdb_id):
    """
    Process the PDB ID input.
    :param pdb_id:
    :return: path_to_fasta_file if the PDB ID is correct, None otherwise
    """
    pdb_id = pdb_id.strip()
    try:
        fasta = requests.get(f'https://www.rcsb.org/fasta/entry/{pdb_id}')
        with open(f"{pdb_id}.fasta", "w") as fasta_file:
            fasta_file.write(fasta.text)
        return fasta_file.name
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
        with open(f"{uniprot_id}.fasta", "w") as fasta_file:
            fasta_file.write(fasta.text)
        return fasta_file.name
    except:
        return None
