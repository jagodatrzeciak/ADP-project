import sys
import requests

from Bio import SeqIO

def run_esmfold(sequence):
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    response = requests.post(url, data=sequence, verify=False)
    if response.status_code == 200:
        pdb_data = response.text
        return pdb_data
    else:
        print(f"Failed to fold sequence: {response.text}")

def main(fasta_file):
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq).upper()
        pdb_data = run_esmfold(sequence)
        if pdb_data:
            pdb_file_name = f"{fasta_file.split('.')[0]}_esm.pdb"

            with open(pdb_file_name, "w") as pdb_file:
                pdb_file.write(pdb_data)
            
            print(f"Saved PDB output for {record.id} to {pdb_file_name}")

        break

if __name__ == "__main__":
    assert len(sys.argv) == 2, "Usage: python esm.py <fasta_file>"
    fasta_file = sys.argv[1]
    main(fasta_file)