import sys
import requests
import os
import gzip
import io

import urllib.request
from Bio import SeqIO

SWISSMODEL_EXPASY_API_TOKEN = os.getenv("SWISSMODEL_EXPASY_API_TOKEN")

# Based on the code from https://swissmodel.expasy.org/docs/help#modelling_api
def run_swissmodel(sequence):
    response = requests.post(
        "https://swissmodel.expasy.org/automodel",
        headers={ "Authorization": f"Token {SWISSMODEL_EXPASY_API_TOKEN}" },
        json={ 
            "target_sequences": 
                [
                    sequence
                ],
    })
    
    # Obtain the project_id from the response created above
    project_id = response.json()["project_id"]

    # And loop until the project completes
    import time
    while True:
        # We wait for some time
        time.sleep(10)

        # Update the status from the server 
        response = requests.get(
            f"https://swissmodel.expasy.org/project/{ project_id }/models/summary/", 
            headers={ "Authorization": f"Token {SWISSMODEL_EXPASY_API_TOKEN}" })

        # Update the status
        status = response.json()["status"]

        print('Job status is now', status)

        if status in ["COMPLETED", "FAILED"]:
            break

    response_object = response.json()
    if response_object['status']=='COMPLETED':
        pdb_url = response_object['models'][0]['coordinates_url']
        pdb_data_gz = urllib.request.urlopen(pdb_url).read()
        pdb_data = gzip.GzipFile(fileobj=io.BytesIO(pdb_data_gz)).read().decode('utf-8')
        return pdb_data
    else:
        print('Model building failed')

def main(fasta_file):
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq)
        pdb_data = run_swissmodel(sequence)
        if pdb_data:
            pdb_file_name = f"{fasta_file.split('.')[0]}_swissmodel.pdb"

            with open(pdb_file_name, "w") as pdb_file:
                pdb_file.write(pdb_data)
            
            print(f"Saved PDB output for {record.id} to {pdb_file_name}")
        break

if __name__ == "__main__":
    assert len(sys.argv) == 2, "Usage: python swissmodel.py <fasta_file>"
    fasta_file = sys.argv[1]

    if not SWISSMODEL_EXPASY_API_TOKEN:
        print("Please set the SWISSMODEL_EXPASY_API_TOKEN environment variable")
        sys.exit(1)

    main(fasta_file)