import sys
import os
import subprocess


def run_omega_fold(fasta_file):
    current_folder = os.getcwd()

    if not os.path.isdir("OmegaFold"):
        raise Exception("OmegaFold not found. Please clone the repository from")

    omega_fold_command = ["python", "OmegaFold/main.py", fasta_file, current_folder]
    subprocess.run(omega_fold_command, capture_output=True, encoding="utf-8")
    result_file = fasta_file.split(".")[0] + ".pdb"

    if not os.path.isfile(result_file):
        raise Exception("OmegaFold failed to generate a PDB file")

    result_file_model = fasta_file.split(".")[0] + "_omegaFold.pdb"

    if os.path.isfile(result_file_model):
        os.remove(result_file_model)

    os.rename(result_file, result_file_model)


def main(fasta_file):
    run_omega_fold(fasta_file)
    print(f"OmegaFold successfully generated a PDB file for {fasta_file}.")


if __name__ == "__main__":
    assert len(sys.argv) == 2, "Usage: python omega_fold.py <fasta_file>"
    fasta_file = sys.argv[1]
    main(fasta_file)
