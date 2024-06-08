import tkinter as tk
from tkinter import messagebox
import subprocess
from utils import check_fasta_correctness, write_fasta_file_to_temp, pdb_correctness, uniprot_correctness


def process_input():
    input_type = input_type_var.get()
    if input_type == "FASTA":
        fasta_sequence = fasta_text.get("1.0", tk.END).strip()
        if not fasta_sequence:
            messagebox.showerror("Input Error", "Please provide a FASTA sequence.")
        else:
            if check_fasta_correctness(fasta_sequence) == -1:
                messagebox.showerror("Input Error", "FASTA sequence is not correctly formatted.")
            elif check_fasta_correctness(fasta_sequence) == 0:
                messagebox.showerror("Input Error", "FASTA sequence contains non-standard amino acids.")
            else:
                path_to_fasta = write_fasta_file_to_temp(fasta_sequence)
                messagebox.showinfo("Success", "FASTA sequence processed successfully!")
                show_model_selection(path_to_fasta)
    elif input_type == "PDB":
        pdb_id = pdb_entry.get().strip()
        if not pdb_id:
            messagebox.showerror("Input Error", "Please provide a PDB ID.")
        else:
            path_to_fasta = pdb_correctness(pdb_id)
            if path_to_fasta is None:
                messagebox.showerror("Input Error", "PDB ID is incorrect.")
            else:
                messagebox.showinfo("Success", "PDB ID processed successfully!")
                show_model_selection(path_to_fasta)
    elif input_type == "UniProt":
        uniprot_id = uniprot_entry.get().strip()
        if not uniprot_id:
            messagebox.showerror("Input Error", "Please provide a UniProt ID.")
        else:
            path_to_fasta = uniprot_correctness(uniprot_id)
            if path_to_fasta is None:
                messagebox.showerror("Input Error", f"UniProt ID ({uniprot_id}) is incorrect.")
            else:
                messagebox.showinfo("Success", "UniProt ID processed successfully!")
                show_model_selection(path_to_fasta)
    else:
        messagebox.showerror("Input Error", "Please select an input type.")


def update_input_fields():
    """
    Update the input fields based on the selected input type to disable/enable the inputting to fields.
    :return:
    """
    input_type = input_type_var.get()
    if input_type == "FASTA":
        fasta_text.config(state=tk.NORMAL)
        pdb_entry.config(state=tk.DISABLED)
        uniprot_entry.config(state=tk.DISABLED)
    elif input_type == "PDB":
        fasta_text.config(state=tk.DISABLED)
        pdb_entry.config(state=tk.NORMAL)
        uniprot_entry.config(state=tk.DISABLED)
    elif input_type == "UniProt":
        fasta_text.config(state=tk.DISABLED)
        pdb_entry.config(state=tk.DISABLED)
        uniprot_entry.config(state=tk.NORMAL)


def show_model_selection(path_to_fasta):
    """
    Show the model selection screen. Run models button is enabled only if at least one model is selected.
    :param path_to_fasta:
    :return:
    """
    for widget in root.winfo_children():
        widget.pack_forget()

    tk.Label(root, text="Select Models to Run:").pack(pady=5)

    model1_checkbutton = tk.Checkbutton(root, text="SwissModel", variable=model1_var)
    model1_checkbutton.pack(anchor=tk.W)

    model2_checkbutton = tk.Checkbutton(root, text="ESMfold", variable=model2_var)
    model2_checkbutton.pack(anchor=tk.W)

    model3_checkbutton = tk.Checkbutton(root, text="OmegaFold", variable=model3_var)
    model3_checkbutton.pack(anchor=tk.W)

    tk.Button(root, text="Run Models", command=lambda: run_models(path_to_fasta)).pack(pady=20)


def run_models(path_to_fasta):
    """
    Run the selected models.
    :param path_to_fasta:
    :return:
    """
    selected_models = []
    if model1_var.get():
        selected_models.append("structure_prediction/swissmodel.py")
    if model2_var.get():
        selected_models.append("structure_prediction/esm.py")
    if model3_var.get():
        selected_models.append("structure_prediction/omega_fold.py")

    if selected_models:
        results = []
        for model in selected_models:
            try:
                result = subprocess.run(["python", model, str(path_to_fasta)], capture_output=True, text=True, check=True)

                result_path = result.stdout.strip().split(" ")[-1]
                if model == "structure_prediction/omega_fold.py":
                    result_path = path_to_fasta.replace(".fasta", "_omegaFold.pdb")
                results.append(result_path)
                print(f"Model {model} executed successfully. Path to pdb: {result_path}")
            except Exception as e:
                messagebox.showerror("Model Execution Error", f"An error occurred while running {model}: {e}")
        messagebox.showinfo("Model Execution", f"Models executed successfully.\nPaths to pdb: {', '.join(results)}")
        ask_iupred_analysis(results, path_to_fasta)
    else:
        messagebox.showerror("Model Selection Error", "Please select at least one model to run.")


def ask_iupred_analysis(results, path_to_fasta):
    """
    Ask the user if they want to run IUPred analysis. I
    :param results: list with pdb file paths
    :param path_to_fasta: path to the fasta file
    :return:
    """
    for widget in root.winfo_children():
        widget.pack_forget()

    tk.Label(root, text="Do you want to run IUPred analysis?").pack(pady=10)
    tk.Button(root, text="Yes", command=lambda: run_iupred(results, path_to_fasta)).pack(side=tk.LEFT, padx=10, pady=20)
    tk.Button(root, text="No", command=lambda: show_final_screen(path_to_fasta, results)).pack(side=tk.RIGHT, padx=10,
                                                                                               pady=20)


def show_final_screen(path_to_fasta, pdb_files, iupred=None):
    """
    Show the final screen with the path to the fasta file.
    :param path_to_fasta:
    :param pdb_files:
    :return:
    """
    for widget in root.winfo_children():
        widget.pack_forget()
    with open(str(path_to_fasta)) as fasta_file:
        sequence = fasta_file.read().strip().split('\n')
        header = sequence[0][1:]
        sequence = "".join(sequence[1:])
    tk.Label(root, text=f"Protein Name: {header}").pack(pady=10)
    tk.Label(root, text=f"Sequence: {sequence}").pack(pady=10)
    tk.Label(root, text=f"Length of sequence: {len(sequence)}").pack(pady=10)
    if iupred is not None:
        tk.Label(root, text="Paths to pymol session with based on disorder colored structures:").pack(pady=10)
        for res in iupred:
            tk.Label(root, text=res).pack(pady=5)
    tk.Label(root, text="PDB files created by models:").pack(pady=10)
    for pdb_file in pdb_files:
        tk.Label(root, text=pdb_file).pack(pady=5)
    tk.Button(root, text="Open PYMOL window with ouput pdbs", command=lambda: open_pymol(pdb_files)).pack(pady=20)
    tk.Button(root, text="Start Over", command=show_main_screen).pack(pady=20)


def open_pymol(pdb_files):
    """
    Open the Pymol windows with the pdb files.
    :param pdb_files:
    :return:
    """
    for pdb_file in pdb_files:
        subprocess.run(["pymol", pdb_file])


def run_iupred(paths_to_pdb, path_to_fasta):
    """
    Run IUPred analysis.
    :param paths_to_pdb:
    :param path_to_fasta:
    :return:
    """
    for widget in root.winfo_children():
        widget.pack_forget()
    tk.Label(root, text="IUPRED running").pack(pady=10)
    iupred_results = []
    for pdb_file in paths_to_pdb:
        try:
            result = subprocess.run(["python", "disorder_analysis/disorder_analysis.py", pdb_file], capture_output=True, text=True, check=True)
            iupred_results.append(result.stdout.strip().split(" ")[-1])
            print(f"IUPred executed successfully. Path to pdb: {result.stdout.strip()}")
        except Exception as e:
            messagebox.showerror("IUPred Execution Error", f"An error occurred while running IUPred: {e}")

    show_final_screen(path_to_fasta, paths_to_pdb, iupred_results)


def show_main_screen():
    """
    Show the main input screen.
    :return:
    """
    for widget in root.winfo_children():
        widget.pack_forget()

    tk.Label(root, text="Select Input Type:").pack(pady=5)
    tk.Radiobutton(root, text="FASTA Sequence", variable=input_type_var, value="FASTA",
                   command=update_input_fields).pack(anchor=tk.W)
    tk.Radiobutton(root, text="PDB ID", variable=input_type_var, value="PDB", command=update_input_fields).pack(
        anchor=tk.W)
    tk.Radiobutton(root, text="UniProt ID", variable=input_type_var, value="UniProt", command=update_input_fields).pack(
        anchor=tk.W)

    tk.Label(root, text="FASTA Sequence:").pack(pady=5)
    fasta_text.pack(padx=10, pady=5, fill=tk.X, expand=True)

    tk.Label(root, text="PDB ID:").pack(pady=5)
    pdb_entry.pack(padx=10, pady=5, fill=tk.X)

    tk.Label(root, text="UniProt ID:").pack(pady=5)
    uniprot_entry.pack(padx=10, pady=5, fill=tk.X)

    process_button.pack(pady=20)

    update_input_fields()


root = tk.Tk()
root.title("Protein Modeling Tool")
root.geometry("800x500")

tk.Label(root, text="3D Structure Prediction and Disorder Analysis Toolkit", font=("Helvetica", 20)).pack(pady=10)
tk.Label(root, text="Version 1.0").pack(pady=5)
tk.Label(root, text="Authors: Michał Rembalski, Jagoda Trzeciak, Bruno Puczko-Szymański, Marta Korpacz").pack(pady=5)
tk.Label(root, text="License: MIT").pack(pady=5)
tk.Label(root, text="Github: https://github.com/jagodatrzeciak/ADP-project").pack(pady=5)
tk.Button(root, text="Start", command=show_main_screen, font=("Helvetica", 12)).pack(pady=20)

input_type_var = tk.StringVar(value="FASTA")

model1_var = tk.BooleanVar()
model2_var = tk.BooleanVar()
model3_var = tk.BooleanVar()

fasta_text = tk.Text(root, wrap=tk.WORD, height=10)
pdb_entry = tk.Entry(root)
uniprot_entry = tk.Entry(root)
process_button = tk.Button(root, text="Process Input", command=process_input)

root.mainloop()
