# 3D Structure Prediction and Disorder Analysis Toolkit

This toolkit provides two primary functionalities: 
- 3D Structure Prediction Using Multiple Methods: TODO
- Disorder Prediction and Visualization: This functionality runs IUPred for a provided UniProt ID or PDB file and colors the structure by disorder using PyMOL.


## 3D Structure Prediction Using Multiple Methods

### Swiss-Model
To use this service, you'll need to register and obtain an API token. Follow the steps below to get started:
- Create account at `https://swissmodel.expasy.org/`
- Go to `https://swissmodel.expasy.org/account`
- Copy your API token 
- Save it as SWISSMODEL_EXPASY_API_TOKEN in your environment

### ESM Metagenomic (ESMFold)
ESMFold is a state-of-the-art method for predicting protein structures from metagenomic data. It leverages deep learning techniques and the ESM model to accurately predict protein structure.

Fortunately, no additional steps are required.

### OmegaFold

To use this model you need follow the instructions below:

```
git clone https://github.com/HeliXonProtein/OmegaFold
cd OmegaFold
python setup.py install
```
This command will download the weight from https://helixon.s3.amazonaws.com/release1.pt to 
```~/.cache/omegafold_ckpt/model.pt``` and load the model.
## Disorder Prediction and Visualization

### Overview
The second component of this toolkit allows users to analyze the intrinsic disorder
of protein structures. It utilizes the IUPred tool to predict disorder regions in protein sequence. The resulting 
predictions are then wizualized by coloring the protein structure in PyMOL.

### Requirements
- biopython
- requests
- pymol

### Example
```bash
# For a UniProt ID
python disorder_analysis/disorder_analysis.py P12345

# For a PDB file
python disorder_analysis/disorder_analysis.py path/to/your.pdb
```