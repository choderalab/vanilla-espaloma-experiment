# Protein-ligand vanilla MD using espaloma-0.3
This repository includes scripts to run vanilla MD using `espaloma-0.3`.


## Description
The robustness and stability of `espaloma-0.3` is tested by running simple vanilla MD simulations of a protein-ligand complex system (Tyk2 protein).
The RMSD profiles of heavy ligand atoms and C-alpha protein atoms are monitored. 
`espaloma-0.3` is used to parametrize both the protein and ligand, and `openff-2.0.0` is used to assign the LJ parameters. 
As a control experiment, the ligand and protein is parametrized with `openff-2.1.0` and Amber `ff14SB`, respectively.
The initial Tyk2 protein and ligand structures are taken from the custom alchemical protein-ligand binding benchmark dataset (https://github.com/kntkb/protein-ligand-benchmark-custom).


## Manifest
- `experiment/`: Stores scripts and directories to run and analyze vanilla MD
    - `script/`
    - `tyk2-lig_ejm_31/`
        - `crd/`
        - `espaloma/`
        - `openff-2.1.0/`
- `envs/`: Stores conda environment files
    - `environment-0.3.0-v3.yaml`: Conda environment to run Perses with `espaloma-0.3` that parameterize both small molecules and proteins


## Citation
If you find this helpful please cite the following:

```
@misc{takaba2023espaloma030,
      title={Espaloma-0.3.0: Machine-learned molecular mechanics force field for the simulation of protein-ligand systems and beyond}, 
      author={Kenichiro Takaba and Iván Pulido and Mike Henry and Hugo MacDermott-Opeskin and John D. Chodera and Yuanqing Wang},
      year={2023},
      eprint={2307.07085},
      archivePrefix={arXiv},
      primaryClass={physics.chem-ph}
}
```