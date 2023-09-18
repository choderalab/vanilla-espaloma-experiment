# Protein-ligand vanilla MD using espaloma-0.3
This repository includes scripts to run vanilla MD using `espaloma-0.3`.

## Description


## Manifest
- `experiment/`: Stores directories and scripts to run vanilla MD
    - `script/`
    - `tyk2-lig_ejm_31/`
- `envs/`: Stores conda environment files
    - `environment-0.3.0-v3.yaml`: Conda environment to run Perses with `espaloma-0.3.0` that parameterize both small molecules and proteins


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