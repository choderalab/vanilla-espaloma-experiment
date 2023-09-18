# Basic Usage
- System preparation using `espaloma`.
    > cd espaloma/prep  
    > bsub < lsf-submit.sh

    Note that `espaloma-0.3.0.pt` and `openff-2.1.0.offxml` are store in a local directory (`/home/[username]/.espaloma/espaloma-0.3.0.pt` and `/home/[username]/.offxml/openff-2.1.0.offxml`), which can be found [here (espaloma)](https://github.com/choderalab/espaloma/releases/tag/0.3.1) and [here (openff)](https://github.com/openforcefield/openff-forcefields), respectively.

- Run vanilla MD
    > cd espaloma/run/sim1  
    > ./submit.sh

- Analyze RMSD
    > cd espaloma/run/analyze  
    > bsub < lsf-submit.sh
