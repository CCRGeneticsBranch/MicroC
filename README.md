Runs the Micro-C pipeline from https://micro-c.readthedocs.io/en/latest/index.html
1. Read alignment
2. Generate pairs file
3. Library QC and complexity
4. Generate Hi-C contact map

Additionally, identifies compartments (FANC), TADs (HiCExplorer and Arrowhead), loops (Mustache), and structural variants (EagleC).

To run the pipeline, clone the repo from github. In an interactive session, initialize the analysis directory with:
```
  bash [path to repo]/MicroC/run_microc --runmode=init --workdir=[directory] --depth=[shallow/deep]
```
Edit the config/samples.tsv and config/config.yaml files (genome version). Perform a dryrun with:
```
  bash [path to repo]/MicroC/run_microc --runmode=dryrun --workdir=[directory]
```
If successful, launch the pipeline with:
```
  bash [path to repo]/MicroC/run_microc --runmode=run --workdir=[directory]
```

Thank you to Hsien-Chao Chou for sharing library merging and topology scripts
