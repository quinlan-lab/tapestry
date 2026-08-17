# PED example

`family.ped` is a minimal six-column PED example. Tapestry reads the PacBio WDL
miniwdl `outputs.json` directly; users do not create a run configuration or
manifest.

The public miniwdl entry point is:

```bash
nextflow run . -profile docker \
  --outputs-json /path/to/miniwdl-run/outputs.json \
  --ped /path/to/family.ped \
  --reference-fasta /path/to/GRCh38.fa \
  --outdir /path/to/tapestry-results \
  -resume
```

Validation always runs first and publishes `resolved-run.json` and
`resolved-manifest.json` under `<outdir>/pipeline_info/`. Add
`-entry validate` to stop after that preflight.
