# Generic run-config example

`family.yaml` and `pacbio-wdl.tapestry.json` demonstrate the schema-v1 public
contracts. Paths are resolved relative to the document containing them.

The `data/` artifacts are intentionally not included: replace those paths with
the localized outputs from a PacBio HiFi human WGS WDL v3.3.0 or v3.3.1 family
run. The canonical manifest is Tapestry-owned JSON, not raw Cromwell, miniwdl,
or Terra metadata.

Run a populated copy with:

```bash
nextflow run . -profile docker \
  --run-config /path/to/family.yaml \
  -resume
```

Validation is always the first stage. To stop after the optional preflight, add
`-entry validate`. The compact YAML contains only user inputs and scientific
choices; the resolved run records the fixed pedigree/GRCh38/`gtg`/model
constraints.
