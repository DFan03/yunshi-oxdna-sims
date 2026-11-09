# Scripts index

Detailed, workflow-ordered guides (download and open in a browser; self-contained, phone-friendly):
[RNA.html](RNA.html) and [DNA.html](DNA.html). Version policy: [VERSIONS.md](VERSIONS.md). Cluster paths: [PORTABILITY.md](PORTABILITY.md).

```mermaid
flowchart LR
  A[setup_run.py] --> B[MC pre-relax<br/>oxDNA CPU] --> C[MD relax<br/>oxDNA CUDA] --> D[Production MD]
  D --> E[Analysis<br/>oat, hb_cal, bend/twist]
  D --> F[Umbrella + WHAM]
  D --> G[FFS - RNA only]
  D --> H[Thermal ramp]
  E --> I[notebooks]
  F --> I
```

| Stage | RNA | DNA |
|---|---|---|
| Setup | `rna/setup/` | `dna/setup/` |
| Relax / production | `templates/`, `jobs/` | `templates/`, `jobs/` |
| Analysis | `rna/analysis/`, `rna/multimer/` | `dna/analysis/`, `dna/multimer/` |
| Umbrella + WHAM | `rna/umbrella/` | `dna/umbrella/` |
| FFS | `rna/ffs/` | none |
| Thermal ramp | `rna/thermal_ramp/` | `dna/thermal_ramp/` |
| Plots | `notebooks/rna/` | `notebooks/dna/` |
