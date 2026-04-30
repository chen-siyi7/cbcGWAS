# `data/`

This directory is intentionally empty in the public repository. The real-data
pipelines (`R/realdata/ad_ldsc.R` and `R/realdata/bmi_t2d_mr.R`) need GWAS
summary statistics and LD reference panels that are not redistributed here.

## What to put in this directory

After cloning the repo, populate this directory with the files below. The
pipeline scripts have a `CONFIG` block at the top that points at these
locations; you can either match this layout or edit `CONFIG` to match yours.

```
data/
├── sumstats/
│   ├── bellenguez2022_ad.tsv.gz          # Alzheimer's disease GWAS
│   ├── yengo2018_bmi.txt.gz              # BMI (continuous exposure)
│   └── mahajan2018_t2d_bmi_adjusted.txt.gz  # BMI-adjusted T2D outcome
├── ldsc_data/
│   ├── eur_w_ld_chr/                     # LDSC EUR LD scores (~5 GB unpacked)
│   │   ├── 1.l2.ldscore.gz
│   │   ├── 1.l2.M
│   │   ├── 1.l2.M_5_50
│   │   ├── ...
│   │   └── 22.l2.ldscore.gz
│   ├── w_hm3.snplist                     # HapMap3 munging reference
│   └── 1kg/
│       ├── EUR.bed                       # 1000G EUR PLINK fileset
│       ├── EUR.bim                       # used by ieugwasr::ld_clump_local
│       └── EUR.fam
└── tools/
    └── plink                             # PLINK 1.9 binary (or use system one)
```

## Where to download

| File / directory                | Source                                                                                  |
| ------------------------------- | --------------------------------------------------------------------------------------- |
| `bellenguez2022_ad.tsv.gz`      | GWAS Catalog accession **GCST90027158** (Bellenguez et al. 2022 Nat Genet 54:412-436)   |
| `yengo2018_bmi.txt.gz`          | GIANT consortium downloads (Yengo et al. 2018 Hum Mol Genet 27:3641-3649)               |
| `mahajan2018_t2d_bmi_adjusted.txt.gz` | DIAGRAM consortium downloads (Mahajan et al. 2018 Nat Genet 50:1505-1513)         |
| `eur_w_ld_chr/` and `w_hm3.snplist` | LDSC reference (Bulik-Sullivan et al. 2015). Zenodo mirror: <https://zenodo.org/records/8182036> |
| `1kg/EUR.{bed,bim,fam}`         | <http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz> (the MR-Base group's 1000G LD reference)  |
| `plink` binary                  | <https://www.cog-genomics.org/plink/1.9/>, or auto-installed by `genetics.binaRies`      |

## What this repository will not redistribute

GWAS summary statistics and reference panels are large, are licensed by
their original authors, and change under the original maintainers'
URLs over time. Distributing copies here would create a stale fork and
license tangle, so the pipelines instead document exactly which file
versions to download and check column names defensively at runtime
(see `R/realdata/helpers.R::sanity_check_sumstats`).

## Verifying your setup before committing to a long run

Both pipeline scripts have a `CONFIG$phase = "PREFLIGHT"` mode that
inspects every path and prints column-name diagnostics without committing
to LDSC munging or MR clumping. Always run preflight first:

```bash
# from repo root, after editing CONFIG paths
Rscript R/realdata/ad_ldsc.R    # CONFIG$phase = "PREFLIGHT" by default
Rscript R/realdata/bmi_t2d_mr.R
```

When preflight prints "all checks passed", switch `CONFIG$phase` to `"RUN"`
and execute the same command again for the full pipeline.
