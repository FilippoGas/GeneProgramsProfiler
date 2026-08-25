# Contributing to GeneProgramsProfiler

Thanks for considering a contribution! This document explains how to report issues, suggest features, and submit pull requests.

## Reporting bugs

Open a [GitHub Issue](https://github.com/FilippoGas/GeneProgramsProfiler/issues/new) with:

- A clear description of the problem
- Steps to reproduce
- Your `config.yaml` (or the relevant section)
- Snakemake version (`snakemake --version`)
- OS and Python version
- Full error output (use `--show-failed-logs`)

## Suggesting features

Open a GitHub Issue with the `enhancement` label. Describe the use case, not just the implementation. If the feature touches a specific module, mention which rule file is involved.

## Development setup

```bash
git clone git@github.com:FilippoGas/GeneProgramsProfiler.git
cd GeneProgramsProfiler

# Create conda environments for the modules you plan to modify
conda env create -f workflow/envs/spectra.yaml
conda env create -f workflow/envs/cNMF.yaml
# ... or let snakemake manage environments with --sdm conda

# Verify the workflow resolves
snakemake --snakefile workflow/Snakefile --dry-run
```

## Code style

| File type | Formatter | Command |
|---|---|---|
| Snakefiles (`.smk`, `Snakefile`) | [snakefmt](https://github.com/snakemake/snakefmt) | `snakefmt workflow/Snakefile workflow/rules/*.smk` |
| YAML (configs, envs, schemas) | [prettier](https://prettier.io/) | `prettier --write config/ workflow/envs/ workflow/schemas/` |
| R scripts | [styler](https://styler.r-lib.org/) (tidyverse style) | `styler::style_file("path/to/script.R")` |
| Python scripts | [ruff](https://docs.astral.sh/ruff/) | `ruff check path/to/script.py` |

CI runs `snakefmt`, `prettier`, and `snakemake --lint` automatically. PRs that fail formatting will be blocked.

## Adding or modifying rules

1. Create or edit a `.smk` file in `workflow/rules/`
2. Include it in `workflow/Snakefile` with `include: "rules/your_rule.smk"`
3. Add a new conda environment in `workflow/envs/` if needed
4. Update `workflow/schemas/config.schema.yaml` with any new config keys
5. Add corresponding defaults in `config/config.yaml` and `.test/config/config.yaml`
6. Update `config/README.md` with documentation for new parameters
7. Verify with `snakemake --snakefile workflow/Snakefile --dry-run`

## Pull requests

1. Create a feature branch from `main`
2. Make your changes following the conventions above
3. PR titles **must** follow [Conventional Commits](https://www.conventionalcommits.org/):
   - `feat(spectra): add configurable HVG selection`
   - `fix(cNMF): handle empty consensus matrix`
   - `docs: update README test instructions`
   - `refactor: standardize output path layout`
4. CI will run formatting, linting, and a full test run before merge

## Testing

CI runs on a self-hosted runner and executes:

1. **Formatting** — snakefmt + prettier check
2. **Linting** — `snakemake --lint`
3. **Test run** — full workflow on the built-in Natri et al. IPF-vs-Control dataset

To test locally:

```bash
snakemake --snakefile workflow/Snakefile --sdm conda --cores 2
```

The test dataset lives in `.test/data/` and is committed to the repo.

## Project structure

```
workflow/
├── Snakefile              # Entrypoint (do not invoke .smk files directly)
├── rules/                 # Rule modules (dependency order)
│   ├── preprocess.smk
│   ├── spectra.smk
│   ├── DE_analysis.smk
│   ├── functional_enrichment.smk
│   ├── cNMF.smk
│   └── collect_results.smk
├── scripts/               # R and Python scripts called by rules
├── envs/                  # Conda environment definitions
└── schemas/               # Config validation schemas
```
