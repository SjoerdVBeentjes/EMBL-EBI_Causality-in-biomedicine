# EMBL-EBI: Causality in Biomedicine

Materials for the Day 1 practical session 'Fundamental problem of causal inference' of the EMBO/EMBL-EBI practical course *"Causality in biomedicine: going beyond associations"* (4–9 October 2026, EMBL-EBI, Hinxton).

## R Environment (renv) {#r-environment-renv}

Package versions for this session are pinned with [`renv`](https://rstudio.github.io/renv/), so everyone - trainer and trainees, runs the exact same set of packages.

### For trainees

From the top level of this repository:

``` bash
Rscript restore_renv.R
```

This installs the exact package versions recorded in [`renv.lock`](renv.lock) into a project-local library — it won't touch or conflict with anything else on your machine. If you're working in RStudio with this repo open as a Project, `renv` activates automatically when the project opens; you can run `renv::restore()` at the console instead.

### Packages included

See [`dependencies.R`](dependencies.R) for a list of packages included.
