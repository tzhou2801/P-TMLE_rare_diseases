# TMLE+SL for covariate adjustment in small rare-disease trials: simulation code

R code for the simulation study (Section 6) in

> Zhou T, Gruber S, Lee H, Lee W, Nie L, van der Laan M (2026). Improving efficiency and power in small rare-disease trials through prespecified, data-adaptive covariate adjustment. arXiv:2611.15071. https://arxiv.org/abs/2611.15071

| Script | Data-generating process | Results in the paper |
|---|---|---|
| [`sim_dgp1.R`](sim_dgp1.R) | DGP 1, eq. (19) | Table 2, Figure 1 |
| [`sim_dgp2.R`](sim_dgp2.R) | DGP 2, eq. (20) | Table 3 |

Each script is self-contained. It simulates 1,000 trials with 50 participants each and applies the five estimators of Table 1 (unadjusted, ANCOVA(W4, W5), ANCOVA(W3, W4, W5), TMLE+SL with the known treatment mechanism and TMLE+SL with an estimated treatment mechanism) in the effect setting (tau = -1.3) and in the null setting (tau = 0).

## Requirements

R and the following CRAN packages:

```r
install.packages(c("glmnet", "future", "future.apply", "ggplot2"))
```

`ggplot2` is needed only for Figure 1 (`sim_dgp1.R`), and `parallel` ships with R. The results in the paper were produced with R 4.2.1 (x86_64, macOS), glmnet 4.1-4, future 1.33.0, future.apply 1.11.0 and ggplot2 3.4.2. Other versions of R or glmnet, or another platform, can change individual replications slightly.

## Running the simulations

```sh
Rscript sim_dgp1.R
Rscript sim_dgp2.R
```

Each script writes to the working directory:

- a CSV file with one row per replication, setting and estimator;
- a CSV file with the rows of Table 2 or Table 3, which are also printed;
- for DGP 1 only, Figure 1 as a PNG file.

The replications run in parallel on all but one of the available cores; set `n_workers` at the top of a script to change this. A full run of one script (1,000 replications in both settings) takes about 8 hours with 6 workers, so run the two scripts one after the other (about 16 hours in all).

Replication i always uses the same random-number stream, derived from the seed 20260622. The results therefore do not depend on the number of workers, and a run of a subset of the replications (for example `reps <- 1:10`) reproduces the corresponding rows of the full run.

## Citation

```bibtex
@misc{zhou2026tmlesl,
  author        = {Zhou, Tianyue and Gruber, Susan and Lee, Hana and Lee, Wonyul and Nie, Lei and van der Laan, Mark},
  title         = {Improving Efficiency and Power in Small Rare-Disease Trials Through Prespecified, Data-Adaptive Covariate Adjustment},
  year          = {2026},
  eprint        = {2611.15071},
  archivePrefix = {arXiv},
  primaryClass  = {stat.ME}
}
```

## License

MIT; see [LICENSE](LICENSE).
