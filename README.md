# Order-based Structure Learning without Score Equivalence

This repository provides an implementation of the paper and demonstrates the reproducibility for the three experiment types in the paper.

# Quick start

## 1. experiment_type == "among" (Table 1, Table 2 and Table 3)

The experiment is to compare the performance of our method ("ORDER") with that of other competing methods ("TD", "LISTEN", "GDS"). The output is log files of performance table. An example of the bash code is below.

```bash
Rscript run.R -x "among" -m "ORDER" -t 'unif' -N 30 -n 1000 -p 40 -s 1 -e 0.1 -l 0.5
```

### description of passing variables.

Rscript run.R

-x "among" : it indicates that the experiment is among competing models.

-m c("ORDER","TD", "LISTEN","GDS") : choosing among 4 models.

-t c("unif", "mis", "gauss", "same") : choosing among 4 options for signals.

-N 30  : Choosing the number of core for parallel computings. (int)

-n 1000  : Choosing the number of samples (int)

-p 40  : Choosing the number of nodes (int)

-s 1  : Choosing the value of error variance (numeric)

-e 0.1  : Choosing the value of edge probability (numeric in (0, 1))

-l 0.3  : Choosing the value of lower bound of signal size (numeric)

-i 0.3  : Choosing the value of degree of misspecification. (numeric in (0, 1); only when t == "mis")

## 2. experiment_type == "mixing" (Fig. 1, Fig. 6 and Fig. 7)

The experiment is to draw a mixing diagram to show chains with three neighborhoods are reasonably mixing well in moderate number of nodes. The output is three mixing plots. An example of the bash code is below.

```bash
Rscript run.R -x "mixing" -m "ORDER"
```

### description of passing variables.

Rscript run.R

-x "mixing" : it indicates that the experiment is for the mixing behavior chapter.

-m c("ORDER", "MINIMAP") : choosing between 2 models.

## 3. experiment_type == "misspecify" (Fig. 2.)

The experiment is to draw boxplots of hamming distance according to the degree of misspecification for three models; "ORDER","TD", "LISTEN". The output is a set boxplots of hamming distance according to the degree of misspecification. An example of the bash code is below.

```bash
Rscript run.R -x "misspecify"
```
