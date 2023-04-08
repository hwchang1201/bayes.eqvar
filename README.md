# Order-based Structure Learning without Score Equivalence

This repository provides an implementation of the paper and demonstrates the reproducibility for the seven experiment types in the paper.

# Quick start

## 1. experiment_type == "among" (Table 1, Table 4 and Table 5)

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

## 4. experiment_type == "highdim" (Table 2.)

The experiment is to show high dimensional experiments.

```bash
Rscript run.R -x "highdim" -k 7
```
### description of passing variables.

Rscript run.R

-x "highdim" : it indicates that the experiment is for high dimensional experiments.

-k c(1,2,3,4,5,6,7) : choosing k = 1, 2, 3, 4, 5, 6, or 7.


## 5. experiment_type == "bias" (Table 3.)

The experiment is to measure bias of the assumption.

```bash
Rscript run.R -x "bias" -i 0.5
```
### description of passing variables.

Rscript run.R

-x "bias" : it indicates that the experiment is for measuring bias of the assumption.

-i 0:9/10 : choosing eps = 0, 0.1, ... , or 0.9.


## 6. experiment_type == "td" 

The experiment is to compare performance between TD and ITD.

```bash
Rscript run.R -x "td" 
```

## 7. experiment_type == "increasing"

The experiment is to show weakly increasing error variance cases.

```bash
Rscript run.R -x "increasing" -i 0.5
```
### description of passing variables.

Rscript run.R

-x "increasing" : it indicates that the experiment is for weakly increasing error variance cases.

-i 0:9/10 : choosing eps = 0, 0.1, ... , or 0.9.
