# bayes.eqvar

We have 5 experiment types in this code. Each type of experiment would give us replication of the results (Figures and Tables) in the paper.

## 1. experiment_type == "among" (Table 1 and Table 2.)

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

## 2. experiment_type == "mixing" (Figure 1.)

The experiment is to draw a mixing diagram to show chains with three neighborhoods are reasonably mixing well in moderate number of nodes. The output is three mixing plots. An example of the bash code is below.

```bash
Rscript run.R -x "mixing"
```

## 3. experiment_type == "misspecify" (Figure 2.)

The experiment is to draw boxplots of hamming distance according to the degree of misspecification for three models; "ORDER","TD", "LISTEN". The output is a set boxplots of hamming distance according to the degree of misspecification. An example of the bash code is below.

```bash
Rscript run.R -x "misspecify"
```

## 3. experiment_type == "within" (Figure 4.)

The experiment is to compare the performance within our method ("ORDER"), within ("RB", "unweighted","ITD"). The output is a set log files of performance table and plots. An example of the bash code is below.

```bash
Rscript run.R -x "within"
```
