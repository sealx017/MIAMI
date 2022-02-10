MIAMI
================
Souvik Seal and Debashis Ghosh

This is an R package implementing the proposed method from the paper,
“MIAMI: Mutual Information-based Analysis of Multiplex Imaging data”.
The package provides a thorough pipeline for performing marker
co-expression analysis of multiplex imaging data, such as Vectra (Huang
et. al. 2013) and MIBI data (Keren et. al. 2019). The package also
provides standalone functions for computing several mutual information
(MI) theoretic measures (Principe 2010), such as EQMI\*, EQMI and CSQMI
for quantifying the dependence between random variables in general
datasets (definitions are available below in the appendix).

## Loading the package

We install and load the developmental version of MIAMI from GitHub.

``` r
suppressMessages(devtools::install_github('sealx017/MIAMI'))
require(MIAMI)
```

## Loading the example datasets

Next, we import the example files named, “Marker_Data.csv” and
“Clinical_Data.csv”. The first one has expression data of five markers
(![p](http://chart.apis.google.com/chart?cht=tx&chl=p "p") = 5), HLA-DR,
CD45RO, H3K27me3, H3K9ac and HLA_Class_1 for 39 subjects. The second one
has data of two clinical outcomes, recurrence and survival and one
covariate, Age for the same set of subjects. Both of these files are
extracted from the triple-negative breast cancer MIBI data first
published in Keren et. al. 2018. The files have a common column named
“ID” denoting subject IDs.

``` r
data("marker_data")
knitr::kable(head(marker_data), format="markdown")
```

|  ID |    HLA.DR |    CD45RO | HLA_Class_1 |  H3K27me3 |    H3K9ac |
|----:|----------:|----------:|------------:|----------:|----------:|
|   1 | 0.0887767 | 0.0029441 |   0.2295286 | 0.1766168 | 0.0704011 |
|   1 | 0.0118106 | 0.1407535 |   0.2150085 | 0.5186082 | 0.2832617 |
|   1 | 0.0000000 | 0.0000000 |   0.0178665 | 0.2305889 | 0.0206786 |
|   1 | 0.0125603 | 0.0000000 |   0.2317846 | 0.4331141 | 0.2314381 |
|   1 | 0.1843811 | 0.0339638 |   0.1492079 | 0.7069783 | 0.4647193 |
|   1 | 0.0832260 | 0.0000000 |   0.1986478 | 0.7122134 | 0.4879173 |

``` r
data("clinical_data")
knitr::kable(head(clinical_data), format="markdown")
```

|  ID | Recurrence | Recurrence_time | Survival | Survival_time | Age |
|----:|-----------:|----------------:|---------:|--------------:|----:|
|   1 |          0 |               9 |        0 |          2612 |  77 |
|   2 |          1 |             745 |        0 |           745 |  67 |
|   3 |          1 |            3130 |        1 |          3130 |  42 |
|   4 |          0 |              31 |        1 |          2523 |  41 |
|   5 |          1 |            1683 |        1 |          1683 |  64 |
|   6 |          1 |            2275 |        1 |          2275 |  53 |

## Compute EQMI\*, EQMI and CSQMI with an arbitrary data matrix

-   We start by showing how to compute the MI-based measures, EQMI\*,
    EQMI, and CSQMI (Principe 2010) between any arbitrary set of random
    variables (r.v.’s) from a general data-frame i.e., when the data
    does not necessarily come from multiplex imaging platforms. We
    create a matrix named Data_matrix with 2500 samples and 5 columns
    corresponding to five r.v.’s and use it in the function named QMI to
    estimate the measures.

-   A novelty of the estimation algorithm is that the quantities are
    computed in a step-wise fashion. It means that along with the EQMI\*
    between the full set of markers, one can easily extract the EQMI\*
    between several smaller sets of markers. For example, with
    ![p](http://chart.apis.google.com/chart?cht=tx&chl=p "p") = 5
    markers, (1, 2, 3, 4, 5), the EQMI\* between all the following sets
    of markers can be extracted in a single estimation procedure,

    -   (1, 2), denoted by EQMI\*\_12
    -   (1, 2, 3), denoted by EQMI\*\_123
    -   (1, 2, 3, 4), denoted by EQMI\*\_1234
    -   (1, 2, 3, 4, 5), denoted by EQMI\*\_12345.

-   The estimation algorithm requires selecting bandwidth parameters for
    each of the r.v.’s. the default option, bandwidth = “Hpi” uses the
    multivariate plug-in bandwidth matrix described in Wand, M. P., &
    Jones, M. C. (1994). However, in larger datasets (especially, for
    large ![p](http://chart.apis.google.com/chart?cht=tx&chl=p "p")) for
    faster computation, bandwidth = “Ind”, which chooses bandwidth by
    Silverman’s rule (Silverman, B.W. (2018)) for every r.v., can be
    used.

``` r
Data_matrix = marker_data[1:2500, -1]
QMIs = QMI(Data_matrix, bandwidth = "HPI", measure = "All", var_names = T)
print(QMIs)
# $EQMI_star
#               Estimate
# EQMI*_12    0.01650183
# EQMI*_123   0.19007064
# EQMI*_1234  0.36213356
# EQMI*_12345 0.56960589
# 
# $EQMI
#               Estimate
# EQMI_12       6.200363
# EQMI_123    283.226354
# EQMI_1234  1574.540741
# EQMI_12345 8456.344853
# 
# $CSQMI
#               Estimate
# CSQMI_12    0.01822507
# CSQMI_123   0.27834795
# CSQMI_1234  0.55462570
# CSQMI_12345 0.96711185
# 
# $Bandwidth_parameters
#             h1           h2           h3           h4           h5
# 1 0.0001129059 3.739611e-05 0.0003045994 0.0005547427 0.0005718942
```

## Computing subject specific EQMI\* values with the imaging data

We return to the analysis of multiplex imaging data. We compute the
EQMI\* of all the subjects using the function QMI_all and store the
values in a matrix whose every row corresponds to a subject.

``` r
EQMI_vector = QMI_all(marker_data, bandwidth = "Ind", measure = "EQMI_star", progress_bar = "False")
knitr::kable(head(EQMI_vector), format="markdown")
```

|  ID | EQMI\*\_12 | EQMI\*\_123 | EQMI\*\_1234 | EQMI\*\_12345 |
|----:|-----------:|------------:|-------------:|--------------:|
|   1 |  0.0122462 |   0.1705994 |    0.3427808 |     0.5413480 |
|   2 |  0.0183228 |   0.1212016 |    0.3560708 |     0.4309078 |
|   3 |  0.0648835 |   0.3260319 |    0.3505626 |     0.5423636 |
|   4 |  0.0076886 |   0.1274159 |    0.2140916 |     0.4323397 |
|   5 |  0.0735815 |   0.3347920 |    0.3927648 |     0.5451000 |
|   6 |  0.0327404 |   0.1719063 |    0.2413455 |     0.5027947 |

## Association testing using CoxPH model

### Association testing with survival outcome

Using the vector of estimated EQMI\* of all the subjects, we perform
association analysis with two clinical outcomes, survival and
recurrence. We have the outcomes, the time to death and the time to
recurrence and the respective censoring indicators (= 0 for an event) as
the columns of the matrix named clinical_data. We have one single
covariate, Age in the same matrix. Below, we create a matrix named
surv_dat with all the subject IDs and their survival outcomes, and
another matrix named covariates with the subject IDs and available
covariates which, in this case, is just Age. We then use the function
Cox_PH to fit the proportional hazard (PH) model outputting a table of
p-values. To add higher order terms in the PH model, change degree to \>
1.

``` r
surv_dat = clinical_data[,c(1,4:5)]
covariates = clinical_data[,c(1,6)]
SurvCox = Cox_PH(surv_dat, covariates, EQMI_vector, degree = 1)
knitr::kable(SurvCox, format="markdown")
```

| Variable      |   p-value |
|:--------------|----------:|
| EQMI\*\_12    | 0.0685692 |
| EQMI\*\_123   | 0.0373026 |
| EQMI\*\_1234  | 0.0467573 |
| EQMI\*\_12345 | 0.0150054 |

### Association testing with recurrence outcome

Below, we create a matrix named recur_dat with all the subject IDs and
their recurrence outcomes, and use the same covariates matrix as
earlier. We again use the Cox_PH function and get the p-values.

``` r
recur_dat = clinical_data[,c(1:3)]
RecurCox = Cox_PH(recur_dat, covariates, EQMI_vector, degree = 1)
knitr::kable(RecurCox, format="markdown")
```

| Variable      |   p-value |
|:--------------|----------:|
| EQMI\*\_12    | 0.0609279 |
| EQMI\*\_123   | 0.0268636 |
| EQMI\*\_1234  | 0.0119667 |
| EQMI\*\_12345 | 0.0128277 |

## Additional tools: Compute and plot univariate and bivariate kernel densities of the markers

#### Univariate marginal density

We provide a few basic functions to estimate and plot the univariate
marginal densities of the random variables. We go back to the matrix
named Data_matrix which has 2500 samples and 5 columns corresponding to
five r.v.’s. We look at first two columns (call them, X_1 and X_2) and
estimate their kernel density estimates as f_1 and f_2 using the
function dens_univ. Next, using the function univ_dens_plot, we plot f_1
and f_2. The default number of grids is, ngrids = 1024 and the bandwidth
parameter is selected using Silverman’s rule.

``` r
X_1 = Data_matrix[,1]
X_2 = Data_matrix[,2]  

f_1 = univ_dens(X_1, ngrids = 1024)
f_2 = univ_dens(X_2, ngrids = 1024)

p1 = univ_dens_plot(f_1)
p2 = univ_dens_plot(f_2)
```

<img src="README_files/figure-gfm/Computing univ density-1.png" width="50%" /><img src="README_files/figure-gfm/Computing univ density-2.png" width="50%" />

#### Bivariate joint density

We estimate the bivariate joint density of X_1 and X_2, as f_12 using
the function biv_dens. Next, using the function biv_dens_plot, we plot
f_12. The default number of grids is, ngrids = 512 and the bandwidth
matrix used is the multivariate plug-in bandwidth matrix described in
Wand, M. P., & Jones, M. C. (1994). It may be easy to interpret the
estimated density focusing only on the smaller values of
![X_1](http://chart.apis.google.com/chart?cht=tx&chl=X_1 "X_1") and
![X_2](http://chart.apis.google.com/chart?cht=tx&chl=X_2 "X_2") since
both the r.v.’s usually have most of the values close to 0. Hence, we
use two thresholds (q1, q2) to display the estimated density only
between 70% quantiles of both the r.v.’s.

``` r

f_12 = biv_dens(cbind(X_1, X_2), ngrids = 512)
q1 = quantile(X_1, 0.7); q2 = quantile(X_2, 0.7)
p = biv_dens_plot(f_12, maxs = c(q1, q2))
```

<img src="README_files/figure-gfm/Computing biv density-1.png" width="100%" />

## References

a\) Huang, W., Hennrick, K., & Drew, S. (2013). A colorful future of
quantitative pathology: validation of Vectra technology using
chromogenic multiplexed immunohistochemistry and prostate tissue
microarrays. Human pathology, 44(1), 29-38.

b\) Keren, L., Bosse, M., Thompson, S., Risom, T., Vijayaragavan, K.,
McCaffrey, E., … & Angelo, M. (2019). MIBI-TOF: A multiplexed imaging
platform relates cellular phenotypes and tissue structure. Science
advances, 5(10), eaax5851.

c\) Principe, J. C. (2010). Information theoretic learning: Renyi’s
entropy and kernel perspectives. Springer Science & Business Media.

d\) Keren, L., Bosse, M., Marquez, D., Angoshtari, R., Jain, S., Varma,
S., … & Angelo, M. (2018). A structured tumor-immune microenvironment in
triple negative breast cancer revealed by multiplexed ion beam imaging.
Cell, 174(6), 1373-1387.

e\) Wand, M. P., & Jones, M. C. (1994). Multivariate plug-in bandwidth
selection. Computational Statistics, 9(2), 97-116.

f\) Silverman, B. W. (2018). Density estimation for statistics and data
analysis. Routledge.

## Appendix

### Definition of the measures

Here, we show the mathematical definitions of the three different
MI-based measures, EQMI\*, EQMI and CSQMI, implemented in this package.

![alt text here](Formulas.png)

![f\_{12 \\ldots p}(x\_{1}, x\_{2}, \\ldots , x\_{p})](http://chart.apis.google.com/chart?cht=tx&chl=f_%7B12%20%5Cldots%20p%7D%28x_%7B1%7D%2C%20x_%7B2%7D%2C%20%5Cldots%20%2C%20x_%7Bp%7D%29 "f_{12 \ldots p}(x_{1}, x_{2}, \ldots , x_{p})")
is the joint PDF and
![f\_{1}(x\_{1}), f\_{2}(x\_{2}), \\ldots f\_{p}(x\_{p})](http://chart.apis.google.com/chart?cht=tx&chl=f_%7B1%7D%28x_%7B1%7D%29%2C%20f_%7B2%7D%28x_%7B2%7D%29%2C%20%5Cldots%20f_%7Bp%7D%28x_%7Bp%7D%29 "f_{1}(x_{1}), f_{2}(x_{2}), \ldots f_{p}(x_{p})")
are the marginal PDFs of the
![p](http://chart.apis.google.com/chart?cht=tx&chl=p "p") variables. The
details of the efficient algorithm used for estimating
![V_J, V_C](http://chart.apis.google.com/chart?cht=tx&chl=V_J%2C%20V_C "V_J, V_C")
and ![V_M](http://chart.apis.google.com/chart?cht=tx&chl=V_M "V_M") can
be found in the manuscript, “MIAMI: Mutual Information-based Analysis of
Multiplex Imaging data”.
