MIAMI
================
Souvik Seal

This is an R package implementing the proposed method from the paper,
“MIAMI: Mutual Information-based Analysis of Multiplex Imaging data.”
The package provides a thorough pipeline for performing marker
co-expression analysis of multiplex imaging datasets, coming from Vectra
or MIBI platforms. The package also provides functions for computing
several mutual information theoretic measures, such as EQMI and CSQMI
for measuring dependence between random variables (Principe 2010) in
general datasets.

## Loading required packages

One can install the developmental version of MIAMI by running the
command: **devtools::install_github(‘sealx017/MIAMI’)**.

``` r
devtools::install_github('sealx017/MIAMI')
require(MIAMI)
#suppressMessages(source("/Users/seals/Documents/GitHub/MIAMI/R/density_and_EQMI.R"))
#require(ks)
#require(survival)
#require(survminer)
#library(formattable)
```

## Loading the dataset

Next, we import the example datasets named, “Marker_Data.csv” and
“Clinical_Data.csv.” The first one is the expression dataset of five
markers, HLA-DR, CD45RO, H3K27me3, H3K9ac and HLA_Class_1, and 39
subjects. The second one is the clinical dataset with recurrence and
survival outcomes. Both of these files are extracted from the
triple-negative breast cancer MIBI data first published in Keren et.
al. 2018. The files have a common column named “ID” denoting individual
subject IDs.

``` r
marker_data = read.csv("Data/Marker_Data.csv")
clinical_data = read.csv("Data/Clinical_Data.csv")
```

``` r
head(marker_data)
#   ID     HLA.DR      CD45RO HLA_Class_1  H3K27me3     H3K9ac
# 1  1 0.08877669 0.002944144  0.22952856 0.1766168 0.07040106
# 2  1 0.01181060 0.140753465  0.21500855 0.5186082 0.28326167
# 3  1 0.00000000 0.000000000  0.01786649 0.2305889 0.02067860
# 4  1 0.01256025 0.000000000  0.23178456 0.4331141 0.23143812
# 5  1 0.18438106 0.033963838  0.14920787 0.7069783 0.46471933
# 6  1 0.08322602 0.000000000  0.19864780 0.7122134 0.48791729
```

``` r
head(clinical_data)
#   ID Recurrence Recurrence_time Survival Survival_time Age
# 1  1          0               9        0          2612  77
# 2  2          1             745        0           745  67
# 3  3          1            3130        1          3130  42
# 4  4          0              31        1          2523  41
# 5  5          1            1683        1          1683  64
# 6  6          1            2275        1          2275  53
```

## Compute EQMI\*, EQMI and CSQMI with an arbitrary data matrix

We start by showing how to compute the measures, EQMI\*, EQMI, and CSQMI
(Principe 2010) between arbitrary random variables (r.v.’s) where the
data does not necessarily have a multiplex imaging data structure. We
create a matrix named Data_matrix with 5000 samples and 5 columns
corresponding to five r.v.’s.

``` r
Data_matrix = marker_data[1:5000,-1]
QMIs = QMI(Data_matrix, bandwidth = "Ind", measure = "All", var_names = T)
print(QMIs)
# $EQMI_star
#              Estimate
# EQMI*_12    0.0127644
# EQMI*_123   0.1749177
# EQMI*_1234  0.3478008
# EQMI*_12345 0.5455553
# 
# $EQMI
#               Estimate
# EQMI_12       1.795842
# EQMI_123     88.907543
# EQMI_1234   497.865526
# EQMI_12345 2479.777068
# 
# $CSQMI
#               Estimate
# CSQMI_12    0.01859524
# CSQMI_123   0.28237864
# CSQMI_1234  0.58020008
# CSQMI_12345 0.98040785
# 
# $Bandwidth_parameters
#             h1           h2           h3           h4           h5
# 1 0.0003397939 0.0001460787 0.0004341861 0.0006256999 0.0006317829
```

## Computing subject specific EQMI values with imaging data

We compute the measure, EQMI\* for every subject. A novelty of the
computation algorithm is that it is computed in a step-wise fashion. It
means that along with the EQMI\* between the full set of markers, one
can easily extract the EQMI\* between several smaller sets of markers.
For example, while computing the EQMI\* between five markers, (1, 2, 3,
4, 5), the EQMI\* between all the following sets of markers can also be
extracted,

a\) (1, 2), denoted by EQMI_12

b\) (1, 2, 3), denoted by EQMI_123

c\) (1, 2, 3, 4), denoted by EQMI_1234

d\) (1, 2, 3, 4, 5), denoted by EQMI_12345.

``` r
EQMI_vector = QMI_all(marker_data[,1:3], bandwidth = "Ind", measure = "EQMI_star", progress_bar = "False")
head(formattable::formattable(EQMI_vector))
```

<table class="table table-condensed">
<thead>
<tr>
<th style="text-align:right;">
ID
</th>
<th style="text-align:right;">
EQMI\*\_12
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.012246224
</td>
</tr>
<tr>
<td style="text-align:right;">
2
</td>
<td style="text-align:right;">
0.018322808
</td>
</tr>
<tr>
<td style="text-align:right;">
3
</td>
<td style="text-align:right;">
0.064883525
</td>
</tr>
<tr>
<td style="text-align:right;">
4
</td>
<td style="text-align:right;">
0.007688614
</td>
</tr>
<tr>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
0.073581489
</td>
</tr>
<tr>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
0.032740427
</td>
</tr>
</tbody>
</table>

## Association testing using CoxPH model

Using the vector of estimated EQMI\* of all the subjects, we perform
association analysis with two clinical outcomes, survival and
recurrence. We have both the survival and recurrence times and
respective censoring indicators ( = 0 for an event) as columns of the
matrix named clinical_data. We have one single covariate “Age” in the
same matrix. Below, we create a matrix named surv_dat with all the
subject IDs and their survival outcomes, and another matrix named
covariates with the subject IDs and available covariates which, in this
case, is just “Age.”

``` r
surv_dat = clinical_data[,c(1,4:5)]
covariates = clinical_data[,c(1,6)]
SurvCox = Cox_PH(surv_dat, covariates, EQMI_vector, degree = 1)
formattable::formattable(SurvCox, align = c('l', 'r'))
```

<table class="table table-condensed">
<thead>
<tr>
<th style="text-align:left;">
Variable
</th>
<th style="text-align:right;">
p-value
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
EQMI\*\_12
</td>
<td style="text-align:right;">
0.06856923
</td>
</tr>
</tbody>
</table>

Below, we create a matrix named recur_dat with all the subject IDs and
their recurrence outcomes, and use the same covariates matrix as
earlier.

``` r
recur_dat = clinical_data[,c(1:3)]
RecurCox = Cox_PH(recur_dat, covariates, EQMI_vector, degree = 1)
formattable::formattable(RecurCox, align = c('l', 'r'))
```

<table class="table table-condensed">
<thead>
<tr>
<th style="text-align:left;">
Variable
</th>
<th style="text-align:right;">
p-value
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
EQMI\*\_12
</td>
<td style="text-align:right;">
0.06092791
</td>
</tr>
</tbody>
</table>

## Few more tools: Compute and plot univariate and bivariate kernel densities of the markers

#### Univariate marginal density

We provide a few basic functions to estimate and plot the univariate
marginal densities of the random variables. We go back to the matrix
named Data_matrix which has 5000 samples and 5 columns corresponding to
five r.v.’s. We look at first two columns (call them, X_1 and X_2) and
estimate their kernel density estimates as f_1 and f_2 using the
function dens_univ. Next, using the function univ_dens_plot, we plot f_1
and f_2. The default number of grids is, ngrids = 1024 and the bandwidth
parameter is selected using Silverman’s ‘rule of thumb.’

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
matrix used is the multivariate plug-in bandwidth matrix Wand, M.P. &
Jones, M.C. (1994).

``` r

f_12 = biv_dens(cbind(X_1, X_2), ngrids = 512)
q1 = quantile(X_1, 0.7)
q2 = quantile(X_2, 0.7)
p = biv_dens_plot(f_12, maxs = c(q1, q2))
```

<img src="README_files/figure-gfm/Computing biv density-1.png" width="100%" />

### References

a\) Keren, L., Bosse, M., Marquez, D., Angoshtari, R., Jain, S., Varma,
S., … & Angelo, M. (2018). A structured tumor-immune microenvironment in
triple negative breast cancer revealed by multiplexed ion beam imaging.
Cell, 174(6), 1373-1387.

b\) Principe, J. C. (2010). Information theoretic learning: Renyi’s
entropy and kernel perspectives. Springer Science & Business Media.

### Appendix

#### Expression of the measures

*x* × *y*
