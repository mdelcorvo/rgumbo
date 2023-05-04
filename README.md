# rgumbo <img src="img/gumbo_logo.png" width="200" align="right" /> 

**<ins>(Under Development...)</ins>**

> *Copyright 2022 [Marcello Del Corvo](https://github.com/mdelcorvo). Licensed under the MIT license.*


[Gumbo](https://www.pillsbury.com/recipes/creole-chicken-gumbo-soup/15686514-d713-4056-bd2e-60cd68b20f9d), like any other soup, is a mixture of different ingredients. This one in particular is made with a dark roux, vegetables, chicken, sausage, shrimp and served over rice. 
Except that instead of tasty food and fresh ingredients, `rgumbo` provides you with R functions.

This package is a result of me constantly breaking the DRY principle
by copy-and-pasting functions from old projects into new ones. Hence, the
functions in `rgumbo` do not have a single common topic, but they are all either
related to manipulating genomic data or general
productivity utilities.

This package  has several functions that can prove to be useful and time-saving if you happen to need
to perform one of the tasks implemented by `rgumbo`.

## Installation

rgumbo is currently only available through GitHub and can be downloaded
easily using devtools.

```
# install.packages("devtools")
devtools::install_github("mdelcorvo/rgumbo")
```

## Getting started

There are many different usecases for rgumbo.  See the
[overview vignette](https://github.com/mdelcorvo/rgumbo/develop/vignettes/overview.md)

```
browseVignettes("rgumbo")
vignette("overview", "rgumbo")
```

Alternatively, see the help file for any specific function for a complete
detailed explanation of the function. For example `?rgumbo::LiftoverVcf`.

Below is a very short introduction to the functions in rgumbo. You need to load 
the package first. `library("rgumbo")`.

### `usePackage()` function: check, install and load all required packages at one
Function to automatically check, install and load all required packages from both CRAN and Bioconductor repositories.

```
usePackage(pkgs = c('data.table','purrr','clusterProfiler'))
```

### `LiftoverVcf()` function: lifts over a VCF file from one genome build to another
It produces a properly headered, sorted and compressed VCF in one go.

```
LiftoverVcf(vcf,output,from="hg38",to="hg19")
```

### `hgLiftOver()` function: lifts over a bed-like files from one genome build to another 
Converts genome coordinates between assemblies in a bed-like dataframe with chromosome number, start and (eventually) end position.

```
hgLiftOver(bed,from="hg38",to="hg19")
```
