# vacamole
`vacaomle` is an R package that features a deterministic age-structured compartmental susceptible-exposed-infectious-recovered (SEIR) model and supporting functions. The model was developed to determine the impacts of different COVID-19 vaccination strategies on disease outcomes; therefore the model is an extended SEIR model to include states for severe disease outcomes and vaccination status.

# Installation
------------

1.  Install [R](http://cran.r-project.org)

2.  Install the development version of serosolver from [GitHub](https://github.com/kylieainslie/vacamole):

``` r
devtools::install_github("kylieainslie/vacamole")
library(serosolver)
```

# Quick start and vignettes
-------------------------

Read the [quick start vignette](https://seroanalytics.github.io/serosolver/articles/serosolver-quick_start_guide.html) to set up and run a simple implementation of the model.

