# Overview
-----------
`vacamole` is an R package that features a deterministic age-structured compartmental susceptible-exposed-infectious-recovered (SEIR) model and supporting functions. The model was developed to determine the impacts of different COVID-19 vaccination strategies on disease outcomes; therefore the model is an extended SEIR model to include states for severe disease outcomes and vaccination status. A brief model description is given below. For full details, the model is described by [Ainslie et al.](https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2022.27.44.2101090)

## Model description
The model is a deterministic compartmental SEIR model. The population is partitioned into nine 10-year age groups (0-9, 10-19, …, 70-79, 80+). Within each age group we further partition the population into those who are unvaccinated, vaccinated with 1 dose, or vaccinated with 2 doses and then finally into disease states: susceptible (S), infected but not yet infectious (E), infectious (I), hospitalized (H), in intensive care (IC), return to the hospital ward after intensive care (HIC), recovered (R), and dead (D) (Figure 1). 

![Figure 1. Basic conceptual model diagram. This diagram does not include the additional states after the second dose of vaccination or the age structure in the model. S = susceptible, E = exposed, I = Infectious, R = Recovered, H = hospitalized, IC = In intensive care, HIC = return to the hospital ward following treatment in IC, Su = vaccinated, but not yet protected, D = dead. States with subscript V indicate individuals who are vaccinated and protected by vaccination. This model assumes the “leaky” vaccine protection, so vaccinated and protected individuals can still be infected, hospitalized, etc. but at a reduced rate.](vignettes/model_diagram.png)

We assume that individuals who recover from infection cannot be re-infected. When a person is vaccinated, they first enter a hold state where they are vaccinated, but not yet (fully) protected (Shold1d or Shold2d). After a delay period, they enter the vaccinated and protected state for the dose they have received (Sv1d or Sv2d). We assume that vaccination only affects susceptible individuals.

The model is designed to incorporate a single vaccine product with a 2-dose regimen that 1) reduces susceptibility to infection, 2) reduces risk of hospitalization if a vaccinated individual is infected, and 3) reduces risk of infecting others (transmission) if a vaccinated person is infected. The vaccine is assumed to provide “leaky” protection, which means that the vaccine reduces the probability of transitioning between states for each vaccinated person. For example, a vaccine with 70\% effectiveness against infection after the first dose would reduce the probability of transitioning from Sv1d to Ev1d by 70\% compared to someone who is unvaccinated (moving from the S compartment to the E compartment). However, since there is more than one vaccine product currently licensed for use against COVID-19, we incorporate different vaccine products by taking the daily weighted average of the number of people with each vaccine product (and dose), the corresponding delay to protection of each vaccine product, and the vaccine effectiveness against each outcome. Rate of vaccination for each day, vaccine product, and dose for each age group are model inputs. 

# Installation
---------------

1.  Install [R](http://cran.r-project.org)

2.  Install the development version of vacamole from [GitHub](https://github.com/kylieainslie/vacamole):

``` r
devtools::install_github("kylieainslie/vacamole")
library(vacamole)
```

# Quick start guide and vignettes
----------------------------

Read the [quick start vignette](https://kylieainslie.github.io/vacamole/articles/quick_start_guide.html) to set up and run a simple implementation of the model.

Additional vignettes are under development.

