This repository contains R functions to derive etimates and generate simulated data according to the methods and models described in the paper "Bias-corrected maximum-likelihood estimation of multiplicity of infection and lineage frequencies". To import/merge molecular data of any type (STR, SNPs, amino acids) and format and apply further analysis please refer to the the R package "MLMOI" (https://cran.r-project.org/web/packages/MLMOI/index.html).

Description of the functions included in the repository:

1. **MLE**: derives the maximum-likelihood estimate (MLE) of the MOI parameter (Poisson parameter) and the lineage (allele) frequencies;
2. **BCMLE**: derives the bias-corrected maximum-likelihood estimate (BCMLE) of the MOI parameter (Poisson parameter) and the lineage (allele) frequencies;
3. **HBCMLE1** derives the 1st version of heuristically bias-corrected maximum-likelihood estimate (HBCMLE1) of the MOI parameter (Poisson parameter) and the lineage (allele) frequencies;
4. **HBCMLE2** derives the 2nd version of heuristically bias-corrected maximum-likelihood estimate (HBCMLE2) of the MOI parameter (Poisson parameter) and the lineage (allele) frequencies;
5. **HBCMLE3** derives the 1st version of heuristically bias-corrected maximum-likelihood estimate (HBCMLE3) of the MOI parameter (Poisson parameter) and the lineage (allele) frequencies;
6. **second_order_bias** derives the approximated second-order bias of the MLE;
7. **prob_pathological** derives the probability of pathological data;
8. **cpoiss** generates conditionally Poisson distributed numbers;
9. **mnom** generates multinomially distributed random vectors;
10. **cnegb** generates negative binomial distributed random numbers.
11. **Cramer-Rao lower bounds**  derives Cramer-Rao lower bounds (CRLB) of the model parameters
12. **simulation_bc** the core function for the simulation study
