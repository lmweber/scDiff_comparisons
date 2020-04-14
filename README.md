# scDiff comparisons

This repository contains scripts to compare `scDiff` (developed by Simone Tiberi) with `diffcyt` methods (Weber et al. 2019), using CyTOF datasets previously used for the `diffcyt` evaluations (Weber et al. 2019) and distributed via the `HDCytoData` package (Weber and Soneson, 2019).


## Contents

- [example](example/): directory containing small example script showing how to run `scDiff`
- [scripts/DS](scripts/DS): scripts to run `diffcyt-DS` methods (i.e. methods for differential state testing) on `Weber_BCR_XL_sim` benchmark datasets from the `HDCytoData` package

    - main simulations
    - null simulations
    - 'less distinct' simulations

- [scripts/DA](scripts/DA): scripts to run `diffcyt-DA` methods (i.e. methods for differential abundance testing) on `Weber_AML_sim` benchmark datasets from the `HDCytoData` package

    - main simulations



## Additional details

See [Weber et al. (2019)](https://www.nature.com/articles/s42003-019-0415-5), [Weber and Soneson (2019)](https://f1000research.com/articles/8-1459), or help files from the [`HDCytoData` package](https://bioconductor.org/packages/HDCytoData) for more details on the benchmark datasets.


