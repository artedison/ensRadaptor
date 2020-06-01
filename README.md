# ensRadaptor

ensRadaptor makes conscutrucing biological ODE model easier. It works as a interface to [ens](https://github.com/artedison/ensRadaptor/blob/master/inst/extdata/enscode/README.md) mehtod<sup>1</sup>, a MCMC based fortran program for rapid ODE system parameterization. Although ens is fast, it has diffcult input and output format, among which considerable information are repeated or not useful for users. Learning the ens grammar oftens takes long time. ensRadaptor is designed to simplify the model construction process and provide easy visualizations. ensRadaptor takes in descriptions of the metabolic systems and experimetnal measurement, and it output input files to the ens program.

For using ensRadaptor, please refer to the script in the folder ./tests/workflowtest_open/

script.input.standalone.R: model contruction

script.output.standalone.R: ens result visualization

To install the package:

```

```

Reference:

1. An ensemble method for identifying regulatory circuits with special reference to the qa gene cluster of Neurospora crassa
D. Battogtokh, D. K. Asch, M. E. Case, J. Arnold, H.-B. Schüttler
Proceedings of the National Academy of Sciences Dec 2002, 99 (26) 16904-16909; DOI: 10.1073/pnas.262658899
