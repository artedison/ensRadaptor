# ensRadaptor

ensRadaptor makes conscutrucing biological ODE model easier. It works as a interface to [ens](https://github.com/artedison/ensRadaptor/blob/master/inst/extdata/enscode/README.md) mehtod<sup>1</sup>, a MCMC based fortran program for rapid ODE system parameterization. ens models biochemcial systems with partial obervation and fit the ODE to experimental measurement. Although ens is fast, it has diffcult input and output format, among which considerable information are repeated or not useful for users. Learning the ens grammar oftens takes long time. ensRadaptor is designed to simplify the model construction process and provide easy visualizations. ensRadaptor takes in descriptions of the metabolic systems and experimetnal measurement, and it output input files to the ens program. ensRadaptor adopted a format based approach ([link](https://github.com/artedison/ensRadaptor/tree/master/inst/extdata/template_format)). For explanation on input and output of ens program, please refer to documents ([link](https://www.dropbox.com/sh/72zd1nxba6xxzvs/AAA4UMF1If-uIq_56iwL1AJ_a?dl=0)) or contact [Bernd Schuttler](https://www.physast.uga.edu/people/heinz_bernd_schuttler). For any questions about ensRadaptor, please contact [Yue Wu](https://mikeaalv.github.io)


## Examples
For using ensRadaptor, please refer to the script in the folder ./tests/workflowtest_open/

script.input.standalone.R: model contruction

script.output.standalone.R: ens result visualization. To produce the input files for this script, the user need to run the ens program based on the files produced by script.input.standalone.R. Specifically, the user can upload the folder 'input' to an HPC environment and run submit.sh. The output ens.o01 and ens.o02 will be used in script.output.standalone.R. Example output figures can be found [here](https://www.dropbox.com/sh/u2qd4llz45400yn/AAAWrqQwEoCTj0jCKOWmBJ7Ka?dl=0).

To install the package:

```
github_install
```

## Reference:

1. An ensemble method for identifying regulatory circuits with special reference to the qa gene cluster of Neurospora crassa
D. Battogtokh, D. K. Asch, M. E. Case, J. Arnold, H.-B. Schüttler
Proceedings of the National Academy of Sciences Dec 2002, 99 (26) 16904-16909; DOI: 10.1073/pnas.262658899
