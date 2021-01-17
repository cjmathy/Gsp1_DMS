# 2020-02-04 Installing bio3d

- install.packages('bio3d') from Rstudio, which is using the base R (if this is a problem, can eventually set up Rstudio to use the R installed by conda, but it will become more complicated with where "install.packages()" puts packages, I think

- install dssp, muscle (and netcdf?) with conda

```
$ conda install -c salilab dssp
$ which mkdssp
> /Users/cjmathy/miniconda3/envs/lab/bin/mkdssp

$ conda install -c bioconda muscle
$ which muscle
> /Users/cjmathy/miniconda3/envs/lab/bin/muscle
```