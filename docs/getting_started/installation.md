# Installation

## Software dependencies
It is recommended to use a [conda](https://docs.conda.io/projects/conda/en/latest/)
environment to handle dependencies. To do so, CrassUS comes with the file
`crassus_env.lock`, enumerating the software needed to run the pipeline. To get
a working conda environment with these software installed:

~~~
$ git clone https://github.com/dcarrillox/CrassUS.git
$ cd CrassUS

# Note the long notation --file flag; -f will not work.
$ conda create -n crassus --file=crassus_env.lock

# Activate it - use the name you gave above, if it is different
$ conda activate crassus

# The (crassus) prefix shows we have activated it
# Check the snakemake version
(crassus) $ snakemake --version
6.6.1
~~~
