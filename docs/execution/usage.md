# Usage

Being in its root directory and with its environment activated, you can run
CrassUS like this:

~~~
(crassus)$ snakemake -j 16 --use-conda --conda-frontend mamba --use-singularity
~~~

CrassUS handles its software dependencies by using the conda environments under
`workflow/envs`. These environments will be created automatically the first time
you execute the pipeline. 
