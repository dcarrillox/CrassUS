# Configuration

## Configuration file

The `config_template.yaml` file provided with this repo has all available
configurable options. Short explanations are provided as commented blocks for
each option.

A configuration for the workflow must be available as a `config.yaml` within the
 `config` directory. You can either copy the `config_template.yaml` and rename
 it to `config.yaml` or make your edits straight on the template and rename it to
`config.yaml`.

**All fields included in the template must be specified**, unless otherwise stated
in the comment.

## Sample sheet

Samples to analyze are specified with a tabular samples sheet. Path to this
file should be in the `sample_sheet` section of the configuration file. The
structure should be as follows:

~~~
$ cat samples.tsv
analysis_id	sample_id	fasta
my_analysis	sample_1	/path/to/sample_1/sample_1.fasta
my_analysis	sample_2	/path/to/sample_2/sample_2.fasta
my_analysis	sample_3	/path/to/sample_3/sample_3.fasta
~~~

Columns names should be **analysis_id**, **sample_id** and **fasta**. The current
workflow accepts multiple samples under `sample_id`, but only one analysis name
under `analysis_id`. Input FASTA files under `fasta` can have any extension but
can not be zipped.
