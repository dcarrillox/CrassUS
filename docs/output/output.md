## Output

## Main output
Final results of the pipeline are in `results/<analysis_id>/crassus_results.tsv`.
This is a tabular file where rows represent contigs identified as candidates
_Crassvirales_. Fields below are provided for each of these contigs

- **crassus_id**: unique CrassUS identifier
- **contig**: user contig identifier
- **sample**: sample of the contig
- **length**: length of the contig
- **len/taxa_len**: relation between contig length and the average genome length
of the deepest predicted _Crassvirales_ taxa.
- **ref_taxa**: deepest predicted _Crassvirales_ taxa for the contig
- **DTR**: Direct Terminal Repeats detected (`True` or `False`)
- **family**: Predicted _Crassvirales_ family, if any
- **subfamily**: Predicted _Crassvirales_ subfamily, if any
- **genus**: Predicted _Crassvirales_ genus, if any
- **species**: Predicted _Crassvirales_ species, if any
- **evidence_family**: Signals supporting the family assignment
- **evidence_genus**: Signals supporting the genus assignment
- **notes**: warnings about the assignments
- **discard**: `True` if there was not enough evidence to classify the contig at
least at the family level.
