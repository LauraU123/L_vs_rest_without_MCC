
# Comparison of L and the rest of the genome

This workflow constructs two parallel phylogenetic analyses with the same sequences.
These can be viewed as a tanglegram on nextstrain.org to check for duplication in RSV.

## The L gene

L is the largest gene in RSV, which encodes RNA-dependent RNA polymerase.
Comparing its evolution to the rest of the genome can check
whether recombination occurs.

## Workflow 

This workflow uses a similar principle to the main nextstrain/rsv workflow, 
but does not include clades or glycosylation, and splits the sequence alignment into two,
constructing two separate analyses based on the same sequences. 


* Input: RSV reference files, metadata and sequences

* Output: 2 annotated phylogenetic trees, one for all genes excluding L, and one for only L



