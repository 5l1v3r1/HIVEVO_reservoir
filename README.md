### Establishment and stability of the latent HIV-1 DNA reservoir

This repository accompanies the manuscript "Establishment and stability of the latent HIV-1 DNA reservoir" by Brodin et al. The code is by Richard Neher, with bits and pieces contributed by Fabio Zanini.

The directory `manuscript` contains the latex source and bibliography of the article, while the directory `src` contains python scripts used to generate the figures in the manuscript. The process data necessary to generate the figures and tables is in the `data` directory.

The source code is split into a number of separate scripts
 * `process_samples.py` generated processed data from the raw reads, which can be obtained from the short read archive
 * `make_alignments.py` splits the reads into hypermutated and normal sequences and aligns them separately to each other and the RNA sequences obtained in a separate study.
 * `make_trees.py` infers approximate phylogenetic trees of DNA and RNA reads from each patients. Annotated figures are produced by the script `draw_trees.py`
 * `root_to_tip.py` analyzes the evolution of RNA and DNA sequences and infers evolutionary rates in both compartments
 * `hypermutation_statistics.py` calculates the fraction of hypermutated reads, the number of stop codons per read, the fraction of obviously deficient reads, and plots the distribution of mutations
 * `clonal_expansion_recapture.py` analyzes the persistence of haplotypes and calculates the fraction of haplotypes that are found in multiple samples.

The scripts require access to deep sequencing data via `HIVEVO_access` and make use
of the `seqan` library with `seqanpy` bindings.