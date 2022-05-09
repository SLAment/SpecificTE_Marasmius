# SpecificTE: Searching for somatic TE insertions in *Marasmius oreades*

A pipeline based on Valentina Peona's "SpecificTE" pipeline to detect transposable element (TE) insertions and empty sites. I modified it to find also sites that might have insertion shared by two or more samples (but at least one with an empty site). The goal is to find novel TE insertions as somatic mutations in a fairy ring of *Marasmius oreades*.

In order to produce the results, there are two Snakemake workflows:

- 1_RepeatsMarasmius
- 2_SpecificTE

You need to run them in order (RepeatsMarasmius first of course).

The results are part of Hiltunen et al. (in prep.).