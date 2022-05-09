# SpecificTE_mar: Searching for somatic TE insertions in *Marasmius oreades*

A pipeline based on Valentina Peona's "SpecificTE" pipeline to detect transposable element (TE) insertions and empty sites. I modified it to find also sites that might have insertion shared by two or more samples (but at least one with an empty site). The goal is to find novel TE insertions as somatic mutations in a fairy ring of *Marasmius oreades*.

I originally developed this pipeline for my PhD thesis working with *Podospora anserina* data, so the strain codes and TE names in the examples of this README do not match the *Marasmius* names (e.g. Podan2).

The pipeline is very stringent in filtering, so it has low sensitivity. However, the final output is a number of alignments with either all samples ("complete") or some ("incomplete"), that are to be inspected manually by the user. So after that, one could say that the accuracy is high. 

See [Valentina's Github](https://github.com/ValentinaBoP/NeurosporaSpecificTE/blob/master/Input_preparation) for the details. Her pipeline is meant to identified only the empty sites themselves.

Unfortunately, at the end of the pipeline I ended up using a trick to get all the alignments in Snakemake, basically producing a dummy file for every repetitive element, which serves as a placeholder for the pipeline to know it's finished. It doesn't check for the alignments themselves.

Naturally, the pipeline needs at least two strains to run. In this version of the pipeline, the input are the genomes and a repeat library for repeat masking with [RepearMasker](http://www.repeatmasker.org/). The assemblies should be in a single folder given in the configuration file. As it is, the pipeline expects the genomes to have a name like "{sample}.fa", where {sample} is the name of each strain/individual. This is easy to modify if you know a bit of Snakemake.


## Building the environment

I ran the pipeline under a [Conda](https://docs.anaconda.com/) environment. Install it first. If you like, you can start by updating it.

    $ conda update -n base conda

To create the environment arbitrarily named `SpecificTE`:

    $ conda create -n SpecificTE -c bioconda snakemake-minimal=5.8.1

To install software, activate the environment.

    $ conda activate SpecificTE

Now install:

    $ conda install -c bioconda samtools=1.10
    $ conda install -c bioconda bedtools=2.29.2
    $ conda install -c bioconda blast=2.10.1
    $ conda install -c bioconda mafft=7.407
    $ conda install numpy=1.19.1
    $ conda install biopython=1.77
    $ conda install -c bioconda repeatmasker=4.0.7

## The configuration file

The configuration file contains the paths to the necessary files to run the pipeline. The pipeline has to be run twice, one for all the assemblies from nucleartype A and one for nucleartype B. Below the configuration file is ready for the type A, but for B you have to comment A and uncomment B. **Run this pipeline in a different folder for each nucleotype!**

    $ cat SpecificTE_mar_config.yaml
```yaml
## Samples (at least two samples needed)
## Nucleotype A (also known as HDa)
SampleIDs: ["E1", "W6", "N11", "S3"]
## Nucleotype B (Hdb)
# SampleIDs: ["E2", "W9", "N18", "S1"]

## Path to assemblies and repeatmasker gffs (without trailing '/'!)
path2assemblies: "path/to/assemblies"
path2repeats: "../1_RepeatsMarasmius/RepeatMasker" # Produced after running this pipeline

## File with the names of all the repeats in the library
outprefix_file: "data/allsamples_repeats_count_min0.txt" # provided in the repository

## Filtering parameters
minRMlen: 120 # Threshold for the final filtering, minimum size of a RepeatMasker feature in the raw data to be considered. Smallest non-satellite element in library is MarorCMC-4, which is 128 bp. 
BUFFER: 5000 # The size of the flanking regions to find empty sites
PERBUF: 0.8 # percentage of BUFFER that should be present in the empty sites
PERID: 99 # Minimal percentage of identity for BLAST searches of the flanking regions of TE insertions
VICINITY: 45000 # Max size of a TE + flanks (the largest normal repeat in the library is maror5-1891#LINE/L2 (~17052 bp))
HOWCLOSE: 20 # How close should the flanks be in the empty sites

```

## Run pipeline locally

Get into the folder with this repo's content, and activate environment:
	
	$ cd path/to/my/analysis
    $ conda activate SpecificTE

First, to get an idea of how the pipeline looks like we can make a rulegraph:

    $ snakemake --snakefile SpecificTE_mar.smk --configfile SpecificTE_mar_config.yaml --rulegraph | dot -Tpng > rulegraph.png

<!-- ![rulegraph](rulegraph.png "rulegraph of SpecificTE_mar.smk") -->

If you want to make sure that all the input files are placed correctly you can check with:
    
    $ conda activate SpecificTE
    $ snakemake --snakefile SpecificTE_mar.smk --configfile SpecificTE_mar_config.yaml -pn

Notice that it will take a while to just calculate the graph (that is, to plan the whole set of steps to perform) if you have many samples and/or many TE to look at. For example, with my configuration file, there were 62790 jobs in total to perform.

If that step above didn't complain and you feel ready, you can finally run the pipeline. I like to make a screen first, then activate the environment, and finally run the pipeline in the background.

    $ screen -R SpecificTE
    $ conda activate SpecificTE
    $ snakemake --snakefile SpecificTE_mar.smk --configfile SpecificTE_mar_config.yaml -p -j 24 --keep-going &> SpecificTE_mar.log &

Notice `-j` stands for the number of threads that you want to give to your pipeline. See [Snakemake](https://snakemake.readthedocs.io/en/stable/) documentation for more information.

----

There is a variant of the pipeline that does that runs RepeatMasker internally instead of taking the path to the RepeatMasker output. In that case, the pipeline is ran like this to test:

    $ conda activate SpecificTE
    $ snakemake --snakefile SpecificTE_rm_mar.smk --configfile SpecificTE_rm_mar_config.yaml -pn

And the full analysis:

    $ screen -R SpecificTE
    $ conda activate SpecificTE
    $ snakemake --snakefile SpecificTE_rm_mar.smk --configfile SpecificTE_rm_mar_config.yaml -p -j 24 --keep-going &> SpecificTE_rm_mar.log &

## Results

The pipeline produces "dummy" files in the end to make sure it ran correctly. Ignore those, but don't erase them unless you want to re-run the last steps again. The interesting output files are:

* "results/shared/{outprefix}_potentialy_shared_loci.txt" -> The BLAST results of potentially shared loci based on pairwise comparisons, so they are redundant. The "outprefix" variable stands for the name of the element.

For example:
    
    Podan2_Rana_LTR__1894   CBS237.71m_Rana_LTR__2200   100.000 1002    0   0   1   1002    1   1002    0.0 1808
    CBS237.71m_Rana_LTR__2200   Podan2_Rana_LTR__1894   100.000 1002    0   0   1   1002    1   1002    0.0 1808

In this case a *solo* of a Rana element (number 1894) in the strain Podan2 matches the location of a Rana LTR (number 2200) in the strain CBS237.71m. The vice-versa is also true (which should be there, otherwise something might have been wrong with the BLAST searches). The columns are: query id, subject id, percentage of identity, alignment length, number of mistmatches, number of gaps, query start, query end, subject start, subject end, e-value, and bit score.

* "results/summary_insertions.tbl" -> The number of unique insertions for each element per strain.

Example:
    
    Genome_TE   Dendrobates Discoglosse Pelobate
    Podan2  0   11  5
    CBS237.71m  0   7   2
    PcWa139m    0   1   1
    CBS112042p  0   3   4
    CBS415.72m  0   8   1
    CBS411.78m  0   2   4
    CBS124.78p  0   8   0

That means that there are 11 Discoglosse unique insertions in the strain Podan2 that are not present in any other strain.

* "results/emptysites/{outprefix}_matrix_emptysites.txt" -> A matrix with pairwise counts of number of confirmed insertions. 

For example:

    Genome  Podan2  CBS237.71m  PcWa139m
    Podan2  0   11  2
    CBS237.71m  20  0   2
    PcWa139m    19  12  0

The target is the rows, and the query is the columns. So there are 20 empty sites in CBS237.71m compared to the genome of Podan2.

* "results/alignments/complete" -> A folder containing all the non-redundant alignments of interesting sites to look at. Most of them are unique insertions that were not caught by the pipeline due to filtering or other issues, and once you are sure it's a unique insertion you can manually add the count to the "results/emptysites/{outprefix}_matrix_emptysites.txt" file. Ideally, you should find a few sites where some strains are empty and some are filled with the TE. Those would be the synapomorphic insertions. Be careful, sometimes you will see two alignments that look the same. Likely, these are not repeated, but in fact there are two different insertions very close to each other (or even nested!). Always check the name of the file to know what the focal TE is.

* "results/alignments/incomplete" -> The same as above, but not all strains made it through the filtering. This happens because sometimes other TEs insert in the flanks and destroy the area, or the locus was deleted altogether. 
