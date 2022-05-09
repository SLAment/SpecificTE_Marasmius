# -*- snakemake -*-

### RepeatsMarasmius: Small pipeline to repeat mask Marasmius oreades genomes
#############################################################################

# The objective is to produce gtf/gffs of repeats in input genomes plus a few
# descriptive plots

#############################################################################
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2020/10/30
# +++++++++++++++++++++++++++++++++++++++++++++++++

from glob import glob
import sys

# -------------------------------------------------
# Data
samples = config["SampleIDs"]
path2data = config["path2data"]
TElib = config["TElib"]

# Scripts
gtfRM2gff = config["gtfRM2gff"]
totalcovergff = config["totalcovergff"]
RepeatsMar = config["RepeatsMar"]

# Environments
plotr = config["plotr"]

# Parameters
MINCOUNT = int(config["mincount"])
outprefix = config["outprefix"]
# -------------------------------------------------


rule all:
	input:
		expand("RepeatMasker/{sample}.repeatmasker.gff3", sample = samples), # Just for vanity
		f"results/{outprefix}_repeats_count_min{MINCOUNT}.txt" # The actual result
	

# -------- Prepare data

## Make a dictionary of the samples (key) and their assemblies (values)
# This way allows me to be more flexible with the names of the files
# assemblies = [(glob(path2data + "/{sample}*.fa".format(sample=sample))) for sample in samples] # the *.fa is hardcoded later
# ASSEMBLIESDIC = dict(zip(samples, assemblies))

# --------

rule RepeatMasker:
	""" Run RepeatMasker on the input genomes """
	input:
		genome = f"{path2data}/" + "{sample}.fa",
		TElib = TElib, # Custom library
		# genome = lambda wildcards: ASSEMBLIESDIC[wildcards.sample][0]
	output:
		"RepeatMasker/{sample}.fa.out.gff"
	params:
		threads = 10,
	shell:
		"RepeatMasker -pa {params.threads} -a -xsmall -gccalc -gff -excln -lib {input.TElib} -dir RepeatMasker {input.genome}; "

rule gtfRM2gff: # This is just for vanity
	""" Transform the gtf of RepeatMasker into a gff3 file """
	input:
		gtfRM2gff = gtfRM2gff,
		gtf = "RepeatMasker/{sample}.fa.out.gff"
	output:
		"RepeatMasker/{sample}.repeatmasker.gff3"
	shell:
		"python {input.gtfRM2gff} {input.gtf} > {output}"

rule plot:
	""" Make histogram of repeat counts and a list of counts"""
	input:
		gtf = "RepeatMasker/{sample}.fa.out.gff",
	output:
		hist = "results/{sample}_TEhist.pdf",
		tbl = "stats/{sample}_repeats_count.txt",
	conda: 
		plotr
	params:
		threads = 1,
	script:
		RepeatsMar	

rule reducelist:
	""" Make list of repeats with a minimum threshold of appearances in all samples """
	input:
		expand("stats/{sample}_repeats_count.txt", sample = samples)
	output:
		f"results/{outprefix}_repeats_count_min{MINCOUNT}.txt"
	run:
		# sampledic = {} #This is for more control but it's not necessary
		chosenrepeats = [] # Here we save all the repeats that meet the criteria

		for file in input:
			samplelist = open(file, 'r')
			# thissample = samples[input.index(file)]
			# sampledic[thissample] = []

			for line in samplelist:
				tab = line.rstrip("\n").split("\t")
				repeat = tab[0]
				countte = int(tab[1])

				if countte >= MINCOUNT:
					chosenrepeats.append(repeat)
					# sampledic[thissample].append(repeat)

		if len(chosenrepeats) == 0:
			print(f"The threshold {MINCOUNT} is too high and no repeat survived :(!")
			sys.exit(1)

		ofile = open(output[0], 'w')

		for repeat in set(chosenrepeats):
			ofile.write(repeat + '\n')

## --- Characterize more the genome

rule totalcovergff: # This is just for vanity
	""" Calculate the total repeat coverage in the genome """
	input:
		totalcovergff = totalcovergff,
		gff = "RepeatMasker/{sample}.repeatmasker.gff3",
		fasta = lambda wildcards: ASSEMBLIESDIC[wildcards.sample][0]
	output:
		"stats/{sample}_repeat_cov.txt"
	shell:
		"python {input.totalcovergff} {input.gff} --fasta {input.fasta} > {output}"


