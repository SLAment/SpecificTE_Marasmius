# -*- snakemake -*-

### SpecificTE_mar: Searching for somatic TE insertions in Marasmius oreades
#############################################################################

# A pipeline based on Valentina Peona's "SpecificTE" pipeline to detect TE
# insertions and empty sites. The goal is to find novel TE insertions as
# somatic mutations in a fairy ring of Marasmius oreades.

# See: https://github.com/ValentinaBoP/NeurosporaSpecificTE

#############################################################################
# ==================================================
# Sandra Lorena Ament Velasquez
# Johannesson Lab, Evolutionary Biology Center, Uppsala University, Sweden
# 2020/10/30
# +++++++++++++++++++++++++++++++++++++++++++++++++

from glob import glob
import numpy as np
import re # For regex
from collections import defaultdict #This tells python that the dictionary contains a list so you can freely append things to the value of each key

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Blast.Applications import * # NCBI BLASTX wrapper from the Bio.Blast.Applications module
from Bio.Align.Applications import MafftCommandline # To Use MAFFT, see http://biopython.org/DIST/docs/api/Bio.Align.Applications._Mafft.MafftCommandline-class.html
import os # To delete files and check if they are empty


# -------------------------------------------------
samples = config["SampleIDs"]
path2assemblies = config["path2assemblies"]
path2repeats = config["path2repeats"]
outprefix_file = config["outprefix_file"]

## Filtering parameters
minRMlen = config["minRMlen"]
BUFFER = config["BUFFER"]
PERBUF = config["PERBUF"]
PERID = config["PERID"]
VICINITY = config["VICINITY"]
HOWCLOSE = config["HOWCLOSE"]
# -------------------------------------------------

## Read names of the repeats
outprefix = [line.rstrip("\n") for line in open(outprefix_file)]

## Make lists to permutate between samples for blasting.
import itertools
combis = list(itertools.permutations(samples, 2))
targets = list(zip(*combis))[0]
queries = list(zip(*combis))[1]

rule all:
	input:
		# For manual inspection
		expand("results/alignments/dummies/{outprefix}.dummy", outprefix = outprefix), # Alignments of the potential shared loci
		expand("results/shared/{outprefix}_potentialy_shared_loci.txt", outprefix = outprefix), # A list with the pairwise homologous elements from different strains

		# Reports
		"results/summary_insertions.tbl", # Count numbers of confirmed insertions unique to one sample
		expand("results/emptysites/{outprefix}_matrix_emptysites.txt", outprefix = outprefix), # A count table of all empty sites of pairwise comparisons

# -------- Prepare data

## Make a dictionary of the samples (key) and their assemblies (values)
# This way allows me to be more flexible with the names of the files
assemblies = [(glob(path2assemblies + "/{sample}*.fa".format(sample=sample))) for sample in samples]
ASSEMBLIESDIC = dict(zip(samples, assemblies))


rule getdatalocally:
	""" Make links to data """
	input:
		lambda wildcards: ASSEMBLIESDIC[wildcards.sample][0]
	output:
		"data/genomes/{sample}.fa"
	shell:
		"ln -s {input} {output}"

rule indexgenome:
	""" Make the index files for each genome """
	input:
		"data/genomes/{sample}.fa"
	output:
		"data/genomes/{sample}.fa.fai"
	shell:
		"samtools faidx {input}"
		# -i, --reverse-complement Reverse complement sequences.

rule gff2bed:
	""" Modify the gff (actually a gtf) from RepeatMasker """ 
	input:
		path2repeats + "/{sample}.fa.out.gff"
		# path2repeats + "/{sample}.repeatmasker.gff"
	output:
		"data/RMSK/{sample}.bed"
	shell:
		"""
		grep -v "#" {input} | grep -v -E ")n|-rich" | cut -f1,4,5,9,6,7 | sed 's;Target \\"Motif:\(.*\)\\" \([0-9 ]*\);\\1;' | awk -v OFS='\\t' '{{print $1,$2,$3,$6,$4,$5}}'  > {output}
		# grep -v "#" {input} | grep -v -E ")n|-rich" | cut -f1,4,5,9,6,7 | sed 's;Target \\"Motif:\([a-zA-Z_-]*\)\\" \([0-9 ]*\);\\1;' | awk -v OFS='\\t' '{{print $1,$2,$3,$6,$4,$5}}'  > {output}
		"""

rule filterbed:
	""" Add an extra column to the bed file with numerated and filter for elements that are of interest """ 
	input:
		"data/RMSK/{sample}.bed"
	output:
		"data/RMSK/{outprefix}_{sample}.bed.filter" # The numbers of the elements match the line numbers of the input file
	# run:
		# # What's the index of the outprefix?
		# outprefixindex = outprefix.index(wildcards.outprefix)

		# # Search string for filtering
		# grepfilter = outprefix[outprefixindex] # This is the TE that we are looking at

		# cmd = """ awk -v OFS='\\t' '{{print $0""" + ',"' + wildcards.sample + """_"$4"__"NR}}' """ + f'{input} | awk "/{wildcards.outprefix}_/" > {output}'
		# shell(cmd)
		
		# grep will return zero only when some string is matched, so it trows an error status 1 if the element is not in the gff. Hence I changed it to awk
		## This bacame so awkward just so I can have the grep string as a variable, but otherwise it would be like:
	shell:
		"""
		awk -v OFS='\\t' '{{print $0,"{wildcards.sample}_"$4"__"NR}}' {input} | grep "{wildcards.outprefix}_" > {output}
		"""

rule filterbed_sense:
	""" Filter out the element in the + sense """ 
	input:
		"data/RMSK/{outprefix}_{sample}.bed.filter"
	output:
		"data/RMSK/{outprefix}_{sample}.bed.filter.sense"
	shell:
		"""
		awk '$6~/+/' {input} > {output}
		"""

rule filterbed_antisense:
	""" Filter out the element in the - sense """ 
	input:
		"data/RMSK/{outprefix}_{sample}.bed.filter"
	output:
		"data/RMSK/{outprefix}_{sample}.bed.filter.antisense"
	shell:
		"""
		awk '$6~/-/' {input} > {output}
		"""

rule get_flanks:
	""" Get the 3' and 5' flanks to the repeats """ 
	# https://bedtools.readthedocs.io/en/latest/content/tools/flank.html
	input:
		index = "data/genomes/{sample}.fa.fai",
		sense = "data/RMSK/{outprefix}_{sample}.bed.filter.sense",
		antisense ="data/RMSK/{outprefix}_{sample}.bed.filter.antisense"
	output:
		sense_flank3 = "data/Flanks/{outprefix}_{sample}.bed.filter.sense.flank3",
		sense_flank5 = "data/Flanks/{outprefix}_{sample}.bed.filter.sense.flank5",
		antisense_flank3 = "data/Flanks/{outprefix}_{sample}.bed.filter.antisense.flank3",
		antisense_flank5 = "data/Flanks/{outprefix}_{sample}.bed.filter.antisense.flank5",
	params:
		buffer = BUFFER
	shell:
		"""
		# antisense flanks
		bedtools flank -g {input.index} -i {input.antisense} -r {params.buffer} -l 0 -s | awk '$7=$7"_flank3"{{print $1,$2,$3,$7,$5,$6}}' OFS="\\t" > {output.antisense_flank3}
		bedtools flank -g {input.index} -i {input.antisense} -r 0 -l {params.buffer} -s | awk '$7=$7"_flank5"{{print $1,$2,$3,$7,$5,$6}}' OFS="\\t" > {output.antisense_flank5}
 		 # sense flanks
		bedtools flank -g {input.index} -i {input.sense} -r {params.buffer} -l 0 -s | awk '$7=$7"_flank3"{{print $1,$2,$3,$7,$5,$6}}' OFS="\\t" > {output.sense_flank3}
		bedtools flank -g {input.index} -i {input.sense} -r 0 -l {params.buffer} -s | awk '$7=$7"_flank5"{{print $1,$2,$3,$7,$5,$6}}' OFS="\\t" > {output.sense_flank5}

		# # antisense flanks
		# bedtools flank -g {input.index} -i {input.antisense} -r {params.buffer} -l 0 -s | awk '$7=$7"_flank3"{{print $0}}' OFS="\\t" > {output.antisense_flank3}
		# bedtools flank -g {input.index} -i {input.antisense} -r 0 -l {params.buffer} -s | awk '$7=$7"_flank5"{{print $0}}' OFS="\\t" > {output.antisense_flank5}
 	# 	 # sense flanks
		# bedtools flank -g {input.index} -i {input.sense} -r {params.buffer} -l 0 -s | awk '$7=$7"_flank3"{{print $0}}' OFS="\\t" > {output.sense_flank3}
		# bedtools flank -g {input.index} -i {input.sense} -r 0 -l {params.buffer} -s | awk '$7=$7"_flank5"{{print $0}}' OFS="\\t" > {output.sense_flank5}
		"""

# -b	Increase the BED/GFF/VCF entry by the same number base pairs in each direction. Integer.
# -l	The number of base pairs to subtract from the start coordinate. Integer.
# -r	The number of base pairs to add to the end coordinate. Integer.
# -s	Define -l and -r based on strand. For example. if used, -l 500 for a negative-stranded feature, it will add 500 bp to the end coordinate.


## --- Get the fasta file for the flanks ---

rule getflankfasta:
	""" Tranform bed into fasta """
	input:
		genome = "data/genomes/{sample}.fa",
		flank = "data/Flanks/{outprefix}_{sample}.bed.filter.{sense}.{flank}"
	output:
		"data/Flanks/{outprefix}_{sample}.bed.filter.{sense}.{flank}.fa"
	shell:
		"bedtools getfasta -fi {input.genome} -bed {input.flank} -nameOnly > {output}"

## --- BLAST the flanks ---

rule blastdbs:
	""" Make a BLAST database of each sample """
	input:
		assembly = "data/genomes/{sample}.fa",
	output:
		# "dbBLAST/{sample}_db/{sample}_db.ndb",
		"dbBLAST/{sample}_db/{sample}_db.nhr",
		"dbBLAST/{sample}_db/{sample}_db.nin",
		"dbBLAST/{sample}_db/{sample}_db.nog",
		# "dbBLAST/{sample}_db/{sample}_db.nos",
		# "dbBLAST/{sample}_db/{sample}_db.not",
		"dbBLAST/{sample}_db/{sample}_db.nsq",
		# "dbBLAST/{sample}_db/{sample}_db.ntf",
		# "dbBLAST/{sample}_db/{sample}_db.nto",
	params:
		threads = 1,
	shell:
		"""
		echo "Making database..."
		makeblastdb -in {input.assembly} -out dbBLAST/{wildcards.sample}_db/{wildcards.sample}_db -dbtype nucl -parse_seqids
		"""

rule BLASTflanks:
	""" BLAST the flanks of every pairwise comparison """
	input:
		# target = "data/genomes/{target}.fa",
		query = "data/Flanks/{outprefix}_{query}.bed.filter.{sense}.{flank}.fa",

		database = "dbBLAST/{target}_db/{target}_db.nhr",

		# Other files, so they get used and then removed by the temp directive
		# dbndb = "dbBLAST/{target}_db/{target}_db.ndb",
		# dbnin = "dbBLAST/{target}_db/{target}_db.nin",
		# dbnog = "dbBLAST/{target}_db/{target}_db.nog",
		# dbnos = "dbBLAST/{target}_db/{target}_db.nos",
		# dbnot = "dbBLAST/{target}_db/{target}_db.not",
		# dbnsq = "dbBLAST/{target}_db/{target}_db.nsq",
		# dbntf = "dbBLAST/{target}_db/{target}_db.ntf",
		# dbnto = "dbBLAST/{target}_db/{target}_db.nto",
	output:
		blast = "BLAST/{target}/{outprefix}_{target}-vs-{query}.{sense}.{flank}.txt",
		# blastfiltered = "BLAST/{target}/{target}-vs-{query}.{sense}.{flank}_stitch.txt",
	params:
		threads = 2, # Vale has 8, but her genomes are much bigger
		perid = PERID
	shell:
		"""
		outputblast=$(echo {output.blast}| sed 's/.txt//')

		database=$(echo {input.database}| sed 's/.nhr//')

		echo "Blasting ..."
		blastn -query {input.query} -out {output.blast} -db $database -perc_identity {params.perid} -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send sstrand length pident evalue' -num_threads {params.threads} -task blastn
		"""

rule BLAST2bed:
	""" Convert blast into BED format """
	input:
		"BLAST/{target}/{outprefix}_{target}-vs-{query}.{sense}.{flank}.txt",
	output:
		"BLAST/{target}/{outprefix}_{target}-vs-{query}.{sense}.{flank}.bed",
	shell:
		"awk '{{print $5, $7, $8, $1, $10, $9, $11, $3, $4}}' OFS='\\t' {input} | sed 's/plus/+/g' | sed 's/minus/-/g' | awk '{{if ($2 > $3) {{print $1, $3, $2, $4, $5, $6, $7, $8, $9, $10}} else {{print $0}} }}' OFS='\t' | bedtools sort -i stdin > {output}"
		# sseqid sstart send qseqid length sstrand pident qstart qend

rule makewindows: # FILTER_ALIGNMENT_putativeEmpty_id70.sh
	""" Make windows with the flank3 and flank5 bed files """
	input:
		flank3 = "BLAST/{target}/{outprefix}_{target}-vs-{query}.{sense}.flank3.bed",
		flank5 = "BLAST/{target}/{outprefix}_{target}-vs-{query}.{sense}.flank5.bed",
	output:
		"BLAST/{target}/{outprefix}_{target}-vs-{query}.{sense}.window",
	params:
		howclose = HOWCLOSE # How close should the flanks be in the empty sites
	shell:
		"bedtools window -sm -w {params.howclose} -a {input.flank5} -b {input.flank3} > {output}"

# -sm	Only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.
# -w	Base pairs added upstream and downstream of each entry in A when searching for overlaps in B. Default is 1000 bp.

rule findemptysites: # filterAlignment.R
	""" Find putative empty sites """
	input:
		windows = "BLAST/{target}/{outprefix}_{target}-vs-{query}.{sense}.window",
	output:
		txt = "Sites/{target}-vs-{query}/{outprefix}_{target}-vs-{query}.{sense}_putativeEmpty.txt",
		bed = "Sites/{target}-vs-{query}/{outprefix}_{target}-vs-{query}.{sense}_putativeEmpty.bed"
	params:
		halfbuffer = BUFFER*PERBUF
	run:
		tabopen = open(input.windows, 'r')

		# Save intermediate file with the coordinates of the putative empty sites
		ofile = open(output.txt, 'w')

		putativeempties = [] # The output list
		for line in tabopen:
			tab = line.rstrip("\n").split("\t")
			# Did both flanks hit the same specific element?
			element_5 = tab[3].replace("_flank5", "")
			element_3 = tab[12].replace("_flank3", "")

			if element_5 == element_3: # Yes!
				# filter for alignment length
				allen_5 = int(tab[4]) # Alignment length for flank 5
				allen_3 = int(tab[10]) # Alignment length for flank 3; since the other one is already large, probably this one will be too?
				qend_5 = int(tab[8]) # for position of the alignment: close to the interface
				qend_3 = int(tab[17])
				
				# filter for position of the alignment: close to the interface
				# if (allen_5 >= 100) and (qend_5 >= params.halfbuffer) and (qend_3 >= params.halfbuffer) : # Vale's original
				if (allen_5 >= params.halfbuffer) and (allen_3 >= params.halfbuffer) and (qend_5 >= params.halfbuffer) and (qend_3 >= params.halfbuffer) : 
					ofile.write(line) # Save it in the new file
					putativeempties.append(tab)
		ofile.close()
		
		# Make BED format of the putative empty sites
		obed = open(output.bed, 'w')
		
		for site in putativeempties:
			chr = site[0]
			start = min( int(site[1]), int(site[2]), int(site[10]), int(site[11]) ) 
			end = max( int(site[1]), int(site[2]), int(site[10]), int(site[11]) )
			name = site[3].replace("_flank5", "")
			score = site[4]
			strand = site[5]

			newline = '\t'.join([chr, str(start), str(end), name, score, strand]) + '\n'
			obed.write(newline)
		obed.close()

rule bedempty3fasta:
	""" Transform bed of putative empty sites into fasta """
	input:
		genome = "data/genomes/{target}.fa",
		bed = "Sites/{target}-vs-{query}/{outprefix}_{target}-vs-{query}.{sense}_putativeEmpty.bed"
	output:
		"Sites/{target}-vs-{query}/{outprefix}_{target}-vs-{query}.{sense}_putativeEmpty.fa"
	shell:
		"bedtools getfasta -fi {input.genome} -bed {input.bed} -nameOnly > {output}"

rule ConfirmSites: # The weakness of this step is a specific case of three LTR elements in a row. If the windows are too short, these nested elements will come back as "empty sites"
	""" BLAST the empty sites into the original genome to see if they do exist """
	# FILTER_ALIGNMENT_confirmEmpty_id70.sh
	input:
		# target = "data/genomes/{target}.fa",
		query = "Sites/{target}-vs-{query}/{outprefix}_{target}-vs-{query}.{sense}_putativeEmpty.fa",

		database = "dbBLAST/{target}_db/{target}_db.nhr",

		# Other files, so they get used and then removed by the temp directive
		# dbndb = "dbBLAST/{target}_db/{target}_db.ndb",
		# dbnin = "dbBLAST/{target}_db/{target}_db.nin",
		# dbnog = "dbBLAST/{target}_db/{target}_db.nog",
		# dbnos = "dbBLAST/{target}_db/{target}_db.nos",
		# dbnot = "dbBLAST/{target}_db/{target}_db.not",
		# dbnsq = "dbBLAST/{target}_db/{target}_db.nsq",
		# dbntf = "dbBLAST/{target}_db/{target}_db.ntf",
		# dbnto = "dbBLAST/{target}_db/{target}_db.nto",
	output:
		blast = "Sites/{target}-vs-{query}/{outprefix}_{target}-vs-{query}.{sense}_confirmEmpty.txt",
	params:
		threads = 2, # Vale has 4, but her genomes are much bigger
		perid = PERID
	shell:
		"""
		outputblast=$(echo {output.blast}| sed 's/.txt//')

		database=$(echo {input.database}| sed 's/.nhr//')

		echo "Blasting ..."
		blastn -query {input.query} -out {output.blast} -db $database -perc_identity {params.perid} -outfmt '6 qseqid qlen qstart qend sseqid slen sstart send sstrand length pident evalue' -num_threads {params.threads} -task blastn
		"""

rule filterConfirmSites:
	""" Filter the confirmed sites """
	input:
		"Sites/{target}-vs-{query}/{outprefix}_{target}-vs-{query}.{sense}_confirmEmpty.txt"
	output:
		"Sites/{target}-vs-{query}/{outprefix}_{target}-vs-{query}.{sense}_confirmEmpty.filter",
	params:
		ratio = PERBUF, # A ratio of the size of both sides of the flanks should be at least 80% for PERBUF = 0.8
		minlen = (BUFFER*2)*PERBUF, # the total window (500 + 500) should be at least 800 bp for PERBUF = 0.8
	shell:
		"awk '($10/$2) >= {params.ratio} {{print $0}}' {input} | awk '$2 >= {params.minlen} {{print $0}}' > {output}"
		
rule finalConfirmSites:
	""" Get the final confirmed sites, which have a single BLAST hit"""
	input:
		"Sites/{target}-vs-{query}/{outprefix}_{target}-vs-{query}.{sense}_confirmEmpty.filter"
	output:
		"Sites/{target}-vs-{query}/{outprefix}_{target}-vs-{query}.{sense}_confirmEmpty.final",
	shell:
		"awk 'NR==FNR {{count[$1]++; next}} count[$1]==1' {input} {input} > {output}"
		# Keep only the sites that had a single good hit

rule filter_fragments:
	""" Retain hits with minimum size and make a nicer table """
	input:
		rm_data = "data/RMSK/{outprefix}_{query}.bed.filter",
		final = "Sites/{target}-vs-{query}/{outprefix}_{target}-vs-{query}.{sense}_confirmEmpty.final",
	output:
		lista = temp("Sites/{target}-vs-{query}/{outprefix}_{target}-vs-{query}.{sense}_confirmEmpty.final.list"),
		listraw = "Sites/{target}-vs-{query}/{outprefix}_{target}-vs-{query}.{sense}_confirmEmpty.final.list.raw",
		finalfilter = "Sites/{target}-vs-{query}/{outprefix}_{target}-vs-{query}.{sense}_confirmEmpty.final.filter",
	params:
		threshold = minRMlen
	shell:
		"""
		# If no empty sites are found one of these commands emits a silent error, triggering Snakemake to fail (probably grep)
		if [ -s {input.final} ]; then
			# Get list of all surviving elements
			cut -f1 {input.final} > {output.lista}
			
			## Recover the original Repeat Masker line for those surviving elements
			# Add a new column with the length of the element
			grep -w -F -f {output.lista} {input.rm_data} | awk '{{print $0, $3 - $2}}' OFS="\\t" > {output.listraw}

			# Filter for sufficiently long elements
			awk -v threshold={params.threshold} '$8 >= threshold {{print $0}}' {output.listraw} > {output.finalfilter}
		else
			echo "File {input.final} is empty!!"
			touch {output.lista} {output.listraw} {output.finalfilter}
		fi

		"""

       # -F, --fixed-strings
       #        Interpret PATTERN as a list of fixed strings (instead of regular expressions), separated by newlines, any of which is to be matched.
       # -f FILE, --file=FILE
       #        Obtain patterns from FILE, one per line.  If this option is used multiple times or is combined with the -e (--regexp) option, search for all patterns given.  The empty file  contains  zero  patterns,  and  therefore  matches
       #        nothing.
       # -w, --word-regexp
       #        Select  only those lines containing matches that form whole words.  The test is that the matching substring must either be at the beginning of the line, or preceded by a non-word constituent character.  Similarly, it must be
       #        either at the end of the line or followed by a non-word constituent character.  Word-constituent characters are letters, digits, and the underscore.  This option has no effect if -x is also specified.

### --- What is next is not in the original SpecificTE pipeline of Valentina as such ---

rule combine_final_empty_sites:
	""" Put the sites of both senses together, sorted """
	input:
		expand("Sites/{{target}}-vs-{{query}}/{{outprefix}}_{{target}}-vs-{{query}}.{sense}_confirmEmpty.final.filter", sense = ["sense", "antisense"])
	output:
		"Sites/{target}-vs-{query}/{outprefix}_{target}-vs-{query}.both_confirmEmpty.final.filter"
	shell:
		"cat {input} | bedtools sort -i stdin > {output}"


rule compare_by_sample:
	""" Get empty sites together but by sample """
	input:
		lambda wildcards: expand("Sites/{target}-vs-{{query}}/{{outprefix}}_{target}-vs-{{query}}.both_confirmEmpty.final.filter", target = [x for x in targets if (x != wildcards.query)])
	output:
		"Empty_sites/{outprefix}_{query}.both_confirmEmpty.final.filter_raw"
		# temp("Empty_sites/{outprefix}_{query}.both_confirmEmpty.final.filter_raw")
	run:
		# Make new file
		ofile = open(output[0], 'w')

		for file in input:
			# Get the name of the target sample
			folderpattern = re.compile("\/([\w.]*)-vs-([\w.]*)/") # Assuming the folder structure of /{target}-vs-{query}/
			matchy = folderpattern.search(file)
			target = matchy.group(1)

			# Write the line of each file but appending the target at the end
			for line in open(file):
				newline = line.rstrip("\n") + f"\t{target}\n" # The target is the strain that is confirmed empty for that site
				ofile.write(newline)

rule sort_compare_by_sample:
	""" Sort and remove duplicates of confirmed empty sites """
	input:
		"Empty_sites/{outprefix}_{sample}.both_confirmEmpty.final.filter_raw"
	output:
		"Empty_sites/{outprefix}_{sample}.both_confirmEmpty.final.filter"
	shell:
		"cat {input} | bedtools sort -i stdin | sort -k1 -k2 -k3 -k9 | uniq > {output}" # I had to sort by the second column too so uniq catches them correctly


## --- Characterizing events ---

rule get_insertions:
	""" Sample specific insertions (everybody else is empty) """
	input:
		"Empty_sites/{outprefix}_{sample}.both_confirmEmpty.final.filter"
	output:
		nonsp = "Nonspecific/{outprefix}_{sample}.both_confirmEmpty.final.filter.non-specific",
		sp = "Insertions/{outprefix}_{sample}.both_confirmEmpty.final.filter.specific"
	run:
		# A dictionary of all elements found
		TEdic_coords = {}
		TEdic_empties = defaultdict(list) #This tells python that the dictionary contains a list so you can freely append things to the value of each key 

		for line in open(input[0]):
			tabs = line.rstrip("\n").split("\t")
			te_chr = tabs[0]
			te_start = tabs[1]
			te_end = tabs[2]
			te_sense = tabs[4]
			element = tabs[6].strip() # Remove a leading white space in the species name (Lore 2020.09.13: I think I fixed this and it's not there anymore, but just in case)
			emptyguy = tabs[8]

			if element not in TEdic_coords.keys(): # Each element is repeated a few times (max the number of strains) so just save this info once
				newline = '\t'.join((te_chr, te_start, te_end, te_sense, element)) + '\n'
				TEdic_coords[element] = newline

			TEdic_empties[element].append(emptyguy)

		# Open files
		spfile = open(output.sp, 'w')
		nonspfile = open(output.nonsp, 'w')

		## Print cases that are lineage specific
		for element in TEdic_empties.keys():
			othersamples = [x for x in samples if (x != wildcards.sample)]

			if set(TEdic_empties[element]) == set(othersamples): # All other samples were confidently detected
				spfile.write(TEdic_coords[element]) # These sites are only filled in one sample and empty in the rest
			else:
				nonspfile.write(TEdic_coords[element]) # Potentially some other individual also has the site filled

rule find_shared_loci:
	""" From the empty sites that were not specific to a single sample, collapse the ones that are homologous into a single one (non-redundant) """
	input:
		genomes = expand("data/genomes/{sample}.fa", sample = samples),
		elements = expand("Nonspecific/{{outprefix}}_{sample}.both_confirmEmpty.final.filter.non-specific", sample = samples),
	output:
		shared = "results/shared/{outprefix}_potentialy_shared_loci.txt", # A file to keep track of potentially shared elements (pairwise)
		nonreds = temp("temp/{outprefix}.nonspecific_nonredundant.txt"), 
		chosenTEs = temp("temp/{outprefix}_chosen_TEs.fa"), # A fasta file with the flanks of the TEs that (potentially) are non-specific 
		chosenTEs_hits = temp("temp/{outprefix}_chosen_TEs_hits.tab"),
		dbchosen = temp(directory("temp/{outprefix}_chosen_TEs_db")),
	params:
		buffer = BUFFER,
		evalue = "5e-10",
		threads = 4,
		perid = PERID,
	run:
		chosen_TEs = {} # A dictionary of all the elements that had confirmed empty sites, but not in all samples
		for file in input.elements:
			if os.stat(file).st_size != 0: # File is not empty
				for line in open(file, 'r'):
					tab = line.rstrip("\n").split("\t")
					element = tab[4]
					chrom = tab[0]
					start = tab[1]
					end = tab[2]

					# Get the genome and element type
					tepattern = re.compile("([a-zA-Z0-9.]*)_([\w-]*)_([\w]*)_([0-9]*)") # eg. PcWa139m_Yeti_LTR__1259 or N11_maror5-28__7513
					# tepattern = re.compile("([a-zA-Z0-9.]*)_([a-zA-Z]*)_([\w]*)_([0-9]*)") # eg. PcWa139m_Yeti_LTR__1259
					matchy = tepattern.search(element)
					te_sample = matchy.group(1)
					typete = matchy.group(2)

					chosen_TEs[element] = (te_sample, typete, chrom, start, end)

		# if the dictionary is not empty
		if chosen_TEs:
			chosen_TEs_seqs = [] # Here we'll keep the flanking sequences of all the chosen elements

			# -----------------------------------
			# Get the fasta sequence of each element
			# -----------------------------------
			for element in chosen_TEs.keys():
				te_sample = chosen_TEs[element][0]
				typete = chosen_TEs[element][1]
				chrom = chosen_TEs[element][2]
				start = int(chosen_TEs[element][3]) # Base 1
				end = int(chosen_TEs[element][4])

				## Find the assembly where this element is
				# --------------------------------
				assembly = f"data/genomes/{te_sample}.fa"
				# --------------------------------

				# Get the right chromosome fasta sequence
				for seq_record in SeqIO.parse(open(assembly, 'r'), "fasta", generic_dna):
					if chrom == seq_record.id:
						thischrseq = seq_record
						break

				## Slice just the flanks, not the TE itself
				# The bed coordinates are base 1
				teseq = thischrseq[start - 1 - params.buffer:start] + thischrseq[end - 1:end + params.buffer]
				teseq.id = element
				teseq.description = ''

				chosen_TEs_seqs.append(teseq)

			# -----------------------------------
			print("BLAST all TEs-flanks against each other to find homologous loci ...")
			# -----------------------------------
			# Write it to a fasta so it can be read
			SeqIO.write(chosen_TEs_seqs, output.chosenTEs, "fasta")

			# Make the mini-database of all sequences (it's temporary so then it will be erased)
			# Define the names of the databases
			databasename = f"{output.dbchosen}/{wildcards.outprefix}_chosen_TEs_db"
			createdb = f"makeblastdb -in {output.chosenTEs} -out {databasename} -dbtype nucl -parse_seqids" # BLAST command
			process = subprocess.Popen(createdb.split(), stdout=subprocess.PIPE) # pipe the command to the shell
			stdout, stderr = process.communicate() # run it

			# BLAST
			blast_command = NcbiblastnCommandline(query=output.chosenTEs, cmd='blastn', out=output.chosenTEs_hits, outfmt=6, db=databasename, evalue=params.evalue, num_threads=params.threads, task='blastn')
			stdout, stderr = blast_command()

			# Make a list of homologous loci
			nonredundant = []
			redundant = []
			selfhits = []

			sharedout = open(output.shared, "w") # Start a file to keep track of potentially shared elements (pairwise)

			for line in open(output.chosenTEs_hits):
				# head = "query_id	subject_id	percent_identity	alignment_length	N_mismatches	N_gaps	query_start	query_end	subject_start	subject_end	evalue	bit_score"
				hit = line.rstrip("\n").split("\t")
				query, ref = hit[0:2]
				pid = float(hit[2])
				allen = float(hit[3])

				# Remove redundancy by keeping track of the queries that already had a hit as references
				if (pid >= params.perid) and (allen >= ((params.buffer)*2)*PERBUF ): # The site might be full, so using the alignment length doesn't help

					if (query != ref):
						sharedout.write(line) # These are actually potentially shared sites, since an element was called for these two strains 
						
						if (query not in nonredundant) and (query not in redundant): 
							nonredundant.append(query)
							redundant.append(ref)
						elif query in nonredundant:
							redundant.append(ref)

					else: # Self-hit, save it for later to look for unique guys
						selfhits.append(query)

			## Add those that had no other homologs
			nonredundant = nonredundant + [hit for hit in selfhits if (hit not in nonredundant) and (hit not in redundant)]

			# # For debugging
			# print("nonredundant", nonredundant)
			# print("redundant", redundant)

			## Save the list of non-redundant elements in a file:
			nonredundants = open(output.nonreds, 'w')
			for ele in nonredundant:
				line = ele + "\n"
				nonredundants.write(line)

		else: # All files are empty so just write the empty output files to satisfy Snakemake
			for file in output[:-1]: # All but the dbchosen directory
				### Print matrix into a file
				ofile = open(file, 'w')
				ofile.close()
			os.mkdir(output.dbchosen)


def remove_overlap(ranges):
	""" Simplify a list of ranges; I got it from https://codereview.stackexchange.com/questions/21307/consolidate-list-of-ranges-that-overlap """
	result = []
	current_start = -1
	current_stop = -1 

	for start, stop in sorted(ranges):
		if start > current_stop:
			# this segment starts after the last segment stops
			# just add a new segment
			result.append( (start, stop) )
			current_start, current_stop = start, stop
		else:
			# current_start already guaranteed to be lower
			current_stop = max(current_stop, stop)
			# segments overlap, replace
			result[-1] = (current_start, current_stop) # SLAV: I modified this to update the stop too.
	return(result)

rule get_alignments:
	""" From the empty sites that were not specific to a single sample, extract the (non-redundant) locus in all samples as an aligned fasta file """
	input:
		nonredundant = "temp/{outprefix}.nonspecific_nonredundant.txt",
		genomes = expand("data/genomes/{sample}.fa", sample = samples),
		elements_nsp = expand("Nonspecific/{{outprefix}}_{sample}.both_confirmEmpty.final.filter.non-specific", sample = samples),
		elements_sp = expand("Insertions/{{outprefix}}_{sample}.both_confirmEmpty.final.filter.specific", sample = samples),
		databases = expand("dbBLAST/{sample}_db/{sample}_db.nhr", sample = samples)
	output:
		# directory("results/alignments")
		dummy = "results/alignments/dummies/{outprefix}.dummy",
	params:
		vicinity = VICINITY, # The length between (and including) flanks that we are willing to tolerate (i.e. maximum size haplotype)
		minhaplo =  0, # Minimum size haplotype
		buffer = BUFFER,
		evalue = "5e-10",
		threads = 4,
		perid = PERID,
	run:
		chosen_TEs = {} # A dictionary of all the elements that had confirmed empty sites, but not in all samples
		### ---  Specific ---
		specific = [] # List of all the insertions specific to one sample
		for file in input.elements_sp:
			if os.stat(file).st_size != 0: # File is not empty
				for line in open(file, 'r'):
					tab = line.rstrip("\n").split("\t")
					element = tab[4]
					chrom = tab[0]
					start = tab[1]
					end = tab[2]

					# Get the genome and element type
					tepattern = re.compile("([a-zA-Z0-9.]*)_([\w-]*)_([\w]*)_([0-9]*)") # eg. PcWa139m_Yeti_LTR__1259 or N11_maror5-28__7513
					# tepattern = re.compile("([a-zA-Z0-9.]*)_([a-zA-Z]*)_([\w]*)_([0-9]*)") # eg. PcWa139m_Yeti_LTR__1259
					matchy = tepattern.search(element)
					te_sample = matchy.group(1)
					typete = matchy.group(2)

					# Add to the big dictionary directly
					chosen_TEs[element] = (te_sample, typete, chrom, start, end)
					specific.append(element) # Keep track that it's specific
		# print(chosen_TEs) 
		### ---  Non-specific ---
		## Read file with the non-redundant elements
		nonreds = []
		for line in open(input.nonredundant, 'r'):
			nonreds.append(line.rstrip("\n"))

		for file in input.elements_nsp:
			if os.stat(file).st_size != 0: # File is not empty
				for line in open(file, 'r'):
					tab = line.rstrip("\n").split("\t")
					element = tab[4]
					chrom = tab[0]
					start = tab[1]
					end = tab[2]

					# Get the genome and element type
					tepattern = re.compile("([a-zA-Z0-9.]*)_([\w-]*)_([\w]*)_([0-9]*)") # eg. PcWa139m_Yeti_LTR__1259 or N11_maror5-28__7513
					# tepattern = re.compile("([a-zA-Z0-9.]*)_([a-zA-Z]*)_([\w]*)_([0-9]*)") # eg. PcWa139m_Yeti_LTR__1259
					matchy = tepattern.search(element)
					te_sample = matchy.group(1)
					typete = matchy.group(2)

					if element in nonreds: ## Filter to avoid redundant alignments
						chosen_TEs[element] = (te_sample, typete, chrom, start, end)

		# -------------------
		
		# Is the dictionary not empty?
		if chosen_TEs:
			# Because I'm not tracking the result of this rule, I have to make the folders myself
			shell("mkdir -p results/alignments/specific")
			shell("mkdir -p results/alignments/specific/complete")
			shell("mkdir -p results/alignments/specific/incomplete")
			shell("mkdir -p results/alignments/nonspecific/")
			shell("mkdir -p results/alignments/nonspecific/complete")
			shell("mkdir -p results/alignments/nonspecific/incomplete")

			# -----------------------------------
			# Get the fasta sequence of each element
			# -----------------------------------
			for element in chosen_TEs.keys():
				print(chosen_TEs[element])
				te_sample = chosen_TEs[element][0]
				typete = chosen_TEs[element][1]
				chrom = chosen_TEs[element][2]
				start = int(chosen_TEs[element][3]) # Base 1
				end = int(chosen_TEs[element][4])

				## Find the assembly where this element is
				# --------------------------------
				assembly = f"data/genomes/{te_sample}.fa"
				# --------------------------------

				# Get the right chromosome fasta sequence
				for seq_record in SeqIO.parse(open(assembly, 'r'), "fasta", generic_dna):
					if chrom == seq_record.id:
						thischrseq = seq_record
						break

				## Slice just the flanks, not the TE itself
				# The bed coordinates are base 1 because they came from a BLAST
				teseq = thischrseq[start - 1 - params.buffer:start] + thischrseq[end - 1:end + params.buffer]
				teseq.id = element
				teseq.description = ''

				# -----------------------------------
				# BLAST the element with its flanks to the sides in every sample
				# -----------------------------------
				slices = [] # Here I'll store the slices (flanks plus whatever is inside) of each sample
				samplecount = {sample:0 for sample in samples} # This is just for me to check that all samples got a hit

				for sample in samples: 
					# --------------------------------
					databasename = f"dbBLAST/{sample}_db/{sample}_db" # This was calculated previously
					outputhits = f"temp/{typete}_{sample}.hits.tab"
					assembly = f"data/genomes/{sample}.fa"
					# --------------------------------

					query_string = '>' + teseq.id + '\n' + str(teseq.seq)
					blast_command = NcbiblastnCommandline(cmd='blastn', out=outputhits, outfmt=6, db=databasename, evalue=params.evalue, num_threads=params.threads, task="blastn")
					stdout, stderr = blast_command(stdin=query_string)

					# -----------------------------------
					# Stitch the resulting BLAST: this is necessary because often the flank blast are divided by the element itself if the site is not empty
					# -----------------------------------
					unsorted_tabs = [line.rstrip("\n").split("\t") for line in open(outputhits, 'r')] 	
			
					stitchedtab = [] # The output

					## Make dic with subjects
					subject_dic = {}
					for tab in unsorted_tabs:
						subject = tab[1]
						if subject not in subject_dic.keys():
							subject_dic[subject] = [tab]
						else:
							subject_dic[subject].append(tab)
					
					# Sort the resulting list for that particular query, check every subject and sort it locally
					for sub in subject_dic.keys(): # For every contig
						currentsubject_sorted = sorted(subject_dic[sub], key = lambda x: int(x[9])) # Sort by the subject_start
						
						# ------------------------------------------------------
						# Stitch pieces together
						# ------------------------------------------------------
						# print("Subject:", sub)

						maxalign = 0
						queryStarts = []
						queryEnds = []
						subjectStarts = [] 
						subjectEnds = [] 

						for i in range(0, len(currentsubject_sorted)): # Find the main piece
							# print(i, len(currentsubject_sorted), currentsubject_sorted[i]) # For debugging
							# query_id, subject_id, percent_identity, alignment_length, N_mismatches, N_gaps, query_start, query_end, subject_start, subject_end, evalue, bit_score = currentsubject_sorted[i]
							query_start = int(currentsubject_sorted[i][6])
							query_end = int(currentsubject_sorted[i][7])	

							alignment_length = int(currentsubject_sorted[i][3])
							subject_start = int(currentsubject_sorted[i][8])
							subject_end = int(currentsubject_sorted[i][9])

							queryStarts.append(query_start)
							queryEnds.append(query_end)
							subjectStarts.append(subject_start)
							subjectEnds.append(subject_end)
				
							# Is this last piece one better than the previous ones?
							if alignment_length > maxalign:
								upper_hit = currentsubject_sorted[i]
								maxalign = alignment_length

							# There is only one clean hit
							if len(currentsubject_sorted) == 1: 
								stitchedtab.append(currentsubject_sorted[i])

							# The hit is broken
							elif i < len(currentsubject_sorted) - 1: 
								next_subject_start = int(currentsubject_sorted[i+1][8])

								if abs(subject_end - next_subject_start) > params.vicinity: # It's probably not part of the same hit # It was subject_start - next_subject_start before
									# print(sub, subject_end, next_subject_start, subject_end - next_subject_start)
									# So write down the previous one
									# -----------------
									# The new values for the subject
									# -----------------
									# query_id, subject_id, percent_identity, alignment_length, N_mismatches, N_gaps, query_start, query_end, subject_start, subject_end, evalue, bit_score
									new_tab = upper_hit[0:4]

									# These are our new values of the "unbroken" hit, but I'm not sure how to retrieve the no. of gaps and mistmatches
									new_tab.extend(['.', '.', min(queryStarts), max(queryEnds)])

									# subject_start > subject_end for the upper_hit
									if int(upper_hit[8]) > int(upper_hit[9]): # hit is reversed
										new_tab.extend([max(subjectStarts), min(subjectEnds)])
									else:
										new_tab.extend([min(subjectStarts), max(subjectEnds)])
									
									# Let's leave the e-val and Bit score the same as the upper hit
									new_tab.extend(upper_hit[10:])

									stitchedtab.append(new_tab) # Write it in the final output

									# -----------------
									# Reset for the next hit
									# -----------------
									queryStarts = []
									queryEnds = []
									subjectStarts = [] 
									subjectEnds = [] 

									maxalign = int(currentsubject_sorted[i+1][3])
									upper_hit = currentsubject_sorted[i+1]

							else:
								# The last hit in that subject
								# -----------------
								# The new values for the subject
								# -----------------
								new_tab = upper_hit[0:4]

								# These are our new values of the "unbroken" hit, but I'm not sure how to retrieve the no. of gaps and mistmatches
								new_tab.extend(['.', '.', min(queryStarts), max(queryEnds)])

								# subject_start > subject_end for the upper_hit
								if int(upper_hit[8]) > int(upper_hit[9]): # hit is reversed
									new_tab.extend([max(subjectStarts), min(subjectEnds)])
								else:
									new_tab.extend([min(subjectStarts), max(subjectEnds)])
								
								# Let's leave the e-val and Bit score the same as the upper hit
								new_tab.extend(upper_hit[10:])

								stitchedtab.append(new_tab) # Write it in the final output

					## Filter
					cleantabs = []
					for tab in stitchedtab:
						print(tab)
						periden = float(tab[2])
						hitlen = float(tab[3])

						if (periden >= params.perid) and (hitlen >= params.buffer*0.75):
							cleantabs.append(tab)				

					# Remove the BLAST file
					os.remove(outputhits)

					# -----------------------------------
					# Get haplotype
					# -----------------------------------
					chunks = {}
					# print(cleantabs)
					for hit5 in cleantabs:
						for hit3 in cleantabs:
							# Is it the same contig
							if hit5[1] == hit3[1]:
								hit5start = int(hit5[8])
								hit5end = int(hit5[9])
								hit3start = int(hit3[8])
								hit3end = int(hit3[9])
								maxcoord = max([hit5start, hit5end, hit3start, hit3end])
								mincoord = min([hit5start, hit5end, hit3start, hit3end])
								
								if ((maxcoord - mincoord) <= params.vicinity) and ((maxcoord - mincoord) >= params.minhaplo):
									# Add it to dictionary and keep the edges
									if hit5[1] not in chunks.keys():
										chunks[hit5[1]] = [(mincoord, maxcoord)]
									else:
										chunks[hit5[1]].append((mincoord, maxcoord))
					
					# Read the assembly of this sample into memory:
					records_dict = SeqIO.to_dict(SeqIO.parse(open(assembly, 'r'), "fasta", generic_dna))

					for ctg in chunks.keys(): # For each contig hit
						hitseq = records_dict[ctg] 	# The actual sequence hit by the query

						for start,end in remove_overlap(chunks[ctg]): # Reduce the list to only the non-overlapping ranges
							if (end - start) >= params.buffer:
								# Slice it 
								slice = hitseq[start - 1:end]
								slice.id = hitseq.id + "_" + str(start + 1) + "-" + str(end)
								slice.description = ''

								# Save it in the list
								slices.append(slice)
								samplecount[sample] += 1
				
				## Print all the sequences recovered for this locus into a fasta for this element
				telocusraw = f"results/alignments/{element}_raw.fa"
				with open(telocusraw, "w") as handle:
					SeqIO.write(slices, telocusraw, "fasta") 

				# -----------------------------------
				## Align with MAFFT
				# -----------------------------------
				print(f"Aligning sequences of locus containing {element} ...")

				mafft_cline = MafftCommandline(input=telocusraw, adjustdirection = True, thread = params.threads)
				stdout, stderr = mafft_cline()


				# Change the name depending on the content
				allone = True
				for key in samplecount.keys():
					if samplecount[key] != 1:
						allone = False

				if allone:
					if element in nonreds:
						telocusalign = f"results/alignments/nonspecific/complete/{element}.fa"
					elif element in specific:
						telocusalign = f"results/alignments/specific/complete/{element}.fa"
				else:
					if element in nonreds:
						telocusalign = f"results/alignments/nonspecific/incomplete/{element}.fa"
					elif element in specific:
						telocusalign = f"results/alignments/specific/incomplete/{element}.fa"

				# Save the file
				with open(telocusalign, "w") as handle:
					handle.write(stdout)

				# Remove the unaligned file 
				os.remove(telocusraw)

		# -----------------------------------
		## Print a dummy to get the folder created
		# -----------------------------------
		ofile = open(output.dummy, 'w')
		ofile.write(f"Number of sample-specific insertions evaluated: {len(specific)}\n")
		ofile.write(f"Number of potential shared loci evaluated: {len(nonreds)}\n")

## --- Reporting ---

rule summary_insertions:
	""" Number of confirmed insertions unique to a given sample """
	input: 
		expand("Insertions/{outprefix}_{sample}.both_confirmEmpty.final.filter.specific", sample = samples, outprefix = outprefix)
	output:
		matrix = "results/summary_insertions.tbl"
	run:
		# Make an empty matrix 
		sppmatrix = np.zeros( (len(samples), len(outprefix)), "int" )

		for sample in samples:
			for prefix in outprefix: # Get how many insertions of that element in this sample

				# Find file
				for file in input:
					if file.find(f"{prefix}_{sample}") > 0: # Index if found and -1 otherwise.
						focalfile = file
						break

				# How many elements were found?
				num_lines = sum(1 for line in open(focalfile))

				indexsample = samples.index(sample)
				indexprefix = outprefix.index(prefix)

				## Add values into the matrix
				sppmatrix[indexsample][indexprefix] = num_lines # [rows][columns]

		### Print matrix into a file
		ofile = open(output.matrix, 'w')
		# Print header
		header = 'Genome_TE\t' + '\t'.join(outprefix) + '\n'
		ofile.write(header)

		for i in range(len(samples)):
			line = f'{samples[i]}\t' + '\t'.join([str(n) for n in sppmatrix[i]]) + '\n'
			ofile.write(line)
			
## The original table in Vale's pipeline, originally called "summary_specific_insertions.tbl"
# I changed the name to avoid confusion with my concept of specific insertions which represent unique insertions compared to all samples (not pairwise)
rule summary_matrix_emptysites:
	""" Make summary of findings """
	input:
		finals = expand("Sites/{target}-vs-{query}/{{outprefix}}_{target}-vs-{query}.both_confirmEmpty.final.filter", zip, target = targets, query = queries),
	output:
		matrix = "results/emptysites/{outprefix}_matrix_emptysites.txt"
	run:
		# Make an empty matrix 
		sppmatrix = np.zeros( (len(samples), len(samples)), "int" )

		### Count how  many empty sites are in the targets
		
		for combi in combis:
			target = combi[0]
			query = combi[1]

			indextarget = samples.index(target)
			indexquery = samples.index(query)

			# How many empty sites were found?
			for file in input.finals:
				if file.find(f"{target}-vs-{query}.both") > 0: # Index if found and -1 otherwise.
					focalfile = file
					break
			
			num_lines = sum(1 for line in open(focalfile))

			## Add values into the matrix
			sppmatrix[indextarget][indexquery] = num_lines # [rows][columns]
			
		### Print matrix into a file
		ofile = open(output.matrix, 'w')
		# Print header
		header = 'Genome\t' + '\t'.join(samples) + '\n'
		ofile.write(header)

		for i in range(len(samples)):
			line = f'{samples[i]}\t' + '\t'.join([str(n) for n in sppmatrix[i]]) + '\n'
			ofile.write(line)

## Example:

    # $ cat results/example_matrix_emptysites.txt
    # Genome  Podan2  CBS237.71m  PcWa139m
    # Podan2  0   11  2
    # CBS237.71m  20  0   2
    # PcWa139m    19  12  0

# The target is the rows, and the query is the columns. So there are 20 empty
# sites in CBS237.71m compared to the genome of Podan2.

