#!/usr/bin/python -u
'''
### 00b.invite_Script_03_in_a_loop.py
- usage: replace
	- GenomeDir [modified] genomes
	- .GTF file's location
	- "RefSeq" tag in line 22
'''

import os, glob
from subprocess import call

GTF_fnp = "/hpc/hub_oudenaarden/gene_models/human_gene_models/hg19_RefSeq_genes_clean.fa"; # see: https://github.com/bow/x_reactivation/issues/190#issuecomment-149318797

GenomeDir = '/hpc/hub_oudenaarden/Abel/genomes/X_React_genomes/hg19_mod/pat'
# GenomeDir = '/hpc/hub_oudenaarden/Abel/genomes/X_React_genomes/hg19_mod/mat'
os.chdir(GenomeDir)

Genomes = glob.glob("*.fa"); print Genomes
for Genome_name in Genomes:
	# print Genome_name
	Out_Ref_name = (Genome_name.split('reformat.fa')[0])+'RefSeq.fa'
	print Out_Ref_name

	bash_command = '/Users/abelvertesy/x_reactivation/analysis/DevCell_analysis/06.5.Reference_generation/03.gtf2fa.pl -in='+GTF_fnp+' -ref='+Genome_name+' > '+Out_Ref_name
	# bash_command = "echo 111"
	print bash_command
	# call(bash_command, shell=True)

print "Sweet Ohio"
