# scseq

This repository contains scripts and tools for single cell RNAseq, as
initially developed by Dominic Grün and Alexander van Oudenaarden at the
[Hubrecht Institute](www.hubrecht.eu/), Utrecht, The Netherlands.

**This code has moved to https://bitbucket.org/princessmaximacenter/scseq/ .**
[Contact me](mailto:p.lijnzaad@gmailcom) if you want access.

It contains (often heavily modified) code from other people,
specifically from Dominic Grün (now at MPI Freiburg), Lennart Kester and
Abel Vertesy. For related code,
see [Abel Vertesy's repository](https://github.com/vertesy/TheCorvinas),
in particular the Python/MapAndGo and Python/Concatenator subdirectories

# Contents

## General stuff
 * env.sh.example - sets PATH en PERLLIB so they don't have to be hard-coded
 * tools.pm - functions used by the perl scripts

## barcodes and external controls

 * holstegelab_EC.fa - kept for historic reasons
 * ERCC92.fa - superset of holstegelab_EC.fa, well-known set of external spike-in controls
 * cel-seq96_barcodes.csv - initial set of barcodes
 * celseq2_bc384.csv - barcodes used for celseq2, 384-well protocols

## preparing the reference transcriptome

 * mask_polyA.pl - masks any /A{10,}/ in reference transcriptome
 * ucsc2gtf.pl - converts UCSC RefSeq.txt table to a gtf file
 * polyXY.pl - finds stretches of e.g. AAAAAAAAGGGGGGGG etc.
 * merge_isoforms_gtf.pl - creates 'supertranscripts', i.e. (virtual) transcripts consisting of all possible exons

## mapping and postprocessing

 * do_mappings_strand_cs2v2.pl - driver script for preprocssing, mapping and bookkeeping. See also Abel Vertesy's Python/MapAndGo/MapAndGo.py in https://github.com/vertesy/TheCorvinas
 * add_bc_to_R2.pl - preprocesses the reads prior to mapping. 
 * process_sam_cel384v2.pl - bookkeeping and statistics per cell. Needs code from https://github.com/plijnzaad/demultiplex for recovering mismatched cell barcodes
