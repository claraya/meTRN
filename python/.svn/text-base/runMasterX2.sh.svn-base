#!/usr/bin/sh


#"fastaFromBed -fi " + genomefile + " -bed " + coordfile + " -fo " + fastafile + " -name"
#python memeCucu.py --path ~/meTRN --organism ce --mode measures --infile in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.bed --parameters pkheradpour/ce_frequency.output --exclude MtDNA


# Generate alternative human ATF3 and MTA3 motifs without intervening base:
#From xx_combos_motif.txt, modified "verts_ATF3_disc_3" into "verts_ATF3_disc_X" in xx_adjust_motif.txt...
#From xx_combos_motif.txt, modified "verts_MTA3_disc_1" into "verts_MTA3_disc_X" in xx_adjust_motif.txt...


# Convert Pouya's motif files into a MEME format file...
#python memeCucu.py --path ~/meTRN --organism ce --mode database --infile pkheradpour/xx_combos_motif.txt --max 200 --parameters "MEME version 4;;ALPHABET= ACGT;;strands: + -;;Background letter frequencies (from);A 0.322165798506 C 0.17769216647 G 0.177616880599 T 0.322525154425;" --background pkheradpour/ce_frequencies.txt
#python memeCucu.py --path ~/meTRN --organism ce --mode database --infile pkheradpour/xx_adjust_motif.txt --max 200 --parameters "MEME version 4;;ALPHABET= ACGT;;strands: + -;;Background letter frequencies (from);A 0.322165798506 C 0.17769216647 G 0.177616880599 T 0.322525154425;" --background pkheradpour/ce_frequencies.txt


# Search for motifs similar to worm motifs:
#python memeCucu.py --path ~/meTRN --organism ce --mode pairwise --peaks ce_selection_reg_cx_hon --name metrn --max 200 --window OFF --infile pkheradpour/xx_combos_motif.meme --parameters worms


# Search for motifs similar to worm motifs (using combined "combos" from P. Kheradpour and motifs from Jolma et al. 2013):
#python memeCucu.py --path ~/meTRN --organism ce --mode pairwise --peaks ce_selection_reg_cx_hon --name jumbo --max 200 --window OFF --infile pkheradpour/xx_jumbos_motif.meme --parameters worms


# Graph all motifs in the combined "combos" from P. Kheradpour and motifs from Jolma et al. 2013):
#python memeCucu.py --path ~/meTRN --organism ce --mode graphing --peaks ce_selection_reg_cx_hon --name jumbo --max 200 --window OFF --infile pkheradpour/xx_jumbos_motif.meme


# Scan binding sites (peaks) for motifs of same factor:
#python memeCucu.py --path ~/meTRN --organism ce --mode scanning --peaks ce_selection_reg_cx_hon --name metrn --max 200 --window OFF --infile pkheradpour/xx_combos_motif.meme --parameters worms --background pkheradpour/ce_frequencies.txt
#python memeCucu.py --path ~/meTRN --organism ce --mode scanning --peaks ce_selection_reg_cx_hon --name metrn --max 400 --window OFF --infile pkheradpour/xx_combos_motif.meme --parameters worms --background pkheradpour/ce_frequencies.txt
#python memeCucu.py --path ~/meTRN --organism ce --mode scanning --peaks ce_selection_reg_cx_hon --name metrn --max 600 --window OFF --infile pkheradpour/xx_combos_motif.meme --parameters worms --background pkheradpour/ce_frequencies.txt
#python memeCucu.py --path ~/meTRN --organism ce --mode scanning --peaks ce_selection_reg_cx_hon --name metrn --max 800 --window OFF --infile pkheradpour/xx_combos_motif.meme --parameters worms --background pkheradpour/ce_frequencies.txt
#python memeCucu.py --path ~/meTRN --organism ce --mode scanning --peaks ce_selection_reg_cx_hon --name metrn --max 1000 --window OFF --infile pkheradpour/xx_combos_motif.meme --parameters worms --background pkheradpour/ce_frequencies.txt


# Search for HOX/PBC motifs in LIN-39 binding sites:
#python memeCucu.py --path ~/meTRN --organism ce --mode scanning --peaks ce_selection_reg_cx_hon --name metrn --max 1000 --window OFF --infile pkheradpour/xx_combos_motif.meme --parameters verts_Hoxa4 --target LIN-39 --background pkheradpour/ce_frequencies.txt
#python memeCucu.py --path ~/meTRN --organism ce --mode scanning --peaks ce_selection_reg_cx_hon --name metrn --max 1000 --window OFF --infile pkheradpour/xx_combos_motif.meme --parameters verts_Hoxa4 --target MAB-5 --background pkheradpour/ce_frequencies.txt


# Search Jolma 2013 motifs in LIN-39 binding sites:
#python memeCucu.py --path ~/meTRN --organism ce --mode scanning --peaks ce_selection_reg_cx_hon --name jolma --max 1000 --window OFF --infile motifs/jolma2013.meme --target LIN-39 --background pkheradpour/ce_frequencies.txt
#python memeCucu.py --path ~/meTRN --organism ce --mode scanning --peaks ce_selection_reg_cx_hon --name jolma --max 1000 --window OFF --infile motifs/jolma2013.meme --target MAB-5 --background pkheradpour/ce_frequencies.txt


# Search Jolma 2013 motifs in exhaustively binding sites:
#python memeCucu.py --path ~/meTRN --organism ce --mode scanning --peaks ce_selection_reg_cx_hon --name jolma --max 200 --window OFF --infile motifs/jolma2013.meme --target ALL --background pkheradpour/ce_frequencies.txt
#python memeCucu.py --path ~/meTRN --organism ce --mode ortholog --peaks ce_selection_reg_cx_hon --name jolma --max 200 --window OFF --infile modencode.family.common.matrix.txt --target ALL --parameters hs --exclude NHR,ZNF


# Build FASTA sequences from HOT-filtered peaks (per context):
#python memeCucu.py --path ~/meTRN --organism hs --mode sequence --peaks hs_standards_reg_cx_hct --max 250 --window 100 --repeatMask OFF
#python memeCucu.py --path ~/meTRN --organism ce --mode sequence --peaks ce_selection_reg_cx_hct --max 250 --window 100 --repeatMask OFF
#python memeCucu.py --path ~/meTRN --organism dm --mode sequence --peaks dm_selection_reg_cx_hct --max 250 --window 100 --repeatMask OFF


# Attempt MEME-based motif discovery:
#python memeCucu.py --path ~/meTRN --organism hs --mode discover --peaks hs_standards_reg_cx_hct --max 250 --window 100 --repeatMask OFF --parameters "-maxsize 10000000 -nmotifs 3 -minsites 50 -maxw 12"
#python memeCucu.py --path ~/meTRN --organism ce --mode discover --peaks ce_selection_reg_cx_hct --max 250 --window 100 --repeatMask OFF --parameters "-maxsize 10000000 -nmotifs 3 -minsites 50 -maxw 12"
#python memeCucu.py --path ~/meTRN --organism dm --mode discover --peaks dm_selection_reg_cx_hct --max 250 --window 100 --repeatMask OFF --parameters "-maxsize 10000000 -nmotifs 3 -minsites 50 -maxw 12"





#python memeCucu.py --path ~/meTRN --organism ce --mode sequence --peaks ce_selection_reg_cx_hon --max 250 --window 100 --repeatMask OFF
#python memeCucu.py --path ~/meTRN --organism ce --mode discover --peaks ce_selection_reg_cx_hon --max 250 --window 100 --repeatMask OFF --parameters "-maxsize 10000000 -nmotifs 3 -minsites 50 -maxw 12"
#python memeCucu.py --path ~/meTRN --organism ce --mode fraction --peaks ce_selection_reg_cx_hon --max 250 --window 100 --repeatMask OFF --sequence GAAAG

#python memeCucu.py --path ~/meTRN --organism ce --mode sequence --peaks ce_selection_reg_cx_hon --max 250 --window 100 --repeatMask OFF
#python memeCucu.py --path ~/meTRN --organism ce --mode fraction --peaks ce_selection_reg_cx_hon --max 250 --window 100 --repeatMask OFF --sequence GAAAG
#python memeCucu.py --path ~/meTRN --organism ce --mode fraction --peaks ce_selection_reg_cx_hon --max 250 --window 100 --repeatMask OFF --sequence TAATT
#python memeCucu.py --path ~/meTRN --organism ce --mode fraction --peaks ce_selection_reg_cx_hon --max 250 --window 100 --repeatMask OFF --sequence ATTGAT
#python memeCucu.py --path ~/meTRN --organism ce --mode fraction --peaks ce_selection_reg_cx_hon --max 250 --window 100 --repeatMask OFF --sequence ATCGAT



#top
#bash runMasterX2.sh