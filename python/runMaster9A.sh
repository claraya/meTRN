#!/usr/bin/sh


# Generate a version of the TSSs without MtDNA annotations:
#grep -v -e "MtDNA" ~/meTRN/data/annotations/in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.bed > ~/meTRN/data/annotations/in2shape_ce_wormbased_TSS_gx_nomtdna_up2000_dn200.bed


# Build transcription regulatory networks from TIP predictions:
#python mapNetwork.py --path ~/meTRN --organism ce --mode build --source extras --infile tip/cel_TIP_target_Q.10.txt --peaks ce_selection_com_cx_raw --name tip --target factor.context --tip ON --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6 
#python mapNetwork.py --path ~/meTRN --organism hs --mode build --source extras --infile tip/hsa_TIP_target_Q.10.txt --peaks hs_selection_com_cx_raw --name tip --target factor.context --tip ON --mapping in2shape_hs_modencode_DAT_xx.txt


# Build transcription regulatory networks from promoter overlap:
#python mapNetwork.py --path ~/meTRN --organism ce --mode build --source annotations --infile in2shape_ce_wormbased_TSS_gx_nomtdna_up2000_dn200.bed --peaks ce_selection_com_cx_raw --name 2kb --target factor.context --tip OFF --overwrite OFF --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6
#python mapNetwork.py --path ~/meTRN --organism dm --mode build --source annotations --infile in2shape_dm_celnikers_TSS_gx_slopbed_up2000_dn200.bed --peaks dm_selection_com_cx_raw --name 2kb --target factor.context --tip OFF --overwrite ON
#python mapNetwork.py --path ~/meTRN --organism hs --mode build --source annotations --infile in2shape_hs_gencode10_TSS_gx_slopbed_up5000_dn500.bed --peaks hs_selection_com_cx_raw --name 5kb --target factor.context --tip OFF --overwrite ON


# Build stage-specific transcription regulatory networks from promoter overlap:
#python mapNetwork.py --path ~/meTRN --organism ce --mode build --source annotations --infile in2shape_ce_wormbased_TSS_gx_nomtdna_up2000_dn200.bed --peaks ce_selection_com_ex_raw --name 2kb --target factor.context --tip OFF --overwrite OFF --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6
#python mapNetwork.py --path ~/meTRN --organism ce --mode build --source annotations --infile in2shape_ce_wormbased_TSS_gx_nomtdna_up2000_dn200.bed --peaks ce_selection_com_l1_raw --name 2kb --target factor.context --tip OFF --overwrite OFF --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6
#python mapNetwork.py --path ~/meTRN --organism ce --mode build --source annotations --infile in2shape_ce_wormbased_TSS_gx_nomtdna_up2000_dn200.bed --peaks ce_selection_com_l2_raw --name 2kb --target factor.context --tip OFF --overwrite OFF --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6
#python mapNetwork.py --path ~/meTRN --organism ce --mode build --source annotations --infile in2shape_ce_wormbased_TSS_gx_nomtdna_up2000_dn200.bed --peaks ce_selection_com_l3_raw --name 2kb --target factor.context --tip OFF --overwrite OFF --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6
#python mapNetwork.py --path ~/meTRN --organism ce --mode build --source annotations --infile in2shape_ce_wormbased_TSS_gx_nomtdna_up2000_dn200.bed --peaks ce_selection_com_l4_raw --name 2kb --target factor.context --tip OFF --overwrite OFF --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6


# Resolve transcription regulatory networks (complete) to protein IDs:
#python mapNetwork.py --path ~/meTRN --organism ce --mode resolve --peaks ce_selection_com_cx_raw --name 2kb --target factor.context
#python mapNetwork.py --path ~/meTRN --organism dm --mode resolve --peaks dm_selection_com_cx_raw --name 2kb --target factor.context
#python mapNetwork.py --path ~/meTRN --organism hs --mode resolve --peaks hs_selection_com_cx_raw --name 5kb --target factor.context --strip ON


# Resolve transcription regulatory networks (contexts) to protein IDs:
#python mapNetwork.py --path ~/meTRN --organism ce --mode resolve --peaks ce_selection_com_ex_raw --name 2kb --target factor.context
#python mapNetwork.py --path ~/meTRN --organism ce --mode resolve --peaks ce_selection_com_l1_raw --name 2kb --target factor.context
#python mapNetwork.py --path ~/meTRN --organism ce --mode resolve --peaks ce_selection_com_l2_raw --name 2kb --target factor.context
#python mapNetwork.py --path ~/meTRN --organism ce --mode resolve --peaks ce_selection_com_l3_raw --name 2kb --target factor.context
#python mapNetwork.py --path ~/meTRN --organism ce --mode resolve --peaks ce_selection_com_l4_raw --name 2kb --target factor.context


# Scan embryo network for cascades of expression:
#python mapNetwork.py --path ~/meTRN --organism ce --mode cascade --peaks ce_selection_com_ex_raw --name 2kb --target factor.context --mapping cetrn_inherit_f0.1 --fraction 0.9


# Compare human and worm networks!
#python mapNetwork.py --path ~/meTRN --organism hs --mode hybrid --infile modencode.merged.system.orth.txt --species hs,ce --a hs_selection_com_cx_raw/5kb --b ce_selection_com_cx_raw/2kb --orthology family





# Scan network for related regulation interactions:
#python mapNetwork.py --path ~/meTRN --organism ce --mode commons --peaks ce_selection_com_cx_raw --name 2kb --target factor.context --mapping worm_interaction_wi2007.txt --a IDA --b IDB --cutoff 1
#python mapNetwork.py --path ~/meTRN --organism ce --mode commons --peaks ce_selection_com_cx_raw --name 2kb --target factor.context --mapping worm_interaction_wi2007.txt --a IDA --b IDB --cutoff 4



#top
#bash runMaster9A.sh