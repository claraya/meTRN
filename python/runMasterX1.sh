#!/usr/bin/sh


# Import worm data:
#python dataImporter.py --path ~/meTRN --mode import.peaks --source extras/brunetlab --peaks ce_brunetlab --rank ON --method OFF --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6 --cutChr ON --idrSource OFF --reformat ON


# Construct peak set for Shuo Han and Anne Brunnet (with the ASH-2 dataset):
#mkdir ~/meTRN/data/peaks/ce_brunetlab_com_cx_raw/
#cp ~/meTRN/data/peaks/ce_selection_com_cx_raw/* ~/meTRN/data/peaks/ce_brunetlab_com_cx_raw/
#cp ~/meTRN/idr/peaks/ce_brunetlab/* ~/meTRN/data/peaks/ce_brunetlab_com_cx_raw/


# Construct context-specific sets:
#mkdir ~/meTRN/data/peaks/ce_brunetlab_com_ex_raw/
#mkdir ~/meTRN/data/peaks/ce_brunetlab_com_l1_raw/
#mkdir ~/meTRN/data/peaks/ce_brunetlab_com_l2_raw/
#mkdir ~/meTRN/data/peaks/ce_brunetlab_com_l3_raw/
#mkdir ~/meTRN/data/peaks/ce_brunetlab_com_l4_raw/

#cp ~/meTRN/data/peaks/ce_brunetlab_com_cx_raw/* ~/meTRN/data/peaks/ce_brunetlab_com_ex_raw/
#mv ~/meTRN/data/peaks/ce_brunetlab_com_ex_raw/*_L1_* ~/meTRN/data/peaks/ce_brunetlab_com_l1_raw/
#mv ~/meTRN/data/peaks/ce_brunetlab_com_ex_raw/*_L2_* ~/meTRN/data/peaks/ce_brunetlab_com_l2_raw/
#mv ~/meTRN/data/peaks/ce_brunetlab_com_ex_raw/*_L3_* ~/meTRN/data/peaks/ce_brunetlab_com_l3_raw/
#mv ~/meTRN/data/peaks/ce_brunetlab_com_ex_raw/*_L4_* ~/meTRN/data/peaks/ce_brunetlab_com_l4_raw/


# Generate HOT-filtered peak sets (occupancy P01 cutoff):
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_brunetlab_com_ex_xot --infile maphot_ce_selection_reg_ex_occP01_hot.bed --source ~/meTRN/data/peaks/ce_brunetlab_com_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_brunetlab_com_l1_xot --infile maphot_ce_selection_reg_l1_occP01_hot.bed --source ~/meTRN/data/peaks/ce_brunetlab_com_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_brunetlab_com_l2_xot --infile maphot_ce_selection_reg_l2_occP01_hot.bed --source ~/meTRN/data/peaks/ce_brunetlab_com_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_brunetlab_com_l3_xot --infile maphot_ce_selection_reg_l3_occP01_hot.bed --source ~/meTRN/data/peaks/ce_brunetlab_com_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_brunetlab_com_l4_xot --infile maphot_ce_selection_reg_l4_occP01_hot.bed --source ~/meTRN/data/peaks/ce_brunetlab_com_l4_raw

#rm -rf ~/meTRN/data/peaks/ce_brunetlab_com_cx_xot/
#mkdir ~/meTRN/data/peaks/ce_brunetlab_com_cx_xot/
#cp ~/meTRN/data/peaks/ce_brunetlab_com_ex_xot/* ~/meTRN/data/peaks/ce_brunetlab_com_cx_xot/
#cp ~/meTRN/data/peaks/ce_brunetlab_com_l1_xot/* ~/meTRN/data/peaks/ce_brunetlab_com_cx_xot/
#cp ~/meTRN/data/peaks/ce_brunetlab_com_l2_xot/* ~/meTRN/data/peaks/ce_brunetlab_com_cx_xot/
#cp ~/meTRN/data/peaks/ce_brunetlab_com_l3_xot/* ~/meTRN/data/peaks/ce_brunetlab_com_cx_xot/
#cp ~/meTRN/data/peaks/ce_brunetlab_com_l4_xot/* ~/meTRN/data/peaks/ce_brunetlab_com_cx_xot/


# Generate complete, collapsed, density and report files:
#python mapPeaks.py --path ~/meTRN --mode build --overwrite OFF


# Launch co-association analyses with IntervalStats:
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_brunetlab_com_cx_xot --source annotations --domain in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub scg3_intervalstats.sh --server ON --job ceBrunet


# Crunch species-specific co-association analyses with IntervalStats:
#python mapCAs.py --path ~/meTRN --mode crunch --peaks ce_brunetlab_com_cx_xot --source annotations --domain in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions


# Report species-specific co-association analyses with IntervalStats:
#python mapCAs.py --path ~/meTRN --mode report --peaks ce_brunetlab_com_cx_xot --source annotations --domain in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions

#top