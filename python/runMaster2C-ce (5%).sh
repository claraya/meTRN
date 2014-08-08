#!/usr/bin/sh


# Map HOT regions (fail) and non-HOT regions (pass); without any actual filtering cutoffs:
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_cx --name basics
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_ex --name basics --tag EX-
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l1 --name basics --tag L1-
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l2 --name basics --tag L2-
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l3 --name basics --tag L3-
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l4 --name basics --tag L4-

	
# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (occupancy 0.05 significance):
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_ex --name occP05 --significance 05 --tag EX- --metric occupancy --cutoff 7
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l1 --name occP05 --significance 05 --tag L1- --metric occupancy --cutoff 12
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l2 --name occP05 --significance 05 --tag L2- --metric occupancy --cutoff 12
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l3 --name occP05 --significance 05 --tag L3- --metric occupancy --cutoff 8
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l4 --name occP05 --significance 05 --tag L4- --metric occupancy --cutoff 11


# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (density 0.05 significance):
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_ex --name denP05 --significance 05 --tag EX- --metric density --cutoff 8.2
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l1 --name denP05 --significance 05 --tag L1- --metric density --cutoff 11.6
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l2 --name denP05 --significance 05 --tag L2- --metric density --cutoff 9.8
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l3 --name denP05 --significance 05 --tag L3- --metric density --cutoff 8.9
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l4 --name denP05 --significance 05 --tag L4- --metric density --cutoff 9.3


# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (combined 0.05 significance):
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_ex --name simP05 --significance 05 --tag EX- --metric combined --cutoff 7,8.2
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l1 --name simP05 --significance 05 --tag L1- --metric combined --cutoff 12,11.6
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l2 --name simP05 --significance 05 --tag L2- --metric combined --cutoff 12,9.8
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l3 --name simP05 --significance 05 --tag L3- --metric combined --cutoff 8,8.9
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l4 --name simP05 --significance 05 --tag L4- --metric combined --cutoff 11,9.3


# Calculate overlap in HOT regions and filter regions against the ubiquitously HOT regions:
#python mapHOT.py --path ~/meTRN --organism ce --mode overlap --name occP05 --target EX:ce_selection_reg_ex,L1:ce_selection_reg_l1,L2:ce_selection_reg_l2,L3:ce_selection_reg_l3,L4:ce_selection_reg_l4 --overlap ce_selection_reg_cx
#python mapHOT.py --path ~/meTRN --organism ce --mode overlap --name denP05 --target EX:ce_selection_reg_ex,L1:ce_selection_reg_l1,L2:ce_selection_reg_l2,L3:ce_selection_reg_l3,L4:ce_selection_reg_l4 --overlap ce_selection_reg_cx
#python mapHOT.py --path ~/meTRN --organism ce --mode overlap --name simP05 --target EX:ce_selection_reg_ex,L1:ce_selection_reg_l1,L2:ce_selection_reg_l2,L3:ce_selection_reg_l3,L4:ce_selection_reg_l4 --overlap ce_selection_reg_cx


# Store HOT and RGB region files for each context and the global context file:
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_ex_occP05 --regions ce_selection_reg_ex_occP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l1_occP05 --regions ce_selection_reg_l1_occP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l2_occP05 --regions ce_selection_reg_l2_occP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l3_occP05 --regions ce_selection_reg_l3_occP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l4_occP05 --regions ce_selection_reg_l4_occP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_cx_occP05 --regions ce_selection_reg_cx_occP05 --source overlap

#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_ex_denP05 --regions ce_selection_reg_ex_denP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l1_denP05 --regions ce_selection_reg_l1_denP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l2_denP05 --regions ce_selection_reg_l2_denP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l3_denP05 --regions ce_selection_reg_l3_denP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l4_denP05 --regions ce_selection_reg_l4_denP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_cx_denP05 --regions ce_selection_reg_cx_denP05 --source overlap

#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_ex_simP05 --regions ce_selection_reg_ex_simP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l1_simP05 --regions ce_selection_reg_l1_simP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l2_simP05 --regions ce_selection_reg_l2_simP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l3_simP05 --regions ce_selection_reg_l3_simP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l4_simP05 --regions ce_selection_reg_l4_simP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_cx_simP05 --regions ce_selection_reg_cx_simP05 --source overlap


# Generate HOT-filtered peak sets (occupancy P05 cutoff):
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_cx_hon --infile maphot_ce_selection_reg_cx_occP05_any.bed --source ~/meTRN/data/peaks/ce_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_ex_hot --infile maphot_ce_selection_reg_ex_occP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l1_hot --infile maphot_ce_selection_reg_l1_occP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l2_hot --infile maphot_ce_selection_reg_l2_occP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l3_hot --infile maphot_ce_selection_reg_l3_occP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l4_hot --infile maphot_ce_selection_reg_l4_occP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_cx_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_ex_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l1_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l2_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l3_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l4_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_cx_hon --infile maphot_ce_selection_reg_cx_occP05_any.bed --source ~/meTRN/data/peaks/ce_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_ex_hot --infile maphot_ce_selection_reg_ex_occP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l1_hot --infile maphot_ce_selection_reg_l1_occP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l2_hot --infile maphot_ce_selection_reg_l2_occP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l3_hot --infile maphot_ce_selection_reg_l3_occP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l4_hot --infile maphot_ce_selection_reg_l4_occP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_cx_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_ex_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l1_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l2_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l3_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l4_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l4_raw

#rm -rf ~/meTRN/data/peaks/ce_selection_reg_cx_hot/
#mkdir ~/meTRN/data/peaks/ce_selection_reg_cx_hot/
#cp ~/meTRN/data/peaks/ce_selection_reg_ex_hot/* ~/meTRN/data/peaks/ce_selection_reg_cx_hot/
#cp ~/meTRN/data/peaks/ce_selection_reg_l1_hot/* ~/meTRN/data/peaks/ce_selection_reg_cx_hot/
#cp ~/meTRN/data/peaks/ce_selection_reg_l2_hot/* ~/meTRN/data/peaks/ce_selection_reg_cx_hot/
#cp ~/meTRN/data/peaks/ce_selection_reg_l3_hot/* ~/meTRN/data/peaks/ce_selection_reg_cx_hot/
#cp ~/meTRN/data/peaks/ce_selection_reg_l4_hot/* ~/meTRN/data/peaks/ce_selection_reg_cx_hot/

#rm -rf ~/meTRN/data/peaks/ce_selection_com_cx_hot/
#mkdir ~/meTRN/data/peaks/ce_selection_com_cx_hot/
#cp ~/meTRN/data/peaks/ce_selection_com_ex_hot/* ~/meTRN/data/peaks/ce_selection_com_cx_hot/
#cp ~/meTRN/data/peaks/ce_selection_com_l1_hot/* ~/meTRN/data/peaks/ce_selection_com_cx_hot/
#cp ~/meTRN/data/peaks/ce_selection_com_l2_hot/* ~/meTRN/data/peaks/ce_selection_com_cx_hot/
#cp ~/meTRN/data/peaks/ce_selection_com_l3_hot/* ~/meTRN/data/peaks/ce_selection_com_cx_hot/
#cp ~/meTRN/data/peaks/ce_selection_com_l4_hot/* ~/meTRN/data/peaks/ce_selection_com_cx_hot/


# Generate HOT-filtered peak sets (density P05 cutoff):
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_cx_hdn --infile maphot_ce_selection_reg_cx_denP05_any.bed --source ~/meTRN/data/peaks/ce_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_ex_hdt --infile maphot_ce_selection_reg_ex_denP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l1_hdt --infile maphot_ce_selection_reg_l1_denP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l2_hdt --infile maphot_ce_selection_reg_l2_denP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l3_hdt --infile maphot_ce_selection_reg_l3_denP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l4_hdt --infile maphot_ce_selection_reg_l4_denP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_cx_hdb --infile maphot_ce_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_ex_hdb --infile maphot_ce_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l1_hdb --infile maphot_ce_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l2_hdb --infile maphot_ce_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l3_hdb --infile maphot_ce_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l4_hdb --infile maphot_ce_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_cx_hdn --infile maphot_ce_selection_reg_cx_denP05_any.bed --source ~/meTRN/data/peaks/ce_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_ex_hdt --infile maphot_ce_selection_reg_ex_denP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l1_hdt --infile maphot_ce_selection_reg_l1_denP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l2_hdt --infile maphot_ce_selection_reg_l2_denP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l3_hdt --infile maphot_ce_selection_reg_l3_denP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l4_hdt --infile maphot_ce_selection_reg_l4_denP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_cx_hdb --infile maphot_ce_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_ex_hdb --infile maphot_ce_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l1_hdb --infile maphot_ce_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l2_hdb --infile maphot_ce_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l3_hdb --infile maphot_ce_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l4_hdb --infile maphot_ce_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l4_raw

#rm -rf ~/meTRN/data/peaks/ce_selection_reg_cx_hdt/
#mkdir ~/meTRN/data/peaks/ce_selection_reg_cx_hdt/
#cp ~/meTRN/data/peaks/ce_selection_reg_ex_hdt/* ~/meTRN/data/peaks/ce_selection_reg_cx_hdt/
#cp ~/meTRN/data/peaks/ce_selection_reg_l1_hdt/* ~/meTRN/data/peaks/ce_selection_reg_cx_hdt/
#cp ~/meTRN/data/peaks/ce_selection_reg_l2_hdt/* ~/meTRN/data/peaks/ce_selection_reg_cx_hdt/
#cp ~/meTRN/data/peaks/ce_selection_reg_l3_hdt/* ~/meTRN/data/peaks/ce_selection_reg_cx_hdt/
#cp ~/meTRN/data/peaks/ce_selection_reg_l4_hdt/* ~/meTRN/data/peaks/ce_selection_reg_cx_hdt/

#rm -rf ~/meTRN/data/peaks/ce_selection_com_cx_hdt/
#mkdir ~/meTRN/data/peaks/ce_selection_com_cx_hdt/
#cp ~/meTRN/data/peaks/ce_selection_com_ex_hdt/* ~/meTRN/data/peaks/ce_selection_com_cx_hdt/
#cp ~/meTRN/data/peaks/ce_selection_com_l1_hdt/* ~/meTRN/data/peaks/ce_selection_com_cx_hdt/
#cp ~/meTRN/data/peaks/ce_selection_com_l2_hdt/* ~/meTRN/data/peaks/ce_selection_com_cx_hdt/
#cp ~/meTRN/data/peaks/ce_selection_com_l3_hdt/* ~/meTRN/data/peaks/ce_selection_com_cx_hdt/
#cp ~/meTRN/data/peaks/ce_selection_com_l4_hdt/* ~/meTRN/data/peaks/ce_selection_com_cx_hdt/


# Generate HOT-filtered peak sets (combined P05 cutoff):
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_cx_hcn --infile maphot_ce_selection_reg_cx_simP05_any.bed --source ~/meTRN/data/peaks/ce_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_ex_hct --infile maphot_ce_selection_reg_ex_simP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l1_hct --infile maphot_ce_selection_reg_l1_simP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l2_hct --infile maphot_ce_selection_reg_l2_simP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l3_hct --infile maphot_ce_selection_reg_l3_simP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l4_hct --infile maphot_ce_selection_reg_l4_simP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_cx_hcb --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_ex_hcb --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l1_hcb --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l2_hcb --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l3_hcb --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l4_hcb --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_cx_hcn --infile maphot_ce_selection_reg_cx_simP05_any.bed --source ~/meTRN/data/peaks/ce_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_ex_hct --infile maphot_ce_selection_reg_ex_simP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l1_hct --infile maphot_ce_selection_reg_l1_simP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l2_hct --infile maphot_ce_selection_reg_l2_simP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l3_hct --infile maphot_ce_selection_reg_l3_simP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l4_hct --infile maphot_ce_selection_reg_l4_simP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_cx_hcb --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_ex_hcb --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l1_hcb --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l2_hcb --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l3_hcb --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l4_hcb --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l4_raw

#rm -rf ~/meTRN/data/peaks/ce_selection_reg_cx_hct/
#mkdir ~/meTRN/data/peaks/ce_selection_reg_cx_hct/
#cp ~/meTRN/data/peaks/ce_selection_reg_ex_hct/* ~/meTRN/data/peaks/ce_selection_reg_cx_hct/
#cp ~/meTRN/data/peaks/ce_selection_reg_l1_hct/* ~/meTRN/data/peaks/ce_selection_reg_cx_hct/
#cp ~/meTRN/data/peaks/ce_selection_reg_l2_hct/* ~/meTRN/data/peaks/ce_selection_reg_cx_hct/
#cp ~/meTRN/data/peaks/ce_selection_reg_l3_hct/* ~/meTRN/data/peaks/ce_selection_reg_cx_hct/
#cp ~/meTRN/data/peaks/ce_selection_reg_l4_hct/* ~/meTRN/data/peaks/ce_selection_reg_cx_hct/

#rm -rf ~/meTRN/data/peaks/ce_selection_com_cx_hct/
#mkdir ~/meTRN/data/peaks/ce_selection_com_cx_hct/
#cp ~/meTRN/data/peaks/ce_selection_com_ex_hct/* ~/meTRN/data/peaks/ce_selection_com_cx_hct/
#cp ~/meTRN/data/peaks/ce_selection_com_l1_hct/* ~/meTRN/data/peaks/ce_selection_com_cx_hct/
#cp ~/meTRN/data/peaks/ce_selection_com_l2_hct/* ~/meTRN/data/peaks/ce_selection_com_cx_hct/
#cp ~/meTRN/data/peaks/ce_selection_com_l3_hct/* ~/meTRN/data/peaks/ce_selection_com_cx_hct/
#cp ~/meTRN/data/peaks/ce_selection_com_l4_hct/* ~/meTRN/data/peaks/ce_selection_com_cx_hct/


# Generate HOT-filtered extended peak sets (occupancy P05 cutoff):
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_cx_hon --infile maphot_ce_selection_reg_cx_occP05_any.bed --source ~/meTRN/data/peaks/ce_extension_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_ex_hot --infile maphot_ce_selection_reg_ex_occP05_hot.bed --source ~/meTRN/data/peaks/ce_extension_com_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_l1_hot --infile maphot_ce_selection_reg_l1_occP05_hot.bed --source ~/meTRN/data/peaks/ce_extension_com_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_l2_hot --infile maphot_ce_selection_reg_l2_occP05_hot.bed --source ~/meTRN/data/peaks/ce_extension_com_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_l3_hot --infile maphot_ce_selection_reg_l3_occP05_hot.bed --source ~/meTRN/data/peaks/ce_extension_com_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_l4_hot --infile maphot_ce_selection_reg_l4_occP05_hot.bed --source ~/meTRN/data/peaks/ce_extension_com_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_cx_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_extension_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_ex_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_extension_com_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_l1_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_extension_com_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_l2_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_extension_com_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_l3_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_extension_com_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_l4_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_extension_com_l4_raw

#rm -rf ~/meTRN/data/peaks/ce_extension_com_cx_hot/
#mkdir ~/meTRN/data/peaks/ce_extension_com_cx_hot/
#cp ~/meTRN/data/peaks/ce_extension_com_ex_hot/* ~/meTRN/data/peaks/ce_extension_com_cx_hot/
#cp ~/meTRN/data/peaks/ce_extension_com_l1_hot/* ~/meTRN/data/peaks/ce_extension_com_cx_hot/
#cp ~/meTRN/data/peaks/ce_extension_com_l2_hot/* ~/meTRN/data/peaks/ce_extension_com_cx_hot/
#cp ~/meTRN/data/peaks/ce_extension_com_l3_hot/* ~/meTRN/data/peaks/ce_extension_com_cx_hot/
#cp ~/meTRN/data/peaks/ce_extension_com_l4_hot/* ~/meTRN/data/peaks/ce_extension_com_cx_hot/


# Generate HOT-filtered corrected peak sets (occupancy P05 cutoff):
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_cx_hon --infile maphot_ce_selection_reg_cx_occP05_any.bed --source ~/meTRN/data/peaks/ce_corrected_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_ex_hot --infile maphot_ce_selection_reg_ex_occP05_hot.bed --source ~/meTRN/data/peaks/ce_corrected_com_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_l1_hot --infile maphot_ce_selection_reg_l1_occP05_hot.bed --source ~/meTRN/data/peaks/ce_corrected_com_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_l2_hot --infile maphot_ce_selection_reg_l2_occP05_hot.bed --source ~/meTRN/data/peaks/ce_corrected_com_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_l3_hot --infile maphot_ce_selection_reg_l3_occP05_hot.bed --source ~/meTRN/data/peaks/ce_corrected_com_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_l4_hot --infile maphot_ce_selection_reg_l4_occP05_hot.bed --source ~/meTRN/data/peaks/ce_corrected_com_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_cx_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_corrected_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_ex_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_corrected_com_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_l1_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_corrected_com_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_l2_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_corrected_com_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_l3_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_corrected_com_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_l4_hob --infile maphot_ce_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/ce_corrected_com_l4_raw

#rm -rf ~/meTRN/data/peaks/ce_corrected_com_cx_hot/
#mkdir ~/meTRN/data/peaks/ce_corrected_com_cx_hot/
#cp ~/meTRN/data/peaks/ce_corrected_com_ex_hot/* ~/meTRN/data/peaks/ce_corrected_com_cx_hot/
#cp ~/meTRN/data/peaks/ce_corrected_com_l1_hot/* ~/meTRN/data/peaks/ce_corrected_com_cx_hot/
#cp ~/meTRN/data/peaks/ce_corrected_com_l2_hot/* ~/meTRN/data/peaks/ce_corrected_com_cx_hot/
#cp ~/meTRN/data/peaks/ce_corrected_com_l3_hot/* ~/meTRN/data/peaks/ce_corrected_com_cx_hot/
#cp ~/meTRN/data/peaks/ce_corrected_com_l4_hot/* ~/meTRN/data/peaks/ce_corrected_com_cx_hot/


# Generate complete, collapsed, density and report files:
#python mapPeaks.py --path ~/meTRN --mode build --overwrite OFF


#coverageBed -a /Volumes/HD1/Users/claraya/meTRN/data/peaks/mappeaks_ce_selection_com_cx_raw_compiled.bed -b /Volumes/HD1/Users/claraya/meTRN/input/ensembl_ws220_nuclear_sizes.bed -hist
#0.2170808 (0.7829192)
#coverageBed -a /Volumes/HD1/Users/claraya/meTRN/data/hot/regions/maphot_ce_selection_reg_cx_occP05_any.bed -b /Volumes/HD1/Users/claraya/meTRN/input/ensembl_ws220_nuclear_sizes.bed -hist
#0.0261704
#coverageBed -a /Volumes/HD1/Users/claraya/meTRN/data/hot/regions/maphot_ce_selection_reg_cx_occP05_all.bed -b /Volumes/HD1/Users/claraya/meTRN/input/ensembl_ws220_nuclear_sizes.bed -hist
#0.005785

#100*0.7829192
#RGB: 100*(0.2170808 - 0.0261704)
#HOT (any): 100*(0.0261704 - 0.0057849)
#HOT (all): 100*0.0057849

#wc -l ~/meTRN/data/peaks/mappeaks_ce_selection_com_cx_raw*
#datasets: 188
#regions: 33833
#peaks: 397539

#wc -l ~/meTRN/data/peaks/mappeaks_ce_selection_com_cx_hob*
#datasets: 188
#regions: 33365
#peaks: 306632
#percent: 306632/397539

#wc -l ~/meTRN/data/peaks/mappeaks_ce_selection_com_cx_hot*
#datasets: 188
#regions: 33823
#peaks: 230325
#percent: 230325/397539

#wc -l ~/meTRN/data/peaks/mappeaks_ce_selection_com_cx_hon*
#datasets: 188
#regions: 31638
#peaks: 186972
#percent: 186972/397539


#top
#bash "runMaster2C-ce (5%).sh"