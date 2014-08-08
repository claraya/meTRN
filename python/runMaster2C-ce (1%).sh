#!/usr/bin/sh


# Map HOT regions (fail) and non-HOT regions (pass); without any actual filtering cutoffs:
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_cx --name basics
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_ex --name basics --tag EX-
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l1 --name basics --tag L1-
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l2 --name basics --tag L2-
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l3 --name basics --tag L3-
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l4 --name basics --tag L4-

	
# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (occupancy 0.01 significance):
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_ex --name occP01 --significance 01 --tag EX- --metric occupancy --cutoff 9
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l1 --name occP01 --significance 01 --tag L1- --metric occupancy --cutoff 18
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l2 --name occP01 --significance 01 --tag L2- --metric occupancy --cutoff 18
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l3 --name occP01 --significance 01 --tag L3- --metric occupancy --cutoff 12
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l4 --name occP01 --significance 01 --tag L4- --metric occupancy --cutoff 16


# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (density 0.01 significance):
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_ex --name denP01 --significance 01 --tag EX- --metric density --cutoff 9.9
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l1 --name denP01 --significance 01 --tag L1- --metric density --cutoff 13.7
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l2 --name denP01 --significance 01 --tag L2- --metric density --cutoff 11.6
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l3 --name denP01 --significance 01 --tag L3- --metric density --cutoff 10.8
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l4 --name denP01 --significance 01 --tag L4- --metric density --cutoff 10.9


# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (combined 0.01 significance):
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_ex --name simP01 --significance 01 --tag EX- --metric combined --cutoff 9,9.9
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l1 --name simP01 --significance 01 --tag L1- --metric combined --cutoff 18,13.7
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l2 --name simP01 --significance 01 --tag L2- --metric combined --cutoff 18,11.6
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l3 --name simP01 --significance 01 --tag L3- --metric combined --cutoff 12,10.8
#python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l4 --name simP01 --significance 01 --tag L4- --metric combined --cutoff 16,10.9


# Calculate overlap in HOT regions and filter regions against the ubiquitously HOT regions:
#python mapHOT.py --path ~/meTRN --organism ce --mode overlap --name occP01 --target EX:ce_selection_reg_ex,L1:ce_selection_reg_l1,L2:ce_selection_reg_l2,L3:ce_selection_reg_l3,L4:ce_selection_reg_l4 --overlap ce_selection_reg_cx
#python mapHOT.py --path ~/meTRN --organism ce --mode overlap --name denP01 --target EX:ce_selection_reg_ex,L1:ce_selection_reg_l1,L2:ce_selection_reg_l2,L3:ce_selection_reg_l3,L4:ce_selection_reg_l4 --overlap ce_selection_reg_cx
#python mapHOT.py --path ~/meTRN --organism ce --mode overlap --name simP01 --target EX:ce_selection_reg_ex,L1:ce_selection_reg_l1,L2:ce_selection_reg_l2,L3:ce_selection_reg_l3,L4:ce_selection_reg_l4 --overlap ce_selection_reg_cx


# Store HOT and RGB region files for each context and the global context file:
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_ex_occP01 --regions ce_selection_reg_ex_occP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l1_occP01 --regions ce_selection_reg_l1_occP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l2_occP01 --regions ce_selection_reg_l2_occP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l3_occP01 --regions ce_selection_reg_l3_occP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l4_occP01 --regions ce_selection_reg_l4_occP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_cx_occP01 --regions ce_selection_reg_cx_occP01 --source overlap

#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_ex_denP01 --regions ce_selection_reg_ex_denP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l1_denP01 --regions ce_selection_reg_l1_denP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l2_denP01 --regions ce_selection_reg_l2_denP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l3_denP01 --regions ce_selection_reg_l3_denP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l4_denP01 --regions ce_selection_reg_l4_denP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_cx_denP01 --regions ce_selection_reg_cx_denP01 --source overlap

#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_ex_simP01 --regions ce_selection_reg_ex_simP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l1_simP01 --regions ce_selection_reg_l1_simP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l2_simP01 --regions ce_selection_reg_l2_simP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l3_simP01 --regions ce_selection_reg_l3_simP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l4_simP01 --regions ce_selection_reg_l4_simP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_cx_simP01 --regions ce_selection_reg_cx_simP01 --source overlap


# Generate HOT-filtered peak sets (occupancy P01 cutoff):
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_cx_xon --infile maphot_ce_selection_reg_cx_occP01_any.bed --source ~/meTRN/data/peaks/ce_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_ex_xot --infile maphot_ce_selection_reg_ex_occP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l1_xot --infile maphot_ce_selection_reg_l1_occP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l2_xot --infile maphot_ce_selection_reg_l2_occP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l3_xot --infile maphot_ce_selection_reg_l3_occP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l4_xot --infile maphot_ce_selection_reg_l4_occP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_cx_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_ex_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l1_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l2_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l3_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l4_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_cx_xon --infile maphot_ce_selection_reg_cx_occP01_any.bed --source ~/meTRN/data/peaks/ce_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_ex_xot --infile maphot_ce_selection_reg_ex_occP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l1_xot --infile maphot_ce_selection_reg_l1_occP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l2_xot --infile maphot_ce_selection_reg_l2_occP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l3_xot --infile maphot_ce_selection_reg_l3_occP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l4_xot --infile maphot_ce_selection_reg_l4_occP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_cx_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_ex_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l1_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l2_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l3_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l4_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l4_raw

#rm -rf ~/meTRN/data/peaks/ce_selection_reg_cx_xot/
#mkdir ~/meTRN/data/peaks/ce_selection_reg_cx_xot/
#cp ~/meTRN/data/peaks/ce_selection_reg_ex_xot/* ~/meTRN/data/peaks/ce_selection_reg_cx_xot/
#cp ~/meTRN/data/peaks/ce_selection_reg_l1_xot/* ~/meTRN/data/peaks/ce_selection_reg_cx_xot/
#cp ~/meTRN/data/peaks/ce_selection_reg_l2_xot/* ~/meTRN/data/peaks/ce_selection_reg_cx_xot/
#cp ~/meTRN/data/peaks/ce_selection_reg_l3_xot/* ~/meTRN/data/peaks/ce_selection_reg_cx_xot/
#cp ~/meTRN/data/peaks/ce_selection_reg_l4_xot/* ~/meTRN/data/peaks/ce_selection_reg_cx_xot/

#rm -rf ~/meTRN/data/peaks/ce_selection_com_cx_xot/
#mkdir ~/meTRN/data/peaks/ce_selection_com_cx_xot/
#cp ~/meTRN/data/peaks/ce_selection_com_ex_xot/* ~/meTRN/data/peaks/ce_selection_com_cx_xot/
#cp ~/meTRN/data/peaks/ce_selection_com_l1_xot/* ~/meTRN/data/peaks/ce_selection_com_cx_xot/
#cp ~/meTRN/data/peaks/ce_selection_com_l2_xot/* ~/meTRN/data/peaks/ce_selection_com_cx_xot/
#cp ~/meTRN/data/peaks/ce_selection_com_l3_xot/* ~/meTRN/data/peaks/ce_selection_com_cx_xot/
#cp ~/meTRN/data/peaks/ce_selection_com_l4_xot/* ~/meTRN/data/peaks/ce_selection_com_cx_xot/


# Generate HOT-filtered peak sets (density P01 cutoff):
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_cx_xdn --infile maphot_ce_selection_reg_cx_denP01_any.bed --source ~/meTRN/data/peaks/ce_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_ex_xdt --infile maphot_ce_selection_reg_ex_denP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l1_xdt --infile maphot_ce_selection_reg_l1_denP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l2_xdt --infile maphot_ce_selection_reg_l2_denP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l3_xdt --infile maphot_ce_selection_reg_l3_denP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l4_xdt --infile maphot_ce_selection_reg_l4_denP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_cx_xdb --infile maphot_ce_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_ex_xdb --infile maphot_ce_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l1_xdb --infile maphot_ce_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l2_xdb --infile maphot_ce_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l3_xdb --infile maphot_ce_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l4_xdb --infile maphot_ce_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_cx_xdn --infile maphot_ce_selection_reg_cx_denP01_any.bed --source ~/meTRN/data/peaks/ce_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_ex_xdt --infile maphot_ce_selection_reg_ex_denP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l1_xdt --infile maphot_ce_selection_reg_l1_denP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l2_xdt --infile maphot_ce_selection_reg_l2_denP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l3_xdt --infile maphot_ce_selection_reg_l3_denP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l4_xdt --infile maphot_ce_selection_reg_l4_denP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_cx_xdb --infile maphot_ce_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_ex_xdb --infile maphot_ce_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l1_xdb --infile maphot_ce_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l2_xdb --infile maphot_ce_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l3_xdb --infile maphot_ce_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l4_xdb --infile maphot_ce_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l4_raw

#rm -rf ~/meTRN/data/peaks/ce_selection_reg_cx_xdt/
#mkdir ~/meTRN/data/peaks/ce_selection_reg_cx_xdt/
#cp ~/meTRN/data/peaks/ce_selection_reg_ex_xdt/* ~/meTRN/data/peaks/ce_selection_reg_cx_xdt/
#cp ~/meTRN/data/peaks/ce_selection_reg_l1_xdt/* ~/meTRN/data/peaks/ce_selection_reg_cx_xdt/
#cp ~/meTRN/data/peaks/ce_selection_reg_l2_xdt/* ~/meTRN/data/peaks/ce_selection_reg_cx_xdt/
#cp ~/meTRN/data/peaks/ce_selection_reg_l3_xdt/* ~/meTRN/data/peaks/ce_selection_reg_cx_xdt/
#cp ~/meTRN/data/peaks/ce_selection_reg_l4_xdt/* ~/meTRN/data/peaks/ce_selection_reg_cx_xdt/

#rm -rf ~/meTRN/data/peaks/ce_selection_com_cx_xdt/
#mkdir ~/meTRN/data/peaks/ce_selection_com_cx_xdt/
#cp ~/meTRN/data/peaks/ce_selection_com_ex_xdt/* ~/meTRN/data/peaks/ce_selection_com_cx_xdt/
#cp ~/meTRN/data/peaks/ce_selection_com_l1_xdt/* ~/meTRN/data/peaks/ce_selection_com_cx_xdt/
#cp ~/meTRN/data/peaks/ce_selection_com_l2_xdt/* ~/meTRN/data/peaks/ce_selection_com_cx_xdt/
#cp ~/meTRN/data/peaks/ce_selection_com_l3_xdt/* ~/meTRN/data/peaks/ce_selection_com_cx_xdt/
#cp ~/meTRN/data/peaks/ce_selection_com_l4_xdt/* ~/meTRN/data/peaks/ce_selection_com_cx_xdt/


# Generate HOT-filtered peak sets (combined P01 cutoff):
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_cx_xcn --infile maphot_ce_selection_reg_cx_simP01_any.bed --source ~/meTRN/data/peaks/ce_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_ex_xct --infile maphot_ce_selection_reg_ex_simP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l1_xct --infile maphot_ce_selection_reg_l1_simP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l2_xct --infile maphot_ce_selection_reg_l2_simP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l3_xct --infile maphot_ce_selection_reg_l3_simP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l4_xct --infile maphot_ce_selection_reg_l4_simP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_cx_xcb --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_ex_xcb --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l1_xcb --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l2_xcb --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l3_xcb --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l4_xcb --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_cx_xcn --infile maphot_ce_selection_reg_cx_simP01_any.bed --source ~/meTRN/data/peaks/ce_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_ex_xct --infile maphot_ce_selection_reg_ex_simP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l1_xct --infile maphot_ce_selection_reg_l1_simP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l2_xct --infile maphot_ce_selection_reg_l2_simP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l3_xct --infile maphot_ce_selection_reg_l3_simP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l4_xct --infile maphot_ce_selection_reg_l4_simP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_cx_xcb --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_ex_xcb --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l1_xcb --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l2_xcb --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l3_xcb --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l4_xcb --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l4_raw

#rm -rf ~/meTRN/data/peaks/ce_selection_reg_cx_xct/
#mkdir ~/meTRN/data/peaks/ce_selection_reg_cx_xct/
#cp ~/meTRN/data/peaks/ce_selection_reg_ex_xct/* ~/meTRN/data/peaks/ce_selection_reg_cx_xct/
#cp ~/meTRN/data/peaks/ce_selection_reg_l1_xct/* ~/meTRN/data/peaks/ce_selection_reg_cx_xct/
#cp ~/meTRN/data/peaks/ce_selection_reg_l2_xct/* ~/meTRN/data/peaks/ce_selection_reg_cx_xct/
#cp ~/meTRN/data/peaks/ce_selection_reg_l3_xct/* ~/meTRN/data/peaks/ce_selection_reg_cx_xct/
#cp ~/meTRN/data/peaks/ce_selection_reg_l4_xct/* ~/meTRN/data/peaks/ce_selection_reg_cx_xct/

#rm -rf ~/meTRN/data/peaks/ce_selection_com_cx_xct/
#mkdir ~/meTRN/data/peaks/ce_selection_com_cx_xct/
#cp ~/meTRN/data/peaks/ce_selection_com_ex_xct/* ~/meTRN/data/peaks/ce_selection_com_cx_xct/
#cp ~/meTRN/data/peaks/ce_selection_com_l1_xct/* ~/meTRN/data/peaks/ce_selection_com_cx_xct/
#cp ~/meTRN/data/peaks/ce_selection_com_l2_xct/* ~/meTRN/data/peaks/ce_selection_com_cx_xct/
#cp ~/meTRN/data/peaks/ce_selection_com_l3_xct/* ~/meTRN/data/peaks/ce_selection_com_cx_xct/
#cp ~/meTRN/data/peaks/ce_selection_com_l4_xct/* ~/meTRN/data/peaks/ce_selection_com_cx_xct/


# Generate HOT-filtered extended peak sets (occupancy P01 cutoff):
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_cx_xon --infile maphot_ce_selection_reg_cx_occP01_any.bed --source ~/meTRN/data/peaks/ce_extension_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_ex_xot --infile maphot_ce_selection_reg_ex_occP01_hot.bed --source ~/meTRN/data/peaks/ce_extension_com_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_l1_xot --infile maphot_ce_selection_reg_l1_occP01_hot.bed --source ~/meTRN/data/peaks/ce_extension_com_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_l2_xot --infile maphot_ce_selection_reg_l2_occP01_hot.bed --source ~/meTRN/data/peaks/ce_extension_com_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_l3_xot --infile maphot_ce_selection_reg_l3_occP01_hot.bed --source ~/meTRN/data/peaks/ce_extension_com_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_l4_xot --infile maphot_ce_selection_reg_l4_occP01_hot.bed --source ~/meTRN/data/peaks/ce_extension_com_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_cx_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_extension_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_ex_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_extension_com_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_l1_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_extension_com_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_l2_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_extension_com_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_l3_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_extension_com_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_extension_com_l4_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_extension_com_l4_raw

#rm -rf ~/meTRN/data/peaks/ce_extension_com_cx_xot/
#mkdir ~/meTRN/data/peaks/ce_extension_com_cx_xot/
#cp ~/meTRN/data/peaks/ce_extension_com_ex_xot/* ~/meTRN/data/peaks/ce_extension_com_cx_xot/
#cp ~/meTRN/data/peaks/ce_extension_com_l1_xot/* ~/meTRN/data/peaks/ce_extension_com_cx_xot/
#cp ~/meTRN/data/peaks/ce_extension_com_l2_xot/* ~/meTRN/data/peaks/ce_extension_com_cx_xot/
#cp ~/meTRN/data/peaks/ce_extension_com_l3_xot/* ~/meTRN/data/peaks/ce_extension_com_cx_xot/
#cp ~/meTRN/data/peaks/ce_extension_com_l4_xot/* ~/meTRN/data/peaks/ce_extension_com_cx_xot/



# Generate HOT-filtered corrected peak sets (occupancy P01 cutoff):
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_cx_xon --infile maphot_ce_selection_reg_cx_occP01_any.bed --source ~/meTRN/data/peaks/ce_corrected_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_ex_xot --infile maphot_ce_selection_reg_ex_occP01_hot.bed --source ~/meTRN/data/peaks/ce_corrected_com_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_l1_xot --infile maphot_ce_selection_reg_l1_occP01_hot.bed --source ~/meTRN/data/peaks/ce_corrected_com_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_l2_xot --infile maphot_ce_selection_reg_l2_occP01_hot.bed --source ~/meTRN/data/peaks/ce_corrected_com_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_l3_xot --infile maphot_ce_selection_reg_l3_occP01_hot.bed --source ~/meTRN/data/peaks/ce_corrected_com_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_l4_xot --infile maphot_ce_selection_reg_l4_occP01_hot.bed --source ~/meTRN/data/peaks/ce_corrected_com_l4_raw

#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_cx_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_corrected_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_ex_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_corrected_com_ex_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_l1_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_corrected_com_l1_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_l2_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_corrected_com_l2_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_l3_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_corrected_com_l3_raw
#python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_corrected_com_l4_xob --infile maphot_ce_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/ce_corrected_com_l4_raw

#rm -rf ~/meTRN/data/peaks/ce_corrected_com_cx_xot/
#mkdir ~/meTRN/data/peaks/ce_corrected_com_cx_xot/
#cp ~/meTRN/data/peaks/ce_corrected_com_ex_xot/* ~/meTRN/data/peaks/ce_corrected_com_cx_xot/
#cp ~/meTRN/data/peaks/ce_corrected_com_l1_xot/* ~/meTRN/data/peaks/ce_corrected_com_cx_xot/
#cp ~/meTRN/data/peaks/ce_corrected_com_l2_xot/* ~/meTRN/data/peaks/ce_corrected_com_cx_xot/
#cp ~/meTRN/data/peaks/ce_corrected_com_l3_xot/* ~/meTRN/data/peaks/ce_corrected_com_cx_xot/
#cp ~/meTRN/data/peaks/ce_corrected_com_l4_xot/* ~/meTRN/data/peaks/ce_corrected_com_cx_xot/




# Generate complete, collapsed, density and report files:
#python mapPeaks.py --path ~/meTRN --mode build --overwrite OFF


#coverageBed -a /Volumes/HD1/Users/claraya/meTRN/data/peaks/mappeaks_ce_selection_com_cx_raw_compiled.bed -b /Volumes/HD1/Users/claraya/meTRN/input/ensembl_ws220_nuclear_sizes.bed -hist
#0.2170808 (0.7829192)
#coverageBed -a /Volumes/HD1/Users/claraya/meTRN/data/hot/regions/maphot_ce_selection_reg_cx_occP01_any.bed -b /Volumes/HD1/Users/claraya/meTRN/input/ensembl_ws220_nuclear_sizes.bed -hist
#0.0152394
#coverageBed -a /Volumes/HD1/Users/claraya/meTRN/data/hot/regions/maphot_ce_selection_reg_cx_occP01_all.bed -b /Volumes/HD1/Users/claraya/meTRN/input/ensembl_ws220_nuclear_sizes.bed -hist
#0.0030759

#100*0.7829192
#RGB: 100*(0.2170808 - 0.0152394)
#XOT (any): 100*(0.0152394 - 0.0030759)
#XOT (all): 100*0.0030759

#wc -l ~/meTRN/data/peaks/mappeaks_ce_selection_com_cx_raw*
#datasets: 188
#regions: 33833
#peaks: 397539

#wc -l ~/meTRN/data/peaks/mappeaks_ce_selection_com_cx_xob*
#datasets: 188
#regions: 33626
#peaks: 342492
#percent: 342492/397539

#wc -l ~/meTRN/data/peaks/mappeaks_ce_selection_com_cx_xot*
#datasets: 188
#regions: 33914
#peaks: 292466
#percent: 292466/397539

#wc -l ~/meTRN/data/peaks/mappeaks_ce_selection_com_cx_xon*
#datasets: 188
#regions: 32720
#peaks: 250953
#percent: 250953/397539


#top
#bash "runMaster2C-ce (1%).sh"