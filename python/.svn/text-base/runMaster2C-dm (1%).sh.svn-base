#!/usr/bin/sh


# Map HOT regions (fail) and non-HOT regions (pass); without any actual filtering cutoffs:
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_cx --name basics
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_ee --name basics --tag EE-
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_le --name basics --tag LE-
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_pp --name basics --tag PP-

	
# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (occupancy 0.01 significance):
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_ee --name occP01 --significance 01 --tag EE- --metric occupancy --cutoff 15
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_le --name occP01 --significance 01 --tag LE- --metric occupancy --cutoff 10
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_pp --name occP01 --significance 01 --tag PP- --metric occupancy --cutoff 12


# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (density 0.01 significance):
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_ee --name denP01 --significance 01 --tag EE- --metric density --cutoff 17.3
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_le --name denP01 --significance 01 --tag LE- --metric density --cutoff 13.7
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_pp --name denP01 --significance 01 --tag PP- --metric density --cutoff 16.4


# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (combined 0.01 significance):
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_ee --name simP01 --significance 01 --tag EE- --metric combined --cutoff 15,17.3
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_le --name simP01 --significance 01 --tag LE- --metric combined --cutoff 10,13.7
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_pp --name simP01 --significance 01 --tag PP- --metric combined --cutoff 12,16.4


# Calculate overlap in HOT regions and filter regions against the ubiquitously HOT regions:
#python mapHOT.py --path ~/meTRN --organism dm --mode overlap --name occP01 --target EE:dm_selection_reg_ee,LE:dm_selection_reg_le,PP:dm_selection_reg_pp --overlap dm_selection_reg_cx
#python mapHOT.py --path ~/meTRN --organism dm --mode overlap --name denP01 --target EE:dm_selection_reg_ee,LE:dm_selection_reg_le,PP:dm_selection_reg_pp --overlap dm_selection_reg_cx
#python mapHOT.py --path ~/meTRN --organism dm --mode overlap --name simP01 --target EE:dm_selection_reg_ee,LE:dm_selection_reg_le,PP:dm_selection_reg_pp --overlap dm_selection_reg_cx


# Store HOT and RGB region files for each context and the global context file:
#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_ee_occP01 --regions dm_selection_reg_ee_occP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_le_occP01 --regions dm_selection_reg_le_occP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_pp_occP01 --regions dm_selection_reg_pp_occP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_cx_occP01 --regions dm_selection_reg_cx_occP01 --source overlap

#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_ee_denP01 --regions dm_selection_reg_ee_denP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_le_denP01 --regions dm_selection_reg_le_denP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_pp_denP01 --regions dm_selection_reg_pp_denP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_cx_denP01 --regions dm_selection_reg_cx_denP01 --source overlap

#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_ee_simP01 --regions dm_selection_reg_ee_simP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_le_simP01 --regions dm_selection_reg_le_simP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_pp_simP01 --regions dm_selection_reg_pp_simP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_cx_simP01 --regions dm_selection_reg_cx_simP01 --source overlap


# Generate HOT-filtered peak sets (occupancy P01 cutoff):
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_cx_xon --infile maphot_dm_selection_reg_cx_occP01_any.bed --source ~/meTRN/data/peaks/dm_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_ee_xot --infile maphot_dm_selection_reg_ee_occP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_le_xot --infile maphot_dm_selection_reg_le_occP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_pp_xot --infile maphot_dm_selection_reg_pp_occP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_pp_raw

#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_cx_xob --infile maphot_dm_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_ee_xob --infile maphot_dm_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_le_xob --infile maphot_dm_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_pp_xob --infile maphot_dm_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_pp_raw

#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_cx_xon --infile maphot_dm_selection_reg_cx_occP01_any.bed --source ~/meTRN/data/peaks/dm_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_ee_xot --infile maphot_dm_selection_reg_ee_occP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_le_xot --infile maphot_dm_selection_reg_le_occP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_pp_xot --infile maphot_dm_selection_reg_pp_occP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_pp_raw

#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_cx_xob --infile maphot_dm_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/dm_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_ee_xob --infile maphot_dm_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/dm_selection_com_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_le_xob --infile maphot_dm_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/dm_selection_com_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_pp_xob --infile maphot_dm_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/dm_selection_com_pp_raw

#rm -rf ~/meTRN/data/peaks/dm_selection_reg_cx_xot/
#mkdir ~/meTRN/data/peaks/dm_selection_reg_cx_xot/
#cp ~/meTRN/data/peaks/dm_selection_reg_ee_xot/* ~/meTRN/data/peaks/dm_selection_reg_cx_xot/
#cp ~/meTRN/data/peaks/dm_selection_reg_le_xot/* ~/meTRN/data/peaks/dm_selection_reg_cx_xot/
#cp ~/meTRN/data/peaks/dm_selection_reg_pp_xot/* ~/meTRN/data/peaks/dm_selection_reg_cx_xot/

#rm -rf ~/meTRN/data/peaks/dm_selection_com_cx_xot/
#mkdir ~/meTRN/data/peaks/dm_selection_com_cx_xot/
#cp ~/meTRN/data/peaks/dm_selection_com_ee_xot/* ~/meTRN/data/peaks/dm_selection_com_cx_xot/
#cp ~/meTRN/data/peaks/dm_selection_com_le_xot/* ~/meTRN/data/peaks/dm_selection_com_cx_xot/
#cp ~/meTRN/data/peaks/dm_selection_com_pp_xot/* ~/meTRN/data/peaks/dm_selection_com_cx_xot/


# Generate HOT-filtered peak sets (density P01 cutoff):
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_cx_xdn --infile maphot_dm_selection_reg_cx_denP01_any.bed --source ~/meTRN/data/peaks/dm_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_ee_xdt --infile maphot_dm_selection_reg_ee_denP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_le_xdt --infile maphot_dm_selection_reg_le_denP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_pp_xdt --infile maphot_dm_selection_reg_pp_denP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_pp_raw

#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_cx_xdb --infile maphot_dm_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_ee_xdb --infile maphot_dm_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_le_xdb --infile maphot_dm_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_pp_xdb --infile maphot_dm_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_pp_raw

#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_cx_xdn --infile maphot_dm_selection_reg_cx_denP01_any.bed --source ~/meTRN/data/peaks/dm_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_ee_xdt --infile maphot_dm_selection_reg_ee_denP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_le_xdt --infile maphot_dm_selection_reg_le_denP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_pp_xdt --infile maphot_dm_selection_reg_pp_denP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_pp_raw

#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_cx_xdb --infile maphot_dm_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/dm_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_ee_xdb --infile maphot_dm_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/dm_selection_com_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_le_xdb --infile maphot_dm_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/dm_selection_com_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_pp_xdb --infile maphot_dm_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/dm_selection_com_pp_raw

#rm -rf ~/meTRN/data/peaks/dm_selection_reg_cx_xdt/
#mkdir ~/meTRN/data/peaks/dm_selection_reg_cx_xdt/
#cp ~/meTRN/data/peaks/dm_selection_reg_ee_xdt/* ~/meTRN/data/peaks/dm_selection_reg_cx_xdt/
#cp ~/meTRN/data/peaks/dm_selection_reg_le_xdt/* ~/meTRN/data/peaks/dm_selection_reg_cx_xdt/
#cp ~/meTRN/data/peaks/dm_selection_reg_pp_xdt/* ~/meTRN/data/peaks/dm_selection_reg_cx_xdt/

#rm -rf ~/meTRN/data/peaks/dm_selection_com_cx_xdt/
#mkdir ~/meTRN/data/peaks/dm_selection_com_cx_xdt/
#cp ~/meTRN/data/peaks/dm_selection_com_ee_xdt/* ~/meTRN/data/peaks/dm_selection_com_cx_xdt/
#cp ~/meTRN/data/peaks/dm_selection_com_le_xdt/* ~/meTRN/data/peaks/dm_selection_com_cx_xdt/
#cp ~/meTRN/data/peaks/dm_selection_com_pp_xdt/* ~/meTRN/data/peaks/dm_selection_com_cx_xdt/


# Generate HOT-filtered peak sets (combined P01 cutoff):
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_cx_xcn --infile maphot_dm_selection_reg_cx_simP01_any.bed --source ~/meTRN/data/peaks/dm_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_ee_xct --infile maphot_dm_selection_reg_ee_simP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_le_xct --infile maphot_dm_selection_reg_le_simP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_pp_xct --infile maphot_dm_selection_reg_pp_simP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_pp_raw

#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_cx_xcb --infile maphot_dm_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_ee_xcb --infile maphot_dm_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_le_xcb --infile maphot_dm_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_pp_xcb --infile maphot_dm_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_pp_raw

#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_cx_xcn --infile maphot_dm_selection_reg_cx_simP01_any.bed --source ~/meTRN/data/peaks/dm_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_ee_xct --infile maphot_dm_selection_reg_ee_simP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_le_xct --infile maphot_dm_selection_reg_le_simP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_pp_xct --infile maphot_dm_selection_reg_pp_simP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_pp_raw

#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_cx_xcb --infile maphot_dm_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/dm_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_ee_xcb --infile maphot_dm_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/dm_selection_com_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_le_xcb --infile maphot_dm_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/dm_selection_com_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_pp_xcb --infile maphot_dm_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/dm_selection_com_pp_raw

#rm -rf ~/meTRN/data/peaks/dm_selection_reg_cx_xct/
#mkdir ~/meTRN/data/peaks/dm_selection_reg_cx_xct/
#cp ~/meTRN/data/peaks/dm_selection_reg_ee_xct/* ~/meTRN/data/peaks/dm_selection_reg_cx_xct/
#cp ~/meTRN/data/peaks/dm_selection_reg_le_xct/* ~/meTRN/data/peaks/dm_selection_reg_cx_xct/
#cp ~/meTRN/data/peaks/dm_selection_reg_pp_xct/* ~/meTRN/data/peaks/dm_selection_reg_cx_xct/

#rm -rf ~/meTRN/data/peaks/dm_selection_com_cx_xct/
#mkdir ~/meTRN/data/peaks/dm_selection_com_cx_xct/
#cp ~/meTRN/data/peaks/dm_selection_com_ee_xct/* ~/meTRN/data/peaks/dm_selection_com_cx_xct/
#cp ~/meTRN/data/peaks/dm_selection_com_le_xct/* ~/meTRN/data/peaks/dm_selection_com_cx_xct/
#cp ~/meTRN/data/peaks/dm_selection_com_pp_xct/* ~/meTRN/data/peaks/dm_selection_com_cx_xct/


# Generate complete, collapsed, density and report files:
#python mapPeaks.py --path ~/meTRN --mode build --overwrite OFF


#coverageBed -a ~/meTRN/data/peaks/mappeaks_dm_selection_com_cx_raw_compiled.bed -b ~/meTRN/input/flybase_dmel_r5-32_nuclear_sizes.bed -hist
#0.199868 (0.800132)
#coverageBed -a ~/meTRN/data/hot/regions/maphot_dm_selection_reg_cx_occP01_any.bed -b ~/meTRN/input/flybase_dmel_r5-32_nuclear_sizes.bed -hist
#0.0029082
#coverageBed -a ~/meTRN/data/hot/regions/maphot_dm_selection_reg_cx_occP01_all.bed -b ~/meTRN/input/flybase_dmel_r5-32_nuclear_sizes.bed -hist
#0.0000692

#100*0.199868
#RGB: 100*(0.199868 - 0.0029082)
#XOT (any): 100*(0.0029082 - 0.0000692)
#XOT (all): 100*0.0000692

#wc -l ~/meTRN/data/peaks/mappeaks_dm_selection_com_cx_raw*
#datasets: 64
#regions: 26104
#peaks: 99358

#wc -l ~/meTRN/data/peaks/mappeaks_dm_selection_com_cx_xob*
#datasets: 64
#regions: 26105
#peaks: 99201
#percent: 99201/99358

#wc -l ~/meTRN/data/peaks/mappeaks_dm_selection_com_cx_xot*
#datasets: 64
#regions: 26192
#peaks: 97351
#percent: 97351/99358

#wc -l ~/meTRN/data/peaks/mappeaks_dm_selection_com_cx_xon*
#datasets: 64
#regions: 26058
#peaks: 95680
#percent: 95680/99358



#top
#bash "runMaster2C-dm (1%).sh"