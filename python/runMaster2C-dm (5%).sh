#!/usr/bin/sh


# Map HOT regions (fail) and non-HOT regions (pass); without any actual filtering cutoffs:
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_cx --name basics
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_ee --name basics --tag EE-
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_le --name basics --tag LE-
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_pp --name basics --tag PP-

	
# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (occupancy 0.05 significance):
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_ee --name occP05 --significance 05 --tag EE- --metric occupancy --cutoff 5
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_le --name occP05 --significance 05 --tag LE- --metric occupancy --cutoff 3
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_pp --name occP05 --significance 05 --tag PP- --metric occupancy --cutoff 4


# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (density 0.05 significance):
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_ee --name denP05 --significance 05 --tag EE- --metric density --cutoff 9.5
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_le --name denP05 --significance 05 --tag LE- --metric density --cutoff 10.8
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_pp --name denP05 --significance 05 --tag PP- --metric density --cutoff 11.9


# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (combined 0.05 significance):
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_ee --name simP05 --significance 05 --tag EE- --metric combined --cutoff 5,9.5
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_le --name simP05 --significance 05 --tag LE- --metric combined --cutoff 3,10.8
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_pp --name simP05 --significance 05 --tag PP- --metric combined --cutoff 4,11.9


# Calculate overlap in HOT regions and filter regions against the ubiquitously HOT regions:
#python mapHOT.py --path ~/meTRN --organism dm --mode overlap --name occP05 --target EE:dm_selection_reg_ee,LE:dm_selection_reg_le,PP:dm_selection_reg_pp --overlap dm_selection_reg_cx
#python mapHOT.py --path ~/meTRN --organism dm --mode overlap --name denP05 --target EE:dm_selection_reg_ee,LE:dm_selection_reg_le,PP:dm_selection_reg_pp --overlap dm_selection_reg_cx
#python mapHOT.py --path ~/meTRN --organism dm --mode overlap --name simP05 --target EE:dm_selection_reg_ee,LE:dm_selection_reg_le,PP:dm_selection_reg_pp --overlap dm_selection_reg_cx


# Store HOT and RGB region files for each context and the global context file:
#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_ee_occP05 --regions dm_selection_reg_ee_occP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_le_occP05 --regions dm_selection_reg_le_occP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_pp_occP05 --regions dm_selection_reg_pp_occP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_cx_occP05 --regions dm_selection_reg_cx_occP05 --source overlap

#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_ee_denP05 --regions dm_selection_reg_ee_denP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_le_denP05 --regions dm_selection_reg_le_denP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_pp_denP05 --regions dm_selection_reg_pp_denP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_cx_denP05 --regions dm_selection_reg_cx_denP05 --source overlap

#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_ee_simP05 --regions dm_selection_reg_ee_simP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_le_simP05 --regions dm_selection_reg_le_simP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_pp_simP05 --regions dm_selection_reg_pp_simP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_cx_simP05 --regions dm_selection_reg_cx_simP05 --source overlap


# Generate HOT-filtered peak sets (occupancy P05 cutoff):
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_cx_hon --infile maphot_dm_selection_reg_cx_occP05_any.bed --source ~/meTRN/data/peaks/dm_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_ee_hot --infile maphot_dm_selection_reg_ee_occP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_le_hot --infile maphot_dm_selection_reg_le_occP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_pp_hot --infile maphot_dm_selection_reg_pp_occP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_pp_raw

#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_cx_hob --infile maphot_dm_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_ee_hob --infile maphot_dm_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_le_hob --infile maphot_dm_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_pp_hob --infile maphot_dm_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_pp_raw

#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_cx_hon --infile maphot_dm_selection_reg_cx_occP05_any.bed --source ~/meTRN/data/peaks/dm_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_ee_hot --infile maphot_dm_selection_reg_ee_occP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_le_hot --infile maphot_dm_selection_reg_le_occP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_pp_hot --infile maphot_dm_selection_reg_pp_occP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_pp_raw

#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_cx_hob --infile maphot_dm_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/dm_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_ee_hob --infile maphot_dm_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/dm_selection_com_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_le_hob --infile maphot_dm_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/dm_selection_com_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_pp_hob --infile maphot_dm_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/dm_selection_com_pp_raw

#rm -rf ~/meTRN/data/peaks/dm_selection_reg_cx_hot/
#mkdir ~/meTRN/data/peaks/dm_selection_reg_cx_hot/
#cp ~/meTRN/data/peaks/dm_selection_reg_ee_hot/* ~/meTRN/data/peaks/dm_selection_reg_cx_hot/
#cp ~/meTRN/data/peaks/dm_selection_reg_le_hot/* ~/meTRN/data/peaks/dm_selection_reg_cx_hot/
#cp ~/meTRN/data/peaks/dm_selection_reg_pp_hot/* ~/meTRN/data/peaks/dm_selection_reg_cx_hot/

#rm -rf ~/meTRN/data/peaks/dm_selection_com_cx_hot/
#mkdir ~/meTRN/data/peaks/dm_selection_com_cx_hot/
#cp ~/meTRN/data/peaks/dm_selection_com_ee_hot/* ~/meTRN/data/peaks/dm_selection_com_cx_hot/
#cp ~/meTRN/data/peaks/dm_selection_com_le_hot/* ~/meTRN/data/peaks/dm_selection_com_cx_hot/
#cp ~/meTRN/data/peaks/dm_selection_com_pp_hot/* ~/meTRN/data/peaks/dm_selection_com_cx_hot/


# Generate HOT-filtered peak sets (density P05 cutoff):
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_cx_hdn --infile maphot_dm_selection_reg_cx_denP05_any.bed --source ~/meTRN/data/peaks/dm_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_ee_hdt --infile maphot_dm_selection_reg_ee_denP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_le_hdt --infile maphot_dm_selection_reg_le_denP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_pp_hdt --infile maphot_dm_selection_reg_pp_denP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_pp_raw

#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_cx_hdb --infile maphot_dm_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_ee_hdb --infile maphot_dm_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_le_hdb --infile maphot_dm_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_pp_hdb --infile maphot_dm_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_pp_raw

#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_cx_hdn --infile maphot_dm_selection_reg_cx_denP05_any.bed --source ~/meTRN/data/peaks/dm_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_ee_hdt --infile maphot_dm_selection_reg_ee_denP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_le_hdt --infile maphot_dm_selection_reg_le_denP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_pp_hdt --infile maphot_dm_selection_reg_pp_denP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_pp_raw

#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_cx_hdb --infile maphot_dm_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/dm_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_ee_hdb --infile maphot_dm_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/dm_selection_com_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_le_hdb --infile maphot_dm_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/dm_selection_com_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_pp_hdb --infile maphot_dm_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/dm_selection_com_pp_raw

#rm -rf ~/meTRN/data/peaks/dm_selection_reg_cx_hdt/
#mkdir ~/meTRN/data/peaks/dm_selection_reg_cx_hdt/
#cp ~/meTRN/data/peaks/dm_selection_reg_ee_hdt/* ~/meTRN/data/peaks/dm_selection_reg_cx_hdt/
#cp ~/meTRN/data/peaks/dm_selection_reg_le_hdt/* ~/meTRN/data/peaks/dm_selection_reg_cx_hdt/
#cp ~/meTRN/data/peaks/dm_selection_reg_pp_hdt/* ~/meTRN/data/peaks/dm_selection_reg_cx_hdt/

#rm -rf ~/meTRN/data/peaks/dm_selection_com_cx_hdt/
#mkdir ~/meTRN/data/peaks/dm_selection_com_cx_hdt/
#cp ~/meTRN/data/peaks/dm_selection_com_ee_hdt/* ~/meTRN/data/peaks/dm_selection_com_cx_hdt/
#cp ~/meTRN/data/peaks/dm_selection_com_le_hdt/* ~/meTRN/data/peaks/dm_selection_com_cx_hdt/
#cp ~/meTRN/data/peaks/dm_selection_com_pp_hdt/* ~/meTRN/data/peaks/dm_selection_com_cx_hdt/


# Generate HOT-filtered peak sets (combined P05 cutoff):
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_cx_hcn --infile maphot_dm_selection_reg_cx_simP05_any.bed --source ~/meTRN/data/peaks/dm_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_ee_hct --infile maphot_dm_selection_reg_ee_simP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_le_hct --infile maphot_dm_selection_reg_le_simP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_pp_hct --infile maphot_dm_selection_reg_pp_simP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_pp_raw

#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_cx_hcb --infile maphot_dm_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_ee_hcb --infile maphot_dm_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_le_hcb --infile maphot_dm_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_pp_hcb --infile maphot_dm_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_pp_raw

#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_cx_hcn --infile maphot_dm_selection_reg_cx_simP05_any.bed --source ~/meTRN/data/peaks/dm_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_ee_hct --infile maphot_dm_selection_reg_ee_simP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_le_hct --infile maphot_dm_selection_reg_le_simP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_pp_hct --infile maphot_dm_selection_reg_pp_simP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_pp_raw

#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_cx_hcb --infile maphot_dm_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/dm_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_ee_hcb --infile maphot_dm_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/dm_selection_com_ee_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_le_hcb --infile maphot_dm_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/dm_selection_com_le_raw
#python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_pp_hcb --infile maphot_dm_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/dm_selection_com_pp_raw

#rm -rf ~/meTRN/data/peaks/dm_selection_reg_cx_hct/
#mkdir ~/meTRN/data/peaks/dm_selection_reg_cx_hct/
#cp ~/meTRN/data/peaks/dm_selection_reg_ee_hct/* ~/meTRN/data/peaks/dm_selection_reg_cx_hct/
#cp ~/meTRN/data/peaks/dm_selection_reg_le_hct/* ~/meTRN/data/peaks/dm_selection_reg_cx_hct/
#cp ~/meTRN/data/peaks/dm_selection_reg_pp_hct/* ~/meTRN/data/peaks/dm_selection_reg_cx_hct/

#rm -rf ~/meTRN/data/peaks/dm_selection_com_cx_hct/
#mkdir ~/meTRN/data/peaks/dm_selection_com_cx_hct/
#cp ~/meTRN/data/peaks/dm_selection_com_ee_hct/* ~/meTRN/data/peaks/dm_selection_com_cx_hct/
#cp ~/meTRN/data/peaks/dm_selection_com_le_hct/* ~/meTRN/data/peaks/dm_selection_com_cx_hct/
#cp ~/meTRN/data/peaks/dm_selection_com_pp_hct/* ~/meTRN/data/peaks/dm_selection_com_cx_hct/


# Generate complete, collapsed, density and report files:
#python mapPeaks.py --path ~/meTRN --mode build --overwrite OFF


#coverageBed -a /Volumes/HD1/Users/claraya/meTRN/data/peaks/mappeaks_dm_selection_com_cx_raw_compiled.bed -b /Volumes/HD1/Users/claraya/meTRN/input/flybase_dmel_r5-32_nuclear_sizes.bed -hist
#0.199868 (0.8001320)
#coverageBed -a /Volumes/HD1/Users/claraya/meTRN/data/hot/regions/maphot_dm_selection_reg_cx_occP05_any.bed -b /Volumes/HD1/Users/claraya/meTRN/input/flybase_dmel_r5-32_nuclear_sizes.bed -hist
#0.0359582
#coverageBed -a /Volumes/HD1/Users/claraya/meTRN/data/hot/regions/maphot_dm_selection_reg_cx_occP05_all.bed -b /Volumes/HD1/Users/claraya/meTRN/input/flybase_dmel_r5-32_nuclear_sizes.bed -hist
#0.0055493

#100*0.199868
#RGB: 100*(0.199868 - 0.0359582)
#HOT (any): 100*(0.0359582 - 0.0055493)
#HOT (all): 100*0.0055493

#wc -l ~/meTRN/data/peaks/mappeaks_dm_selection_com_cx_raw*
#datasets: 64
#regions: 26104
#peaks: 99358

#wc -l ~/meTRN/data/peaks/mappeaks_dm_selection_com_cx_hob*
#datasets: 64
#regions: 25976
#peaks: 90684
#percent: 90684/99358

#wc -l ~/meTRN/data/peaks/mappeaks_dm_selection_com_cx_hot*
#datasets: 64
#regions: 26862
#peaks: 73192
#percent: 73192/99358

#wc -l ~/meTRN/data/peaks/mappeaks_dm_selection_com_cx_hon*
#datasets: 64
#regions: 24338
#peaks: 61540
#percent: 61540/99358


#top
#bash "runMaster2C-dm (5%).sh"