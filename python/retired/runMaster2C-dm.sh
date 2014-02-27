#!/usr/bin/sh


# Map HOT regions (fail) and non-HOT regions (pass); without any actual filtering cutoffs:
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_cx --name basics
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_ee --name basics --tag EE-
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_le --name basics --tag LE-
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_pp --name basics --tag PP-

	
# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (0.01 significance):
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_ee --name simP01 --cutoff 15 --tag EE-
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_le --name simP01 --cutoff 10 --tag LE-
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_pp --name simP01 --cutoff 12 --tag PP-


# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (0.05 significance):
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_ee --name simP05 --cutoff 5 --tag EE-
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_le --name simP05 --cutoff 3 --tag LE-
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_pp --name simP05 --cutoff 4 --tag PP-


# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (0.10 significance):
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_ee --name simP10 --cutoff 4 --tag EE-
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_le --name simP10 --cutoff 2 --tag LE-
#python mapHOT.py --path ~/meTRN --organism dm --mode scan --target peaks --peaks dm_selection_reg_pp --name simP10 --cutoff 3 --tag PP-


# Calculate overlap in HOT regions and filter regions against the ubiquitously HOT regions:
python mapHOT.py --path ~/meTRN --organism dm --mode overlap --name simP01 --target EE:dm_selection_reg_ee,LE:dm_selection_reg_le,PP:dm_selection_reg_pp --overlap dm_selection_reg_cx
python mapHOT.py --path ~/meTRN --organism dm --mode overlap --name simP05 --target EE:dm_selection_reg_ee,LE:dm_selection_reg_le,PP:dm_selection_reg_pp --overlap dm_selection_reg_cx


# Store HOT and RGB region files for each context and the global context file:
python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_ee_simP01 --regions dm_selection_reg_ee_simP01 --source analysis
python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_le_simP01 --regions dm_selection_reg_le_simP01 --source analysis
python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_pp_simP01 --regions dm_selection_reg_pp_simP01 --source analysis
python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_cx_simP01 --regions dm_selection_reg_cx_simP01 --source overlap

python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_ee_simP05 --regions dm_selection_reg_ee_simP05 --source analysis
python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_le_simP05 --regions dm_selection_reg_le_simP05 --source analysis
python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_pp_simP05 --regions dm_selection_reg_pp_simP05 --source analysis
python mapHOT.py --path ~/meTRN --organism dm --mode regions --name dm_selection_reg_cx_simP05 --regions dm_selection_reg_cx_simP05 --source overlap


# Generate HOT-filtered peak sets (P01 cutoff):
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_cx_xno --infile maphot_dm_selection_reg_cx_simP01_any.bed --source ~/meTRN/data/peaks/dm_selection_reg_cx_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_ee_xct --infile maphot_dm_selection_reg_ee_simP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_ee_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_le_xct --infile maphot_dm_selection_reg_le_simP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_le_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_pp_xct --infile maphot_dm_selection_reg_pp_simP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_pp_raw

python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_cx_xub --infile maphot_dm_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_cx_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_ee_xub --infile maphot_dm_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_ee_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_le_xub --infile maphot_dm_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_le_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_pp_xub --infile maphot_dm_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_pp_raw

python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_cx_xno --infile maphot_dm_selection_reg_cx_simP01_any.bed --source ~/meTRN/data/peaks/dm_selection_com_cx_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_ee_xct --infile maphot_dm_selection_reg_ee_simP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_ee_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_le_xct --infile maphot_dm_selection_reg_le_simP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_le_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_pp_xct --infile maphot_dm_selection_reg_pp_simP01_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_pp_raw

python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_cx_xub --infile maphot_dm_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/dm_selection_com_cx_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_ee_xub --infile maphot_dm_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/dm_selection_com_ee_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_le_xub --infile maphot_dm_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/dm_selection_com_le_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_pp_xub --infile maphot_dm_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/dm_selection_com_pp_raw


rm -rf ~/meTRN/data/peaks/dm_selection_reg_cx_xct/
mkdir ~/meTRN/data/peaks/dm_selection_reg_cx_xct/
cp ~/meTRN/data/peaks/dm_selection_reg_ee_xct/* ~/meTRN/data/peaks/dm_selection_reg_cx_xct/
cp ~/meTRN/data/peaks/dm_selection_reg_le_xct/* ~/meTRN/data/peaks/dm_selection_reg_cx_xct/
cp ~/meTRN/data/peaks/dm_selection_reg_pp_xct/* ~/meTRN/data/peaks/dm_selection_reg_cx_xct/

rm -rf ~/meTRN/data/peaks/dm_selection_com_cx_xct/
mkdir ~/meTRN/data/peaks/dm_selection_com_cx_xct/
cp ~/meTRN/data/peaks/dm_selection_com_ee_xct/* ~/meTRN/data/peaks/dm_selection_com_cx_xct/
cp ~/meTRN/data/peaks/dm_selection_com_le_xct/* ~/meTRN/data/peaks/dm_selection_com_cx_xct/
cp ~/meTRN/data/peaks/dm_selection_com_pp_xct/* ~/meTRN/data/peaks/dm_selection_com_cx_xct/


# Generate HOT-filtered peak sets (P05 cutoff):
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_cx_hno --infile maphot_dm_selection_reg_cx_simP05_any.bed --source ~/meTRN/data/peaks/dm_selection_reg_cx_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_ee_hct --infile maphot_dm_selection_reg_ee_simP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_ee_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_le_hct --infile maphot_dm_selection_reg_le_simP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_le_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_pp_hct --infile maphot_dm_selection_reg_pp_simP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_reg_pp_raw

python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_cx_hub --infile maphot_dm_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_cx_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_ee_hub --infile maphot_dm_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_ee_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_le_hub --infile maphot_dm_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_le_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_reg_pp_hub --infile maphot_dm_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/dm_selection_reg_pp_raw

python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_cx_hno --infile maphot_dm_selection_reg_cx_simP05_any.bed --source ~/meTRN/data/peaks/dm_selection_com_cx_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_ee_hct --infile maphot_dm_selection_reg_ee_simP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_ee_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_le_hct --infile maphot_dm_selection_reg_le_simP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_le_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_pp_hct --infile maphot_dm_selection_reg_pp_simP05_hot.bed --source ~/meTRN/data/peaks/dm_selection_com_pp_raw

python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_cx_hub --infile maphot_dm_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/dm_selection_com_cx_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_ee_hub --infile maphot_dm_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/dm_selection_com_ee_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_le_hub --infile maphot_dm_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/dm_selection_com_le_raw
python mapHOT.py --path ~/meTRN --organism dm --mode filter:remove --peaks dm_selection_com_pp_hub --infile maphot_dm_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/dm_selection_com_pp_raw

rm -rf ~/meTRN/data/peaks/dm_selection_reg_cx_hct/
mkdir ~/meTRN/data/peaks/dm_selection_reg_cx_hct/
cp ~/meTRN/data/peaks/dm_selection_reg_ee_hct/* ~/meTRN/data/peaks/dm_selection_reg_cx_hct/
cp ~/meTRN/data/peaks/dm_selection_reg_le_hct/* ~/meTRN/data/peaks/dm_selection_reg_cx_hct/
cp ~/meTRN/data/peaks/dm_selection_reg_pp_hct/* ~/meTRN/data/peaks/dm_selection_reg_cx_hct/

rm -rf ~/meTRN/data/peaks/dm_selection_com_cx_hct/
mkdir ~/meTRN/data/peaks/dm_selection_com_cx_hct/
cp ~/meTRN/data/peaks/dm_selection_com_ee_hct/* ~/meTRN/data/peaks/dm_selection_com_cx_hct/
cp ~/meTRN/data/peaks/dm_selection_com_le_hct/* ~/meTRN/data/peaks/dm_selection_com_cx_hct/
cp ~/meTRN/data/peaks/dm_selection_com_pp_hct/* ~/meTRN/data/peaks/dm_selection_com_cx_hct/


# Generate complete, collapsed, density and report files:
#python mapPeaks.py --path ~/meTRN --mode build --overwrite OFF

#coverageBed -a ~/meTRN/data/peaks/mappeaks_dm_selection_com_cx_raw_compiled.bed -b ~/meTRN/input/ensembl_ws220_nuclear_sizes.bed -hist
#0.0310155 (0.9689844)
#coverageBed -a ~/meTRN/data/hot/regions/maphot_dm_selection_reg_cx_simP05_any.bed -b ~/meTRN/input/ensembl_ws220_nuclear_sizes.bed -hist
#0.0075204
#coverageBed -a ~/meTRN/data/hot/regions/maphot_dm_selection_reg_cx_simP05_all.bed -b ~/meTRN/input/ensembl_ws220_nuclear_sizes.bed -hist
#0.0015460
#coverageBed -a ~/meTRN/data/hot/regions/maphot_dm_selection_reg_cx_simP01_any.bed -b ~/meTRN/input/ensembl_ws220_nuclear_sizes.bed -hist
#0.0010604
#coverageBed -a ~/meTRN/data/hot/regions/maphot_dm_selection_reg_cx_simP01_all.bed -b ~/meTRN/input/ensembl_ws220_nuclear_sizes.bed -hist
#0.0000895

#100*0.9689844
#RGB: 100*(0.0310155 - 0.0075204)
#HOT (any): 100*(0.0075204 - 0.0015460)
#HOT (all): 100*0.0015460
#XOT (any): 100*(0.0010604 - 0.0000895)
#XOT (all): 100*0.0000895


wc -l ~/meTRN/data/peaks/mappeaks_dm_selection_com_cx_raw*
#datasets: 62
#regions: 24758
#peaks: 85598

wc -l ~/meTRN/data/peaks/mappeaks_dm_selection_com_cx_hub*
#datasets: 62
#regions: 24600
#peaks: 77579
#percent: 77579/85598

wc -l ~/meTRN/data/peaks/mappeaks_dm_selection_com_cx_hct*
#datasets: 62
#regions: 25518
#peaks: 61062
#percent: 61062/85598

wc -l ~/meTRN/data/peaks/mappeaks_dm_selection_com_cx_hno*
#datasets: 62
#regions: 22764
#peaks: 51338
#percent: 51338/85598

wc -l ~/meTRN/data/peaks/mappeaks_dm_selection_com_cx_xub*
#datasets: 62
#regions: 24759
#peaks: 85449
#percent: 85449/85598

wc -l ~/meTRN/data/peaks/mappeaks_dm_selection_com_cx_xct*
#datasets: 62
#regions: 24844
#peaks: 83689
#percent: 83689/85598

wc -l ~/meTRN/data/peaks/mappeaks_dm_selection_com_cx_xno*
#datasets: 62
#regions: 24697
#peaks: 82160
#percent: 82160/85598



# Generate HOT-region Circos plots:
#python mapCircos.py --path ~/meTRN --mode karyotype --organism c.elegans --source ce1
#python mapCircos.py --path ~/meTRN --mode import --organism c.elegans --source ce1 --infile ~/meTRN/data/hot/analysis/maphot_peaks_simple_gffkde2_dm_selection_reg_sx_simple_cut38_conXX_bw300_cs1_cp00001_pl30_peaks_eetended_fail.bed --name hot_regions_simple
#python mapCircos.py --path ~/meTRN --mode import --organism c.elegans --source ce1 --infile ~/meTRN/data/hot/analysis/maphot_filter_select_sx_shared_peaks_basics_gffkde2_dm_selection_reg_sx_basics_cutXX_conXX_bw300_cs1_cp00001_pl30_peaks_eetended_pass.bed --name hot_regions_filter
#python mapCircos.py --path ~/meTRN --mode import --organism c.elegans --source ce1 --infile ~/meTRN/data/hot/overlap/maphot_overlap_simple_dm_selection_reg_sx_shared --name hot_regions_filter

#circos -conf ./circos.conf 


#top