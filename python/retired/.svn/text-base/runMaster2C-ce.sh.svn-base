#!/usr/bin/sh


# Map HOT regions (fail) and non-HOT regions (pass); without any actual filtering cutoffs:
python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_cx --name basics
python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_ex --name basics --tag EX-
python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l1 --name basics --tag L1-
python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l2 --name basics --tag L2-
python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l3 --name basics --tag L3-
python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l4 --name basics --tag L4-

	
# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (0.01 significance):
python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_ex --name simP01 --cutoff 9 --tag EX-
python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l1 --name simP01 --cutoff 18 --tag L1-
python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l2 --name simP01 --cutoff 18 --tag L2-
python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l3 --name simP01 --cutoff 12 --tag L3-
python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l4 --name simP01 --cutoff 16 --tag L4-


# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (0.05 significance):
python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_ex --name simP05 --cutoff 7 --tag EX-
python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l1 --name simP05 --cutoff 12 --tag L1-
python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l2 --name simP05 --cutoff 12 --tag L2-
python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l3 --name simP05 --cutoff 8 --tag L3-
python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l4 --name simP05 --cutoff 11 --tag L4-


# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (0.10 significance):
python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_ex --name simP10 --cutoff 5 --tag EX-
python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l1 --name simP10 --cutoff 10 --tag L1-
python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l2 --name simP10 --cutoff 9 --tag L2-
python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l3 --name simP10 --cutoff 6 --tag L3-
python mapHOT.py --path ~/meTRN --organism ce --mode scan --target peaks --peaks ce_selection_reg_l4 --name simP10 --cutoff 8 --tag L4-


# Calculate overlap in HOT regions and filter regions against the ubiquitously HOT regions:
python mapHOT.py --path ~/meTRN --organism ce --mode overlap --name simP01 --target EX:ce_selection_reg_ex,L1:ce_selection_reg_l1,L2:ce_selection_reg_l2,L3:ce_selection_reg_l3,L4:ce_selection_reg_l4 --overlap ce_selection_reg_cx

python mapHOT.py --path ~/meTRN --organism ce --mode overlap --name simP05 --target EX:ce_selection_reg_ex,L1:ce_selection_reg_l1,L2:ce_selection_reg_l2,L3:ce_selection_reg_l3,L4:ce_selection_reg_l4 --overlap ce_selection_reg_cx


# Store HOT and RGB region files for each context and the global context file:
python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_ex_simP01 --regions ce_selection_reg_ex_simP01 --source analysis
python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l1_simP01 --regions ce_selection_reg_l1_simP01 --source analysis
python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l2_simP01 --regions ce_selection_reg_l2_simP01 --source analysis
python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l3_simP01 --regions ce_selection_reg_l3_simP01 --source analysis
python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l4_simP01 --regions ce_selection_reg_l4_simP01 --source analysis
python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_cx_simP01 --regions ce_selection_reg_cx_simP01 --source overlap

python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_ex_simP05 --regions ce_selection_reg_ex_simP05 --source analysis
python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l1_simP05 --regions ce_selection_reg_l1_simP05 --source analysis
python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l2_simP05 --regions ce_selection_reg_l2_simP05 --source analysis
python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l3_simP05 --regions ce_selection_reg_l3_simP05 --source analysis
python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_l4_simP05 --regions ce_selection_reg_l4_simP05 --source analysis
python mapHOT.py --path ~/meTRN --organism ce --mode regions --name ce_selection_reg_cx_simP05 --regions ce_selection_reg_cx_simP05 --source overlap


# Generate HOT-filtered peak sets (P01 cutoff):
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_cx_xno --infile maphot_ce_selection_reg_cx_simP01_any.bed --source ~/meTRN/data/peaks/ce_selection_reg_cx_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_ex_xct --infile maphot_ce_selection_reg_ex_simP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_ex_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l1_xct --infile maphot_ce_selection_reg_l1_simP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l1_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l2_xct --infile maphot_ce_selection_reg_l2_simP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l2_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l3_xct --infile maphot_ce_selection_reg_l3_simP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l3_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l4_xct --infile maphot_ce_selection_reg_l4_simP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l4_raw

python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_cx_xub --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_cx_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_ex_xub --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_ex_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l1_xub --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l1_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l2_xub --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l2_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l3_xub --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l3_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l4_xub --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l4_raw

python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_cx_xno --infile maphot_ce_selection_reg_cx_simP01_any.bed --source ~/meTRN/data/peaks/ce_selection_com_cx_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_ex_xct --infile maphot_ce_selection_reg_ex_simP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_ex_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l1_xct --infile maphot_ce_selection_reg_l1_simP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l1_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l2_xct --infile maphot_ce_selection_reg_l2_simP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l2_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l3_xct --infile maphot_ce_selection_reg_l3_simP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l3_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l4_xct --infile maphot_ce_selection_reg_l4_simP01_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l4_raw

python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_cx_xub --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_cx_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_ex_xub --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_ex_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l1_xub --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l1_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l2_xub --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l2_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l3_xub --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l3_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l4_xub --infile maphot_ce_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l4_raw

rm -rf ~/meTRN/data/peaks/ce_selection_reg_cx_xct/
mkdir ~/meTRN/data/peaks/ce_selection_reg_cx_xct/
cp ~/meTRN/data/peaks/ce_selection_reg_ex_xct/* ~/meTRN/data/peaks/ce_selection_reg_cx_xct/
cp ~/meTRN/data/peaks/ce_selection_reg_l1_xct/* ~/meTRN/data/peaks/ce_selection_reg_cx_xct/
cp ~/meTRN/data/peaks/ce_selection_reg_l2_xct/* ~/meTRN/data/peaks/ce_selection_reg_cx_xct/
cp ~/meTRN/data/peaks/ce_selection_reg_l3_xct/* ~/meTRN/data/peaks/ce_selection_reg_cx_xct/
cp ~/meTRN/data/peaks/ce_selection_reg_l4_xct/* ~/meTRN/data/peaks/ce_selection_reg_cx_xct/

rm -rf ~/meTRN/data/peaks/ce_selection_com_cx_xct/
mkdir ~/meTRN/data/peaks/ce_selection_com_cx_xct/
cp ~/meTRN/data/peaks/ce_selection_com_ex_xct/* ~/meTRN/data/peaks/ce_selection_com_cx_xct/
cp ~/meTRN/data/peaks/ce_selection_com_l1_xct/* ~/meTRN/data/peaks/ce_selection_com_cx_xct/
cp ~/meTRN/data/peaks/ce_selection_com_l2_xct/* ~/meTRN/data/peaks/ce_selection_com_cx_xct/
cp ~/meTRN/data/peaks/ce_selection_com_l3_xct/* ~/meTRN/data/peaks/ce_selection_com_cx_xct/
cp ~/meTRN/data/peaks/ce_selection_com_l4_xct/* ~/meTRN/data/peaks/ce_selection_com_cx_xct/


# Generate HOT-filtered peak sets (P05 cutoff):
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_cx_hno --infile maphot_ce_selection_reg_cx_simP05_any.bed --source ~/meTRN/data/peaks/ce_selection_reg_cx_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_ex_hct --infile maphot_ce_selection_reg_ex_simP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_ex_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l1_hct --infile maphot_ce_selection_reg_l1_simP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l1_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l2_hct --infile maphot_ce_selection_reg_l2_simP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l2_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l3_hct --infile maphot_ce_selection_reg_l3_simP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l3_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l4_hct --infile maphot_ce_selection_reg_l4_simP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_reg_l4_raw

python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_cx_hub --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_cx_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_ex_hub --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_ex_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l1_hub --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l1_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l2_hub --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l2_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l3_hub --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l3_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_reg_l4_hub --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_reg_l4_raw

python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_cx_hno --infile maphot_ce_selection_reg_cx_simP05_any.bed --source ~/meTRN/data/peaks/ce_selection_com_cx_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_ex_hct --infile maphot_ce_selection_reg_ex_simP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_ex_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l1_hct --infile maphot_ce_selection_reg_l1_simP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l1_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l2_hct --infile maphot_ce_selection_reg_l2_simP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l2_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l3_hct --infile maphot_ce_selection_reg_l3_simP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l3_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l4_hct --infile maphot_ce_selection_reg_l4_simP05_hot.bed --source ~/meTRN/data/peaks/ce_selection_com_l4_raw

python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_cx_hub --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_cx_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_ex_hub --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_ex_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l1_hub --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l1_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l2_hub --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l2_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l3_hub --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l3_raw
python mapHOT.py --path ~/meTRN --organism ce --mode filter:remove --peaks ce_selection_com_l4_hub --infile maphot_ce_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/ce_selection_com_l4_raw

rm -rf ~/meTRN/data/peaks/ce_selection_reg_cx_hct/
mkdir ~/meTRN/data/peaks/ce_selection_reg_cx_hct/
cp ~/meTRN/data/peaks/ce_selection_reg_ex_hct/* ~/meTRN/data/peaks/ce_selection_reg_cx_hct/
cp ~/meTRN/data/peaks/ce_selection_reg_l1_hct/* ~/meTRN/data/peaks/ce_selection_reg_cx_hct/
cp ~/meTRN/data/peaks/ce_selection_reg_l2_hct/* ~/meTRN/data/peaks/ce_selection_reg_cx_hct/
cp ~/meTRN/data/peaks/ce_selection_reg_l3_hct/* ~/meTRN/data/peaks/ce_selection_reg_cx_hct/
cp ~/meTRN/data/peaks/ce_selection_reg_l4_hct/* ~/meTRN/data/peaks/ce_selection_reg_cx_hct/

rm -rf ~/meTRN/data/peaks/ce_selection_com_cx_hct/
mkdir ~/meTRN/data/peaks/ce_selection_com_cx_hct/
cp ~/meTRN/data/peaks/ce_selection_com_ex_hct/* ~/meTRN/data/peaks/ce_selection_com_cx_hct/
cp ~/meTRN/data/peaks/ce_selection_com_l1_hct/* ~/meTRN/data/peaks/ce_selection_com_cx_hct/
cp ~/meTRN/data/peaks/ce_selection_com_l2_hct/* ~/meTRN/data/peaks/ce_selection_com_cx_hct/
cp ~/meTRN/data/peaks/ce_selection_com_l3_hct/* ~/meTRN/data/peaks/ce_selection_com_cx_hct/
cp ~/meTRN/data/peaks/ce_selection_com_l4_hct/* ~/meTRN/data/peaks/ce_selection_com_cx_hct/


# Generate complete, collapsed, density and report files:
#python mapPeaks.py --path ~/meTRN --mode build --overwrite OFF

#coverageBed -a ~/meTRN/data/peaks/mappeaks_ce_selection_com_cx_raw_compiled.bed -b ~/meTRN/input/ensembl_ws220_nuclear_sizes.bed -hist
#0.2170808 (0.7829192)
#coverageBed -a ~/meTRN/data/hot/regions/maphot_ce_selection_reg_cx_simP05_any.bed -b ~/meTRN/input/ensembl_ws220_nuclear_sizes.bed -hist
#0.0261704
#coverageBed -a ~/meTRN/data/hot/regions/maphot_ce_selection_reg_cx_simP05_all.bed -b ~/meTRN/input/ensembl_ws220_nuclear_sizes.bed -hist
#0.0057849
#coverageBed -a ~/meTRN/data/hot/regions/maphot_ce_selection_reg_cx_simP01_any.bed -b ~/meTRN/input/ensembl_ws220_nuclear_sizes.bed -hist
#0.0152395
#coverageBed -a ~/meTRN/data/hot/regions/maphot_ce_selection_reg_cx_simP01_all.bed -b ~/meTRN/input/ensembl_ws220_nuclear_sizes.bed -hist
#0.0030759

#100*0.7829192
#RGB: 100*(0.2170808 - 0.0261704)
#HOT (any): 100*(0.0261704 - 0.0057849)
#HOT (all): 100*0.0057849
#XOT (any): 100*(0.0152395 - 0.0030759)
#XOT (all): 100*0.0030759


wc -l ~/meTRN/data/peaks/mappeaks_ce_selection_com_cx_raw*
#datasets: 188
#regions: 33833
#peaks: 397539

wc -l ~/meTRN/data/peaks/mappeaks_ce_selection_com_cx_hub*
#datasets: 188
#regions: 33365
#peaks: 306632
#percent: 306632/397539

wc -l ~/meTRN/data/peaks/mappeaks_ce_selection_com_cx_hct*
#datasets: 188
#regions: 33823
#peaks: 230325
#percent: 230325/397539

wc -l ~/meTRN/data/peaks/mappeaks_ce_selection_com_cx_hno*
#datasets: 188
#regions: 31638
#peaks: 186972
#percent: 186972/397539

wc -l ~/meTRN/data/peaks/mappeaks_ce_selection_com_cx_xub*
#datasets: 188
#regions: 33626
#peaks: 342492
#percent: 342492/397539

wc -l ~/meTRN/data/peaks/mappeaks_ce_selection_com_cx_xct*
#datasets: 188
#regions: 33914
#peaks: 292466
#percent: 292466/397539

wc -l ~/meTRN/data/peaks/mappeaks_ce_selection_com_cx_xno*
#datasets: 188
#regions: 32720
#peaks: 250953
#percent: 250953/397539



# Generate HOT-region Circos plots:
#python mapCircos.py --path ~/meTRN --mode karyotype --organism c.elegans --source ce1
#python mapCircos.py --path ~/meTRN --mode import --organism c.elegans --source ce1 --infile ~/meTRN/data/hot/analysis/maphot_peaks_simple_gffkde2_ce_selection_reg_sx_simple_cut38_conXX_bw300_cs1_cp00001_pl30_peaks_extended_fail.bed --name hot_regions_simple
#python mapCircos.py --path ~/meTRN --mode import --organism c.elegans --source ce1 --infile ~/meTRN/data/hot/analysis/maphot_filter_select_sx_shared_peaks_basics_gffkde2_ce_selection_reg_sx_basics_cutXX_conXX_bw300_cs1_cp00001_pl30_peaks_extended_pass.bed --name hot_regions_filter
#python mapCircos.py --path ~/meTRN --mode import --organism c.elegans --source ce1 --infile ~/meTRN/data/hot/overlap/maphot_overlap_simple_ce_selection_reg_sx_shared --name hot_regions_filter

#circos -conf ./circos.conf 


#top