#!/usr/bin/sh


# Launch TF binding signal analyses near TSSs (WS220):
#python mapProfile.py --path ~/meTRN --organism ce --mode signal --group regulation --peaks ce_selection_com_cx_raw --name ce_wormbased_TSS_gx --window 1000 --source annotations --infile in2shape_ce_wormbased_TSS_gx.bed  --strand ON --cutChr ON --chunks 1 --overwrite OFF
#python mapProfile.py --path ~/meTRN --organism ce --mode matrix --group regulation --peaks ce_selection_com_cx_raw --name ce_wormbased_TSS_gx --window 1000 --target factor.context --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6
#python mapProfile.py --path ~/meTRN --organism ce --mode deltas --group regulation --peaks ce_selection_com_cx_raw --name ce_wormbased_TSS_gx --window 1000 --target factor.context --metric signal.rng


# Launch TF binding signal analyses near EM and L3 TSSs (Waterston Lab):
#python mapProfile.py --path ~/meTRN --organism ce --mode signal --group regulation --peaks ce_selection_pol_ex_raw --name ce_modencode_TSS_ex --window 1000 --source annotations --infile in2shape_ce_modencode_TSS_ee.bed  --strand ON --cutChr ON --chunks 1 --overwrite OFF --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6
#python mapProfile.py --path ~/meTRN --organism ce --mode signal --group regulation --peaks ce_selection_pol_l3_raw --name ce_modencode_TSS_l3 --window 1000 --source annotations --infile in2shape_ce_modencode_TSS_l3.bed  --strand ON --cutChr ON --chunks 1 --overwrite OFF --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6

#python mapProfile.py --path ~/meTRN --organism ce --mode matrix --group regulation --peaks ce_selection_pol_ex_raw --name ce_modencode_TSS_ex --window 1000 --target factor.context --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6 
#python mapProfile.py --path ~/meTRN --organism ce --mode matrix --group regulation --peaks ce_selection_pol_l3_raw --name ce_modencode_TSS_l3 --window 1000 --target factor.context --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6 


# Launch TF binding signal analyses near CapSeq/CIP-TAP TSSs (Gu et al. 2012, Cell):
#python mapProfile.py --path ~/meTRN --organism ce --mode signal --group regulation --peaks ce_selection_com_cx_raw --name ce_gu2012pro_TSS_gx --window 1000 --source annotations --infile in2shape_ce_gu2012pro_TSS_gx.bed  --strand ON --cutChr ON --chunks 1 --overwrite OFF --precision 3 --multiarray OFF
#python mapProfile.py --path ~/meTRN --organism ce --mode matrix --group regulation --peaks ce_selection_com_cx_raw --name ce_gu2012pro_TSS_gx --window 1000 --target factor.context --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6
#python mapProfile.py --path ~/meTRN --organism ce --mode deltas --group regulation --peaks ce_selection_com_cx_raw --name ce_gu2012pro_TSS_gx --window 1000 --target factor.context --metric signal.rng
#python mapProfile.py --path ~/meTRN --organism ce --mode ratios --group regulation --peaks ce_selection_com_cx_raw --name ce_gu2012pro_TSS_gx --window 1000 --target factor.context --metric signal.rng


# Import chromatin signal files to bedGraph:
#python mapProfile.py --path ~/meTRN --organism ce --mode import --group chromatin --source extras/chromatin/ce/ee/fold --overwrite OFF --rename DPL1:DPL-1,DPY27:DPY-27,EFL1:EFL-1,EPC1:EPC-1,HDA1:HDA-1,HPL2:HPL-2,LET418:LET-418,LIN35:LIN-35,LIN37:LIN-37,LIN52:LIN-52,LIN54:LIN-54,LIN61:LIN-61,LIN9:LIN-9,NURF1:NURF-1,MRG1:MRG-1,RPC1:RPC-1,PolII:AMA-1

#python mapProfile.py --path ~/meTRN --organism ce --mode import --group chromatin --source extras/chromatin/ce/l3/fold --overwrite OFF --rename DPL1:DPL-1,DPY27:DPY-27,EFL1:EFL-1,EPC1:EPC-1,HDA1:HDA-1,HPL2:HPL-2,LET418:LET-418,LIN35:LIN-35,LIN37:LIN-37,LIN52:LIN-52,LIN54:LIN-54,LIN61:LIN-61,LIN9:LIN-9,NURF1:NURF-1,MRG1:MRG-1,RPC1:RPC-1,PolII:AMA-1


# Launch chromatin signal analyses near TSSs:
#python mapProfile.py --path ~/meTRN --organism ce --mode signal --group chromatin --name ce_chromatin_TSS_gx --window 1000 --source annotations --infile in2shape_ce_wormbased_TSS_gx.bed --strand ON --cutChr ON --chunks 1 --overwrite OFF
#python mapProfile.py --path ~/meTRN --organism ce --mode matrix --group chromatin --name ce_chromatin_TSS_gx --window 1000 --target factor.context --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6



# Cluster TF binding signals in select examples:
#python mapProfile.py --path ~/meTRN --organism ce --mode cluster --name ce_wormbased_TSS_gx --infile mapprofile_signal_N2_POL2_L2_yale_stn.fc.signal.txt --target factor.context




#top
#bash runMaster8A.sh