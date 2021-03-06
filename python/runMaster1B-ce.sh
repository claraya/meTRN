#!/usr/bin/sh

# Download worm data:
#python dataImporter.py --path ~/meTRN --mode download.peaks --source ce/blacklisted/v1/ --server claraya@snively.stanford.edu --parameters peaks/complete/blacklist/ce/v2/*


# Add new datasets (n=3):
#scp claraya@snively.stanford.edu:/srv/gs1/projects/snyder/modENCODE/data/peaks/complete/blacklist/ce/*bed ~/meTRN/idr/final/ce/blacklisted/v2/


# Import worm data:
#python dataImporter.py --path ~/meTRN --mode import.peaks --source ce/blacklisted/v1 --peaks ce_optimal --rank ON --method OFF --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6,RPC-15:RPC-1 --cutChr ON
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_optimal --peaks ce_selection_stn --include _stn_ --exclude ce_XE,ce_YL

#python dataImporter.py --path ~/meTRN --mode import.peaks --source ce/blacklisted/v2 --peaks ce_updated --rank ON --method OFF --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6,RPC-15:RPC-1 --cutChr ON
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_updated --peaks ce_extension_stn --include _stn_ --exclude ce_XE,ce_YL

#python dataImporter.py --path ~/meTRN --mode import.peaks --source ce/blacklisted/v2 --peaks ce_correct --rank ON --method OFF --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6,RPC-15:RPC-1,DCC4640_LIN-35:DCC4640_EFL-1,DCC_HPL-2:SDQ2340_HPL-2,DCC_HDA-1:SDQ2354_HDA-1 --cutChr ON

#python dataImporter.py --path ~/meTRN --mode import.peaks --source ce/blacklisted/v2 --peaks ce_release --rank ON --method OFF --rename POL2:AMA-1,POLIII:RPC-1,C01312.1:C01G12.1,R0GF6.6:R06F6.6,RPC-15:RPC-1,DCC4640_LIN-35:DCC4640_EFL-1,DCC_HPL-2:SDQ2340_HPL-2,DCC_HDA-1:SDQ2354_HDA-1,HAM-1:MEP-1,CEH-28:CEH-38 --cutChr ON


# These are datasets that won't be included for release because we don't have the full data:
#DCC_HPL-2:SDQ2340_HPL-2
#DCC_HDA-1:SDQ2354_HDA-1

# Generate quality-approved datasets for reporting:
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_release --peaks ce_reporting_com --nonredundant OFF --label factor.context.method --exclude DCC3217_W03F9.2_YA_modencode_stn,DCC4033_PHA-4_L4_modencode_stn,DCC3160_EOR-1_L3_modencode_stn,DCC4027_FKH-10_L3_modencode_stn,DCC4028_FKH-10_L4_modencode_stn,OP177_EGL-27_L1_yale,OP184_LIN-15B_L3_stanford,CEH-39_EM,NHR-2_EM,PAX-1_EM
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_reporting_com --peaks ce_reporting_com_xx_raw


# Generate non-redundant datasets for analysis:
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_release --peaks ce_finalized_com --nonredundant ON --label factor.context.method --exclude DCC3217_W03F9.2_YA_modencode_stn,DCC4033_PHA-4_L4_modencode_stn,DCC3160_EOR-1_L3_modencode_stn,DCC4027_FKH-10_L3_modencode_stn,DCC4028_FKH-10_L4_modencode_stn,OP177_EGL-27_L1_yale,OP184_LIN-15B_L3_yale,CEH-39_EM,NHR-2_EM,PAX-1_EM,ce_XE,ce_YL,DCC_,TBP-1,SDQ2340_HPL-2,SDQ2354_HDA-1 --include _stn_ --allow YL409

#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_finalized_com --peaks ce_finalized_com_cx --include _EM_,_EE_,_LE_,_L1_,_L2_,_L3_,_L4_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_finalized_com_cx --peaks ce_finalized_com_ex --include _EM_,_EE_,_LE_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_finalized_com_cx --peaks ce_finalized_com_l1 --include _L1_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_finalized_com_cx --peaks ce_finalized_com_l2 --include _L2_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_finalized_com_cx --peaks ce_finalized_com_l3 --include _L3_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_finalized_com_cx --peaks ce_finalized_com_l4 --include _L4_

#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_finalized_com_cx --peaks ce_finalized_com_cx_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_finalized_com_ex --peaks ce_finalized_com_ex_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_finalized_com_l1 --peaks ce_finalized_com_l1_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_finalized_com_l2 --peaks ce_finalized_com_l2_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_finalized_com_l3 --peaks ce_finalized_com_l3_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_finalized_com_l4 --peaks ce_finalized_com_l4_raw


# Generate non-redundant datasets for analysis:
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_correct --peaks ce_corrected_com --nonredundant ON --label factor.context.method --exclude DCC3217_W03F9.2_YA_modencode_stn,DCC4033_PHA-4_L4_modencode_stn,DCC3160_EOR-1_L3_modencode_stn,DCC4027_FKH-10_L3_modencode_stn,DCC4028_FKH-10_L4_modencode_stn,OP177_EGL-27_L1_yale,OP184_LIN-15B_L3_yale,CEH-39_EM,NHR-2_EM,PAX-1_EM,ce_XE,ce_YL,DCC_,TBP-1,SDQ2340_HPL-2,SDQ2354_HDA-1 --include _stn_ --allow YL409

#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_corrected_com --peaks ce_corrected_com_cx --include _EM_,_EE_,_LE_,_L1_,_L2_,_L3_,_L4_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_corrected_com_cx --peaks ce_corrected_com_ex --include _EM_,_EE_,_LE_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_corrected_com_cx --peaks ce_corrected_com_l1 --include _L1_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_corrected_com_cx --peaks ce_corrected_com_l2 --include _L2_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_corrected_com_cx --peaks ce_corrected_com_l3 --include _L3_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_corrected_com_cx --peaks ce_corrected_com_l4 --include _L4_

#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_corrected_com_cx --peaks ce_corrected_com_cx_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_corrected_com_ex --peaks ce_corrected_com_ex_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_corrected_com_l1 --peaks ce_corrected_com_l1_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_corrected_com_l2 --peaks ce_corrected_com_l2_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_corrected_com_l3 --peaks ce_corrected_com_l3_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_corrected_com_l4 --peaks ce_corrected_com_l4_raw


# Generate non-redundant datasets for analysis:
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_selection_stn --peaks ce_selection_com --include _stn_ --nonredundant ON --label factor.context --exclude DCC3217_W03F9.2_YA_modencode_stn,DCC4033_PHA-4_L4_modencode_stn,DCC3160_EOR-1_L3_modencode_stn,CEH-39_EM,NHR-2_EM,PAX-1_EM
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_selection_com --peaks ce_selection_reg --exclude AMA-1,RPC-1_,RPC-15
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_selection_com --peaks ce_selection_pol --include AMA-1,RPC-1_,RPC-15

#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_selection_com --peaks ce_selection_com_cx --include _EM_,_EE_,_LE_,_L1_,_L2_,_L3_,_L4_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_selection_com_cx --peaks ce_selection_com_ex --include _EM_,_EE_,_LE_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_selection_com_cx --peaks ce_selection_com_l1 --include _L1_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_selection_com_cx --peaks ce_selection_com_l2 --include _L2_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_selection_com_cx --peaks ce_selection_com_l3 --include _L3_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_selection_com_cx --peaks ce_selection_com_l4 --include _L4_

#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_selection_reg --peaks ce_selection_reg_cx --include _EM_,_EE_,_LE_,_L1_,_L2_,_L3_,_L4_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_selection_reg_cx --peaks ce_selection_reg_ex --include _EM_,_EE_,_LE_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_selection_reg_cx --peaks ce_selection_reg_l1 --include _L1_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_selection_reg_cx --peaks ce_selection_reg_l2 --include _L2_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_selection_reg_cx --peaks ce_selection_reg_l3 --include _L3_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_selection_reg_cx --peaks ce_selection_reg_l4 --include _L4_

#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_selection_pol --peaks ce_selection_pol_cx --include _EM_,_EE_,_LE_,_L1_,_L2_,_L3_,_L4_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_selection_pol_cx --peaks ce_selection_pol_ex --include _EM_,_EE_,_LE_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_selection_pol_cx --peaks ce_selection_pol_l1 --include _L1_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_selection_pol_cx --peaks ce_selection_pol_l2 --include _L2_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_selection_pol_cx --peaks ce_selection_pol_l3 --include _L3_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_selection_pol_cx --peaks ce_selection_pol_l4 --include _L4_


#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_extension_stn --peaks ce_extension_com --include _stn_ --nonredundant ON --label factor.context --exclude DCC3217_W03F9.2_YA_modencode_stn,DCC4033_PHA-4_L4_modencode_stn,DCC3160_EOR-1_L3_modencode_stn,CEH-39_EM,NHR-2_EM,PAX-1_EM
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_extension_com --peaks ce_extension_reg --exclude AMA-1,RPC-1_,RPC-15
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_extension_com --peaks ce_extension_pol --include AMA-1,RPC-1_,RPC-15

#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_extension_com --peaks ce_extension_com_cx --include _EM_,_EE_,_LE_,_L1_,_L2_,_L3_,_L4_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_extension_com_cx --peaks ce_extension_com_ex --include _EM_,_EE_,_LE_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_extension_com_cx --peaks ce_extension_com_l1 --include _L1_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_extension_com_cx --peaks ce_extension_com_l2 --include _L2_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_extension_com_cx --peaks ce_extension_com_l3 --include _L3_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_extension_com_cx --peaks ce_extension_com_l4 --include _L4_

#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_extension_reg --peaks ce_extension_reg_cx --include _EM_,_EE_,_LE_,_L1_,_L2_,_L3_,_L4_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_extension_reg_cx --peaks ce_extension_reg_ex --include _EM_,_EE_,_LE_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_extension_reg_cx --peaks ce_extension_reg_l1 --include _L1_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_extension_reg_cx --peaks ce_extension_reg_l2 --include _L2_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_extension_reg_cx --peaks ce_extension_reg_l3 --include _L3_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source ce_extension_reg_cx --peaks ce_extension_reg_l4 --include _L4_


# Transfer worm data:
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_selection_com --peaks ce_selection_com_xx_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_selection_com_cx --peaks ce_selection_com_cx_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_selection_com_ex --peaks ce_selection_com_ex_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_selection_com_l1 --peaks ce_selection_com_l1_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_selection_com_l2 --peaks ce_selection_com_l2_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_selection_com_l3 --peaks ce_selection_com_l3_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_selection_com_l4 --peaks ce_selection_com_l4_raw

#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_selection_reg --peaks ce_selection_reg_xx_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_selection_reg_cx --peaks ce_selection_reg_cx_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_selection_reg_ex --peaks ce_selection_reg_ex_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_selection_reg_l1 --peaks ce_selection_reg_l1_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_selection_reg_l2 --peaks ce_selection_reg_l2_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_selection_reg_l3 --peaks ce_selection_reg_l3_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_selection_reg_l4 --peaks ce_selection_reg_l4_raw

#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_selection_pol --peaks ce_selection_pol_xx_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_selection_pol_cx --peaks ce_selection_pol_cx_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_selection_pol_ex --peaks ce_selection_pol_ex_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_selection_pol_l1 --peaks ce_selection_pol_l1_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_selection_pol_l2 --peaks ce_selection_pol_l2_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_selection_pol_l3 --peaks ce_selection_pol_l3_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_selection_pol_l4 --peaks ce_selection_pol_l4_raw


#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_extension_com --peaks ce_extension_com_xx_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_extension_com_cx --peaks ce_extension_com_cx_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_extension_com_ex --peaks ce_extension_com_ex_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_extension_com_l1 --peaks ce_extension_com_l1_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_extension_com_l2 --peaks ce_extension_com_l2_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_extension_com_l3 --peaks ce_extension_com_l3_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_extension_com_l4 --peaks ce_extension_com_l4_raw

#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_extension_reg --peaks ce_extension_reg_xx_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_extension_reg_cx --peaks ce_extension_reg_cx_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_extension_reg_ex --peaks ce_extension_reg_ex_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_extension_reg_l1 --peaks ce_extension_reg_l1_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_extension_reg_l2 --peaks ce_extension_reg_l2_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_extension_reg_l3 --peaks ce_extension_reg_l3_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source ce_extension_reg_l4 --peaks ce_extension_reg_l4_raw


#bash runMaster1B-ce.sh