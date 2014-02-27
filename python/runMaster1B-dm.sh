#!/usr/bin/sh

# Download worm data:
#python dataImporter.py --path ~/meTRN --mode download.peaks --source dm/blacklisted --server claraya@snively.stanford.edu --parameters peaks/complete/blacklist/dm/v2/*
#python dataImporter.py --path ~/meTRN --mode download.peaks --source dm/joined --server claraya@snively.stanford.edu --parameters peaks/joined/blacklist/dm/v2/EE/*
#python dataImporter.py --path ~/meTRN --mode download.peaks --source dm/joined --server claraya@snively.stanford.edu --parameters peaks/joined/blacklist/dm/v2/LE/*
#python dataImporter.py --path ~/meTRN --mode download.peaks --source dm/joined --server claraya@snively.stanford.edu --parameters peaks/joined/blacklist/dm/v2/PP/*

# Transfer the new PolII data to the joined folders (locally):
#cp dm_unk_RpII215_Embryos-16-24-hr_Harvard_stn_peaks.bed ../joined/
#cp dm_unk_RpII215_W3L_Harvard_stn_peaks.bed  ../joined/
#cd ../joined/
#mv dm_unk_mv dm_unk_RpII215_Embryos-16-24-hr_Harvard_stn_peaks.bed dm_na_RpII215_LE_Harvard_stn_peaks.bed
#mv dm_unk_RpII215_W3L_Harvard_stn_peaks.bed dm_na_RpII215_PP_Harvard_stn_peaks.bed

# Import fly data:
#python dataImporter.py --path ~/meTRN --mode import.peaks --source dm/joined --peaks dm_joined --rank ON --method OFF --cutChr ON --fixed ON --rename .bed.fixed:.bed --filterChr U,Uextra
#python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_joined --peaks dm_selection_stn --include _stn_ --exclude HDAC3_EE,HP1b_LE,suHw_PP


# Generate non-redundant datasets for analysis:
#python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_selection_stn --peaks dm_selection_com --include _stn_ --nonredundant ON --label factor.context
#python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_selection_com --peaks dm_selection_reg --exclude lola,RpII215
#python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_selection_com --peaks dm_selection_pol --include lola,RpII215

python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_selection_com --peaks dm_selection_com_cx --include _EE_,_LE_,_PP_
python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_selection_com_cx --peaks dm_selection_com_ee --include _EE_
python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_selection_com_cx --peaks dm_selection_com_le --include _LE_
python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_selection_com_cx --peaks dm_selection_com_pp --include _PP_

python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_selection_reg --peaks dm_selection_reg_cx --include _EE_,_LE_,_PP_
python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_selection_reg_cx --peaks dm_selection_reg_ee --include _EE_
python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_selection_reg_cx --peaks dm_selection_reg_le --include _LE_
python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_selection_reg_cx --peaks dm_selection_reg_pp --include _PP_

python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_selection_pol --peaks dm_selection_pol_cx --include _EE_,_LE_,_PP_
python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_selection_pol_cx --peaks dm_selection_pol_ee --include _EE_
python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_selection_pol_cx --peaks dm_selection_pol_le --include _LE_
python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_selection_pol_cx --peaks dm_selection_pol_pp --include _PP_


# Transfer fly data:
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source dm_selection_com --peaks dm_selection_com_xx_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source dm_selection_com_cx --peaks dm_selection_com_cx_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source dm_selection_com_ee --peaks dm_selection_com_ee_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source dm_selection_com_le --peaks dm_selection_com_le_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source dm_selection_com_pp --peaks dm_selection_com_pp_raw

python dataImporter.py --path ~/meTRN --mode transfer.peaks --source dm_selection_reg --peaks dm_selection_reg_xx_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source dm_selection_reg_cx --peaks dm_selection_reg_cx_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source dm_selection_reg_ee --peaks dm_selection_reg_ee_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source dm_selection_reg_le --peaks dm_selection_reg_le_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source dm_selection_reg_pp --peaks dm_selection_reg_pp_raw


# Generate extended peak sets:
#python dataImporter.py --path ~/meTRN --mode import.peaks --source dm/blacklisted --peaks dm_extension --rank ON --method OFF --cutChr ON --fixed ON --rename .bed.fixed:.bed --filterChr U,Uextra
#python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_extension --peaks dm_extension_stn --include _stn_

#python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_extension_stn --peaks dm_extension_com --include _stn_ --nonredundant ON --label factor.context
#python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_extension_com --peaks dm_extension_reg --exclude lola,RpII215
#python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_extension_com --peaks dm_extension_pol --include lola,RpII215

#python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_extension_com --peaks dm_extension_com_kc --include _Kc167_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_extension_com --peaks dm_extension_com_s2 --include _S2_

#python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_extension_reg --peaks dm_extension_reg_kc --include _Kc167_
#python dataImporter.py --path ~/meTRN --mode select.peaks --source dm_extension_reg --peaks dm_extension_reg_s2 --include _S2_


# Transfer fly data:
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source dm_extension_com --peaks dm_extension_com_xx_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source dm_extension_com_kc --peaks dm_extension_com_kc_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source dm_extension_com_s2 --peaks dm_extension_com_s2_raw

python dataImporter.py --path ~/meTRN --mode transfer.peaks --source dm_extension_reg --peaks dm_extension_reg_xx_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source dm_extension_reg_kc --peaks dm_extension_reg_kc_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source dm_extension_reg_s2 --peaks dm_extension_reg_s2_raw