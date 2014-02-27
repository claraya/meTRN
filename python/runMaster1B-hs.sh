#!/usr/bin/sh


# Download human data:
#python dataImporter.py --path ~/meTRN --mode download.peaks --source hs/blacklisted --server claraya@snively.stanford.edu --parameters peaks/complete/blacklist/hs/v2/*

# Import human data:
#python dataImporter.py --path ~/meTRN --mode import.peaks --source hs/blacklisted --peaks hs_optimal --rank ON --method OFF
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_optimal --peaks hs_standards_stn --include .none.stn
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_optimal --peaks hs_alternate_stn --include .none


# Generate non-redundant datasets for analysis (selection):
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_alternate_stn --peaks hs_alternate_com --include .none --nonredundant ON --label factor.context
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_alternate_com --peaks hs_alternate_reg --exclude POLR
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_alternate_com --peaks hs_alternate_pol --include POLR

#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_alternate_com --peaks hs_alternate_com_cx --include GM12878,H1-hESC,HepG2,HeLa-S3,K562
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_alternate_com_cx --peaks hs_alternate_com_gm --include GM12878
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_alternate_com_cx --peaks hs_alternate_com_h1 --include H1-hESC
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_alternate_com_cx --peaks hs_alternate_com_hg --include HepG2
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_alternate_com_cx --peaks hs_alternate_com_hl --include HeLa-S3
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_alternate_com_cx --peaks hs_alternate_com_k5 --include K562

#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_alternate_reg --peaks hs_alternate_reg_cx --include GM12878,H1-hESC,HepG2,HeLa-S3,K562
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_alternate_reg_cx --peaks hs_alternate_reg_gm --include GM12878
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_alternate_reg_cx --peaks hs_alternate_reg_h1 --include H1-hESC
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_alternate_reg_cx --peaks hs_alternate_reg_hg --include HepG2
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_alternate_reg_cx --peaks hs_alternate_reg_hl --include HeLa-S3
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_alternate_reg_cx --peaks hs_alternate_reg_k5 --include K562

#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_alternate_pol --peaks hs_alternate_pol_cx --include GM12878,H1-hESC,HepG2,HeLa-S3,K562
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_alternate_pol_cx --peaks hs_alternate_pol_gm --include GM12878
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_alternate_pol_cx --peaks hs_alternate_pol_h1 --include H1-hESC
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_alternate_pol_cx --peaks hs_alternate_pol_hg --include HepG2
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_alternate_pol_cx --peaks hs_alternate_pol_hl --include HeLa-S3
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_alternate_pol_cx --peaks hs_alternate_pol_k5 --include K562


# Transfer human data (selection):
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_alternate_com --peaks hs_alternate_com_xx_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_alternate_com_cx --peaks hs_alternate_com_cx_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_alternate_com_gm --peaks hs_alternate_com_gm_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_alternate_com_h1 --peaks hs_alternate_com_h1_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_alternate_com_hg --peaks hs_alternate_com_hg_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_alternate_com_hl --peaks hs_alternate_com_hl_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_alternate_com_k5 --peaks hs_alternate_com_k5_raw

python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_alternate_reg --peaks hs_alternate_reg_xx_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_alternate_reg_cx --peaks hs_alternate_reg_cx_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_alternate_reg_gm --peaks hs_alternate_reg_gm_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_alternate_reg_h1 --peaks hs_alternate_reg_h1_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_alternate_reg_hg --peaks hs_alternate_reg_hg_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_alternate_reg_hl --peaks hs_alternate_reg_hl_raw
python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_alternate_reg_k5 --peaks hs_alternate_reg_k5_raw



# Generate non-redundant datasets for analysis (standard):
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_stn --peaks hs_standards_com --include .none.stn --nonredundant ON --label factor.context
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_com --peaks hs_standards_reg --exclude POLR
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_com --peaks hs_standards_pol --include POLR

#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_com --peaks hs_standards_com_cx --include GM12878,H1-hESC,HepG2,HeLa-S3,K562
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_com_cx --peaks hs_standards_com_gm --include GM12878
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_com_cx --peaks hs_standards_com_h1 --include H1-hESC
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_com_cx --peaks hs_standards_com_hg --include HepG2
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_com_cx --peaks hs_standards_com_hl --include HeLa-S3
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_com_cx --peaks hs_standards_com_k5 --include K562

#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_reg --peaks hs_standards_reg_cx --include GM12878,H1-hESC,HepG2,HeLa-S3,K562
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_reg_cx --peaks hs_standards_reg_gm --include GM12878
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_reg_cx --peaks hs_standards_reg_h1 --include H1-hESC
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_reg_cx --peaks hs_standards_reg_hg --include HepG2
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_reg_cx --peaks hs_standards_reg_hl --include HeLa-S3
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_reg_cx --peaks hs_standards_reg_k5 --include K562

#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_pol --peaks hs_standards_pol_cx --include GM12878,H1-hESC,HepG2,HeLa-S3,K562
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_pol_cx --peaks hs_standards_pol_gm --include GM12878
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_pol_cx --peaks hs_standards_pol_h1 --include H1-hESC
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_pol_cx --peaks hs_standards_pol_hg --include HepG2
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_pol_cx --peaks hs_standards_pol_hl --include HeLa-S3
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_pol_cx --peaks hs_standards_pol_k5 --include K562


# Generate redundant datasets for testing (standard):
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_standards_stn --peaks hs_redundant_com --include .none.stn --label factor.context
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_redundant_com --peaks hs_redundant_reg --exclude POLR
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_redundant_com --peaks hs_redundant_pol --include POLR

#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_redundant_reg --peaks hs_redundant_reg_cx --include GM12878,H1-hESC,HepG2,HeLa-S3,K562
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_redundant_reg_cx --peaks hs_redundant_reg_gm --include GM12878
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_redundant_reg_cx --peaks hs_redundant_reg_h1 --include H1-hESC
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_redundant_reg_cx --peaks hs_redundant_reg_hg --include HepG2
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_redundant_reg_cx --peaks hs_redundant_reg_hl --include HeLa-S3
#python dataImporter.py --path ~/meTRN --mode select.peaks --source hs_redundant_reg_cx --peaks hs_redundant_reg_k5 --include K562


# Transfer human data (standard):
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_standards_com --peaks hs_standards_com_xx_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_standards_com_cx --peaks hs_standards_com_cx_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_standards_com_gm --peaks hs_standards_com_gm_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_standards_com_h1 --peaks hs_standards_com_h1_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_standards_com_hg --peaks hs_standards_com_hg_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_standards_com_hl --peaks hs_standards_com_hl_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_standards_com_k5 --peaks hs_standards_com_k5_raw

#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_standards_reg --peaks hs_standards_reg_xx_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_standards_reg_cx --peaks hs_standards_reg_cx_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_standards_reg_gm --peaks hs_standards_reg_gm_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_standards_reg_h1 --peaks hs_standards_reg_h1_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_standards_reg_hg --peaks hs_standards_reg_hg_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_standards_reg_hl --peaks hs_standards_reg_hl_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_standards_reg_k5 --peaks hs_standards_reg_k5_raw

#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_redundant_reg --peaks hs_redundant_reg_xx_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_redundant_reg_cx --peaks hs_redundant_reg_cx_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_redundant_reg_gm --peaks hs_redundant_reg_gm_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_redundant_reg_h1 --peaks hs_redundant_reg_h1_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_redundant_reg_hg --peaks hs_redundant_reg_hg_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_redundant_reg_hl --peaks hs_redundant_reg_hl_raw
#python dataImporter.py --path ~/meTRN --mode transfer.peaks --source hs_redundant_reg_k5 --peaks hs_redundant_reg_k5_raw


