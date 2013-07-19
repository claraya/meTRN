#!/usr/bin/sh

# Map expressed and repressed genes per stage:
#python mapExpression.py --path ~/meTRN --organism ce --mode parse --infile in2shape_ce_expressed_RNA_ee.bed --folder annotations --cutoff 0 --genes in2shape_ce_wormbased_GEN_gx.bed --transcripts in2shape_ce_wormbased_RNA_gx.bed --coord TSS --name in2shape_ce_RENAME_TSS_ee.bed
#python mapExpression.py --path ~/meTRN --organism ce --mode parse --infile in2shape_ce_expressed_RNA_le.bed --folder annotations --cutoff 0 --genes in2shape_ce_wormbased_GEN_gx.bed --transcripts in2shape_ce_wormbased_RNA_gx.bed --coord TSS --name in2shape_ce_RENAME_TSS_le.bed
#python mapExpression.py --path ~/meTRN --organism ce --mode parse --infile in2shape_ce_expressed_RNA_l1.bed --folder annotations --cutoff 0 --genes in2shape_ce_wormbased_GEN_gx.bed --transcripts in2shape_ce_wormbased_RNA_gx.bed --coord TSS --name in2shape_ce_RENAME_TSS_l1.bed
#python mapExpression.py --path ~/meTRN --organism ce --mode parse --infile in2shape_ce_expressed_RNA_l2.bed --folder annotations --cutoff 0 --genes in2shape_ce_wormbased_GEN_gx.bed --transcripts in2shape_ce_wormbased_RNA_gx.bed --coord TSS --name in2shape_ce_RENAME_TSS_l2.bed
#python mapExpression.py --path ~/meTRN --organism ce --mode parse --infile in2shape_ce_expressed_RNA_l3.bed --folder annotations --cutoff 0 --genes in2shape_ce_wormbased_GEN_gx.bed --transcripts in2shape_ce_wormbased_RNA_gx.bed --coord TSS --name in2shape_ce_RENAME_TSS_l3.bed


# Launch co-association analysis near TSSs with IntervalStats (SCG3):
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_selection_com_ex_xot --source annotations --domain in2shape_ce_expressed_TSS_ee_slopbed_up2000_dn200.bed --name promoter_express --qsub scg3_intervallight.sh --server ON --job expressionStatsEEE
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_selection_com_ex_xot --source annotations --domain in2shape_ce_repressed_TSS_ee_slopbed_up2000_dn200.bed --name promoter_repress --qsub scg3_intervallight.sh --server ON --job expressionStatsEER
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_selection_com_ex_xot --source annotations --domain in2shape_ce_reference_TSS_ee_slopbed_up2000_dn200.bed --name promoter_observe --qsub scg3_intervallight.sh --server ON --job expressionStatsEET

#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_selection_com_ex_xot --source annotations --domain in2shape_ce_expressed_TSS_le_slopbed_up2000_dn200.bed --name promoter_express --qsub scg3_intervallight.sh --server ON --job expressionStatsLEE
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_selection_com_ex_xot --source annotations --domain in2shape_ce_repressed_TSS_le_slopbed_up2000_dn200.bed --name promoter_repress --qsub scg3_intervallight.sh --server ON --job expressionStatsLER
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_selection_com_ex_xot --source annotations --domain in2shape_ce_reference_TSS_le_slopbed_up2000_dn200.bed --name promoter_observe --qsub scg3_intervallight.sh --server ON --job expressionStatsLET

#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_selection_com_l1_xot --source annotations --domain in2shape_ce_expressed_TSS_l1_slopbed_up2000_dn200.bed --name promoter_express --qsub scg3_intervallight.sh --server ON --job expressionStatsL1E
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_selection_com_l1_xot --source annotations --domain in2shape_ce_repressed_TSS_l1_slopbed_up2000_dn200.bed --name promoter_repress --qsub scg3_intervallight.sh --server ON --job expressionStatsL1R
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_selection_com_l1_xot --source annotations --domain in2shape_ce_reference_TSS_l1_slopbed_up2000_dn200.bed --name promoter_observe --qsub scg3_intervallight.sh --server ON --job expressionStatsL1T

#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_selection_com_l2_xot --source annotations --domain in2shape_ce_expressed_TSS_l3_slopbed_up2000_dn200.bed --name promoter_express --qsub scg3_intervallight.sh --server ON --job expressionStatsL2E
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_selection_com_l2_xot --source annotations --domain in2shape_ce_repressed_TSS_l3_slopbed_up2000_dn200.bed --name promoter_repress --qsub scg3_intervallight.sh --server ON --job expressionStatsL2R
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_selection_com_l2_xot --source annotations --domain in2shape_ce_reference_TSS_l3_slopbed_up2000_dn200.bed --name promoter_observe --qsub scg3_intervallight.sh --server ON --job expressionStatsL2T

#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_selection_com_l3_xot --source annotations --domain in2shape_ce_expressed_TSS_l3_slopbed_up2000_dn200.bed --name promoter_express --qsub scg3_intervallight.sh --server ON --job expressionStatsL3E
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_selection_com_l3_xot --source annotations --domain in2shape_ce_repressed_TSS_l3_slopbed_up2000_dn200.bed --name promoter_repress --qsub scg3_intervallight.sh --server ON --job expressionStatsL3R
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_selection_com_l3_xot --source annotations --domain in2shape_ce_reference_TSS_l3_slopbed_up2000_dn200.bed --name promoter_observe --qsub scg3_intervallight.sh --server ON --job expressionStatsL3T


# Crunch co-association results:
python mapCAs.py --path ~/meTRN --mode crunch --peaks ce_selection_com_ex_xot --name promoter_express --cutoff 0.01
python mapCAs.py --path ~/meTRN --mode crunch --peaks ce_selection_com_ex_xot --name promoter_repress --cutoff 0.01
python mapCAs.py --path ~/meTRN --mode crunch --peaks ce_selection_com_ex_xot --name promoter_observe --cutoff 0.01

python mapCAs.py --path ~/meTRN --mode report --peaks ce_selection_com_ex_xot --name promoter_express --cutoff 0.01
python mapCAs.py --path ~/meTRN --mode report --peaks ce_selection_com_ex_xot --name promoter_repress --cutoff 0.01
python mapCAs.py --path ~/meTRN --mode report --peaks ce_selection_com_ex_xot --name promoter_observe --cutoff 0.01


python mapCAs.py --path ~/meTRN --mode crunch --peaks ce_selection_com_ex_xot --name promoter_express --cutoff 0.01
python mapCAs.py --path ~/meTRN --mode crunch --peaks ce_selection_com_ex_xot --name promoter_repress --cutoff 0.01
python mapCAs.py --path ~/meTRN --mode crunch --peaks ce_selection_com_ex_xot --name promoter_observe --cutoff 0.01

python mapCAs.py --path ~/meTRN --mode report --peaks ce_selection_com_ex_xot --name promoter_express --cutoff 0.01
python mapCAs.py --path ~/meTRN --mode report --peaks ce_selection_com_ex_xot --name promoter_repress --cutoff 0.01
python mapCAs.py --path ~/meTRN --mode report --peaks ce_selection_com_ex_xot --name promoter_observe --cutoff 0.01


python mapCAs.py --path ~/meTRN --mode crunch --peaks ce_selection_com_l1_xot --name promoter_express --cutoff 0.01
python mapCAs.py --path ~/meTRN --mode crunch --peaks ce_selection_com_l1_xot --name promoter_repress --cutoff 0.01
python mapCAs.py --path ~/meTRN --mode crunch --peaks ce_selection_com_l1_xot --name promoter_observe --cutoff 0.01

python mapCAs.py --path ~/meTRN --mode report --peaks ce_selection_com_l1_xot --name promoter_express --cutoff 0.01
python mapCAs.py --path ~/meTRN --mode report --peaks ce_selection_com_l1_xot --name promoter_repress --cutoff 0.01
python mapCAs.py --path ~/meTRN --mode report --peaks ce_selection_com_l1_xot --name promoter_observe --cutoff 0.01


python mapCAs.py --path ~/meTRN --mode crunch --peaks ce_selection_com_l2_xot --name promoter_express --cutoff 0.01
python mapCAs.py --path ~/meTRN --mode crunch --peaks ce_selection_com_l2_xot --name promoter_repress --cutoff 0.01
python mapCAs.py --path ~/meTRN --mode crunch --peaks ce_selection_com_l2_xot --name promoter_observe --cutoff 0.01

python mapCAs.py --path ~/meTRN --mode report --peaks ce_selection_com_l2_xot --name promoter_express --cutoff 0.01
python mapCAs.py --path ~/meTRN --mode report --peaks ce_selection_com_l2_xot --name promoter_repress --cutoff 0.01
python mapCAs.py --path ~/meTRN --mode report --peaks ce_selection_com_l2_xot --name promoter_observe --cutoff 0.01


python mapCAs.py --path ~/meTRN --mode crunch --peaks ce_selection_com_l3_xot --name promoter_express --cutoff 0.01
python mapCAs.py --path ~/meTRN --mode crunch --peaks ce_selection_com_l3_xot --name promoter_repress --cutoff 0.01
python mapCAs.py --path ~/meTRN --mode crunch --peaks ce_selection_com_l3_xot --name promoter_observe --cutoff 0.01

python mapCAs.py --path ~/meTRN --mode report --peaks ce_selection_com_l3_xot --name promoter_express --cutoff 0.01
python mapCAs.py --path ~/meTRN --mode report --peaks ce_selection_com_l3_xot --name promoter_repress --cutoff 0.01
python mapCAs.py --path ~/meTRN --mode report --peaks ce_selection_com_l3_xot --name promoter_observe --cutoff 0.01



#top