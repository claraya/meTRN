#!/usr/bin/sh


# Launch species-specific co-association analyses with IntervalStats:
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks hs_selection_com_cx_xot --source annotations --domain in2shape_hs_modencode_TSS_gx_slopbed_up5000_dn500.nh.bed --name promoter_regions --qsub scg3_intervalstats.sh --server ON --job hsStats
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_selection_com_cx_xot --source annotations --domain in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub scg3_intervalstats.sh --server ON --job ceStats
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks dm_selection_com_cx_xot --source annotations --domain in2shape_dm_modencode_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub scg3_intervalstats.sh --server ON --job dmStats


# Crunch species-specific co-association analyses with IntervalStats:
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode crunch --peaks hs_selection_com_cx_xot --source annotations --domain in2shape_hs_modencode_TSS_gx_slopbed_up5000_dn500.nh.bed --name promoter_regions --qsub scg3_intervalstats.sh --server ON --job hsStats
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode crunch --peaks ce_selection_com_cx_xot --source annotations --domain in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub scg3_intervalstats.sh --server ON --job ceStats
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode crunch --peaks dm_selection_com_cx_xot --source annotations --domain in2shape_dm_modencode_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub scg3_intervalstats.sh --server ON --job dmStats


# Crunch species-specific co-association analyses with IntervalStats:
#python mapCAs.py --path ~/meTRN --mode crunch --peaks hs_selection_com_cx_xot --source annotations --domain in2shape_hs_modencode_TSS_gx_slopbed_up5000_dn500.nh.bed --name promoter_regions
#python mapCAs.py --path ~/meTRN --mode crunch --peaks ce_selection_com_cx_xot --source annotations --domain in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions
#python mapCAs.py --path ~/meTRN --mode crunch --peaks dm_selection_com_cx_xot --source annotations --domain in2shape_dm_modencode_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions


# Report species-specific co-association analyses with IntervalStats:
#python mapCAs.py --path ~/meTRN --mode report --peaks hs_selection_com_cx_xot --source annotations --domain in2shape_hs_modencode_TSS_gx_slopbed_up5000_dn500.nh.bed --name promoter_regions
#python mapCAs.py --path ~/meTRN --mode report --peaks ce_selection_com_cx_xot --source annotations --domain in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions
#python mapCAs.py --path ~/meTRN --mode report --peaks dm_selection_com_cx_xot --source annotations --domain in2shape_dm_modencode_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions


# Classify paralogous and non-paralogous species-specific co-association analyses:
#python mapCAs.py --path ~/meTRN --mode family --peaks hs_selection_com_cx_xot --source annotations --domain in2shape_hs_modencode_TSS_gx_slopbed_up5000_dn500.nh.bed --name promoter_regions
#python mapCAs.py --path ~/meTRN --mode family --peaks ce_selection_com_cx_xot --source annotations --domain in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions
#python mapCAs.py --path ~/meTRN --mode family --peaks dm_selection_com_cx_xot --source annotations --domain in2shape_dm_modencode_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions


# Carry-out procedure for the worm extension datasets:
#cp -r ~/meTRN/data/coassociations/ce_selection_com_cx_xot ~/meTRN/data/coassociations/ce_extension_com_cx_xot
#python mapCAs.py --path ~/meTRN --mode launch --peaks ce_extension_com_cx_xot --source annotations --domain in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --overwrite OFF


#top
