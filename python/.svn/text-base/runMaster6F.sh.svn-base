#!/usr/bin/sh



# Intersect genomic co-associations with cellular expression overlap:
#python mapCells.py --path ~/meTRN --mode cell.hybrid --organism ce --A ce_selection_com_cx_xot/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --B cetrn_tracked_f0.1/mapcells_cetrn_tracked_f0.1_matrix_overlap --indexes i.factor,j.factor --values i.context,j.context --contexts embryonic
#python mapCells.py --path ~/meTRN --mode cell.hybrid --organism ce --A ce_selection_com_cx_xot/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --B cetrn_focused_f0.1/mapcells_cetrn_focused_f0.1_matrix_overlap --indexes i.factor,j.factor --values i.context,j.context --contexts embryonic
#python mapCells.py --path ~/meTRN --mode cell.hybrid --organism ce --A ce_selection_com_cx_xot/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --B cetrn_fime244_f0.1/mapcells_cetrn_fime244_f0.1_matrix_overlap --indexes i.factor,j.factor --values i.context,j.context --contexts embryonic

#python mapCells.py --path ~/meTRN --mode cell.hybrid --organism ce --A ce_selection_com_cx_xot/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --B cetrn_inherit_f0.1/mapcells_cetrn_inherit_f0.1_matrix_overlap --indexes i.factor,j.factor --values i.context,j.context --contexts embryonic
#python mapCells.py --path ~/meTRN --mode cell.hybrid --organism ce --A ce_selection_com_cx_xot/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --B cetrn_inleafs_f0.1/mapcells_cetrn_inleafs_f0.1_matrix_overlap --indexes i.factor,j.factor --values i.context,j.context --contexts embryonic
#python mapCells.py --path ~/meTRN --mode cell.hybrid --organism ce --A ce_selection_com_cx_xot/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --B cetrn_maximal_f0.1/mapcells_cetrn_maximal_f0.1_matrix_overlap --indexes i.factor,j.factor --values i.context,j.context --contexts embryonic
#python mapCells.py --path ~/meTRN --mode cell.hybrid --organism ce --A ce_selection_com_cx_xot/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --B cetrn_mxleafs_f0.1/mapcells_cetrn_mxleafs_f0.1_matrix_overlap --indexes i.factor,j.factor --values i.context,j.context --contexts embryonic

#python mapCells.py --path ~/meTRN --mode cell.hybrid --organism ce --A ce_selection_com_cx_xot/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --B cetrn_inherit_f0.1/mapcells_cetrn_inherit_f0.1_matrix_overlap --indexes i.factor,j.factor --values i.context,j.context --contexts order.condensed
#python mapCells.py --path ~/meTRN --mode cell.hybrid --organism ce --A ce_selection_com_cx_xot/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --B cetrn_inleafs_f0.1/mapcells_cetrn_inleafs_f0.1_matrix_overlap --indexes i.factor,j.factor --values i.context,j.context --contexts order.condensed
#python mapCells.py --path ~/meTRN --mode cell.hybrid --organism ce --A ce_selection_com_cx_xot/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --B cetrn_maximal_f0.1/mapcells_cetrn_maximal_f0.1_matrix_overlap --indexes i.factor,j.factor --values i.context,j.context --contexts order.condensed
#python mapCells.py --path ~/meTRN --mode cell.hybrid --organism ce --A ce_selection_com_cx_xot/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --B cetrn_mxleafs_f0.1/mapcells_cetrn_mxleafs_f0.1_matrix_overlap --indexes i.factor,j.factor --values i.context,j.context --contexts order.condensed


# Intersect genomic co-associations with cellular expression overlap (across developmental time-points):
#python mapCells.py --path ~/meTRN --mode master:cell.hybrid --organism ce --A ce_selection_com_cx_xot/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --indexes i.factor,j.factor --values i.context,j.context --contexts embryonic --collection cetrn_fime --nametag _f0.1 --start 50 --stop 250 --total 350 --step 5 --chunks 5 --threads 8


# Integrate genomic co-associations with cellular expression overlap as a function of developmental time-points:
#python mapCells.py --path ~/meTRN --mode cell.dynamics --organism ce --collection cetrn_fime --nametag _f0.1 --start 50 --stop 250 --total 350 --step 5 --name cetrn_focused_f0.1
#python mapCells.py --path ~/meTRN --mode cell.dynamics --organism ce --collection cetrn_fime --nametag _f0.1 --start 50 --stop 250 --total 350 --step 5 --name cetrn_focused_f0.1 --peaks ce_selection_com_cx_xot --domain promoter_regions --contexts embryonic


# Launch cellular-resolution co-association analyses with IntervalStats:
#python mapCAs.py --path ~/meTRN --mode launch --peaks ce_xpfocused_reg_ex_all --cells ON --source annotations --domain mapcells_cetrn_focused_f0.1_in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub uwgs_intervalstats.sh --server GS --job focusAll --chunks 20
#python mapCAs.py --path ~/meTRN --mode launch --peaks ce_xpinherit_reg_ex_all --cells ON --source annotations --domain mapcells_cetrn_inherit_f0.1_in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub uwgs_intervalstats.sh --server GS --job inherAll --chunks 20

#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_xpfocused_reg_ex_all --cells ON --source annotations --domain mapcells_cetrn_focused_f0.1_in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub scg3_intervalextra.sh --server ON --job focusAll --chunks 20
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_xpfocused_reg_ex_240 --cells ON --source annotations --domain mapcells_cetrn_time240_f0.1_in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub scg3_intervalextra.sh --server ON --job focus240 --chunks 20


# Launch cellular-resolution co-association analyses with IntervalStats:
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_xptracked_reg_ex_050 --cells ON --source annotations --domain mapcells_cetrn_time050_f0.1_in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub scg3_intervalextra.sh --server ON --job time050
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_xptracked_reg_ex_070 --cells ON --source annotations --domain mapcells_cetrn_time070_f0.1_in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub scg3_intervalextra.sh --server ON --job time070 --chunks 10
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_xptracked_reg_ex_090 --cells ON --source annotations --domain mapcells_cetrn_time090_f0.1_in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub scg3_intervalextra.sh --server ON --job time090 --chunks 10
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_xptracked_reg_ex_110 --cells ON --source annotations --domain mapcells_cetrn_time110_f0.1_in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub scg3_intervalextra.sh --server ON --job time110 --chunks 10
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_xptracked_reg_ex_130 --cells ON --source annotations --domain mapcells_cetrn_time130_f0.1_in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub scg3_intervalextra.sh --server ON --job time130 --chunks 10
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_xptracked_reg_ex_150 --cells ON --source annotations --domain mapcells_cetrn_time150_f0.1_in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub scg3_intervalextra.sh --server ON --job time150 --chunks 10
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_xptracked_reg_ex_170 --cells ON --source annotations --domain mapcells_cetrn_time170_f0.1_in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub scg3_intervalextra.sh --server ON --job time170 --chunks 10
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_xptracked_reg_ex_190 --cells ON --source annotations --domain mapcells_cetrn_time190_f0.1_in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub scg3_intervalextra.sh --server ON --job time190 --chunks 10


# Crunch species-specific co-association analyses with IntervalStats:
#python mapCAs.py --path ~/meTRN --mode crunch --peaks ce_xptracked_reg_ex_170 --source annotations --domain in2shape_hs_modencode_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions
#python mapCAs.py --path ~/meTRN --mode crunch --peaks ce_xptracked_reg_ex_190 --source annotations --domain in2shape_hs_modencode_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions


# Report species-specific co-association analyses with IntervalStats:
#python mapCAs.py --path ~/meTRN --mode report --peaks ce_xptracked_reg_ex_170 --source annotations --domain in2shape_hs_modencode_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions
#python mapCAs.py --path ~/meTRN --mode report --peaks ce_xptracked_reg_ex_190 --source annotations --domain in2shape_hs_modencode_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions


# Re-run this shit (greetings from Peet's)!
#python mapHybrid.py --path ~/meTRN --mode merge.direct --organism ce --species ce,ce --source coassociations --A ce_xptracked_reg_ex_170/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --B ce_selection_com_cx_xot/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --indexes i:j,i:j --values mirror.passing,mirror.passing --name xptracked_reg_ex_170


#top
#bash runMaster6F.sh