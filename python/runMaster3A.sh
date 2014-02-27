#!/usr/bin/sh

# Grab new ortholog lists:
#scp claraya@snively.stanford.edu:/srv/gs1/projects/snyder/modENCODE/data/orthologs/family_group/formatted_v3/* ~/meTRN/data/orthologs/groups/


# Generate orthologous peak sets (raw):
#python mapPeaks.py --path ~/meTRN --mode orthologs --rename selection --peaks hs_selection_com_xx_raw --orthology groups  --nametag ortho --species hs,ce
#python mapPeaks.py --path ~/meTRN --mode orthologs --rename extension --peaks ce_extension_com_xx_raw --orthology groups  --nametag ortho --species hs,ce

#python mapPeaks.py --path ~/meTRN --mode orthologs --rename selection --peaks hs_selection_com_xx_raw --orthology groups  --nametag ortho --species hs,dm
#python mapPeaks.py --path ~/meTRN --mode orthologs --rename selection --peaks dm_selection_com_xx_raw --orthology groups  --nametag ortho --species hs,dm

#python mapPeaks.py --path ~/meTRN --mode orthologs --rename extension --peaks ce_extension_com_xx_raw --orthology groups  --nametag ortho --species ce,dm
#python mapPeaks.py --path ~/meTRN --mode orthologs --rename selection --peaks dm_selection_com_xx_raw --orthology groups  --nametag ortho --species ce,dm


# Generate orthologous peak sets (raw):
#python mapPeaks.py --path ~/meTRN --mode orthologs --rename selection --peaks hs_selection_com_cx_raw --orthology groups  --nametag ortho --species hs,ce
#python mapPeaks.py --path ~/meTRN --mode orthologs --rename extension --peaks ce_extension_com_cx_raw --orthology groups  --nametag ortho --species hs,ce

#python mapPeaks.py --path ~/meTRN --mode orthologs --rename selection --peaks hs_selection_com_cx_raw --orthology groups  --nametag ortho --species hs,dm
#python mapPeaks.py --path ~/meTRN --mode orthologs --rename selection --peaks dm_selection_com_cx_raw --orthology groups  --nametag ortho --species hs,dm

#python mapPeaks.py --path ~/meTRN --mode orthologs --rename extension --peaks ce_extension_com_cx_raw --orthology groups  --nametag ortho --species ce,dm
#python mapPeaks.py --path ~/meTRN --mode orthologs --rename selection --peaks dm_selection_com_cx_raw --orthology groups  --nametag ortho --species ce,dm


# Generate orthologous peak sets (hot):
#python mapPeaks.py --path ~/meTRN --mode orthologs --rename selection --peaks hs_selection_com_cx_hot --orthology groups  --nametag ortho --species hs,ce
#python mapPeaks.py --path ~/meTRN --mode orthologs --rename extension --peaks ce_extension_com_cx_hot --orthology groups  --nametag ortho --species hs,ce

#python mapPeaks.py --path ~/meTRN --mode orthologs --rename selection --peaks hs_selection_com_cx_hot --orthology groups  --nametag ortho --species hs,dm
#python mapPeaks.py --path ~/meTRN --mode orthologs --rename selection --peaks dm_selection_com_cx_hot --orthology groups  --nametag ortho --species hs,dm

#python mapPeaks.py --path ~/meTRN --mode orthologs --rename extension --peaks ce_extension_com_cx_hot --orthology groups  --nametag ortho --species ce,dm
#python mapPeaks.py --path ~/meTRN --mode orthologs --rename selection --peaks dm_selection_com_cx_hot --orthology groups  --nametag ortho --species ce,dm


# Generate orthologous peak sets (xot):
#python mapPeaks.py --path ~/meTRN --mode orthologs --rename selection --peaks hs_selection_com_cx_xot --orthology groups  --nametag ortho --species hs,ce
#python mapPeaks.py --path ~/meTRN --mode orthologs --rename extension --peaks ce_extension_com_cx_xot --orthology groups  --nametag ortho --species hs,ce

#python mapPeaks.py --path ~/meTRN --mode orthologs --rename selection --peaks hs_selection_com_cx_xot --orthology groups  --nametag ortho --species hs,dm
#python mapPeaks.py --path ~/meTRN --mode orthologs --rename selection --peaks dm_selection_com_cx_xot --orthology groups  --nametag ortho --species hs,dm

#python mapPeaks.py --path ~/meTRN --mode orthologs --rename extension --peaks ce_extension_com_cx_xot --orthology groups  --nametag ortho --species ce,dm
#python mapPeaks.py --path ~/meTRN --mode orthologs --rename selection --peaks dm_selection_com_cx_xot --orthology groups  --nametag ortho --species ce,dm


# Generate missing peak-report files:
#python mapPeaks.py --path ~/meTRN --mode build --overwrite OFF


# Launch comparative co-association analyses with IntervalStats (promoter regions, on SCG3 at Stanford):
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks hs_orthoHsCe_com_cx_xot --source annotations --domain in2shape_hs_modencode_TSS_gx_slopbed_up5000_dn500.nh.bed --name promoter_regions --qsub scg3_intervalstats.sh --server ON --job hs_orthoHsCe_promoter
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_orthoHsCe_com_cx_xot --source annotations --domain in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub scg3_intervalstats.sh --server ON --job ce_orthoHsCe_promoter

#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks hs_orthoHsDm_com_cx_xot --source annotations --domain in2shape_hs_modencode_TSS_gx_slopbed_up5000_dn500.nh.bed --name promoter_regions --qsub scg3_intervalstats.sh --server ON --job hs_orthoHsDm_promoter
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks dm_orthoHsDm_com_cx_xot --source annotations --domain in2shape_dm_modencode_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub scg3_intervalstats.sh --server ON --job dm_orthoHsDm_promoter

#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_orthoCeDm_com_cx_xot --source annotations --domain in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub scg3_intervalstats.sh --server ON --job ce_orthoCeDm_promoter
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks dm_orthoCeDm_com_cx_xot --source annotations --domain in2shape_dm_modencode_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub scg3_intervalstats.sh --server ON --job dm_orthoCeDm_promoter


# Launch comparative co-association analyses with IntervalStats (promoter regions; on GS clusters at UW):
#python mapCAs.py --path ~/meTRN --mode launch --peaks hs_orthoHsCe_com_cx_xot --source annotations --domain in2shape_hs_modencode_TSS_gx_slopbed_up5000_dn500.nh.bed --name promoter_regions --qsub uwgs_intervalstats.sh --server GS --job hs_orthoHsCe_promoter
#python mapCAs.py --path ~/meTRN --mode launch --peaks ce_orthoHsCe_com_cx_xot --source annotations --domain in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub uwgs_intervalstats.sh --server GS --job ce_orthoHsCe_promoter

#python mapCAs.py --path ~/meTRN --mode launch --peaks hs_orthoHsDm_com_cx_xot --source annotations --domain in2shape_hs_modencode_TSS_gx_slopbed_up5000_dn500.nh.bed --name promoter_regions --qsub uwgs_intervalstats.sh --server GS --job hs_orthoHsDm_promoter
#python mapCAs.py --path ~/meTRN --mode launch --peaks dm_orthoHsDm_com_cx_xot --source annotations --domain in2shape_dm_modencode_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub uwgs_intervalstats.sh --server GS --job dm_orthoHsDm_promoter

#python mapCAs.py --path ~/meTRN --mode launch --peaks ce_orthoCeDm_com_cx_xot --source annotations --domain in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub uwgs_intervalstats.sh --server GS --job ce_orthoCeDm_promoter
#python mapCAs.py --path ~/meTRN --mode launch --peaks dm_orthoCeDm_com_cx_xot --source annotations --domain in2shape_dm_modencode_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --qsub uwgs_intervalstats.sh --server GS --job dm_orthoCeDm_promoter


# Crunch co-association results from IntervalStats (promoter regions):
#python mapCAs.py --path ~/meTRN --mode crunch --peaks hs_orthoHsCe_com_cx_xot --name promoter_regions --threads 5 --cutoff 0.05
#python mapCAs.py --path ~/meTRN --mode crunch --peaks ce_orthoHsCe_com_cx_xot --name promoter_regions --threads 5 --cutoff 0.05

#python mapCAs.py --path ~/meTRN --mode crunch --peaks hs_orthoHsDm_com_cx_xot --name promoter_regions --threads 5 --cutoff 0.05
#python mapCAs.py --path ~/meTRN --mode crunch --peaks dm_orthoHsDm_com_cx_xot --name promoter_regions --threads 5 --cutoff 0.05

#python mapCAs.py --path ~/meTRN --mode crunch --peaks ce_orthoCeDm_com_cx_xot --name promoter_regions --threads 5 --cutoff 0.05
#python mapCAs.py --path ~/meTRN --mode crunch --peaks dm_orthoCeDm_com_cx_xot --name promoter_regions --threads 5 --cutoff 0.05


# Generate co-association reports from IntervalStats (promoter regions):
#python mapCAs.py --path ~/meTRN --mode report --peaks hs_orthoHsCe_com_cx_xot --name promoter_regions --threads 5 --cutoff 0.05 --label 'factor(context)'
#python mapCAs.py --path ~/meTRN --mode report --peaks ce_orthoHsCe_com_cx_xot --name promoter_regions --threads 5 --cutoff 0.05 --label 'factor(context)'

#python mapCAs.py --path ~/meTRN --mode report --peaks hs_orthoHsDm_com_cx_xot --name promoter_regions --threads 5 --cutoff 0.05 --label 'factor(context)'
#python mapCAs.py --path ~/meTRN --mode report --peaks dm_orthoHsDm_com_cx_xot --name promoter_regions --threads 5 --cutoff 0.05 --label 'factor(context)'

#python mapCAs.py --path ~/meTRN --mode report --peaks ce_orthoCeDm_com_cx_xot --name promoter_regions --threads 5 --cutoff 0.05 --label 'factor(context)'
#python mapCAs.py --path ~/meTRN --mode report --peaks dm_orthoCeDm_com_cx_xot --name promoter_regions --threads 5 --cutoff 0.05 --label 'factor(context)'


# Compare co-association results between species orthologs (promoter regions):
#python mapHybrid.py --path ~/meTRN --mode merge.matrix --organism hs --species hs,ce --orthology groups --source coassociations --A hs_orthoHsCe_com_cx_xot/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --B ce_orthoHsCe_com_cx_xot/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --indexes i:j,i:j --values mirror.passing,mirror.passing --name orthoHsCe_com_cx_xot

#python mapHybrid.py --path ~/meTRN --mode merge.matrix --organism hs --species hs,dm --orthology groups --source coassociations --A hs_orthoHsDm_com_cx_xot/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --B dm_orthoHsDm_com_cx_xot/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --indexes i:j,i:j --values mirror.passing,mirror.passing --name orthoHsDm_com_cx_xot

#python mapHybrid.py --path ~/meTRN --mode merge.matrix --organism ce --species ce,dm --orthology groups --source coassociations --A ce_orthoCeDm_com_cx_xot/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --B dm_orthoCeDm_com_cx_xot/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --indexes i:j,i:j --values mirror.passing,mirror.passing --name orthoCeDm_com_cx_xot


# Launch comparative co-association analyses with IntervalStats (combined regions):
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks hs_orthoHsCe_com_cx_xot --source peaks --domain mappeaks_hs_orthoHsCe_com_cx_xot_compiled.bed --name compiled_regions --qsub scg3_intervalstats.sh --server ON --job hs_orthoHsCe_compiled
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_orthoHsCe_com_cx_xot --source peaks --domain mappeaks_ce_orthoHsCe_com_cx_xot_compiled.bed --name compiled_regions --qsub scg3_intervalstats.sh --server ON --job ce_orthoHsCe_compiled

#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks hs_orthoHsDm_com_cx_xot --source peaks --domain mappeaks_hs_orthoHsDm_com_cx_xot_compiled.bed --name compiled_regions --qsub scg3_intervalstats.sh --server ON --job hs_orthoHsDm_compiled
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks dm_orthoHsDm_com_cx_xot --source peaks --domain mappeaks_dm_orthoHsDm_com_cx_xot_compiled.bed --name compiled_regions --qsub scg3_intervalstats.sh --server ON --job dm_orthoHsDm_compiled

#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks ce_orthoCeDm_com_cx_xot --source peaks --domain mappeaks_ce_orthoCeDm_com_cx_xot_compiled.bed --name compiled_regions --qsub scg3_intervalstats.sh --server ON --job ce_orthoCeDm_compiled
#python mapCAs.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peaks dm_orthoCeDm_com_cx_xot --source peaks --domain mappeaks_dm_orthoCeDm_com_cx_xot_compiled.bed --name compiled_regions --qsub scg3_intervalstats.sh --server ON --job dm_orthoCeDm_compiled








#top
#bash runMaster3A.sh