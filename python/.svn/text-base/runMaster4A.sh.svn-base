#!/usr/bin/sh


# Launch orthologs GO analysis:
#r --slave < ~/Desktop/Dropbox/meTRN/Code/R/mapGO_ortho-hs.r
#r --slave < ~/Desktop/Dropbox/meTRN/Code/R/mapGO_ortho-ce.r
#r --slave < ~/Desktop/Dropbox/meTRN/Code/R/mapGO_ortho-dm.r


# Launch selection GO analysis:
#r --slave < ~/Desktop/Dropbox/meTRN/Code/R/mapGO_analysis-hs.r
#r --slave < ~/Desktop/Dropbox/meTRN/Code/R/mapGO_analysis-ce.r
#r --slave < ~/Desktop/Dropbox/meTRN/Code/R/mapGO_analysis-dm.r


# Parse GO analysis from ortholog binding data (at 0.01 cutoff):
#python mapGO.py --path ~/meTRN --organism hs --mode report --peaks hs_orthoHsCe_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.01
#python mapGO.py --path ~/meTRN --organism ce --mode report --peaks ce_orthoHsCe_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.01

#python mapGO.py --path ~/meTRN --organism hs --mode report --peaks hs_orthoHsDm_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.01
#python mapGO.py --path ~/meTRN --organism dm --mode report --peaks dm_orthoHsDm_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.01

#python mapGO.py --path ~/meTRN --organism ce --mode report --peaks ce_orthoCeDm_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.01
#python mapGO.py --path ~/meTRN --organism dm --mode report --peaks dm_orthoCeDm_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.01


# Perform between-species comparisons (at 0.01 cutoff):
#python mapHybrid.py --path ~/meTRN --mode merge.overlap --organism hs --species hs,ce --orthology groups --source go --A hs_orthoHsCe_com_cx_xot/p5e-1/summary/mapgo_complete_hs_orthoHsCe_com_cx_xot_p5_hc1_hp1e-02_summary --B  ce_orthoHsCe_com_cx_xot/p5e-1/summary/mapgo_complete_ce_orthoHsCe_com_cx_xot_p5_hc1_hp1e-02_summary --values id,id --label rebuild --target 'factor(context)' --name orthoHsCe_com_cx_xot_hp1e-02 --fraction 20

#python mapHybrid.py --path ~/meTRN --mode merge.overlap --organism hs --species hs,dm --orthology groups --source go --A hs_orthoHsDm_com_cx_xot/p5e-1/summary/mapgo_complete_hs_orthoHsDm_com_cx_xot_p5_hc1_hp1e-02_summary --B  dm_orthoHsDm_com_cx_xot/p5e-1/summary/mapgo_complete_dm_orthoHsDm_com_cx_xot_p5_hc1_hp1e-02_summary --values id,id --label rebuild --target 'factor(context)' --name orthoHsDm_com_cx_xot_hp1e-02 --fraction 20

#python mapHybrid.py --path ~/meTRN --mode merge.overlap --organism ce --species ce,dm --orthology groups --source go --A ce_orthoCeDm_com_cx_xot/p5e-1/summary/mapgo_complete_ce_orthoCeDm_com_cx_xot_p5_hc1_hp1e-02_summary --B  dm_orthoCeDm_com_cx_xot/p5e-1/summary/mapgo_complete_dm_orthoCeDm_com_cx_xot_p5_hc1_hp1e-02_summary --values id,id --label rebuild --target 'factor(context)' --name orthoCeDm_com_cx_xot_hp1e-02 --fraction 20


# Parse GO analysis from ortholog binding data (at 0.05 cutoff):
#python mapGO.py --path ~/meTRN --organism hs --mode report --peaks hs_orthoHsCe_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05
#python mapGO.py --path ~/meTRN --organism ce --mode report --peaks ce_orthoHsCe_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05

#python mapGO.py --path ~/meTRN --organism hs --mode report --peaks hs_orthoHsDm_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05
#python mapGO.py --path ~/meTRN --organism dm --mode report --peaks dm_orthoHsDm_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05

#python mapGO.py --path ~/meTRN --organism ce --mode report --peaks ce_orthoCeDm_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05
#python mapGO.py --path ~/meTRN --organism dm --mode report --peaks dm_orthoCeDm_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05


# Perform between-species comparisons (at 0.05 cutoff):
#python mapHybrid.py --path ~/meTRN --mode merge.overlap --organism hs --species hs,ce --orthology groups --source go --A hs_orthoHsCe_com_cx_xot/p5e-1/summary/mapgo_complete_hs_orthoHsCe_com_cx_xot_p5_hc1_hp5e-02_summary --B  ce_orthoHsCe_com_cx_xot/p5e-1/summary/mapgo_complete_ce_orthoHsCe_com_cx_xot_p5_hc1_hp5e-02_summary --values id,id --label rebuild --target 'factor(context)' --name orthoHsCe_com_cx_xot_hp5e-02 --fraction 20

#python mapHybrid.py --path ~/meTRN --mode merge.overlap --organism hs --species hs,dm --orthology groups --source go --A hs_orthoHsDm_com_cx_xot/p5e-1/summary/mapgo_complete_hs_orthoHsDm_com_cx_xot_p5_hc1_hp5e-02_summary --B  dm_orthoHsDm_com_cx_xot/p5e-1/summary/mapgo_complete_dm_orthoHsDm_com_cx_xot_p5_hc1_hp5e-02_summary --values id,id --label rebuild --target 'factor(context)' --name orthoHsDm_com_cx_xot_hp5e-02 --fraction 20

#python mapHybrid.py --path ~/meTRN --mode merge.overlap --organism ce --species ce,dm --orthology groups --source go --A ce_orthoCeDm_com_cx_xot/p5e-1/summary/mapgo_complete_ce_orthoCeDm_com_cx_xot_p5_hc1_hp5e-02_summary --B  dm_orthoCeDm_com_cx_xot/p5e-1/summary/mapgo_complete_dm_orthoCeDm_com_cx_xot_p5_hc1_hp5e-02_summary --values id,id --label rebuild --target 'factor(context)' --name orthoCeDm_com_cx_xot_hp5e-02 --fraction 20


# Parse GO analysis from ortholog binding data (at 0.1 cutoff):
#python mapGO.py --path ~/meTRN --organism hs --mode report --peaks hs_orthoHsCe_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.1
#python mapGO.py --path ~/meTRN --organism ce --mode report --peaks ce_orthoHsCe_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.1

#python mapGO.py --path ~/meTRN --organism hs --mode report --peaks hs_orthoHsDm_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.1
#python mapGO.py --path ~/meTRN --organism dm --mode report --peaks dm_orthoHsDm_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.1

#python mapGO.py --path ~/meTRN --organism ce --mode report --peaks ce_orthoCeDm_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.1
#python mapGO.py --path ~/meTRN --organism dm --mode report --peaks dm_orthoCeDm_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.1


# Perform between-species comparisons (at 0.1 cutoff):
#python mapHybrid.py --path ~/meTRN --mode merge.overlap --organism hs --species hs,ce --orthology groups --source go --A hs_orthoHsCe_com_cx_xot/p5e-1/summary/mapgo_complete_hs_orthoHsCe_com_cx_xot_p5_hc1_hp1e-01_summary --B  ce_orthoHsCe_com_cx_xot/p5e-1/summary/mapgo_complete_ce_orthoHsCe_com_cx_xot_p5_hc1_hp1e-01_summary --values id,id --label rebuild --target 'factor(context)' --name orthoHsCe_com_cx_xot_hp1e-01 --fraction 20

#python mapHybrid.py --path ~/meTRN --mode merge.overlap --organism hs --species hs,dm --orthology groups --source go --A hs_orthoHsDm_com_cx_xot/p5e-1/summary/mapgo_complete_hs_orthoHsDm_com_cx_xot_p5_hc1_hp1e-01_summary --B  dm_orthoHsDm_com_cx_xot/p5e-1/summary/mapgo_complete_dm_orthoHsDm_com_cx_xot_p5_hc1_hp1e-01_summary --values id,id --label rebuild --target 'factor(context)' --name orthoHsDm_com_cx_xot_hp1e-01 --fraction 20

#python mapHybrid.py --path ~/meTRN --mode merge.overlap --organism ce --species ce,dm --orthology groups --source go --A ce_orthoCeDm_com_cx_xot/p5e-1/summary/mapgo_complete_ce_orthoCeDm_com_cx_xot_p5_hc1_hp1e-01_summary --B  dm_orthoCeDm_com_cx_xot/p5e-1/summary/mapgo_complete_dm_orthoCeDm_com_cx_xot_p5_hc1_hp1e-01_summary --values id,id --label rebuild --target 'factor(context)' --name orthoCeDm_com_cx_xot_hp1e-01 --fraction 20



# Parse GO analysis from selection binding data:
#python mapGO.py --path ~/meTRN --organism hs --mode report --peaks hs_selection_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05
#python mapGO.py --path ~/meTRN --organism ce --mode report --peaks ce_selection_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05
#python mapGO.py --path ~/meTRN --organism dm --mode report --peaks dm_selection_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05

#python mapGO.py --path ~/meTRN --organism hs --mode report --peaks hs_selection_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 2 --hitPvalue 0.05
#python mapGO.py --path ~/meTRN --organism ce --mode report --peaks ce_selection_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 2 --hitPvalue 0.05
#python mapGO.py --path ~/meTRN --organism dm --mode report --peaks dm_selection_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 2 --hitPvalue 0.05

#python mapGO.py --path ~/meTRN --organism hs --mode report --peaks hs_selection_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 3 --hitPvalue 0.05
#python mapGO.py --path ~/meTRN --organism ce --mode report --peaks ce_selection_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 3 --hitPvalue 0.05
#python mapGO.py --path ~/meTRN --organism dm --mode report --peaks dm_selection_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 3 --hitPvalue 0.05


# Parse GO analysis for worm-analysis specifically:
#python mapGO.py --path ~/meTRN --organism ce --mode report --peaks ce_selection_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 1 --hitPvalue 0.01
#python mapGO.py --path ~/meTRN --organism ce --mode report --peaks ce_selection_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 2 --hitPvalue 0.01
#python mapGO.py --path ~/meTRN --organism ce --mode report --peaks ce_selection_com_cx_xot --analysis p5e-1 --target factor.context --minPvalue 0.5 --hitCount 3 --hitPvalue 0.01


#top
#bash runMaster4A.sh