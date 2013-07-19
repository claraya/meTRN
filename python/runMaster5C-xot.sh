#!/usr/bin/sh


# Generate matrix region reference files:
#python mapNeurons.py --path ~/meTRN --organism ce --mode matrix.regions --peaks ce_selection_com_cx_xot
#python mapNeurons.py --path ~/meTRN --organism ce --mode matrix.regions --peaks ce_selection_com_ex_xot
#python mapNeurons.py --path ~/meTRN --organism ce --mode matrix.regions --peaks ce_selection_com_l1_xot
#python mapNeurons.py --path ~/meTRN --organism ce --mode matrix.regions --peaks ce_selection_com_l2_xot
#python mapNeurons.py --path ~/meTRN --organism ce --mode matrix.regions --peaks ce_selection_com_l3_xot
#python mapNeurons.py --path ~/meTRN --organism ce --mode matrix.regions --peaks ce_selection_com_l4_xot


# Generate BED files from SOM regions:
#python mapNeurons.py --path ~/meTRN --organism ce --mode convert.regions --peaks ce_selection_com_ex_xot --neurons any.ex.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode convert.regions --peaks ce_selection_com_l1_xot --neurons any.l1.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode convert.regions --peaks ce_selection_com_l2_xot --neurons any.l2.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode convert.regions --peaks ce_selection_com_l3_xot --neurons any.l3.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode convert.regions --peaks ce_selection_com_l4_xot --neurons any.l4.som

#python mapNeurons.py --path ~/meTRN --organism ce --mode convert.regions --peaks ce_selection_com_cx_xot --neurons all.elc.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode convert.regions --peaks ce_selection_com_cx_xot --neurons all.e1c.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode convert.regions --peaks ce_selection_com_cx_xot --neurons all.e2c.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode convert.regions --peaks ce_selection_com_cx_xot --neurons all.e3c.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode convert.regions --peaks ce_selection_com_cx_xot --neurons all.e4c.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode convert.regions --peaks ce_selection_com_cx_xot --neurons all.1v2.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode convert.regions --peaks ce_selection_com_cx_xot --neurons all.2v3.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode convert.regions --peaks ce_selection_com_cx_xot --neurons all.3v4.som


# Generate neuron reports from SOM codes and regions:
#python mapNeurons.py --path ~/meTRN --organism ce --mode neuron.summary --peaks ce_selection_com_ex_xot --neurons any.ex.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode neuron.summary --peaks ce_selection_com_l1_xot --neurons any.l1.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode neuron.summary --peaks ce_selection_com_l2_xot --neurons any.l2.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode neuron.summary --peaks ce_selection_com_l3_xot --neurons any.l3.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode neuron.summary --peaks ce_selection_com_l4_xot --neurons any.l4.som

#python mapNeurons.py --path ~/meTRN --organism ce --mode neuron.summary --peaks ce_selection_com_cx_xot --neurons all.elc.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode neuron.summary --peaks ce_selection_com_cx_xot --neurons all.e1c.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode neuron.summary --peaks ce_selection_com_cx_xot --neurons all.e2c.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode neuron.summary --peaks ce_selection_com_cx_xot --neurons all.e3c.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode neuron.summary --peaks ce_selection_com_cx_xot --neurons all.e4c.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode neuron.summary --peaks ce_selection_com_cx_xot --neurons all.1v2.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode neuron.summary --peaks ce_selection_com_cx_xot --neurons all.2v3.som
#python mapNeurons.py --path ~/meTRN --organism ce --mode neuron.summary --peaks ce_selection_com_cx_xot --neurons all.3v4.som


# Annotate neuron regions:
#python mapNeurons.py --path ~/meTRN --organism ce --mode matrix.expression --infile mapneurons_matrix_context_any_ORC_FAC_ce_selection_com_ex_xot.bed --peaks ce_selection_com_ex_xot --source annotations --input in2shape_ce_expressed_RNA_ee.bed --tss in2shape_ce_modencode_TSS_gx_slopbed_up5000_dn500.nh.bed --contexts EE
#python mapNeurons.py --path ~/meTRN --organism ce --mode matrix.expression --infile mapneurons_matrix_context_any_ORC_FAC_ce_selection_com_ex_xot.bed --peaks ce_selection_com_ex_xot --source annotations --input in2shape_ce_expressed_RNA_le.bed --tss in2shape_ce_modencode_TSS_gx_slopbed_up5000_dn500.nh.bed --contexts LE
#python mapNeurons.py --path ~/meTRN --organism ce --mode matrix.expression --infile mapneurons_matrix_context_any_ORC_FAC_ce_selection_com_ex_xot.bed --peaks ce_selection_com_ex_xot --source annotations --input in2shape_ce_expressed_RNA_em.bed --tss in2shape_ce_modencode_TSS_gx_slopbed_up5000_dn500.nh.bed --contexts embryonic --name EX

#python mapNeurons.py --path ~/meTRN --organism ce --mode matrix.expression --infile mapneurons_matrix_context_any_ORC_FAC_ce_selection_com_l1_xot.bed --peaks ce_selection_com_l1_xot --source annotations --input in2shape_ce_expressed_RNA_l1.bed --tss in2shape_ce_modencode_TSS_gx_slopbed_up5000_dn500.nh.bed --contexts L1
#python mapNeurons.py --path ~/meTRN --organism ce --mode matrix.expression --infile mapneurons_matrix_context_any_ORC_FAC_ce_selection_com_l2_xot.bed --peaks ce_selection_com_l2_xot --source annotations --input in2shape_ce_expressed_RNA_l2.bed --tss in2shape_ce_modencode_TSS_gx_slopbed_up5000_dn500.nh.bed --contexts L2
#python mapNeurons.py --path ~/meTRN --organism ce --mode matrix.expression --infile mapneurons_matrix_context_any_ORC_FAC_ce_selection_com_l3_xot.bed --peaks ce_selection_com_l3_xot --source annotations --input in2shape_ce_expressed_RNA_l3.bed --tss in2shape_ce_modencode_TSS_gx_slopbed_up5000_dn500.nh.bed --contexts L3


# Collect cells from neurons:
#python mapCells.py --path ~/meTRN --organism ce/ --mode cell.collection --neurons all.3v4.som --collection all.3v4.som



#top