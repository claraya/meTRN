#!/usr/bin/sh


# Prepare inputs for neuron region analyses:
#rm -rf ~/meTRN/data/binding/ce_selection_com_ex_som/
#mkdir ~/meTRN/data/binding/ce_selection_com_ex_som/
#mkdir ~/meTRN/data/binding/ce_selection_com_ex_som/input/
#cp ~/meTRN/data/neurons/ce_selection_com_ex_raw/binary/results/any.ex.som/regions/bed/*bed ~/meTRN/data/binding/ce_selection_com_ex_som/input/

#rm -rf ~/meTRN/data/binding/ce_selection_com_l3_som/
#mkdir ~/meTRN/data/binding/ce_selection_com_l3_som/
#mkdir ~/meTRN/data/binding/ce_selection_com_l3_som/input/
#cp ~/meTRN/data/neurons/ce_selection_com_l3_raw/binary/results/any.l3.som/regions/bed/*bed ~/meTRN/data/binding/ce_selection_com_l3_som/input/


# Determine overlap with promoter/enhancer distances for neuron regions:
#python mapBinding.py --path ~/meTRN --organism ce --mode map:regions --peaks ce_selection_com_ex_som --infile in2shape_ce_modencode_MIX_ee.bed --name enhancer.overlaps.ee --headerDict auto --target window --label factor --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer --reference 0:500

#python mapBinding.py --path ~/meTRN --organism ce --mode map:regions --peaks ce_selection_com_l3_som --infile in2shape_ce_modencode_MIX_l3.bed --name enhancer.overlaps.l3 --headerDict auto --target window --label factor --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer --reference 0:500


#top
#bash runMaster2G.sh