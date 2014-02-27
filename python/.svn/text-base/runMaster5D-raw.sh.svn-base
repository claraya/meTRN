#!/usr/bin/sh


# Prepare inputs for neuron GO enrichment analyses:
#mkdir ~/meTRN/data/go/any.ex.som
#mkdir ~/meTRN/data/go/any.l1.som
#mkdir ~/meTRN/data/go/any.l2.som
#mkdir ~/meTRN/data/go/any.l3.som
#mkdir ~/meTRN/data/go/any.l4.som

#mkdir ~/meTRN/data/go/any.ex.som/input
#mkdir ~/meTRN/data/go/any.l1.som/input
#mkdir ~/meTRN/data/go/any.l2.som/input
#mkdir ~/meTRN/data/go/any.l3.som/input
#mkdir ~/meTRN/data/go/any.l4.som/input

#cp ~/meTRN/data/neurons/ce_selection_com_ex_raw/binary/results/any.ex.som/regions/bed/*bed ~/meTRN/data/go/any.ex.som/input/
#cp ~/meTRN/data/neurons/ce_selection_com_l1_raw/binary/results/any.l1.som/regions/bed/*bed ~/meTRN/data/go/any.l1.som/input/
#cp ~/meTRN/data/neurons/ce_selection_com_l2_raw/binary/results/any.l2.som/regions/bed/*bed ~/meTRN/data/go/any.l2.som/input/
#cp ~/meTRN/data/neurons/ce_selection_com_l3_raw/binary/results/any.l3.som/regions/bed/*bed ~/meTRN/data/go/any.l3.som/input/
#cp ~/meTRN/data/neurons/ce_selection_com_l4_raw/binary/results/any.l4.som/regions/bed/*bed ~/meTRN/data/go/any.l4.som/input/


# Launch SOM map GO analysis:
#r --slave < ~/Desktop/Dropbox/meTRN/Code/R/mapGO_neurons.r


# Parse GO analysis from SOM neurons:
#python mapGO.py --path ~/meTRN --organism ce --peaks any.ex.som --analysis p5e-1 --target som.neurons --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05
#python mapGO.py --path ~/meTRN --organism ce --peaks any.l1.som --analysis p5e-1 --target som.neurons --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05
#python mapGO.py --path ~/meTRN --organism ce --peaks any.l2.som --analysis p5e-1 --target som.neurons --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05
#python mapGO.py --path ~/meTRN --organism ce --peaks any.l3.som --analysis p5e-1 --target som.neurons --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05
#python mapGO.py --path ~/meTRN --organism ce --peaks any.l4.som --analysis p5e-1 --target som.neurons --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05


# Build sunburst-visualization:
#python mapGO.py --path ~/meTRN --organism ce --mode sunburst --peaks any.ex.som,any.l1.som,any.l2.som,any.l3.som,any.l4.som --analysis p5e-1 --target som.neurons --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05 --minCount 5 --maxCount 25 --maxColor 6


# Collect neuron data from multiple SOMs:
#python mapNeurons.py --path ~/meTRN --organism ce --mode neuron.compile --peaks ce_selection_com_ex_raw,ce_selection_com_l1_raw,ce_selection_com_l2_raw,ce_selection_com_l3_raw,ce_selection_com_l4_raw --neurons any.ex.som,any.l1.som,any.l2.som,any.l3.som,any.l4.som --id binary.signature --rename neuron:N --fraction 0.75 --name any.xx.som

#python mapGO.py --path ~/meTRN --organism ce --mode compile --peaks any.ex.som,any.l1.som,any.l2.som,any.l3.som,any.l4.som --analysis p5e-1 --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05 --name any.xx.som


# Compare GO results from individual datasets and co-associations:
#python mapGO.py --path ~/meTRN --organism ce --mode compare --peaks ce_selection_com_cx_xot,any.ex.som --analysis p5e-1 --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05 --target embryonic --include BP


# Gather example targets from GO Terms
python mapGO.py --path ~/meTRN --organism ce --mode example --peaks ce_selection_com_cx_xot --analysis p5e-1 --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05 --target synaptic --include HAM-1,DPL-1,MEP-1,MAB-5,CES-1,FKH-10,ZAG-1



#top
#bash runMaster5D-raw.sh