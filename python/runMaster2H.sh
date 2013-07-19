#!/usr/bin/sh


# Prepare inputs for RGB region analyses:
#mkdir ~/meTRN/data/binding/hs_selection_com_cx_raw/
#mkdir ~/meTRN/data/binding/hs_selection_com_cx_raw/input/
#cp ~/meTRN/data/peaks/mappeaks_hs_selection_com_cx_raw_compiled.bed ~/meTRN/data/binding/hs_selection_com_cx_raw/input/


# Determine overlap with chromatin states for the aggregate/compiled binding regions:
#python mapBinding.py --path ~/meTRN --organism hs --mode map:regions --peaks hs_selection_com_cx_raw --infile in2shape_hs_modencode_HMM_gm.bed --name genomic.states.gm --queries auto --target feature --label factor.context --reference 1_Pro --exclude 17_Unmap --rename mappeaks_:,_compiled:
#python mapBinding.py --path ~/meTRN --organism hs --mode map:regions --peaks hs_selection_com_cx_raw --infile in2shape_hs_modencode_HMM_h1.bed --name genomic.states.h1 --queries auto --target feature --label factor.context --reference 1_Pro --exclude 17_Unmap --rename mappeaks_:,_compiled:

#python mapBinding.py --path ~/meTRN --organism hs --mode map:regions --peaks hs_selection_com_cx_raw --infile in2shape_hs_modencode_XMM_gm.bed --name combine.states.gm --queries auto --target feature --label factor.context --reference TSS --rename mappeaks_:,_compiled:
#python mapBinding.py --path ~/meTRN --organism hs --mode map:regions --peaks hs_selection_com_cx_raw --infile in2shape_hs_modencode_XMM_h1.bed --name combine.states.h1 --queries auto --target feature --label factor.context --reference TSS --rename mappeaks_:,_compiled:
#python mapBinding.py --path ~/meTRN --organism hs --mode map:regions --peaks hs_selection_com_cx_raw --infile in2shape_hs_modencode_XMM_hg.bed --name combine.states.hg --queries auto --target feature --label factor.context --reference TSS --rename mappeaks_:,_compiled:
#python mapBinding.py --path ~/meTRN --organism hs --mode map:regions --peaks hs_selection_com_cx_raw --infile in2shape_hs_modencode_XMM_hl.bed --name combine.states.hl --queries auto --target feature --label factor.context --reference TSS --rename mappeaks_:,_compiled:
#python mapBinding.py --path ~/meTRN --organism hs --mode map:regions --peaks hs_selection_com_cx_raw --infile in2shape_hs_modencode_XMM_k5.bed --name combine.states.k5 --queries auto --target feature --label factor.context --reference TSS --rename mappeaks_:,_compiled:

#python mapBinding.py --path ~/meTRN --organism hs --mode map:regions --peaks hs_selection_com_cx_raw --infile in2shape_hs_modencode_KMM_gm.bed --name kundaje.states.gm --queries auto --target feature --label factor.context --reference Tss --rename mappeaks_:,_compiled:
#python mapBinding.py --path ~/meTRN --organism hs --mode map:regions --peaks hs_selection_com_cx_raw --infile in2shape_hs_modencode_KMM_h1.bed --name kundaje.states.h1 --queries auto --target feature --label factor.context --reference Tss --rename mappeaks_:,_compiled:
#python mapBinding.py --path ~/meTRN --organism hs --mode map:regions --peaks hs_selection_com_cx_raw --infile in2shape_hs_modencode_KMM_hg.bed --name kundaje.states.hg --queries auto --target feature --label factor.context --reference Tss --rename mappeaks_:,_compiled:
#python mapBinding.py --path ~/meTRN --organism hs --mode map:regions --peaks hs_selection_com_cx_raw --infile in2shape_hs_modencode_KMM_hl.bed --name kundaje.states.hl --queries auto --target feature --label factor.context --reference Tss --rename mappeaks_:,_compiled:
#python mapBinding.py --path ~/meTRN --organism hs --mode map:regions --peaks hs_selection_com_cx_raw --infile in2shape_hs_modencode_KMM_k5.bed --name kundaje.states.k5 --queries auto --target feature --label factor.context --reference Tss --rename mappeaks_:,_compiled:


# Aggregate binding into matrix, recording for each region/feature the classifications across diverse comparisons:
#python mapBinding.py --path ~/meTRN --organism hs --mode map:compile --peaks hs_selection_com_cx_raw --name genomic.states.xx --queries auto --target feature --label factor.context --reference 1_Pro --exclude 17_Unmap --rename mappeaks_:,_compiled: --source genomic.states.gm,genomic.states.h1 --index regions --group P:1_Pro-E:2_Enh1,3_Enh2 --min 0.5 --policy score

#python mapBinding.py --path ~/meTRN --organism hs --mode map:compile --peaks hs_selection_com_cx_raw --name combine.states.xx --queries auto --target feature --label factor.context --reference TSS --rename mappeaks_:,_compiled: --source combine.states.gm,combine.states.h1,combine.states.hg,combine.states.hl,combine.states.k5 --index regions

python mapBinding.py --path ~/meTRN --organism hs --mode map:compile --peaks hs_selection_com_cx_raw --name kundaje.states.xx --queries auto --target feature --label factor.context --reference Tss --rename mappeaks_:,_compiled: --source kundaje.states.gm,kundaje.states.h1,kundaje.states.hg,kundaje.states.hl,kundaje.states.k5 --index regions --group P:Tss,TssF,PromP,PromF-E:Enh,EnhF,EnhW,EnhWF --min 0.5 --policy score

#python mapBinding.py --path ~/meTRN --organism hs --mode map:compile --peaks hs_selection_com_cx_raw --name kundaje.states.xx --queries auto --target feature --label factor.context --reference Tss --rename mappeaks_:,_compiled: --source kundaje.states.gm --index regions --group P:Tss,TssF,PromP,PromF-E:Enh,EnhF,EnhW,EnhWF --min 0.5 --policy score


#top
#bash runMaster2H.sh