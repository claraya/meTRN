#!/usr/bin/sh



# Calibrate adjustment sizes:
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:montecarlo --peaks hs_selection_reg_gm_raw --start 1 --stop 1 --contexts total.extended --regions hs_selection_reg_gm --shuffle ON --adjust 197 --name complete --qsub scg3_simulation.sh --server ON --job gmMCsim --chunks 10 --overwrite ON
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:montecarlo --peaks hs_selection_reg_h1_raw --start 1 --stop 1 --contexts total.extended --regions hs_selection_reg_h1 --shuffle ON --adjust 240 --name complete --qsub scg3_simulation.sh --server ON --job h1MCsim --chunks 10 --overwrite ON
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:montecarlo --peaks hs_selection_reg_hg_raw --start 1 --stop 1 --contexts total.extended --regions hs_selection_reg_hg --shuffle ON --adjust 221 --name complete --qsub scg3_simulation.sh --server ON --job hgMCsim --chunks 10 --overwrite ON
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:montecarlo --peaks hs_selection_reg_hl_raw --start 1 --stop 1 --contexts total.extended --regions hs_selection_reg_hl --shuffle ON --adjust 220 --name complete --qsub scg3_simulation.sh --server ON --job hlMCsim --chunks 10 --overwrite ON
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:montecarlo --peaks hs_selection_reg_k5_raw --start 1 --stop 1 --contexts total.extended --regions hs_selection_reg_k5 --shuffle ON --adjust 210 --name complete --qsub scg3_simulation.sh --server ON --job k5MCsim --chunks 10 --overwrite ON


# Launch peak simulations:
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:montecarlo --peaks hs_selection_reg_gm_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_gm --shuffle ON --adjust 197 --name complete --qsub scg3_simulation.sh --server ON --job gmMCsim --chunks 10
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:montecarlo --peaks hs_selection_reg_h1_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_h1 --shuffle ON --adjust 240 --name complete --qsub scg3_simulation.sh --server ON --job h1MCsim --chunks 10
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:montecarlo --peaks hs_selection_reg_hg_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_hg --shuffle ON --adjust 221 --name complete --qsub scg3_simulation.sh --server ON --job hgMCsim --chunks 10
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:montecarlo --peaks hs_selection_reg_hl_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_hl --shuffle ON --adjust 220 --name complete --qsub scg3_simulation.sh --server ON --job hlMCsim --chunks 10
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:montecarlo --peaks hs_selection_reg_k5_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_k5 --shuffle ON --adjust 210 --name complete --qsub scg3_simulation.sh --server ON --job k5MCsim --chunks 10


# Analyze simulation densities (SCG3):
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:simulation --peaks hs_selection_reg_gm_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_gm --shuffle ON --adjust 197 --name complete --qsub scg3_simulation.sh --server ON --job gmMCreg --chunks 20
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:simulation --peaks hs_selection_reg_h1_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_h1 --shuffle ON --adjust 240 --name complete --qsub scg3_simulation.sh --server ON --job h1MCreg --chunks 20
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:simulation --peaks hs_selection_reg_hg_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_hg --shuffle ON --adjust 221 --name complete --qsub scg3_simulation.sh --server ON --job hgMCreg --chunks 20
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:simulation --peaks hs_selection_reg_hl_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_hl --shuffle ON --adjust 220 --name complete --qsub scg3_simulation.sh --server ON --job hlMCreg --chunks 20
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:simulation --peaks hs_selection_reg_k5_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_k5 --shuffle ON --adjust 210 --name complete --qsub scg3_simulation.sh --server ON --job k5MCreg --chunks 20


# Collect occupancy simulation reports (SCG3):
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:collecting --peaks hs_selection_reg_gm_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_gm --shuffle ON --adjust 197 --name complete --qsub scg3_simulation.sh --server ON --job gmMCcol
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:collecting --peaks hs_selection_reg_h1_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_h1 --shuffle ON --adjust 240 --name complete --qsub scg3_simulation.sh --server ON --job h1MCcol
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:collecting --peaks hs_selection_reg_hg_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_hg --shuffle ON --adjust 221 --name complete --qsub scg3_simulation.sh --server ON --job hgMCcol
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:collecting --peaks hs_selection_reg_hl_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_hl --shuffle ON --adjust 220 --name complete --qsub scg3_simulation.sh --server ON --job hlMCcol
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:collecting --peaks hs_selection_reg_k5_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_k5 --shuffle ON --adjust 210 --name complete --qsub scg3_simulation.sh --server ON --job k5MCcol


# Collect density simulation reports (SCG3):
python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:collecting --peaks hs_selection_reg_gm_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_gm --shuffle ON --adjust 197 --name complete --qsub scg3_simulation.sh --server ON --job gmMCden --metric density
python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:collecting --peaks hs_selection_reg_h1_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_h1 --shuffle ON --adjust 240 --name complete --qsub scg3_simulation.sh --server ON --job h1MCden --metric density
python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:collecting --peaks hs_selection_reg_hg_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_hg --shuffle ON --adjust 221 --name complete --qsub scg3_simulation.sh --server ON --job hgMCden --metric density
python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:collecting --peaks hs_selection_reg_hl_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_hl --shuffle ON --adjust 220 --name complete --qsub scg3_simulation.sh --server ON --job hlMCden --metric density
python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism hs --mode master:collecting --peaks hs_selection_reg_k5_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_k5 --shuffle ON --adjust 210 --name complete --qsub scg3_simulation.sh --server ON --job k5MCden --metric density


# Download simulation reports from cluster (local):
#python mapHOT.py --path ~/meTRN --organism hs --mode downloader --peaks hs_selection_reg_gm_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_gm --shuffle ON --adjust 225 --name complete --parameters claraya@snively.stanford.edu
#python mapHOT.py --path ~/meTRN --organism hs --mode downloader --peaks hs_selection_reg_h1_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_h1 --shuffle ON --adjust 225 --name complete --parameters claraya@snively.stanford.edu
#python mapHOT.py --path ~/meTRN --organism hs --mode downloader --peaks hs_selection_reg_hg_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_hg --shuffle ON --adjust 225 --name complete --parameters claraya@snively.stanford.edu
#python mapHOT.py --path ~/meTRN --organism hs --mode downloader --peaks hs_selection_reg_hl_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_hl --shuffle ON --adjust 225 --name complete --parameters claraya@snively.stanford.edu
#python mapHOT.py --path ~/meTRN --organism hs --mode downloader --peaks hs_selection_reg_k5_raw --start 1 --stop 1000 --contexts total.extended --regions hs_selection_reg_k5 --shuffle ON --adjust 225 --name complete --parameters claraya@snively.stanford.edu



#top