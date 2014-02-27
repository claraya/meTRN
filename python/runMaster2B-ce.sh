#!/usr/bin/sh

# Calibrate adjustment sizes:
#python mapHOT.py --path ~/meTRN --organism ce --mode montecarlo --peaks ce_selection_reg_ex_raw --start 1 --stop 1 --contexts total.extended --regions ce_selection_reg_ex --shuffle ON --adjust 345 --name complete --overwrite ON

#python mapHOT.py --path ~/meTRN --organism ce --mode montecarlo --peaks ce_selection_reg_l1_raw --start 1 --stop 1 --contexts total.extended --regions ce_selection_reg_l1 --shuffle ON --adjust 300 --name complete --overwrite ON

#python mapHOT.py --path ~/meTRN --organism ce --mode montecarlo --peaks ce_selection_reg_l2_raw --start 1 --stop 1 --contexts total.extended --regions ce_selection_reg_l2 --shuffle ON --adjust 345 --name complete --overwrite ON

#python mapHOT.py --path ~/meTRN --organism ce --mode montecarlo --peaks ce_selection_reg_l3_raw --start 1 --stop 1 --contexts total.extended --regions ce_selection_reg_l3 --shuffle ON --adjust 335 --name complete --overwrite ON

#python mapHOT.py --path ~/meTRN --organism ce --mode montecarlo --peaks ce_selection_reg_l4_raw --start 1 --stop 1 --contexts total.extended --regions ce_selection_reg_l4 --shuffle ON --adjust 335 --name complete --overwrite ON


# Launch peak simulations:
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism ce --mode master:montecarlo --peaks ce_selection_reg_ex_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_ex --shuffle ON --adjust 345 --name complete --qsub scg3_simulation.sh --server ON --job exMCsim --chunks 10
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism ce --mode master:montecarlo --peaks ce_selection_reg_l1_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_l1 --shuffle ON --adjust 300 --name complete --qsub scg3_simulation.sh --server ON --job l1MCsim --chunks 10
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism ce --mode master:montecarlo --peaks ce_selection_reg_l2_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_l2 --shuffle ON --adjust 345 --name complete --qsub scg3_simulation.sh --server ON --job l2MCsim --chunks 10
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism ce --mode master:montecarlo --peaks ce_selection_reg_l3_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_l3 --shuffle ON --adjust 335 --name complete --qsub scg3_simulation.sh --server ON --job l3MCsim --chunks 10
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism ce --mode master:montecarlo --peaks ce_selection_reg_l4_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_l4 --shuffle ON --adjust 335 --name complete --qsub scg3_simulation.sh --server ON --job l4MCsim --chunks 10


# Analyze simulation densities (SCG3):
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism ce --mode master:simulation --peaks ce_selection_reg_ex_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_ex --shuffle ON --adjust 345 --name complete --qsub scg3_simulation.sh --server ON --job exMCreg --chunks 20
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism ce --mode master:simulation --peaks ce_selection_reg_l1_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_l1 --shuffle ON --adjust 300 --name complete --qsub scg3_simulation.sh --server ON --job l1MCreg --chunks 20
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism ce --mode master:simulation --peaks ce_selection_reg_l2_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_l2 --shuffle ON --adjust 345 --name complete --qsub scg3_simulation.sh --server ON --job l2MCreg --chunks 20
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism ce --mode master:simulation --peaks ce_selection_reg_l3_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_l3 --shuffle ON --adjust 335 --name complete --qsub scg3_simulation.sh --server ON --job l3MCreg --chunks 20
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism ce --mode master:simulation --peaks ce_selection_reg_l4_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_l4 --shuffle ON --adjust 335 --name complete --qsub scg3_simulation.sh --server ON --job l4MCreg --chunks 20


# Collect occupancy simulation reports (SCG3):
python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism ce --mode master:collecting --peaks ce_selection_reg_ex_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_ex --shuffle ON --adjust 345 --name complete --qsub scg3_simulation.sh --server ON --job exMCcol
python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism ce --mode master:collecting --peaks ce_selection_reg_l1_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_l1 --shuffle ON --adjust 300 --name complete --qsub scg3_simulation.sh --server ON --job l1MCcol
python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism ce --mode master:collecting --peaks ce_selection_reg_l2_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_l2 --shuffle ON --adjust 345 --name complete --qsub scg3_simulation.sh --server ON --job l2MCcol
python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism ce --mode master:collecting --peaks ce_selection_reg_l3_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_l3 --shuffle ON --adjust 335 --name complete --qsub scg3_simulation.sh --server ON --job l3MCcol
python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism ce --mode master:collecting --peaks ce_selection_reg_l4_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_l4 --shuffle ON --adjust 335 --name complete --qsub scg3_simulation.sh --server ON --job l4MCcol


# Collect density simulation reports (SCG3):
python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism ce --mode master:collecting --peaks ce_selection_reg_ex_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_ex --shuffle ON --adjust 345 --name complete --qsub scg3_simulation.sh --server ON --job exMCden --metric density
python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism ce --mode master:collecting --peaks ce_selection_reg_l1_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_l1 --shuffle ON --adjust 300 --name complete --qsub scg3_simulation.sh --server ON --job l1MCden --metric density
python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism ce --mode master:collecting --peaks ce_selection_reg_l2_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_l2 --shuffle ON --adjust 345 --name complete --qsub scg3_simulation.sh --server ON --job l2MCden --metric density
python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism ce --mode master:collecting --peaks ce_selection_reg_l3_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_l3 --shuffle ON --adjust 335 --name complete --qsub scg3_simulation.sh --server ON --job l3MCden --metric density
python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism ce --mode master:collecting --peaks ce_selection_reg_l4_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_l4 --shuffle ON --adjust 335 --name complete --qsub scg3_simulation.sh --server ON --job l4MCden --metric density


# Download simulation reports from cluster (local):
#python mapHOT.py --path ~/meTRN --organism ce --mode downloader --peaks ce_selection_reg_ex_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_ex --shuffle ON --adjust 345 --name complete --parameters claraya@snively.stanford.edu
#python mapHOT.py --path ~/meTRN --organism ce --mode downloader --peaks ce_selection_reg_l1_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_l1 --shuffle ON --adjust 300 --name complete --parameters claraya@snively.stanford.edu
#python mapHOT.py --path ~/meTRN --organism ce --mode downloader --peaks ce_selection_reg_l2_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_l2 --shuffle ON --adjust 345 --name complete --parameters claraya@snively.stanford.edu
#python mapHOT.py --path ~/meTRN --organism ce --mode downloader --peaks ce_selection_reg_l3_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_l3 --shuffle ON --adjust 335 --name complete --parameters claraya@snively.stanford.edu
#python mapHOT.py --path ~/meTRN --organism ce --mode downloader --peaks ce_selection_reg_l4_raw --start 1 --stop 1000 --contexts total.extended --regions ce_selection_reg_l4 --shuffle ON --adjust 335 --name complete --parameters claraya@snively.stanford.edu


#top