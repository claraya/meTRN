#!/usr/bin/sh

# Calibrate adjustment sizes:
#python mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism dm --mode montecarlo --peaks dm_selection_reg_ee_raw --start 1 --stop 1 --contexts total.extended --regions dm_selection_reg_ee --shuffle ON --adjust 700 --name complete --overwrite ON --server ON
#python mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism dm --mode montecarlo --peaks dm_selection_reg_le_raw --start 1 --stop 1 --contexts total.extended --regions dm_selection_reg_le --shuffle ON --adjust 2100 --name complete --overwrite ON --server ON
#python mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism dm --mode montecarlo --peaks dm_selection_reg_pp_raw --start 1 --stop 1 --contexts total.extended --regions dm_selection_reg_pp --shuffle ON --adjust 700 --name complete --overwrite ON --server ON


# Launch peak simulations:
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism dm --mode master:montecarlo --peaks dm_selection_reg_ee_raw --start 1 --stop 1000 --contexts total.extended --regions dm_selection_reg_ee --shuffle ON --adjust 700 --name complete --qsub scg3_simulation.sh --server ON --job eeMCsim --chunks 10
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism dm --mode master:montecarlo --peaks dm_selection_reg_le_raw --start 1 --stop 1000 --contexts total.extended --regions dm_selection_reg_le --shuffle ON --adjust 2100 --name complete --qsub scg3_simulation.sh --server ON --job leMCsim --chunks 10
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism dm --mode master:montecarlo --peaks dm_selection_reg_pp_raw --start 1 --stop 1000 --contexts total.extended --regions dm_selection_reg_pp --shuffle ON --adjust 700 --name complete --qsub scg3_simulation.sh --server ON --job ppMCsim --chunks 10


# Analyze simulation densities (SCG3):
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism dm --mode master:simulation --peaks dm_selection_reg_ee_raw --start 1 --stop 1000 --contexts total.extended --regions dm_selection_reg_ee --shuffle ON --adjust 700 --name complete --qsub scg3_simulation.sh --server ON --job eeMCreg --chunks 20
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism dm --mode master:simulation --peaks dm_selection_reg_le_raw --start 1 --stop 1000 --contexts total.extended --regions dm_selection_reg_le --shuffle ON --adjust 2100 --name complete --qsub scg3_simulation.sh --server ON --job leMCreg --chunks 20
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism dm --mode master:simulation --peaks dm_selection_reg_pp_raw --start 1 --stop 1000 --contexts total.extended --regions dm_selection_reg_pp --shuffle ON --adjust 700 --name complete --qsub scg3_simulation.sh --server ON --job ppMCreg --chunks 20


# Collect occupancy simulation reports (SCG3):
python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism dm --mode master:collecting --peaks dm_selection_reg_ee_raw --start 1 --stop 1000 --contexts total.extended --regions dm_selection_reg_ee --shuffle ON --adjust 700 --name complete --qsub scg3_simulation.sh --server ON --job eeMCcol
python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism dm --mode master:collecting --peaks dm_selection_reg_le_raw --start 1 --stop 1000 --contexts total.extended --regions dm_selection_reg_le --shuffle ON --adjust 2100 --name complete --qsub scg3_simulation.sh --server ON --job leMCcol
python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism dm --mode master:collecting --peaks dm_selection_reg_pp_raw --start 1 --stop 1000 --contexts total.extended --regions dm_selection_reg_pp --shuffle ON --adjust 700 --name complete --qsub scg3_simulation.sh --server ON --job ppMCcol


# Collect density simulation reports (SCG3):
python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism dm --mode master:collecting --peaks dm_selection_reg_ee_raw --start 1 --stop 1000 --contexts total.extended --regions dm_selection_reg_ee --shuffle ON --adjust 700 --name complete --qsub scg3_simulation.sh --server ON --job eeMCden --metric density
python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism dm --mode master:collecting --peaks dm_selection_reg_le_raw --start 1 --stop 1000 --contexts total.extended --regions dm_selection_reg_le --shuffle ON --adjust 2100 --name complete --qsub scg3_simulation.sh --server ON --job leMCden --metric density
python /srv/gs1/projects/snyder/claraya/meTRN/python/mapHOT.py --path /srv/gs1/projects/snyder/claraya/meTRN --organism dm --mode master:collecting --peaks dm_selection_reg_pp_raw --start 1 --stop 1000 --contexts total.extended --regions dm_selection_reg_pp --shuffle ON --adjust 700 --name complete --qsub scg3_simulation.sh --server ON --job ppMCden --metric density


# Download simulation reports from cluster (local):
#python mapHOT.py --path ~/meTRN --organism dm --mode downloader --peaks dm_selection_reg_ee_raw --start 1 --stop 1000 --contexts total.extended --regions dm_selection_reg_ee --shuffle ON --adjust 700 --name complete --parameters claraya@snively.stanford.edu
#python mapHOT.py --path ~/meTRN --organism dm --mode downloader --peaks dm_selection_reg_le_raw --start 1 --stop 1000 --contexts total.extended --regions dm_selection_reg_le --shuffle ON --adjust 2100 --name complete --parameters claraya@snively.stanford.edu
#python mapHOT.py --path ~/meTRN --organism dm --mode downloader --peaks dm_selection_reg_pp_raw --start 1 --stop 1000 --contexts total.extended --regions dm_selection_reg_pp --shuffle ON --adjust 700 --name complete --parameters claraya@snively.stanford.edu


#top