#!/usr/bin/sh


# Launch SOM analyses on worm binding matrixes (Stanford):
#python /srv/gs1/projects/snyder/claraya/meTRN/python/mapSOM.py --path /srv/gs1/projects/snyder/claraya/meTRN --mode launch --peakflag raw --technique binary --iterations 100 --qsub scg3_neurons.sh --server ON --job rawBinary


# Launch SOM analyses on worm binding matrixes (GS; second run includes selected samples):
#python /net/fields/vol1/home/araya/meTRN/python/mapSOM.py --path /net/fields/vol1/home/araya/meTRN --mode launch --peakflag raw --technique binary --iterations 100 --qsub uwgs_neurons.sh --server GS --job rawBinary
#python /net/fields/vol1/home/araya/meTRN/python/mapSOM.py --path /net/fields/vol1/home/araya/meTRN --mode launch --peakflag raw --technique binary --iterations 100 --qsub uwgs_neurons.sh --server GS --job rawSample


# Recover SOM analyses on worm binding matrixes:
#python mapSOM.py --path ~/meTRN --mode recover
#python mapSOM.py --path ~/meTRN --mode examine


# Explore subset SOMs (where factors have been removed):
python mapSOM.py --path ~/meTRN --mode filter --peakflag raw --technique binary --input E1C:fil.e1c.XXX --iterations 10 --threads 4
python mapSOM.py --path ~/meTRN --mode filter --peakflag raw --technique binary --input 1V2:fil.1v2.XXX --iterations 10 --threads 4
python mapSOM.py --path ~/meTRN --mode filter --peakflag raw --technique binary --input 2V3:fil.2v3.XXX --iterations 10 --threads 4
python mapSOM.py --path ~/meTRN --mode filter --peakflag raw --technique binary --input 3V4:fil.3v4.XXX --iterations 10 --threads 4

python mapSOM.py --path ~/meTRN --mode remove --peakflag raw --technique binary --input E1C:rem.e1c.XXX --iterations 15 --threads 4
python mapSOM.py --path ~/meTRN --mode remove --peakflag raw --technique binary --input 1V2:rem.1v2.XXX --iterations 17 --threads 4
python mapSOM.py --path ~/meTRN --mode remove --peakflag raw --technique binary --input 2V3:rem.2v3.XXX --iterations 22 --threads 4
python mapSOM.py --path ~/meTRN --mode remove --peakflag raw --technique binary --input 3V4:rem.3v4.XXX --iterations 26 --threads 4



#####ELC:100.elc.som,E1C:100.e1c.som,1V2:100.1v2.som,2V3:100.2v3.som,3V4:100.3v4.som



#top
#bash runMaster7A.sh