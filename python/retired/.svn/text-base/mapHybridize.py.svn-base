#!/usr/bin/env python
# organize orthologs and files!

import sys
import time
import optparse
import general
import numpy
import pickle
import pdb
import cetrn
import modencode
import multiprocessing
import itertools
import os

print "Command:", " ".join(sys.argv)
print "Timestamp:", time.asctime(time.localtime())

def findOrganisms(infile, organismTags):
	organisms = list()
	for organismTag in organismTags:
		if organismTag in infile.split("_"):
			organisms.append(organismTag)
	return organisms

def buildOrthologs(infile):
	orthologs_dict = dict()
	inlines = open(infile).readlines()
	for inline in inlines:
		species, ortholog, factor =  inline.strip().split(":")
		if not species in orthologs_dict:
			orthologs_dict[species] = dict()
		orthologs_dict[species][ortholog] = factor
	return orthologs_dict

def loader(xfile, xhd, target="factor", header=True):
	xdict = dict()
	xlines = open(xfile).readlines()
	if header:
		xlines.pop(0)
	for xline in xlines:
		xitems = xline.strip().split("\t")
		branch = xitems[xhd[target]]
		if not branch in xdict:
			xdict[branch] = list()
		xdict[branch].append(xline.strip())
	return xdict

def main():
	
	parser = optparse.OptionParser()
	parser.add_option("--path", action = "store", type = "string", dest = "path", help = "path from script to files")
	parser.add_option("--mode", action = "store", type = "string", dest = "mode", help = "expand", default="OFF")
	parser.add_option("--analysis", action = "store", type = "string", dest = "analysis", help = "family, orthologs, paralogs")
	parser.add_option("--organism", action = "store", type = "string", dest = "organism", help = "which organism is this for?", default="OFF")
	parser.add_option("--source", action = "store", type = "string", dest = "source", help = "organism-version codename", default="OFF")
	parser.add_option("--peaks", action = "store", type = "string", dest = "peaks", help = "idr, blacklist, filtered", default="OFF")
	parser.add_option("--name", action = "store", type = "string", dest = "name", help = "output name?", default="OFF")
	parser.add_option("--A", action = "store", type = "string", dest = "a", help = "paths to files of interest", default="OFF")
	parser.add_option("--B", action = "store", type = "string", dest = "b", help = "files to be hybridized", default="OFF")
	parser.add_option("--comparison", action = "store", type = "string", dest = "comparison", help = "which organisms are to be compared?", default="OFF")
	parser.add_option("--cutoff", action = "store", type = "float", dest = "cutoff", help = "p-value cutoff for fraction calculation", default=0.05)
	(option, args) = parser.parse_args()
	
	# define organisms:
	organismTags = ["hs","mm","ce","dm"]
	
	# define organism parameters:
	if option.organism == "h.sapiens" or option.organism == "human" or option.organism == "hs":
		organismTag = "hs"
		contextTag = "cells"
		idColumns = ["name", "code", "hgcn","ensembl"]
		idComplexList = list()
	elif option.organism == "m.musculus" or option.organism == "mouse" or option.organism == "mm":
		organismTag = "mm"
		contextTag = "cells"
		idColumns = ["name", "code", "hgcn","ensembl"]
		idComplexList = list()
	elif option.organism == "c.elegans" or option.organism == "worm" or option.organism == "ce":
		organismTag = "ce"
		contextTag = "stage"
		idColumns = ["name", "code", "wormbase","ensembl"]
		idComplexList = list()
	elif option.organism == "d.melanogaster" or option.organism == "fly" or option.organism == "dm":
		organismTag = "dm"
		contextTag = "stage"
		idColumns = ["name", "code", "flybase","ensembl"]
		idComplexList = ["dataset",":","url"]
	
	# define version tag:
	if "," in option.source:
		source1, source2 = option.source.split(",")
		versionTag1 = "v" + source1[2:]
		versionTag2 = "v" + source2[2:]
	
	# import paths:
	path_dict = modencode.configBuild(option.path + "/input/" + "configure_path.txt")
	
	# specify input and output paths:
	inpath = path_dict["input"]
	extraspath = path_dict["extras"]
	pythonpath = path_dict["python"]
	downloadpath = path_dict["download"]
	fastqpath = path_dict["fastq"]
	bowtiepath = path_dict["bowtie"]
	bwapath = path_dict["bwa"]
	macspath = path_dict["macs"]
	memepath = path_dict["meme"]
	idrpath = path_dict["idr"]
	igvpath = path_dict["igv"]
	testpath = path_dict["test"]
	processingpath = path_dict["processing"]
	annotationspath = path_dict["annotations"]
	coassociationspath = path_dict["coassociations"]
	neuronspath = path_dict["neurons"]
	orthologspath = path_dict["orthologs"]
	peakspath = path_dict["peaks"]
	gopath = path_dict["go"]
	hotpath = path_dict["hot"]
	
	# update analysis path:
	if option.analysis == "families":
		analysispath = orthologspath
		outheader = ["family.id", "species.a", "species.b", "gene.a", "gene.b"]
		matchTag = "family.txt"
	elif option.analysis == "orthologs":
		analysispath = orthologspath + "orthologs/"
		outheader = ["family.id", "species.a", "species.b", "gene.a", "gene.b", "count.a", "count.b"]
		matchTag = "orthologs.txt"
	elif option.analysis == "paralogs":
		analysispath = orthologspath + "paralogous/"
		outheader = ["family.id", "species", "gene.a", "gene.b"]
		matchTag = "paralog.txt"
	
	# define P-value cutoff handle:
	pvaluecutoff_handle = "%.0e" % (float(option.cutoff))
	
	# grab files (other) mode:
	if option.mode == "go":
			
		# import orthologs dictionary:
		orthologs_dict = buildOrthologs(inpath + "configure_orthologs_" + option.analysis + ".txt")
		
		# find species-comparison orthologs:
		aspecies, bspecies = option.comparison.split(",")
		comparison_files = os.listdir(orthologspath + "expanded/" + option.analysis)
		for comparison_file in comparison_files:
			comparisonOrganisms = findOrganisms(comparison_file, organismTags)
			if sorted([aspecies, bspecies]) == sorted(comparisonOrganisms):
				ortholog_file = comparison_file
		
		# make output folders:
		if not "comparative" in os.listdir(gopath):
			command = "mkdir " + gopath + "comparative/"
			os.system(command)
			
		if not aspecies in os.listdir(gopath + "comparative/"):
			command = "mkdir " + gopath + "comparative/" + aspecies
			os.system(command)
		
		if not bspecies in os.listdir(gopath + "comparative/" + aspecies):
			command = "mkdir " + gopath + "comparative/" + aspecies + "/" + bspecies + "/"
			os.system(command)
		
		# make output file:
		f_output = open(gopath + "comparative/" + aspecies + "/" + bspecies + "/maphybridize_go_comparison_p" + pvaluecutoff_handle + ".txt", "w")
		
		# load GO lines from input files:
		orthologs = list()
		asublines, bsublines = list(), list()
		ahd = general.build_header_dict(option.a)
		bhd = general.build_header_dict(option.b)
		adict = loader(option.a, ahd)
		bdict = loader(option.b, bhd)
		
		# generate a-file and b-file headers, as well as output header:
		aheader, bheader = list(), list()
		for header in general.valuesort(ahd):
			aheader.append(header + ".a")
		for header in general.valuesort(bhd):
			bheader.append(header + ".b")
		outheader = ["i", "j", "items.a", "items.b", "overlap.a", "overlap.b", "overlap.avg", "overlap.sum", "overlap.max", "overlap.count", "items.count"] #"a.only.goids", "b.only.goids", "overlap.goids"]
		print >>f_output, "\t".join(outheader)
		
		# prefilter goids:
		print
		print "Finding shared GO ids..."
		gxids, axids, bxids = list(), list(), list()
		ghits, ahits, bhits = list(), list(), list()
		for afactor in adict:
			for aline in adict[afactor]:
				aitems = aline.strip().split("\t")
				dataset, strain, factor, stage, institute, method = aitems[ahd["dataset"]], aitems[ahd["strain"]], aitems[ahd["factor"]], aitems[ahd["stage"]], aitems[ahd["institute"]], aitems[ahd["method"]]
				goid, goterm, gocount, pvalue = aitems[ahd["go.id"]], aitems[ahd["go.term"]], aitems[ahd["go.count"]], aitems[ahd["adj.pvalue"]]
				if float(pvalue) < option.cutoff: # and int(gocount) > 50 and int(gocount) < 500 :
					ahits.append(goid)
					ghits.append(goid)
				axids.append(goid)
				gxids.append(goid)
				
		for bfactor in bdict:
			for bline in bdict[bfactor]:
				bitems = bline.strip().split("\t")
				dataset, strain, factor, stage, institute, method = bitems[bhd["dataset"]], bitems[bhd["strain"]], bitems[bhd["factor"]], bitems[bhd["stage"]], bitems[bhd["institute"]], bitems[bhd["method"]]
				goid, goterm, gocount, pvalue = bitems[bhd["go.id"]], bitems[bhd["go.term"]], bitems[bhd["go.count"]], bitems[bhd["adj.pvalue"]]
				if float(pvalue) < option.cutoff: # and int(gocount) > 50 and int(gocount) < 500 :
					bhits.append(goid)
					ghits.append(goid)
				bxids.append(goid)
				gxids.append(goid)
					
		gxids, axids, bxids = list(set(gxids)), list(set(axids)), list(set(bxids))
		ghits, ahits, bhits = list(set(ghits)), list(set(ahits)), list(set(ahits))
		#dbhits = list(set(axids).intersection(set(bxids)))
		dbhits = list(set(ahits).intersection(set(bhits)))
		#dbhits = ghits
		
		# load GO info for orthologs:
		print "Loading GO information..."
		ainfo, binfo, comparisons, goids, aoids, boids = dict(), dict(), list(), list(), list(), list()
		ihd = general.build_header_dict(orthologspath + "expanded/" + option.analysis + "/" + ortholog_file)
		inlines = open(orthologspath + "expanded/" + option.analysis + "/" + ortholog_file).readlines()
		inlines.pop(0)
		for inline in inlines:
			initems = inline.strip().split("\t")
			agene, bgene = initems[ihd["gene.a"]], initems[ihd["gene.b"]]
			if agene in orthologs_dict[aspecies] and bgene in orthologs_dict[bspecies]:
				afactor, bfactor = orthologs_dict[aspecies][agene], orthologs_dict[bspecies][bgene]
				
				if afactor in adict and bfactor in bdict:
					
					print afactor, bfactor
					comparisons.append([afactor, bfactor])
					
					for aline in adict[afactor]:
						aitems = aline.strip().split("\t")
						dataset, strain, factor, stage, institute, method = aitems[ahd["dataset"]], aitems[ahd["strain"]], aitems[ahd["factor"]], aitems[ahd["stage"]], aitems[ahd["institute"]], aitems[ahd["method"]]
						goid, goterm, gocount, pvalue = aitems[ahd["go.id"]], aitems[ahd["go.term"]], aitems[ahd["go.count"]], aitems[ahd["adj.pvalue"]]
						if float(pvalue) < option.cutoff: # and int(gocount) > 50 and int(gocount) < 500 :
							if not afactor in ainfo:
								ainfo[afactor] = dict()
							if not stage in ainfo[afactor]:
								ainfo[afactor][stage] = dict()
							ainfo[afactor][stage][goid] = [dataset, strain, factor, stage, institute, method, goid, goterm, pvalue]
							aoids.append(goid)
							goids.append(goid)
						
					for bline in bdict[bfactor]:
						bitems = bline.strip().split("\t")
						dataset, strain, factor, stage, institute, method = bitems[bhd["dataset"]], bitems[bhd["strain"]], bitems[bhd["factor"]], bitems[bhd["stage"]], bitems[bhd["institute"]], bitems[bhd["method"]]
						goid, goterm, gocount, pvalue = bitems[bhd["go.id"]], bitems[bhd["go.term"]], bitems[bhd["go.count"]], bitems[bhd["adj.pvalue"]]
						if float(pvalue) < option.cutoff: # and int(gocount) > 50 and int(gocount) < 500 :
							if not bfactor in binfo:
								binfo[bfactor] = dict()
							if not stage in binfo[bfactor]:
								binfo[bfactor][stage] = dict()
							binfo[bfactor][stage][goid] = [dataset, strain, factor, stage, institute, method, goid, goterm, pvalue]
							boids.append(goid)
							goids.append(goid)
		
		apass, bpass = list(), list()
		for comparison in comparisons:
			aquery, bquery = comparison
			if aquery in ainfo and bquery in binfo:
				apass.append(aquery)
				bpass.append(bquery)
				if bquery == "Hr78":
					print aquery, bquery
		
		#goids = list(set(dbhits))
		goids, aoids, boids = list(set(gxids)), list(set(axids)), list(set(bxids))
		goids = set(aoids).intersection(set(boids))
		universe = len(goids)
		matrix = dict()
		print "Calculating GO-term overlaps..."
		for afactor in sorted(ainfo.keys()):
			for astage in ainfo[afactor]:
				adataset = afactor + "." + astage
				if not adataset in matrix:
					matrix[adataset] = dict()
				
				for bfactor in sorted(binfo.keys()):
					for bstage in binfo[bfactor]:
						bdataset = bfactor + "." + bstage		
						
						if afactor in apass and bfactor in bpass:
						
							aids = set(ainfo[afactor][astage].keys()).intersection(set(goids))
							bids = set(binfo[bfactor][bstage].keys()).intersection(set(goids))
							aonly = set(aids).difference(set(bids))
							bonly = set(bids).difference(set(aids))
							overlap = set(aids).intersection(set(bids))
							total = set(aids).union(set(bids))
							if len(overlap) == 0:
								aoverlap, boverlap, overlap_avg, overlap_max, overlap_sum = 0, 0, 0, 0, 0
							else:
								aoverlap = float(len(overlap))/len(aids)
								boverlap = float(len(overlap))/len(bids)
								overlap_avg = numpy.mean([aoverlap, boverlap])
								overlap_max = max([aoverlap, boverlap])
								overlap_sum = float(len(overlap))/len(total)
							output = [adataset, bdataset, len(aids), len(bids), aoverlap, boverlap, overlap_avg, overlap_sum, overlap_max, len(overlap), universe] #",".join(sorted(list(aonly))), ",".join(sorted(list(bonly))), ",".join(sorted(list(overlap)))]
							matrix[adataset][bdataset] = output
							
							
							#print afactor, aonly
							#print bfactor, bonly
							#print overlap
							#for goid in aids:
							#	print goid
							#	command = "grep " + '"' + goid + '"' + " " + option.a + " | grep " + '"' + afactor + '"' + " | grep " + '"' + astage + '"'
							#	os.system(command)
							#	command = "grep " + '"' + goid + '"' + " " + option.b + " | grep " + '"' + bfactor + '"' + " | grep " + '"' + bstage + '"'
							#	os.system(command)
							#	print
							#pdb.set_trace()
		
		for adataset in matrix:
			for bdataset in matrix[adataset]:
				print >>f_output, "\t".join(map(str, matrix[adataset][bdataset]))
		
		print
		# close output file:
		f_output.close()
		
		
			
if __name__ == "__main__":
	main()
	
	print "Completed:", time.asctime(time.localtime())

#python mapHybridize.py --path ~/meTRN --mode go --analysis families --source ce1,dm1 --peaks idr --comparison ce,dm --A ~/ceTRN/data/go/optimal_standard_totals_sx_filter/p5e-1/summary/mapgo_complete_optimal_standard_totals_sx_filter_p5_hc1_hp1e-02_matrix --B ~/dmTRN/data/go/optimal_standard_totals_sx_filter/p5e-1/summary/mapgo_complete_optimal_standard_totals_sx_filter_p5_hc1_hp1e-02_matrix --cutoff 0.01

#python mapHybridize.py --path ~/meTRN --mode go --analysis families --source ce1,dm1 --peaks idr --comparison ce,dm --A ~/ceTRN/data/go/optimal_standard_totals_sx_filter/p5e-1/summary/mapgo_complete_optimal_standard_totals_sx_filter_p5_hc1_hp1e-02_matrix --B ~/dmTRN/data/go/optimal_standard_totals_sx_filter/p5e-1/summary/mapgo_complete_optimal_standard_totals_sx_filter_p5_hc1_hp1e-02_matrix --cutoff 0.05

#python mapHybridize.py --path ~/meTRN --mode go --analysis families --source ce1,dm1 --peaks idr --comparison ce,dm --A ~/ceTRN/data/go/optimal_standard_totals_sx_filter/p5e-1/summary/mapgo_complete_optimal_standard_totals_sx_filter_p5_hc1_hp1e-02_matrix --B ~/dmTRN/data/go/optimal_standard_totals_sx_filter/p5e-1/summary/mapgo_complete_optimal_standard_totals_sx_filter_p5_hc1_hp1e-02_matrix --cutoff 0.1
