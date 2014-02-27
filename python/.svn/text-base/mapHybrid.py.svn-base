#!/usr/bin/env python
# generate peak set complete files, binding region files, and report files!

import sys
import time
import optparse
import general
import numpy
import hyper
import pickle
import pdb
import metrn
import modencode
import os

from scipy import stats


print "Command:", " ".join(sys.argv)
print "Timestamp:", time.asctime(time.localtime())

	
def main():
	
	parser = optparse.OptionParser()
	parser.add_option("--path", action = "store", type = "string", dest = "path", help = "Path from script to files")
	parser.add_option("--peaks", action = "store", type = "string", dest = "peaks", help = "Peaks to be used for analysis")
	parser.add_option("--mode", action = "store", type = "string", dest = "mode", help = "Operations to be performed")
	parser.add_option("--organism", action = "store", type = "string", dest = "organism", help = "Target organism for operations...", default="OFF")
	parser.add_option("--species", action = "store", type = "string", dest = "species", help = "Species to compare", default="OFF")
	parser.add_option("--orthology", action = "store", type = "string", dest = "orthology", help = "Use 'direct', 'family' (Yong's), or 'group' (Pouya's) orthologs?", default="direct")
	parser.add_option("--nametag", action = "store", type = "string", dest = "nametag", help = "Orthology nametag: nametagHsCe", default="ortho")
	parser.add_option("--commonNames", action = "store", type = "string", dest = "commonNames", help = "Grab common names file?", default="ON")
	parser.add_option("--familyFiles", action = "store", type = "string", dest = "familyFiles", help = "Grab cleaned files?", default="formatted")
	parser.add_option("--target", action = "store", type = "string", dest = "target", help = "Target identification", default="OFF")
	parser.add_option("--label", action = "store", type = "string", dest = "label", help = "How should labels be generated?", default="rebuild")
	parser.add_option("--indexes", action = "store", type = "string", dest = "indexes", help = "Indexes for matrix construction...", default="OFF")
	parser.add_option("--values", action = "store", type = "string", dest = "values", help = "Values for matrix construction...", default="OFF")
	parser.add_option("--contexts", action = "store", type = "string", dest = "contexts", help = "What contexts of development should I track?", default="OFF")
	parser.add_option("--source", action = "store", type = "string", dest = "source", help = "Path to source files...", default="OFF")
	parser.add_option("--server", action = "store", type = "string", dest = "server", help = "Are we on the server?", default="OFF")
	parser.add_option("--name", action = "store", type = "string", dest = "name", help = "Output name?", default="OFF")
	parser.add_option("--A", action = "store", type = "string", dest = "a", help = "Paths to files of interest", default="OFF")
	parser.add_option("--B", action = "store", type = "string", dest = "b", help = "Files to be hybridized", default="OFF")
	parser.add_option("--cutoff", action = "store", type = "float", dest = "cutoff", help = "P-value cutoff for fraction calculation", default=0.05)
	parser.add_option("--fraction", action = "store", type = "int", dest = "fraction", help = "Tolerated fractional representation of GO terms", default=10)
	parser.add_option("--rename", action = "store", type = "string", dest = "rename", help = "Targets to rename. Comma-separated list of 'target:replacement' pairs to search and replace.", default="OFF")
	(option, args) = parser.parse_args()
	
	# import paths:
	if option.server == "OFF":
		path_dict = modencode.configBuild(option.path + "/input/" + "configure_path.txt")
	elif option.server == "ON":
		path_dict = modencode.configBuild(option.path + "/input/" + "configure_server.txt")
	
	# specify input and output paths:
	inpath = path_dict["input"]
	extraspath = path_dict["extras"]
	pythonpath = path_dict["python"]
	scriptspath = path_dict["scripts"]
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
	orthologspath = path_dict["orthologs"]
	coassociationspath = path_dict["coassociations"]
	cellspath = path_dict["cells"]
	peakspath = path_dict["peaks"]
	gopath = path_dict["go"]
	hotpath = path_dict["hot"]
	qsubpath = path_dict["qsub"]
	
	# standardize paths for analysis:
	alignerpath = bwapath
	indexpath = alignerpath + "index/"
	alignmentpath = alignerpath + "alignment/"
	qcfilterpath = alignerpath + "qcfilter/"
	qcmergepath = alignerpath + "qcmerge/"
	
	# import configuration dictionaries:
	source_dict = modencode.configBuild(inpath + "configure_source.txt")
	method_dict = modencode.configBuild(inpath + "configure_method.txt")
	context_dict = modencode.configBuild(inpath + "configure_context.txt")
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
	
	# update analysis path:
	#if option.analysis == "families":
	#	analysispath = orthologspath
	#	outheader = ["family.id", "species.a", "species.b", "gene.a", "gene.b"]
	#	matchTag = "family.txt"
	#elif option.analysis == "orthologs":
	#	analysispath = orthologspath + "orthologs/"
	#	outheader = ["family.id", "species.a", "species.b", "gene.a", "gene.b", "count.a", "count.b"]
	#	matchTag = "orthologs.txt"
	#elif option.analysis == "paralogs":
	#	analysispath = orthologspath + "paralogous/"
	#	outheader = ["family.id", "species", "gene.a", "gene.b"]
	#	matchTag = "paralog.txt"
	
	# define P-value cutoff handle:
	pvaluecutoff_handle = "%.0e" % (float(option.cutoff))
	
	# merge (datatypes) matrix mode:
	if option.mode == "merge.direct":
		
		# find species-comparison orthologs:
		speciesTags = option.species.split(",")
		aspecies, bspecies = speciesTags
		
		# make comparison output folders:
		comparisonpath = path_dict[option.source] + "comparison/" + aspecies + "/" + bspecies + "/"
		general.pathGenerator(comparisonpath)
		
		# define input files:
		ainfile = str(path_dict[option.source] + "/" + option.a).replace("//","/")
		binfile = str(path_dict[option.source] + "/" + option.b).replace("//","/")
		
		# find target matrix indexes (keys):
		aindexes, bindexes = option.indexes.split(",")
		ai, aj = aindexes.split(":")
		bi, bj = aindexes.split(":")
		
		# find target matrix values:
		ax, bx = option.values.split(",")
		
		# load input matrixes:
		amatrix = general.matrixBuilder(i=ai, j=aj, x=ax, infile=ainfile, datatype="float")
		bmatrix = general.matrixBuilder(i=bi, j=bj, x=bx, infile=binfile, datatype="float")
		
		#print amatrix.keys()
		#print bmatrix.keys()
		#print set(amatrix.keys()).intersection(set(bmatrix.keys()))
		#pdb.set_trace()
		
		# define output file:
		f_outfile = comparisonpath + "maphybrid_" + option.source + "_" + option.name + "_combined.txt"
		f_output = open(f_outfile, "w")
		print >>f_output, "\t".join(["i", "j", "a.value", "b.value", "difference", "log2.ratio"])
		
		# merge matrixes:
		for ai in amatrix:
			for aj in amatrix:
				if ai in bmatrix and aj in bmatrix:
					output = [ai, aj, amatrix[ai][aj], bmatrix[ai][aj], amatrix[ai][aj]-bmatrix[ai][aj], numpy.log2(amatrix[ai][aj]/bmatrix[ai][aj])]
					print >>f_output, "\t".join(map(str, output))
		f_output.close()
		
	
	# merge ortholog matrix mode:
	if option.mode == "merge.matrix":
			
		# import orthologs dictionary:
		#orthologs_dict = buildOrthologs(inpath + "configure_orthologs_" + option.analysis + ".txt")
		
		# find species-comparison orthologs:
		speciesTags = option.species.split(",")
		aspecies, bspecies = speciesTags
		
		# generate output peaks name:
		orthologTag = option.nametag + metrn.orthologLabel(aspecies, speciesTags)
		
		# define orthology path:
		if option.orthology == "direct":
			orthologypath = orthologspath + "orthologs/"
		elif option.orthology == "family":
			orthologypath = orthologspath + "families/"
		elif option.orthology == "groups":
			orthologypath = orthologspath + "groups/"
			
		# generate orthology dictionary:
		ortholog_dict = metrn.orthologBuilder(speciesTags, path=orthologypath, orthology=option.orthology, commonNames=option.commonNames, familyFiles=option.familyFiles, verbose="OFF")
		
		# target specie orthologs:
		aorthologs = metrn.orthologFinder(aspecies, speciesTags, path=orthologypath, orthology=option.orthology, commonNames=option.commonNames, familyFiles=option.familyFiles, verbose="OFF")
		borthologs = metrn.orthologFinder(bspecies, speciesTags, path=orthologypath, orthology=option.orthology, commonNames=option.commonNames, familyFiles=option.familyFiles, verbose="OFF")
		
		# define input files:
		ainfile = str(path_dict[option.source] + "/" + option.a).replace("//","/")
		binfile = str(path_dict[option.source] + "/" + option.b).replace("//","/")
		
		# find target matrix indexes (keys):
		aindexes, bindexes = option.indexes.split(",")
		ai, aj = aindexes.split(":")
		bi, bj = aindexes.split(":")
		
		# find target matrix values:
		ax, bx = option.values.split(",")
		
		# load organism matrixes:
		amatrix = general.matrixBuilder(i=ai, j=aj, x=ax, infile=ainfile, datatype="float")
		bmatrix = general.matrixBuilder(i=bi, j=bj, x=bx, infile=binfile, datatype="float")
		
		# load organism matrixes:
		aicontext = general.matrixBuilder(i=ai, j=aj, x="i.context", infile=ainfile)
		ajcontext = general.matrixBuilder(i=ai, j=aj, x="j.context", infile=ainfile)
		bicontext = general.matrixBuilder(i=bi, j=bj, x="i.context", infile=binfile)
		bjcontext = general.matrixBuilder(i=bi, j=bj, x="j.context", infile=binfile)
		
		# build expanded matrixes:
		aexpand, bexpand = dict(), dict()
		
		# generate comparison matrix:
		ak, acombined, acomplete = 0, dict(), dict()
		for aifactor in ortholog_dict[aspecies]:
			for ajfactor in ortholog_dict[aspecies]:
				for ai in amatrix:
					for aj in amatrix[ai]:
						if aifactor in ai and ajfactor in aj:
							#print ai, aj, aifactor, ajfactor, aicontext[ai][aj], ajcontext[ai][aj]
							ak += 1
							
							if not aifactor in acombined:
								acombined[aifactor] = dict()
							if not ajfactor in acombined[aifactor]:
								acombined[aifactor][ajfactor] = list()
							acombined[aifactor][ajfactor].append(amatrix[ai][aj])
							
							if not aifactor in acomplete:
								acomplete[aifactor] = dict()
							if not ai in acomplete[aifactor]:
								acomplete[aifactor][ai] = dict()
							if not ajfactor in acomplete[aifactor][ai]:
								acomplete[aifactor][ai][ajfactor] = dict()
							acomplete[aifactor][ai][ajfactor][aj] = amatrix[ai][aj]
		
		# generate comparison matrix:
		bk, bcombined, bcomplete = 0, dict(), dict()
		for aifactor in ortholog_dict[aspecies]:
			for ajfactor in ortholog_dict[aspecies]:
				for bifactor in ortholog_dict[aspecies][aifactor][bspecies]:
					for bjfactor in ortholog_dict[aspecies][ajfactor][bspecies]:
						
						for bi in bmatrix:
							for bj in bmatrix[bi]:
								if bifactor in bi and bjfactor in bj:
									#print bi, bj, bifactor, bjfactor, bicontext[bi][bj], bjcontext[bi][bj]
									bk += 1
							
									if not aifactor in bcombined:
										bcombined[aifactor] = dict()
									if not ajfactor in bcombined[aifactor]:
										bcombined[aifactor][ajfactor] = list()
									bcombined[aifactor][ajfactor].append(bmatrix[bi][bj])
									
									if not aifactor in bcomplete:
										bcomplete[aifactor] = dict()
									if not bi in bcomplete[aifactor]:
										bcomplete[aifactor][bi] = dict()
									if not ajfactor in bcomplete[aifactor][bi]:
										bcomplete[aifactor][bi][ajfactor] = dict()
									bcomplete[aifactor][bi][ajfactor][bj] = bmatrix[bi][bj]
		
								
				
		# make comparison output folders:
		comparisonpath = path_dict[option.source] + "comparison/" + aspecies + "/" + bspecies + "/"
		general.pathGenerator(comparisonpath)
		
		# make compined output file:
		x = 0
		processed = list()
		f_outfile = comparisonpath + "maphybrid_" + option.source + "_" + option.name + "_combined.txt"
		f_output = open(f_outfile, "w")
		print >>f_output, "\t".join(["i", "j", "a.species", "b.species", "label", "a.mean", "b.mean", "a.max", "b.max", "a.std", "b.std"])
		for aifactor in acombined:
			for ajfactor in acombined:
				if aifactor in bcombined and ajfactor in bcombined:
					label = ":".join(sorted([aifactor, ajfactor]))
					print >>f_output, "\t".join(map(str, [aifactor, ajfactor, aspecies, bspecies, label, numpy.mean(acombined[aifactor][ajfactor]), numpy.mean(bcombined[aifactor][ajfactor]), max(acombined[aifactor][ajfactor]), max(bcombined[aifactor][ajfactor]), numpy.std(acombined[aifactor][ajfactor]), numpy.std(bcombined[aifactor][ajfactor])]))
					x += 1
		f_output.close()
		
		# make complete output file:
		y = 0
		f_outfile = comparisonpath + "maphybrid_" + option.source + "_" + option.name + "_complete.txt"
		f_output = open(f_outfile, "w")
		print >>f_output, "\t".join(["i.ortholog", "j.ortholog", "a.species", "b.species", "a.i", "a.j", "b.i", "b.j", "a.value", "b.value", "a.comparison", "b.comparison", "i", "j"])
		for aifactor in acomplete:
			for ai in acomplete[aifactor]:
				for ajfactor in acomplete[aifactor][ai]:
					for aj in acomplete[aifactor][ai][ajfactor]:
						
						process = False
						if aifactor in bcomplete:
							for bi in bcomplete[aifactor]:
								if ajfactor in bcomplete[aifactor][bi]:
									for bj in bcomplete[aifactor][bi][ajfactor]:
										process = True
										
						if process:
							print >>f_output, "\t".join(map(str, [aifactor, ajfactor, aspecies, bspecies, ai, aj, bi, bj, acomplete[aifactor][ai][ajfactor][aj], bcomplete[aifactor][bi][ajfactor][bj], ":".join([ai, aj]), ":".join([bi, bj]), ":".join([aifactor, ai, aj]), ":".join([ajfactor, bi, bj])]))
							y += 1
		f_output.close()
		
		print ak, bk
		print x, y
		print
		
		
		# load GO lines from input files:
		#orthologs = list()
		#asublines, bsublines = list(), list()
		#ahd = general.build_header_dict(option.a)
		#bhd = general.build_header_dict(option.b)
		#adict = loader(option.a, ahd)
		#bdict = loader(option.b, bhd)
		
		# generate a-file and b-file headers, as well as output header:
		#aheader, bheader = list(), list()
		#for header in general.valuesort(ahd):
		#	aheader.append(header + ".a")
		#for header in general.valuesort(bhd):
		#	bheader.append(header + ".b")
		#outheader = ["i", "j", "items.a", "items.b", "overlap.a", "overlap.b", "overlap.avg", "overlap.sum", "overlap.max", "overlap.count", "items.count"] #"a.only.goids", "b.only.goids", "overlap.goids"]
		#print >>f_output, "\t".join(outheader)
		
		
		"""
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
		
		print
		"""
		
	# merge ortholog values mode:
	if option.mode == "merge.overlap":
			
		# import orthologs dictionary:
		#orthologs_dict = buildOrthologs(inpath + "configure_orthologs_" + option.analysis + ".txt")
		
		# find species-comparison orthologs:
		speciesTags = option.species.split(",")
		aspecies, bspecies = speciesTags
		
		# generate output peaks name:
		orthologTag = option.nametag + metrn.orthologLabel(aspecies, speciesTags)
		
		# define orthology path:
		if option.orthology == "direct":
			orthologypath = orthologspath + "orthologs/"
		elif option.orthology == "family":
			orthologypath = orthologspath + "families/"
		elif option.orthology == "groups":
			orthologypath = orthologspath + "groups/"
		
		# generate orthology dictionary:
		ortholog_dict = metrn.orthologBuilder(speciesTags, path=orthologypath, orthology=option.orthology, commonNames=option.commonNames, familyFiles=option.familyFiles, verbose="OFF")
		
		print
		print "Evaluating orthologs:"
		for targetFactor in ortholog_dict[option.organism]:
			for specieTag in ortholog_dict[option.organism][targetFactor]:
				print option.organism, targetFactor, specieTag, ":", ",".join(ortholog_dict[option.organism][targetFactor][specieTag])
		print
		#pdb.set_trace()
		
		# target specie orthologs:
		aorthologs = metrn.orthologFinder(aspecies, speciesTags, path=orthologypath, orthology=option.orthology, commonNames=option.commonNames, familyFiles=option.familyFiles, verbose="OFF")
		borthologs = metrn.orthologFinder(bspecies, speciesTags, path=orthologypath, orthology=option.orthology, commonNames=option.commonNames, familyFiles=option.familyFiles, verbose="OFF")
		
		# define input files:
		ainfile = str(path_dict[option.source] + "/" + option.a).replace("//","/")
		binfile = str(path_dict[option.source] + "/" + option.b).replace("//","/")
		
		# find target matrix indexes:
		#aindexes, bindexes = option.indexes.split(",")
		#ai, aj = aindexes.split(":")
		#bi, bj = aindexes.split(":")
		
		# find target matrix values:
		ax, bx = option.values.split(",")
		
		# load header dictionaries:
		aHeader = general.build_header_dict(ainfile)
		bHeader = general.build_header_dict(binfile)
		
		# capture universe of values:
		universe = list()
		
		# load input dictionaries:
		adict, acomplex = dict(), dict()
		for inline in open(ainfile).readlines()[1:]:
			initems = inline.strip().split("\t")
			invalue = initems[aHeader[ax]]
			inlabel = metrn.labelExtractor(initems, target="dataset", mode=option.label, headerDict=aHeader)
			if int(initems[aHeader["genome.count"]]) < int(initems[aHeader["genome.total"]])/option.fraction:
				universe.append(invalue)
				if not inlabel in adict:
					adict[inlabel] = list()
				adict[inlabel].append(invalue)
			
			if option.source == "go":
				if not inlabel in acomplex:
					acomplex[inlabel] = dict()
				acomplex[inlabel][invalue] = [initems[aHeader["dataset.count"]], initems[aHeader["genome.count"]], initems[aHeader["adjusted.pvalue"]]]
				
		bdict, bcomplex = dict(), dict()
		for inline in open(binfile).readlines()[1:]:
			initems = inline.strip().split("\t")
			invalue = initems[bHeader[bx]]
			inlabel = metrn.labelExtractor(initems, target="dataset", mode=option.label, headerDict=bHeader)
			if int(initems[bHeader["genome.count"]]) < int(initems[bHeader["genome.total"]])/option.fraction:
				universe.append(invalue)
				if not inlabel in bdict:
					bdict[inlabel] = list()
				bdict[inlabel].append(invalue)
		
			if option.source == "go":
				if not inlabel in bcomplex:
					bcomplex[inlabel] = dict()
				bcomplex[inlabel][invalue] = [initems[bHeader["dataset.count"]], initems[bHeader["genome.count"]], initems[bHeader["adjusted.pvalue"]]]
		
		# reduce universe of values to set:
		universe = set(universe)
		
		# make output folders:
		comparisonpath = path_dict[option.source] + "comparison/" + aspecies + "/" + bspecies + "/"
		general.pathGenerator(comparisonpath)
		
		# setup output file:
		f_outfile = comparisonpath + "maphybrid_" + option.source + "_" + option.name + "_combined.txt"
		f_output = open(f_outfile, "w")
		print >>f_output, "\t".join(["i", "j", "match", "i.values", "j.values", "overlap", "total", "i.fraction", "j.fraction", "overlap.avg", "overlap.sum", "overlap.max", "pvalue", "adjusted.pvalue", "overlap.values"])
		
		# count number of tests:
		adjust = 0
		for alabel in adict:
			for blabel in bdict:
				adjust += 1
				
		# generate matrix:
		matrix = dict()
		for alabel in adict:
			for blabel in bdict:
				
				# extract label info:
				aorganism, astrain, afactor, acontext, ainstitute, amethod = alabel.split("_")[:6]
				borganism, bstrain, bfactor, bcontext, binstitute, bmethod = blabel.split("_")[:6]
				
				# determine orthology:
				if bfactor in ortholog_dict[aorganism][afactor][borganism]:
					match = "+"
				else:
					match = ""
				
				# regenerate labels:
				i = metrn.labelGenerator(target=option.target, mode="label", dataset=alabel)
				j = metrn.labelGenerator(target=option.target, mode="label", dataset=blabel)
				
				if not alabel in matrix:
					matrix[alabel] = dict()
				if not blabel in matrix[alabel]:
					matrix[alabel][blabel] = dict()
				
				avalues = set(adict[alabel]).intersection(universe)
				bvalues = set(bdict[blabel]).intersection(universe)
				aonly = set(avalues).difference(set(bvalues))
				bonly = set(bvalues).difference(set(avalues))
				overlap = set(avalues).intersection(set(bvalues))
				total = set(avalues).union(set(bvalues))
				
				#if afactor in ["MXI1", "MDL-1"] and bfactor in ["MXI1", "MDL-1"]:
				#	print len(avalues)
				#	print len(bvalues)
				#	print len(overlap)
				#	print overlap
				#	pdb.set_trace()
				
				if len(overlap) == 0:
					aoverlap, boverlap, overlap_avg, overlap_max, overlap_sum = 0, 0, 0, 0, 0
					pvalue, adjPvalue = 1, 1
				else:
					aoverlap = float(len(overlap))/len(avalues)
					boverlap = float(len(overlap))/len(bvalues)
					overlap_avg = numpy.mean([aoverlap, boverlap])
					overlap_max = max([aoverlap, boverlap])
					overlap_sum = float(len(overlap))/len(total)
				
					# Hypergeometric paramters:
					m = len(avalues) # number of white balls in urn
					n = len(universe) - len(avalues) # number of black balls in urn
					N = len(bvalues) # number of balls drawn from urn
					x = len(overlap) # number of white balls in drawn
				
					# If I pull out all balls with elephant tatoos (N), is the draw enriched in white balls?:
					pvalue = hyper.fishers(x, m+n, m, N, method="right")
					adjPvalue = hyper.limit(pvalue*adjust)
					
				i = i.replace("-S3", "S3").replace("-hESC", "hesc")
				j = j.replace("-S3", "S3").replace("-hESC", "hesc")
				output = [i, j, match, len(avalues), len(bvalues), len(overlap), len(universe), aoverlap, boverlap, overlap_avg, overlap_sum, overlap_max, pvalue, adjPvalue, ",".join(sorted(list(overlap)))]
				matrix[alabel][blabel] = output
				print >>f_output, "\t".join(map(str, output))
				
		# close output file:
		f_output.close()
		
	
	# merge ortholog binding frequencies mode:
	if option.mode == "merge.binding":
			
		# find species-comparison orthologs:
		speciesTags = option.species.split(",")
		aspecies, bspecies = speciesTags
		
		# generate output peaks name:
		orthologTag = option.nametag + metrn.orthologLabel(aspecies, speciesTags)
		
		# define orthology path:
		if option.orthology == "direct":
			orthologypath = orthologspath + "orthologs/"
		elif option.orthology == "family":
			orthologypath = orthologspath + "families/"
		elif option.orthology == "groups":
			orthologypath = orthologspath + "groups/"
		
		# generate orthology dictionary:
		ortholog_dict = metrn.orthologBuilder(speciesTags, path=orthologypath, orthology=option.orthology, commonNames=option.commonNames, familyFiles=option.familyFiles, verbose="OFF")
		
		#for targetFactor in ortholog_dict[option.organism]:
		#	for specieTag in ortholog_dict[option.organism][targetFactor]:
		#		print option.organism, targetFactor, specieTag, ":", ",".join(ortholog_dict[option.organism][targetFactor][specieTag])
		
		# target specie orthologs:
		aorthologs = metrn.orthologFinder(aspecies, speciesTags, path=orthologypath, orthology=option.orthology, commonNames=option.commonNames, familyFiles=option.familyFiles, verbose="OFF")
		borthologs = metrn.orthologFinder(bspecies, speciesTags, path=orthologypath, orthology=option.orthology, commonNames=option.commonNames, familyFiles=option.familyFiles, verbose="OFF")
		
		# define input files:
		ainfile = str(path_dict[option.source] + "/" + option.a).replace("//","/")
		binfile = str(path_dict[option.source] + "/" + option.b).replace("//","/")
		
		# load header dictionaries:
		aHeader = general.build_header_dict(ainfile)
		bHeader = general.build_header_dict(binfile)
		
		# load binding data:
		print
		print "Loading binding data..."
		aDict = general.build2(ainfile, id_column="dataset")
		bDict = general.build2(binfile, id_column="dataset")
		
		# make comparison path	
		comparisonpath = path_dict[option.source] + "comparison/" + aspecies + "/" + bspecies + "/" + option.name + "/"
		general.pathGenerator(comparisonpath)
		
		# explicit lables for the chromatin states or promoter regions:
		if option.indexes == "iHMM":
			inlabels = ["1_Pro", "2_Enh1", "3_Enh2", "4_Egn1", "5_Egn2", "6_Egn3", "7_Egn4", "8_Egn5", "9_Egn6", "10_Rep1", "11_Rep2", "12_Het1", "13_Het2", "14_Low1", "15_Low2", "16_Low3"]
		elif option.indexes == "125kb":
			inlabels = ["0:1000", "1001:2000", "2001:3000", "3001:4000", "4001:5000", "others"]
		elif option.indexes == "125EN":
			inlabels = ["0:500", "501:1000", "1001:2000", "2001:10000", "others", "enhancer"]
		
		# export fractions:
		print "Exporting binding ratios..."
		index = 1
		f_output = open(comparisonpath + "maphybrid_binding_" + option.a.split("/")[0] + "_vs_" + option.b.split("/")[0] + "_summary.txt", "w")
		print >>f_output, "\t".join(["i", "j", "index", "label", "type", "color", "i.value", "j.value", "i.fraction", "j.fraction"])
		for afactor in sorted(list(set(aDict.keys()).intersection(set(aorthologs)))):
			for bfactor in sorted(list(set(ortholog_dict[aspecies][afactor][bspecies]).intersection(set(bDict)))):
				label = ":".join([afactor, bfactor])
				avalues, bvalues, xratios = list(), list(), list()
				color = 1
				for inlabel in inlabels:
					#print aspecies, bspecies, afactor, bfactor, aDict[afactor][inlabel], bDict[bfactor][inlabel]
					avalue = float(aDict[afactor][inlabel])
					bvalue = float(bDict[bfactor][inlabel])
					if (avalue + bvalue) > 0:
						aratio = float(avalue)/(avalue + bvalue)
						bratio = float(bvalue)/(avalue + bvalue)
					else:
						aratio, bratio = 0, 0
					avalues.append(avalue)
					bvalues.append(bvalue)
					xratios.append(aratio)
					output = [afactor, bfactor, index, label, inlabel, color, avalue, bvalue, aratio, bratio]
					print >>f_output, "\t".join(map(str, output))
					color += 1
				#print afactor, bfactor
				#print avalues
				#print bvalues
				#print xratios
				#print
				index += 1
		f_output.close()
		print
			
if __name__ == "__main__":
	main()
	
	print "Completed:", time.asctime(time.localtime())

#python mapHybrid.py --path ~/meTRN --mode merge.matrix --organism hs --species hs,ce --orthology family --source coassociations --A hs_orthoHsCe_com_cx_xot/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --B ce_orthoHsCe_com_cx_xot/promoter_regions/summary/mapcas_report_promoter_regions_p5e-02_matrix.txt --indexes i:j,i:j --values mirror.passing,mirror.passing --name orthoHsCe_com_cx_xot

#python mapHybrid.py --path ~/meTRN --mode merge.overlap --organism hs --species hs,ce --orthology family --source go --A hs_orthoHsCe_com_cx_xot/p5e-1/summary/mapgo_complete_hs_orthoHsCe_com_cx_xot_p5_hc1_hp5e-02_summary --B  ce_orthoHsCe_com_cx_xot/p5e-1/summary/mapgo_complete_ce_orthoHsCe_com_cx_xot_p5_hc1_hp5e-02_summary --values id,id --label rebuild --target 'factor(context)' --name orthoHsCe_com_cx_xot