#!/usr/bin/env python
# generate peak set complete files, binding region files, and report files!

import sys
import time
import optparse
import general
import numpy
import pickle
import pdb
import metrn
import modencode
import os

import itertools
import simplejson as json
import palette

print "Command:", " ".join(sys.argv)
print "Timestamp:", time.asctime(time.localtime())

def main():
	
	parser = optparse.OptionParser()
	parser.add_option("--path", action = "store", type = "string", dest = "path", help = "path from script to files")
	parser.add_option("--mode", action = "store", type = "string", dest = "mode", help = "operations to be performed")
	parser.add_option("--infile", action = "store", type = "string", dest = "infile", help = "Input file", default="OFF")
	parser.add_option("--source", action = "store", type = "string", dest = "source", help = "Input path", default="OFF")
	parser.add_option("--peaks", action = "store", type = "string", dest = "peaks", help = "peaks to be used for analysis", default="OFF")
	parser.add_option("--organism", action = "store", type = "string", dest = "organism", help = "Target organism for operations...", default="OFF")
	parser.add_option("--library", action = "store", type = "string", dest = "library", help = "library configuration final report card", default="OFF")
	parser.add_option("--rmdup", action = "store", type = "string", dest = "rmdup", help = "Remove duplicates library", default="OFF")
	parser.add_option("--overwrite", action = "store", type = "string", dest = "overwrite", help = "Overwrite files?", default="OFF")
	parser.add_option("--rename", action = "store", type = "string", dest = "rename", help = "Targets to rename. Comma-separated list of 'target:replacement' pairs to search and replace.", default="OFF")
	parser.add_option("--include", action = "store", type = "string", dest = "include", help = "Targets to include", default="OFF")
	parser.add_option("--exclude", action = "store", type = "string", dest = "exclude", help = "Targets to exclude", default="OFF")
	parser.add_option("--contexts", action = "store", type = "string", dest = "contexts", help = "Contexts to compare", default="OFF")
	parser.add_option("--species", action = "store", type = "string", dest = "species", help = "Species to compare", default="OFF")
	parser.add_option("--orthology", action = "store", type = "string", dest = "orthology", help = "Use 'direct' or 'family' orthologs?", default="direct")
	parser.add_option("--nametag", action = "store", type = "string", dest = "nametag", help = "Orthology nametag: nametagHsCe", default="ortho")
	parser.add_option("--commonNames", action = "store", type = "string", dest = "commonNames", help = "Grab common names file?", default="ON")
	parser.add_option("--familyFiles", action = "store", type = "string", dest = "familyFiles", help = "Grab cleaned files?", default="formatted")
	parser.add_option("--server", action = "store", type = "string", dest = "server", help = "Are we on the server?", default="OFF")
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
	peakspath = path_dict["peaks"]
	reportspath = path_dict["reports"]
	gopath = path_dict["go"]
	hotpath = path_dict["hot"]
	qsubpath = path_dict["qsub"]
	cellspath = path_dict["cells"]
	
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
	
	# define organism parameters:
	if option.organism == "hs" or option.organism == "h.sapiens":
		organismTag = "hs"
		#organismIGV = "ce6"
	elif option.organism == "mm" or option.organism == "m.musculus":
		organismTag = "mm"
		#organismIGV = "ce6"
	elif option.organism == "ce" or option.organism == "c.elegans":
		organismTag = "ce"
		#organismIGV = "ce6"
	elif option.organism == "dm" or option.organism == "d.melanogaster":
		organismTag = "dm"
		#organismIGV = "dm5"
	
	# define organisms:
	organismTags = ["hs","mm","ce","dm"]
	
	# sunburst visualization mode:
	if option.mode == "sunburst":
	
		# define reports paths:
		reportspath = reportspath + organismTag + "/"
		general.pathGenerator(reportspath)
	
		# load report data:
		headerDict = { "dataset":0, "factor":1, "context.main":2, "context.sub":3, "factor.class":4, "peaks":5 }
		reportDict = general.build2(path_dict[option.source] + option.infile, id_column="dataset", header_dict=headerDict, header=False)
		outfile = option.infile.replace(".tab", "")
		
		#k = 0
		#for dataset in reportDict:
		#	if reportDict[dataset]["factor.class"] == "POL":
		#		k += 1
		#print "POL:", k
		
		# generate sunburst matrix:
		matrix, sizing = dict(), dict()
		for dataset in sorted(reportDict.keys()):
			factor, context, subtext, factorClass, size = reportDict[dataset]["factor"], reportDict[dataset]["context.main"], reportDict[dataset]["context.sub"], reportDict[dataset]["factor.class"], reportDict[dataset]["peaks"]
			if option.peaks == "OFF" or dataset in os.listdir(peakspath + option.peaks):
				if not context in matrix:
					matrix[context] = dict()
					sizing[context] = dict()
				if not factor in matrix[context]:
					matrix[context][factor] = dict()
				matrix[context][factor][dataset] = numpy.log10(float(size))
				sizing[context][factor] = sum(matrix[context][factor].values())
		
		#print set(os.listdir(peakspath + option.peaks)).difference(set(reportDict.keys()))
		#pdb.set_trace()
		
		# define color function:
		def colorMaker(organismTag, context, colors, exception="#B8BEC2"): 
			if organismTag + "." + context in colors:
				return colors[organismTag + "." + context]
			else:
				return exception
		
		# define color spaces:
		colors = {
			"ce.EE" : "#F3C520",
			"ce.LE" : "#F3D020",
			"ce.EM" : "#F3DB20",
			"ce.EX" : "#F3DB20",
			"ce.L1" : "#20F263",
			"ce.L2" : "#20F2AC",
			"ce.L3" : "#20E4F2",
			"ce.L4" : "#209EF2",
			"dm.EE" : "#F3C520",
			"dm.LE" : "#F3D020",
			"dm.PP" : "#209EF2",
			"hs.H1-hESC" : "#F3DB20",
			"hs.GM12878" : "#20F263",
			"hs.HepG2" : "#20F2AC",
			"hs.HeLa-S3" : "#20E4F2",
			"hs.K562" : "#209EF2",
			"hs.H1" : "#F3DB20",
			"hs.GM" : "#20F263",
			"hs.HG" : "#20F2AC",
			"hs.HL" : "#20E4F2",
			"hs.K5" : "#209EF2"
			}
		
		# define context order:
		order = {
			#"ce": ["EE", "LE", "EM", "L1","L2", "L3", "L4", "LY", "YA", "D4", "S1"],
			"ce": ["EX", "L1","L2", "L3", "L4"],
			"dm": ["EE", "LE", "PP"],
			"hs": ["H1", "GM", "HG", "HL","K5"] 
			#"hs": ["H1-hESC", "GM12878", "HepG2", "HeLa-S3","K562"] 
			} 
		
		# order contexts:
		contextOrder = order[organismTag]
		for context in sorted(matrix.keys()):
			if not context in contextOrder:
				contextOrder.append(context)
			
		# generate sunburst matrix:
		sunburst = { "name":"sunburst", "color":"#B8BEC2", "children":list() }
		contexts = list()
		finalDatasets, finalFactors, finalPOL, finalTF, finalCH = list(), list(), list(), list(), list()
		for context in contextOrder:
			factors = list()
			colorBase = palette.Color(colorMaker(organismTag, context, colors))
			colorContext = colorBase
			colorLimit = float(0.30)/len(matrix[context].keys())
			if colorLimit > 0.05:
				colorLimit = 0.05
			colorEdit = colorLimit
			factorOrder = general.valuesort(sizing[context])
			factorOrder.reverse()
			for factor in factorOrder:
				finalFactors.append(factor)
				colorFactor = colorBase.darker(amt=0.08).lighter(amt=colorEdit)
				colorEdit += colorLimit
				datasets = list()
				for dataset in general.valuesort(matrix[context][factor]):
					finalDatasets.append(dataset)
					factorClass = reportDict[dataset]["factor.class"]
					if factorClass == "POL":
						finalPOL.append(factor)
						colorDataset = "#FF5757"
					elif factorClass == "TF":
						finalTF.append(factor)
						colorDataset = colorContext.lighter(amt=0.30).hex
						#colorDataset = colorContext.darker(amt=0.05).hex
						if context == "OT":
							colorDataset = "#D6D6D6"
					elif factorClass == "CH":
						finalCH.append(factor)
						colorDataset = colorContext.darker(amt=0.05).hex
					else:
						colorDataset = "#FFFFFF"
					#colorDataset = colorContext.darker(amt=0.05).hex
					#colorDataset = colorContext.lighter(amt=0.30).hex
					datasets.append({"name":dataset, "color":colorDataset, "size":matrix[context][factor][dataset]})
				factors.append({"name":factor, "color":colorFactor.hex, "children":datasets})
			contexts.append({"name":context, "color":colorContext.hex, "children":factors})
		sunburst["children"] = contexts
		
		print
		print "Inputs:", len(reportDict)
		print "Datasets:", len(finalDatasets)
		print "Factors:", len(set(finalFactors))
		print
		print "Class:Polymerase:", len(finalPOL)
		print "Class:Transcription Factor:", len(finalTF)
		print "Class:Chromatin-associated:", len(finalCH)
		print
		print "Factors per context:"
		for context in contextOrder:
			print context + ":", len(matrix[context].keys())
		print
		
		# export sunburst json file (to report path):
		f_output = open(reportspath + "mapdata_" + outfile + ".json", "w")
		json.dump(sunburst, f_output)
		f_output.close()
		
		# export sunburst json file (to web directory):
		reportspath = "/Library/WebServer/Documents/metrn/data/reports/"
		f_output = open(reportspath + "mapdata_" + outfile + ".json", "w")
		json.dump(sunburst, f_output)
		f_output.close()
		
	
	# paralogs file copy mode:
	if option.mode == "paralogs":
	
		# indicate target species:
		organismKeys = {
			"Cel-" : "ce",
			"dmel_" : "dm",
			"human_" : "hs"
			}
	
		# indicate gene names:
		#id2name_dict, name2id_dict = modencode.idBuild()
		print
		print "Loading gene IDs and common names..."
		organismDict = dict()
		for organismTag in ["ce", "hs"]:
			if organismTag == "ce":
				nameUpper, idUpper = True, False
			if organismTag == "dm":
				nameUpper, idUpper = False, False
			if organismTag == "hs":
				nameUpper, idUpper = False, False
			organismDict[organismTag] = modencode.idBuild(inpath + metrn.reference[organismTag]["gene_ids"], metrn.reference[organismTag]["gene.link"], metrn.reference[organismTag]["symbol.link"], mode="label", header=True, nameUpper=nameUpper, idUpper=idUpper)
					
		# parse family file:
		paralogs = dict()
		inlines = open(path_dict[option.source] + option.infile).readlines()
		print "Parsing family-level orthologs..."
		for inline in inlines:
			genes = inline.strip().split("\t")
			familyID = int(genes.pop(0))
			paralogs[familyID] = dict()
			#print familyID, len(genes)
			for gene in genes:
				for key in organismKeys:
					if gene[:len(key)] == key:
						species, gene = organismKeys[key], gene[len(key):]
						if species == "ce" and not "WBGene" in gene:
							gene = gene.upper()
						if species == "dm" and "-P" in gene:
							gene = gene.split("-P")[0]
						if not species in paralogs[familyID]:
							paralogs[familyID][species] = list()
						paralogs[familyID][species].append(gene)
						#print species, gene
			
		# export paralogs per species:
		print "Exporting paralogs per family and species..."
		general.pathGenerator(orthologspath + "paralogs/")
		f_output = open(orthologspath + "paralogs/mapdata_paralogs_xx_family.txt", "w")
		c_output = open(orthologspath + "paralogs/mapdata_paralogs_ce_matrix.txt", "w")
		d_output = open(orthologspath + "paralogs/mapdata_paralogs_dm_matrix.txt", "w")
		h_output = open(orthologspath + "paralogs/mapdata_paralogs_hs_matrix.txt", "w")
		print >>f_output, "\t".join(["family.id", "species", "paralogs"])
		print >>c_output, "\t".join(["i", "j"])
		print >>d_output, "\t".join(["i", "j"])
		print >>h_output, "\t".join(["i", "j"])
		for familyID in sorted(paralogs.keys()):
			for species in sorted(paralogs[familyID].keys()):
				print >>f_output, "\t".join(map(str, [familyID, species, ",".join(sorted(paralogs[familyID][species]))]))
				for permutation in itertools.permutations(sorted(paralogs[familyID][species]), 2):
					if species == "ce":
						print >>c_output, "\t".join(permutation)
					if species == "dm":
						print >>d_output, "\t".join(permutation)
					if species == "hs":
						print >>h_output, "\t".join(permutation)
		f_output.close()
		c_output.close()
		d_output.close()
		h_output.close()
		print
			
	
	# scan numbers of peaks, binding regions, HOT regions, XOT regions, HOT cutoffs, and XOT cutoffs:
	elif option.mode == "scanner":
	
		# grab stage and parts:
		organism, selection, factors, stage, subset = option.peaks.split("_")
		if stage == "cx":
			hotFlag = "any.bed"
			ubiFlag = "all.bed"
		else:
			hotFlag = "hot.bed"
			ubiFlag = False
		# define input files
		completefile = peakspath + "mappeaks_" + option.peaks + "_complete.bed"
		compiledfile = peakspath + "mappeaks_" + option.peaks + "_compiled.bed"
		hotregionfile = hotpath + "regions/maphot_" + "_".join([organism, selection, "reg", stage, "occP05", hotFlag])
		xotregionfile = hotpath + "regions/maphot_" + "_".join([organism, selection, "reg", stage, "occP01", hotFlag])
		if ubiFlag:
			hotubiquitousfile = hotpath + "regions/maphot_" + "_".join([organism, selection, "reg", stage, "occP05", ubiFlag])
			xotubiquitousfile = hotpath + "regions/maphot_" + "_".join([organism, selection, "reg", stage, "occP01", ubiFlag])
		
		# count binding sites and binding regions:
		sites = len(open(completefile).readlines())
		regions = len(open(compiledfile).readlines())
		
		# get HOT region and XOT region cutoffs:
		hotCutoff, xotCutoff = list(), list()
		for inline in open(hotregionfile).readlines():
			chrm, start, stop, feature, score, strand = inline.strip().split("\t")[:6]
			hotCutoff.append(int(score))
		for inline in open(xotregionfile).readlines():
			chrm, start, stop, feature, score, strand = inline.strip().split("\t")[:6]
			xotCutoff.append(int(score))
		hotRegions, xotRegions = len(hotCutoff), len(xotCutoff)
		hotMean, xotMean = round(numpy.mean(hotCutoff), 1), round(numpy.mean(xotCutoff), 1)
		hotCutoff, xotCutoff = min(hotCutoff), min(xotCutoff)
		print
		print "Stage Region Report:", stage
		print sites, regions, str(hotRegions) + " (" + str(round(100*float(hotRegions)/regions, 1)) + "%)", str(xotRegions) + " (" + str(round(100*float(xotRegions)/regions, 1)) + "%)", hotCutoff, xotCutoff, hotMean, xotMean
		print
		
		# get ubiquitously HOT region and XOT regions:
		if ubiFlag:
			hotCutoff, xotCutoff = list(), list()
			for inline in open(hotubiquitousfile).readlines():
				chrm, start, stop, feature, score, strand = inline.strip().split("\t")[:6]
				hotCutoff.append(int(score))
			for inline in open(xotubiquitousfile).readlines():
				chrm, start, stop, feature, score, strand = inline.strip().split("\t")[:6]
				xotCutoff.append(int(score))
			hotRegions, xotRegions = len(hotCutoff), len(xotCutoff)
			print "Ubiquitously HOT regions:", str(hotRegions) + " (" + str(round(100*float(hotRegions)/regions, 1)) + "%)"
			print "Ubiquitously XOT regions:", str(xotRegions) + " (" + str(round(100*float(xotRegions)/regions, 1)) + "%)"
			print
		
	
	# file survey mode, counts numbers of files and reads:
	elif option.mode == "survey":
	
		print
		print "Loading configuration dictionary..."
		configDict = general.build2(path_dict[option.source] + option.infile, id_column="url")
		process = False
		
		if "mapdata_revise_" + option.peaks + "_library_report.txt" in os.listdir(reportspath + option.organism):
			datasetReport = reportspath + option.organism + "/" + "mapdata_revise_" + option.peaks + "_dataset_report.txt"
			libraryReport = reportspath + option.organism + "/" + "mapdata_revise_" + option.peaks + "_library_report.txt"
			process = True
		
		if process:
			print "Loading revised configurations..."
			datasetDict = general.build2(datasetReport, id_column="dataset")
			libraryDict = general.build2(libraryReport, id_column="filename")
		
		print "Loading peak set files of interest..."
		peakfiles = os.listdir(peakspath + option.peaks)
		
		print "Collecting configuration keys for peaks..."
		inputKeys = list()
		for peakfile in peakfiles:
			for inputKey in sorted(configDict.keys()):
				
				if option.rename != "OFF":
					for scheme in option.rename.split(","):
						target, replacement = scheme.split(":")
						searchKey = inputKey.replace(target, replacement)
				else:
					searchKey = str(inputKey)
					
				if peakfile.replace("_peaks.bed", "") in option.organism + "_" + searchKey:
					inputKeys.append(inputKey)
		
		print
		print "Input files (input FASTQs):", len(configDict)
		print "Peaks files (ChIP-seq experiments):", len(peakfiles)
		print "Match files (ChIP-seq FASTQs):", len(inputKeys)
		
		# collect input reads: 
		sumReads = 0
		for inputKey in inputKeys:
			sumReads += int(configDict[inputKey]["input.reads"])
		
		# collect aligned reads:
		if process:
			sumSignal = 0
			for library in libraryDict:
				sumSignal += int(libraryDict[library]["distinct.reads"])
		
		print "Match reads (Millions)", round(float(sumReads)/(1000000), 2)
		if process:
			print "Signal reads (Millions)", round(float(sumSignal)/(1000000), 2)	
		print
	
	
	# review file mode:
	elif option.mode == "review":
	
		print
		print "Loading configuration dictionary..."
		configDict = general.build2(path_dict[option.source] + option.infile, id_column="url")
		
		print "Loading peak set files of interest..."
		peakfiles = os.listdir(peakspath + option.peaks)
		
		print "Examining experiments attempted (ChIP-seq)..."
		comboDict = dict()
		for inputFile in configDict:
		
			# obtain factor and context for experiment:
			dataset = str(inputFile)
			if option.rename != "OFF":
				for scheme in option.rename.split(","):
					target, replacement = scheme.split(":")
					dataset = dataset.replace(target, replacement)
			
			# reload the corresponding dataset information:
			organism, strain, factor, context, institute, method = metrn.labelComponents(option.organism + "_" + "_".join(dataset.split("_")[:5]))
			configDict[inputFile]["strain"] = strain
			configDict[inputFile]["factor"] = factor
			configDict[inputFile]["stage"] = context
			
			# record factor and context combination:
			if not factor in comboDict:
				comboDict[factor] = dict()
			if not context in comboDict[factor]:
				comboDict[factor][context] = list()
			comboDict[factor][context].append(inputFile)
				
		print "Scanning quality-filtered experiments..."
		successDict, failureDict = dict(), dict()
		successCount, failureCount = 0, 0
		for factor in comboDict:
			for context in comboDict[factor]:
				success, searchKey = False, factor + "_" + context
				for peakfile in peakfiles:
					if searchKey in peakfile:
						success = True
				if success:
					if not factor in successDict:
						successDict[factor] = list()
					successDict[factor].append(context)
					successCount += 1
				else:
					if not factor in failureDict:
						failureDict[factor] = list()
					failureDict[factor].append(context)
					failureCount += 1
		print
		print "Input factor-stage combinations (ChIP-seq experiments):", successCount + failureCount
		print "Successful factor-stage combinations:", successCount, "(" + str(round(100*float(successCount)/(successCount+failureCount), 2)) + "%)"
		print "Failed factor-stage combinations:", failureCount, "(" + str(round(100*float(failureCount)/(successCount+failureCount), 2)) + "%)"
		print

		f_outfile = reportspath + option.organism + "/" + "mapdata_review_" + option.peaks + "_failed_experiments.txt"
		f_output = open(f_outfile, "w")
		print >>f_output, "\t".join(["factor", "context"])
		for factor in sorted(failureDict.keys()):
			for context in sorted(failureDict[factor]):
				print >>f_output, "\t".join([factor, context])
		f_output.close()
		
		f_outfile = reportspath + option.organism + "/" + "mapdata_review_" + option.peaks + "_success_experiments.txt"
		f_output = open(f_outfile, "w")
		print >>f_output, "\t".join(["factor", "context"])
		for factor in sorted(successDict.keys()):
			for context in sorted(successDict[factor]):
				print >>f_output, "\t".join([factor, context])
		f_output.close()
	
	
	# revise quality mode:
	elif option.mode == "revise":
	
		# parse input report dataset and library files:
		datasetReport, libraryReport = option.infile.replace("summary", "report-pool"), option.infile.replace("summary", "report-rmdup")
	
		print
		print "Loading configuration data..."
		configsDict = general.build2(path_dict[option.source] + option.infile, id_column="url")
		datasetDict = general.build2(path_dict[option.source] + datasetReport, id_column="dataset")
		libraryDict = general.build2(path_dict[option.source] + libraryReport, id_column="filename")
		
		print "Loading target peak set files..."
		peakfiles = os.listdir(peakspath + option.peaks)
		
		print "Detecting duplicate library files..."
		combinations, exclusionLibs, exclusionSets = list(), list(), list()
		for library in sorted(libraryDict.keys()):
			for compare in sorted(libraryDict.keys()):
				if library != compare and libraryDict[library]["distinct.reads"] == libraryDict[compare]["distinct.reads"]:
				
					# apply renaming schemes:
					if option.rename != "OFF":
						for scheme in option.rename.split(","):
							target, replacement = scheme.split(":")
							library = library.replace(target, replacement)
							compare = compare.replace(target, replacement)
					
					# examine combination:
					combination = ",".join(sorted([library, compare]))
					if not combination in combinations:
						combinations.append(combination)
						
						# find duplicated datasets in peaks:
						libraryHit, compareHit = False, False
						librarySet, compareSet = "_".join(library.split("_")[:5]), "_".join(compare.split("_")[:5])
						for peakfile in peakfiles:
							if librarySet in peakfile:
								libraryHit = True
							if compareSet in peakfile:
								compareHit = True
						
						# remove duplicated datasets in peaks:
						if compareHit and libraryHit:
							print "Error: Duplicate libraries found!", librarySet, compareSet
							pdb.set_trace()
						
						# determine which dataset/library to exclude:
						elif libraryHit:
							exclusionLibs.append(compare)
							exclusionSets.append(compareSet)
						elif compareHit:
							exclusionLibs.append(library)
							exclusionSets.append(librarySet)
							
		# store the exclusion libraries and datasets:
		exclusionLibs, exclusionSets = list(set(exclusionLibs)), list(set(exclusionSets))
		
		# load dataset and library columns:
		datasetColumns = open(path_dict[option.source] + datasetReport).readline().strip().split("\t")
		libraryColumns = open(path_dict[option.source] + libraryReport).readline().strip().split("\t")
		
		# define output files:
		d_output = open(reportspath + option.organism + "/" + "mapdata_revise_" + option.peaks + "_dataset_report.txt", "w")
		l_output = open(reportspath + option.organism + "/" + "mapdata_revise_" + option.peaks + "_library_report.txt", "w")
		p_output = open(reportspath + option.organism + "/" + "mapdata_revise_" + option.peaks + "_replica_paired.txt", "w")
		r_output = open(reportspath + option.organism + "/" + "mapdata_revise_" + option.peaks + "_replica_report.txt", "w")
		#s_output = open(reportspath + option.organism + "/" + "mapdata_revise_" + option.peaks + "_replica_sample.txt", "w")
		print >>d_output, "\t".join(datasetColumns)
		print >>l_output, "\t".join(libraryColumns)
		print >>p_output, "\t".join(["factor", "context", "datasets", "max.overlap", "method", "strain", "i", "j"])
		print >>r_output, "\t".join(["factor", "context", "datasets", "max.overlap", "method", "strain", "i", "j"])
		#print >>s_output, "\t".join(["i", "j", "overlap"])
		
		print "Exporting pruned dataset report..."
		for dataset in sorted(datasetDict.keys()):
		
			# apply renaming schemes:
			original = str(dataset)
			if option.rename != "OFF":
				for scheme in option.rename.split(","):
					target, replacement = scheme.split(":")
					dataset = dataset.replace(target, replacement)
			
			# reload the corresponding dataset information:
			organism, strain, factor, context, institute, method = metrn.labelComponents(option.organism + "_" + "_".join(dataset.split("_")[:6]))
			datasetDict[dataset] = datasetDict[original]
			datasetDict[dataset]["strain"] = strain
			datasetDict[dataset]["factor"] = factor
			datasetDict[dataset]["stage"] = context
			datasetDict[dataset]["institute"] = institute
			datasetDict[dataset]["method"] = method
			if dataset != original:
				del datasetDict[original]
			
			# export the dataset information:
			if not dataset in exclusionSets:
				output = list()
				for column in datasetColumns:
					if column == "dataset":
						output.append(dataset)
					else:
						output.append(datasetDict[dataset][column])
				print >>d_output, "\t".join(output)
			
		print "Exporting pruned library report..."
		missing = list()
		for library in sorted(libraryDict.keys()):
			dataset = "_".join(library.split("_")[:5])
			
			# apply renaming schemes:
			original = str(library)
			if option.rename != "OFF":
				for scheme in option.rename.split(","):
					target, replacement = scheme.split(":")
					dataset = dataset.replace(target, replacement)
					library = library.replace(target, replacement)
			
			# reload the corresponding dataset information:
			organism, strain, factor, context, institute, method = metrn.labelComponents(option.organism + "_" + "_".join(dataset.split("_")[:6]))
			libraryDict[library] = libraryDict[original]
			libraryDict[library]["strain"] = strain
			libraryDict[library]["factor"] = factor
			libraryDict[library]["stage"] = context
			libraryDict[library]["institute"] = institute
			libraryDict[library]["method"] = method
			if library != original:
				del libraryDict[original]
			
			if not dataset in exclusionSets:
				if not dataset in datasetDict:
					print "Dataset missing:", dataset, "(" + library + ")"
					missing.append(dataset)
				
				else:
					output = list()
					for column in libraryColumns:
						if column == "filename":
							output.append(library)
						else:
							output.append(libraryDict[library][column])
					print >>l_output, "\t".join(output)
		
		print "Missing datasets:", len(set(missing))
		
		print
		print "Identifying replicate datasets..." 
		experimentDict, k, r = dict(), 0, 0
		for dataset in datasetDict:
			if not dataset in exclusionSets:
				strain, factor, context, institute, method = dataset.split("_")
				combination = ":".join([factor, context])
				if not combination in experimentDict:
					experimentDict[combination] = list()
				experimentDict[combination].append(dataset)
				k += 1
		for combination in sorted(experimentDict.keys()):
			if len(experimentDict[combination]) > 1:
				print combination, ":", ", ".join(sorted(experimentDict[combination]))
				r += 1
		
		print
		print "Combinations:", len(experimentDict)
		print "Replicates:", r
		print "Datasets:", k
		print
		
		if option.library != "OFF":
			print "Scoring replicate overlaps..."
			import itertools
			import random
			replicateDict = dict()
			librarypath = idrpath + "final/" + option.library + "/"
			libraryfiles = os.listdir(librarypath)
			for combination in sorted(experimentDict.keys()):
				if len(experimentDict[combination]) > 1:
					for xombination in itertools.combinations(experimentDict[combination], 2):
						i, j = xombination
						ifile = option.peaks[:2] + "_" + i + "_peaks.bed"
						jfile = option.peaks[:2] + "_" + j + "_peaks.bed"
						istrain, ifactor, icontext, iinstitute, imethod = i.split("_")
						jstrain, jfactor, jcontext, jinstitute, jmethod = j.split("_")
						if imethod == jmethod:
							method = str(imethod)
						else:
							method = "mix"
						if istrain == jstrain:
							strain = "match"
						else:
							strain = "mix"
						
						process = True
						if option.exclude != "OFF":
							for exclusion in option.exclude.split(","):
								if exclusion in ifile or exclusion in jfile:
									process = False
						
						if process and datasetDict[i]['quality'] in ['0','1'] and datasetDict[j]['quality'] in ['0','1'] and ifile in libraryfiles and jfile in libraryfiles:
							ipeaks = librarypath + ifile
							jpeaks = librarypath + jfile
							xpeaks = reportspath + "mapdata_revise_" + option.peaks + "_overlap.bed"
							command = "intersectBed -a " + ipeaks + " -b " + jpeaks + " > " + xpeaks
							os.system(command)
							
							isize, jsize, xsize = len(open(ipeaks).readlines()), len(open(jpeaks).readlines()), len(open(xpeaks).readlines())
							key = ":".join(sorted([i,j]))
							if not combination in replicateDict:
								replicateDict[combination] = dict()
							replicateDict[combination][key] = min(float(xsize)/min(isize, jsize), 1)
							command = "rm -rf " + xpeaks
							os.system(command)
							
							output = [ifactor, icontext, 2, replicateDict[combination][key], method, strain, i, j]
							print >>p_output, "\t".join(map(str, output))
				
					if combination in replicateDict:
						replicateRanks = general.valuesort(replicateDict[combination])
						replicateRanks.reverse()
						i, j = replicateRanks[0].split(":")
						istrain, ifactor, icontext, iinstitute, imethod = i.split("_")
						jstrain, jfactor, jcontext, jinstitute, jmethod = j.split("_")
						if imethod == jmethod:
							method = str(imethod)
						else:
							method = "mix"
						if istrain == jstrain:
							strain = "match"
						else:
							strain = "mix"
						
						factor, context = combination.split(":")
						output = [factor, context, len(experimentDict[combination]), replicateDict[combination][replicateRanks[0]], method, strain, i, j]
						print >>r_output, "\t".join(map(str, output))
			
			"""
			iterations = 500
			samplingDict = dict()
			print "Performing random sampling..."
			for iteration in range(1, iterations):
				ifile = random.choice(libraryfiles)
				jfile = random.choice(libraryfiles)
				iorganism, istrain, ifactor, icontext, iinstitute, imethod = ifile.split("_")[:6]
				jorganism, jstrain, jfactor, jcontext, jinstitute, jmethod = jfile.split("_")[:6]
				while ifactor == jfactor:
					jfile = random.choice(libraryfiles)
					jorganism, jstrain, jfactor, jcontext, jinstitute, jmethod = jfile.split("_")[:6]
				ipeaks = librarypath + ifile
				jpeaks = librarypath + jfile
				xpeaks = reportspath + "mapdata_revise_" + option.peaks + "_overlap.bed"
				command = "intersectBed -a " + ipeaks + " -b " + jpeaks + " > " + xpeaks
				os.system(command)
				
				isize, jsize, xsize = len(open(ipeaks).readlines()), len(open(jpeaks).readlines()), len(open(xpeaks).readlines())
				output = [ifile, jfile, min(float(xsize)/min(isize, jsize), 1)]
				print >>s_output, "\t".join(map(str, output))
				command = "rm -rf " + xpeaks
				os.system(command)
			"""
			
		# close output files:
		d_output.close()
		l_output.close()
		p_output.close()
		r_output.close()
		#s_output.close()
		
	
	# report quality mode:
	elif option.mode == "report":
	
		print
		print "Loading revised configurations..."
		datasetDict = general.build2(reportspath + option.organism + "/" + "mapdata_revise_" + option.peaks + "_dataset_report.txt", id_column="dataset")
		libraryDict = general.build2(reportspath + option.organism + "/" + "mapdata_revise_" + option.peaks + "_library_report.txt", id_column="filename")
		
		# define output files:
		datasetColumns = open(reportspath + option.organism + "/" + "mapdata_revise_" + option.peaks + "_dataset_report.txt").readline().strip().split("\t")
		libraryColumns = open(reportspath + option.organism + "/" + "mapdata_revise_" + option.peaks + "_library_report.txt").readline().strip().split("\t")
		d_output = open(reportspath + option.organism + "/" + "mapdata_report_" + option.peaks + "_dataset_report.txt", "w")
		l_output = open(reportspath + option.organism + "/" + "mapdata_report_" + option.peaks + "_library_report.txt", "w")
		print >>d_output, "\t".join(datasetColumns)
		print >>l_output, "\t".join(libraryColumns)
		
		print "Loading peak set files of interest..."
		experimentDict = dict()
		peakfiles = os.listdir(peakspath + option.peaks)
		for peakfile in peakfiles:
			organism, strain, factor, context, institute, method = peakfile.split("_")[:6]
			experiment = factor + ":" + context
			if not experiment in experimentDict:
				experimentDict[experiment] = list()
			experimentDict[experiment].append(peakfile)
		
		print "Tallying peak set qualities..."
		matchedDict, missing = dict(), list()
		qualityDict = { "1":0, "0":0, "-1":0 }
		for dataset in datasetDict:
			
			# obtain factor and context for experiment:
			original = str(dataset)
			if option.rename != "OFF":
				for scheme in option.rename.split(","):
					target, replacement = scheme.split(":")
					dataset = dataset.replace(target, replacement)
			
			# update the relevant information if necessary...
			organism, strain, factor, context, institute, method = metrn.labelComponents(option.organism + "_" + "_".join(dataset.split("_")[:6]))
			datasetDict[dataset] = datasetDict[original]
			datasetDict[dataset]["strain"] = strain
			datasetDict[dataset]["factor"] = factor
			datasetDict[dataset]["stage"] = context
			datasetDict[dataset]["institute"] = institute
			datasetDict[dataset]["method"] = method
			if dataset != original:
				del datasetDict[original]
			
			# reconstruct dataset ID:
			datasetID = metrn.labelGenerator(target="dataset", mode="label", organism=option.organism, strain=datasetDict[dataset]["strain"], factor=datasetDict[dataset]["factor"], context=datasetDict[dataset]["stage"], institute=datasetDict[dataset]["institute"], method=datasetDict[dataset]["method"], separator="_")
			
			# check datasets of interest:
			if datasetID + "_peaks.bed" in peakfiles:
				qualityDict[datasetDict[dataset]["quality"]] += 1
				matchedDict[dataset] = datasetID + "_peaks.bed"
			else:
				experiment = factor + ":" + context
				if not experiment in experimentDict:
					missing.append(dataset)
				#if datasetDict[dataset]["quality"] != "-1":
				#	print experiment
				#	print organism, strain, factor, context, institute, method
				#	print original
				#	print dataset
				#	print datasetID
				#	pdb.set_trace()
					
		# export review information:
		print "Examining approved datasets not found..."
		for dataset in missing:
			if datasetDict[dataset]["quality"] != "-1":
				print dataset
		print
			
		print "Experiments:", len(experimentDict)
		print "Datasets (input):", len(datasetDict)
		print "Datasets (found):", len(matchedDict), "(" + str(round(float(100*len(matchedDict))/len(datasetDict), 0)) + "%)"
		print "Datasets (match):", len(matchedDict), "(" + str(round(float(100*len(matchedDict))/len(peakfiles), 0)) + "%)"
		print "High:", qualityDict["1"], "(" + str(round(float(100*qualityDict["1"])/sum(qualityDict.values()), 0)) + "%)"
		print "Mid:", qualityDict["0"], "(" + str(round(float(100*qualityDict["0"])/sum(qualityDict.values()), 0)) + "%)"
		print "Low:", qualityDict["-1"], "(" + str(round(float(100*qualityDict["-1"])/sum(qualityDict.values()), 0)) + "%)"
		print
		if len(matchedDict) != len(peakfiles):
			print "Note: There are datasets in the folder for which no information is available."
			print "Missing:", ", ".join(list(set(peakfiles).difference(set(matchedDict.values()))))
			print
		
		print "Exporting pruned dataset report..."
		for dataset in sorted(datasetDict.keys()):
			if dataset in matchedDict:
				output = list()
				for column in datasetColumns:
					if column == "dataset":
						output.append(dataset)
					else:
						output.append(datasetDict[dataset][column])
				print >>d_output, "\t".join(output)
				
		print "Exporting pruned library report..."
		for library in sorted(libraryDict.keys()):
			dataset = "_".join(library.split("_")[:5])
			if dataset in matchedDict:
				output = list()
				for column in libraryColumns:
					if column == "filename":
						output.append(library)
					else:
						output.append(libraryDict[library][column])
				print >>l_output, "\t".join(output)
		print
		
		# close output files:
		d_output.close()
		l_output.close()
	
		
	# waterston new genes mode:
	elif option.mode == "waterston":
	
		print
		print "Loading cellular expression data..."
		dataDict = general.build2(cellspath + "expression/" + option.infile, i="gene", j="cell", x="time.series", mode="matrix")
		
		print "Loading published expression data..."
		publishedData = list()
		for source in option.source.split(","):
			publishedData.extend(general.clean(open(extraspath + source).read().replace("\r","\n").split("\n")))
		publishedData = sorted(list(set(publishedData)))
		#print publishedData
		#print len(publishedData)
		
		oldGenes, oldTraces, xxxGenes, xxxTraces = list(), list(), list(), list()
		for gene in dataDict:
			for cell in dataDict[gene]:
				for trace in dataDict[gene][cell].split(","):
					if not trace in publishedData:
						xxxGenes.append(gene)
						xxxTraces.append(trace)
					else:
						oldGenes.append(gene)
						oldTraces.append(trace)
		newGenes = set(xxxGenes).difference(set(oldGenes))
		newTraces = set(xxxTraces).difference(set(oldTraces))
		oldGenes, oldTraces, newGenes, newTraces = list(set(oldGenes)), list(set(oldTraces)), list(set(newGenes)), list(set(newTraces))
		
		print "Total Genes:", len(dataDict), "(" + str(len(set(xxxTraces + oldTraces))) + " series)"
		print "Old Genes:", len(oldGenes), "(" + str(len(oldTraces)) + " series)"
		print "New Genes:", len(newGenes), "(" + str(len(newTraces)) + " series)"
		print
	
								
if __name__ == "__main__":
	main()
	
	print "Completed:", time.asctime(time.localtime())


#python mapData.py --path ~/meTRN/ --mode paralogs --organism ce --source extras --infile modencode.family.common.report.txt
#python mapData.py --path ~/meTRN/ --mode sunburst --organism ce --source extras --infile ce_modencode_report.tab --peaks ce_reporting_com_xx_raw
#python mapData.py --path ~/meTRN/ --mode sunburst --organism dm --source extras --infile dm_modencode_report.tab 
#python mapData.py --path ~/meTRN/ --mode sunburst --organism hs --source extras --infile hs_modencode_report.tab 
#python mapData.py --path ~/meTRN/ --mode reinke --source extras --infile reinke_missing_files.txt
#python mapData.py --path ~/meTRN/ --mode waterston --infile mapcells_avgExp_waterston_expression_assayed --source waterston_mace_published.txt,waterston_murray_published.txt
#python mapData.py --path ~/meTRN/ --mode scanner --peaks ce_selection_com_ex_raw
#python mapData.py --path ~/meTRN/ --mode scanner --peaks ce_selection_com_l1_raw
#python mapData.py --path ~/meTRN/ --mode scanner --peaks ce_selection_com_l2_raw
#python mapData.py --path ~/meTRN/ --mode scanner --peaks ce_selection_com_l3_raw
#python mapData.py --path ~/meTRN/ --mode scanner --peaks ce_selection_com_l4_raw
#python mapData.py --path ~/meTRN/ --mode scanner --peaks ce_selection_com_cx_raw



