#!/usr/bin/env python
# extract GO analysis from R script outputs!

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

import simplejson as json
import palette

print "Command:", " ".join(sys.argv)
print "Timestamp:", time.asctime(time.localtime())

# define a function to extract shared factors:
def signatureSharing(sign1, sign2, factors):
	vector1, vector2, factors = sign1.split(","), sign2.split(","), factors.split(",")
	shared, total = list(), len(vector1)
	if len(vector1) == len(vector2):
		for i in range(0, total):
			if vector1[i] == "1" and vector2[i] == "1":
				shared.append(factors[i])
	return shared
	
# define a function to extract union of factors:
def signatureUnion(sign1, sign2, factors):
	vector1, vector2, factors = sign1.split(","), sign2.split(","), factors.split(",")
	union, total = list(), len(vector1)
	if len(vector1) == len(vector2):
		for i in range(0, total):
			if vector1[i] == "1" or vector2[i] == "1":
				union.append(factors[i])
	return union

		
def main():
	
	parser = optparse.OptionParser()
	parser.add_option("--path", action = "store", type = "string", dest = "path", help = "Path from script to files")
	parser.add_option("--mode", action = "store", type = "string", dest = "mode", help = "Operation mode", default="report")
	parser.add_option("--peaks", action = "store", type = "string", dest = "peaks", help = "Basename for target peaks", default="OFF")
	parser.add_option("--organism", action = "store", type = "string", dest = "organism", help = "Target organism for operations...", default="OFF")
	parser.add_option("--analysis", action = "store", type = "string", dest = "analysis", help = "What GO analysis results should be used?")
	parser.add_option("--target", action = "store", type = "string", dest = "target", help = "How should the parsing be performed?", default="dataset")
	parser.add_option("--minPvalue", action = "store", type = "float", dest = "minPvalue", help = "Minimum P-value", default=0.01)
	parser.add_option("--hitCount", action = "store", type = "int", dest = "hitCount", help = "Minimum number of significance hits for each GO ID", default=1)
	parser.add_option("--hitPvalue", action = "store", type = "float", dest = "hitPvalue", help = "GO IDs must have at least one hit with a P-value below this cutoff", default=0.01)
	parser.add_option("--include", action = "store", type = "string", dest = "include", help = "Include GO-results files with this in name; comma-separated", default=False)
	parser.add_option("--exclude", action = "store", type = "string", dest = "exclude", help = "Exclude GO-results files with this in name; comma-separated", default=False)
	parser.add_option("--minCount", action = "store", type = "int", dest = "minCount", help = "Minimum number of GO term instances", default=0)
	parser.add_option("--maxCount", action = "store", type = "int", dest = "maxCount", help = "Maximum number of GO term instances", default=1000000000000)
	parser.add_option("--maxColor", action = "store", type = "float", dest = "maxColor", help = "Maximum value for coloring", default=100)
	parser.add_option("--name", action = "store", type = "string", dest = "name", help = "Output file name-tag.", default="OFF")
	parser.add_option("--threads", action = "store", type = "int", dest = "threads", help = "Parallel processing threads", default=1)
	parser.add_option("--chunks", action = "store", type = "int", dest = "chunks", help = "", default=100)
	parser.add_option("--module", action = "store", type = "string", dest = "module", help = "", default="md1")
	parser.add_option("--qsub", action = "store", type = "string", dest = "qsub", help = "Qsub configuration header", default="OFF")
	parser.add_option("--server", action = "store", type = "string", dest = "server", help = "Are we on the server?", default="OFF")
	parser.add_option("--job", action = "store", type = "string", dest = "job", help = "Job name for cluster", default="OFF")
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
	peakspath = path_dict["peaks"]
	gopath = path_dict["go"]
	hotpath = path_dict["hot"]
	qsubpath = path_dict["qsub"]
	neuronspath = path_dict["neurons"]
	
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

	# report mode:
	if option.mode == "report":
	
		# define report path:
		gopath = gopath + option.peaks + "/" + option.analysis + "/"
		
		# setup output key:
		hitPvalue_handle = "%.0e" % (float(option.hitPvalue))
		configuration_handle = "p" + str(option.minPvalue).split(".")[1] + "_hc" + str(option.hitCount) + "_hp" + hitPvalue_handle
		s_outfile = "mapgo_complete_" + option.peaks + "_" + configuration_handle + "_summary"
		m_outfile = "mapgo_complete_" + option.peaks + "_" + configuration_handle + "_matrix"
		b_outfile = "mapgo_complete_" + option.peaks + "_" + configuration_handle + "_matrix_bp"
		c_outfile = "mapgo_complete_" + option.peaks + "_" + configuration_handle + "_matrix_cc"
		f_outfile = "mapgo_complete_" + option.peaks + "_" + configuration_handle + "_matrix_mf"
		
		# generate input and output paths:
		resultspath = gopath + "results/"
		summarypath = gopath + "summary/"
		sunburstpath = gopath + "sunburst/"
		general.pathGenerator(resultspath)
		general.pathGenerator(summarypath)
		general.pathGenerator(sunburstpath)
		
		# open input and output files:
		s_output = open(summarypath + s_outfile, "w")
		m_output = open(summarypath + m_outfile, "w")
		b_output = open(summarypath + b_outfile, "w")
		c_output = open(summarypath + c_outfile, "w")
		f_output = open(summarypath + f_outfile, "w")
		
		# define and print output headers:
		go_headers = ["id","term","ontology","pvalue","dataset.count","genome.count","dataset.total","genome.total","adjusted.pvalue","entrezID.collapsed","entrezID.count"]
		print >>s_output, "\t".join(["organism", "strain","factor","context","institute","method"] + go_headers)
		
		# scan GO analysis outputs:
		go_dict, ontology_dict, experiment_dict  = dict(), dict(), dict()
		infiles = os.listdir(resultspath)
		mx = 1
		print
		print "Processing:", len(infiles), "files..."
		for infile in sorted(infiles):
			
			include = True
			if option.exclude:
				for exlcusionString in option.exclude.split(","):
					if exlcusionString in infile:
						include = False
		
			if "_enrichedGO" in infile and not "_enrichedGO_p" in infile and include:
				
				# get experiment:
				if "_enrichedGO_bp" in infile:
					experiment = infile.replace("mapgo_","").replace("_enrichedGO_bp","")
				elif "_enrichedGO_cc" in infile:
					experiment = infile.replace("mapgo_","").replace("_enrichedGO_cc","")
				elif "_enrichedGO_mf" in infile:
					experiment = infile.replace("mapgo_","").replace("_enrichedGO_mf","")
				else:
					print infile
					print "Error: Experiment could not be determined!"
					pdb.set_trace()
				
				print "...File #" + str(mx) + ":", experiment
				mx += 1
				
				# set dataset identifier:
				if option.target == "hot.regions":
					label, script, overlap, organism, selection, types, context, method, region, subset = experiment.split("_")
					datasetID = ":".join([context, subset])
				elif option.target == "som.neurons":
					datasetID = experiment.replace("input_neuron", "N")
				else:
					organism, strain, factor, context, institute, method = metrn.labelComponents(experiment.split(option.peaks)[1].strip("_"))
					datasetID = metrn.labelGenerator(target=option.target, mode="label", organism=organism, strain=strain, factor=factor, context=context, institute=institute, method=method)
					
				# load dataset into GO results:
				if not datasetID in go_dict:
					go_dict[datasetID] = dict()
				
				# store experiment information:
				experiment_dict[datasetID] = experiment
				
				# build header dictionary:
				header_dict = general.build_header_dict(resultspath + infile)
				
				# process lines:
				inlines = open(resultspath + infile, "U").readlines()
				inlines.pop(0)
				for inline in inlines:
					initems  = inline.strip().split("\t")
					adjPvalue = float(initems[header_dict["BH.adjusted.p.value"]])
					if adjPvalue <= option.minPvalue:	
						goID = initems[header_dict["go.id"]]
						if not goID in go_dict[datasetID]:
							go_dict[datasetID][goID] = dict()
							go_dict[datasetID][goID]["entrezID.collapsed"] = list()
							go_dict[datasetID][goID]["experiment"] = list()
						go_dict[datasetID][goID]["id"] = initems[header_dict["go.id"]]
						go_dict[datasetID][goID]["term"] = initems[header_dict["go.term"]]
						go_dict[datasetID][goID]["ontology"] = initems[header_dict["Ontology"]]
						go_dict[datasetID][goID]["pvalue"] = float(initems[header_dict["pvalue"]])
						go_dict[datasetID][goID]["dataset.count"] = int(initems[header_dict["count.InDataset"]])
						go_dict[datasetID][goID]["genome.count"] = int(initems[header_dict["count.InGenome"]])
						go_dict[datasetID][goID]["dataset.total"] = int(initems[header_dict["totaltermInDataset"]])
						go_dict[datasetID][goID]["genome.total"] = int(initems[header_dict["totaltermInGenome"]])
						go_dict[datasetID][goID]["adjusted.pvalue"] = float(initems[header_dict["BH.adjusted.p.value"]])
						go_dict[datasetID][goID]["entrezID.collapsed"].append(initems[header_dict["EntrezID"]])
						go_dict[datasetID][goID]["entrezID.count"] = len(go_dict[datasetID][goID]["entrezID.collapsed"])
						go_dict[datasetID][goID]["experiment"].append(experiment)
						
						# store gene ontology information:
						ontology = initems[header_dict["Ontology"]]
						if not ontology in ontology_dict:
							ontology_dict[ontology] = dict()
						if not goID in ontology_dict[ontology]:
							ontology_dict[ontology][goID] = dict()
						ontology_dict[ontology][goID]["term"] = initems[header_dict["go.term"]]
						ontology_dict[ontology][goID]["genome.count"] = initems[header_dict["count.InGenome"]]
					
		# order goIDs in ontologies and by occurrences:
		goIDs = list()
		print
		for ontology in sorted(ontology_dict):
			go_counts_dict = dict()
			for goID in sorted(ontology_dict[ontology]):
				go_counts_dict[goID] = int(ontology_dict[ontology][goID]["genome.count"])
			goIDs.extend(general.valuesort(go_counts_dict))
			print "Ontology", ontology + ":", len(go_counts_dict)
		print "Total (list):", len(goIDs)
		print "Total (set):", len(set(goIDs))
		print
		
		# filter goIDs based on the number of enriched datasets:
		filteredGOids = list()
		for goID in goIDs:
			hitCounts, hitFilter = 0, False
			for datasetID in sorted(go_dict.keys()):
				if goID in go_dict[datasetID]:
					hitCounts += 1
					if go_dict[datasetID][goID]["adjusted.pvalue"] < option.hitPvalue:
						hitFilter = True
			if hitCounts >= int(option.hitCount) and hitFilter:
				filteredGOids.append(goID)
		
		# filter enrichment data:
		filtered_dict = dict()
		for datasetID in go_dict:
			for goID in go_dict[datasetID]:
				if goID in filteredGOids:
					if not datasetID in filtered_dict:
						filtered_dict[datasetID] = dict()
					filtered_dict[datasetID][goID] = float(go_dict[datasetID][goID]["adjusted.pvalue"])
		
		# print out filtered enrichments to a matrix:
		matrixHeader = ["dataset", "organism", "strain", "factor", "context", "institute", "method", "go.id", "go.term", "ontology", "go.label", "go.count", "pvalue", "adj.pvalue"]
		print >>m_output, "\t".join(matrixHeader)	
		print >>b_output, "\t".join(matrixHeader)	
		print >>c_output, "\t".join(matrixHeader)	
		print >>f_output, "\t".join(matrixHeader)	
		for datasetID in sorted(filtered_dict.keys()):
			experiment = experiment_dict[datasetID]
			
			if option.target == "hot.regions":
				label, script, overlap, organism, institute, factor, context, method, region, strain = experiment.split("_")
			elif option.target == "som.neurons":
				organism, strain, factor, context, institute, method = datasetID, "na", datasetID, "na", "na", "na"
			else:
				organism, strain, factor, context, institute, method = metrn.labelComponents(experiment.split(option.peaks)[1].strip("_"))
			
			for goID in goIDs:
				if goID in filteredGOids:
					
					# gather basic go term information:
					go_onto, go_term, go_count = False, False, False
					for ontology in ontology_dict:
						if goID in ontology_dict[ontology]:
							go_onto = ontology
							go_term = ontology_dict[ontology][goID]["term"]
							go_count = ontology_dict[ontology][goID]["genome.count"]
					go_label = '"' + go_term + ' ' + go_onto.lower() + '"'
					#go_label = '"' + go_term + ' ' + '[' + go_onto + ']"'
					
					# gather p-values where enrichments have occurred. otherwise default to 1:
					if goID in go_dict[datasetID]:
						pvalue = go_dict[datasetID][goID]["pvalue"]
						adjusted_pvalue = go_dict[datasetID][goID]["adjusted.pvalue"]
					else:
						pvalue = "1"
						adjusted_pvalue = "1"
				
					# export data to matrix:
					output = map(str, [datasetID, organism, strain, factor, context, institute, method, goID, '"' + go_term + '"', go_onto, go_label, go_count, pvalue, adjusted_pvalue])
					if go_onto == "BP":
						print >>b_output, "\t".join(output)
					elif go_onto == "CC":
						print >>c_output, "\t".join(output)
					elif go_onto == "MF":
						print >>f_output, "\t".join(output)
					print >>m_output, "\t".join(output)
					
		# print out filtered enrichments to summary report:
		processed = list()
		for datasetID in sorted(filtered_dict.keys()):
			experiment = experiment_dict[datasetID]
			if option.target == "hot.regions":
				label, script, overlap, organism, institute, factor, context, method, region, strain = experiment.split("_")
			elif option.target == "som.neurons":
				organism, strain, factor, context, institute, method = datasetID, "na", datasetID, "na", "na", "na"
			else:
				organism, strain, factor, context, institute, method = metrn.labelComponents(experiment.split(option.peaks)[1].strip("_"))
			
			e_outfile = "mapgo_distinct_" + experiment + "_enrichedGO_" + configuration_handle
			e_output = open(summarypath + e_outfile, "w")
			print >>e_output, "\t".join(["organism", "strain","factor","context","institute","method"] + go_headers)
			
			processedID = "_".join([organism, strain, factor, context, institute, method])
			if not processedID in processed:
				processed.append(processedID)
				sortedGOids = general.valuesort(filtered_dict[datasetID])
				for goID in sortedGOids:
					output = [organism, strain, factor, context, institute, method]
					for header in go_headers:
						if "collapsed" in header:
						 	#output.append(",".join(go_dict[datasetID][goID][header]))
						 	output.append(go_dict[datasetID][goID]["entrezID.count"])
						else:
						 	output.append(go_dict[datasetID][goID][header])
					if go_dict[datasetID][goID]["adjusted.pvalue"] < option.hitPvalue:
						print >>s_output, "\t".join(map(str, output))
						print >>e_output, "\t".join(map(str, output))
			e_output.close()
		
		# close output file:
		s_output.close()
		m_output.close()
		b_output.close()
		c_output.close()
		f_output.close()
		
		
	# sunburst mode:
	if option.mode == "sunburst":
	
		from pylab import *
		import color
		
		startcolor = '#586323'  # a dark olive 
		midcolor = '#fcffc9'    # a bright yellow
		endcolor = '#bd2309'    # medium dark red
		oliveColors = [startcolor,midcolor,endcolor]
		colorPalette = color.gradient("oliveColors", oliveColors)
		
		#print colorPalette(float(50)/N)
		#pdb.set_trace()
		
		# preload GO results:
		print
		print "Loading GO summary results..."
		goIDs, goMaster, goScores, goCounts, goOntologies = list(), dict(), dict(), dict(), dict()
		for inputSet in option.peaks.split(","):
	
			# define report path:
			gopath = path_dict["go"]
			gopath = gopath + inputSet + "/" + option.analysis + "/"
			
			# setup output key:
			hitPvalue_handle = "%.0e" % (float(option.hitPvalue))
			configuration_handle = "p" + str(option.minPvalue).split(".")[1] + "_hc" + str(option.hitCount) + "_hp" + hitPvalue_handle
			s_outfile = "mapgo_complete_" + inputSet + "_" + configuration_handle + "_summary"
			m_outfile = "mapgo_complete_" + inputSet + "_" + configuration_handle + "_matrix"
			b_outfile = "mapgo_complete_" + inputSet + "_" + configuration_handle + "_matrix_bp"
			c_outfile = "mapgo_complete_" + inputSet + "_" + configuration_handle + "_matrix_cc"
			f_outfile = "mapgo_complete_" + inputSet + "_" + configuration_handle + "_matrix_mf"
			
			# generate input and output paths:
			resultspath = gopath + "results/"
			summarypath = gopath + "summary/"
			sunburstpath = gopath + "sunburst/"
			general.pathGenerator(resultspath)
			general.pathGenerator(summarypath)
			general.pathGenerator(sunburstpath)
			
			# load GO summary data:
			goMatrix = general.build2(summarypath + m_outfile, i="go.term", j="dataset", x="adj.pvalue", mode="matrix")
			goTables = general.build2(summarypath + m_outfile, i="go.term", j="ontology", mode="matrix", counter=True, skip=True, verbose=False)
			goLimits = general.build2(summarypath + m_outfile, i="go.term", j="go.count", mode="matrix", counter=True, skip=True, verbose=False)
			print inputSet, ":", len(goMatrix)
			
			# add GO results to matrix:
			for goID in goMatrix:
				if not goID in goMaster:
					goMaster[goID] = dict()
					goIDs.append(goID)
				if not inputSet in goMaster[goID]:
					goMaster[goID][inputSet] = dict()
				for dataset in goMatrix[goID]:
					goMaster[goID][inputSet][dataset] = float(goMatrix[goID][dataset])
			
			# store ontology information:
			for goID in goTables:
				goOntologies[goID] = goTables[goID].keys()[0]
			
			# store counts information:
			for goID in goLimits:
				goCounts[goID] = int(goLimits[goID].keys()[0])
				
			# store scores:
			for goID in goMatrix:
				if not goID in goScores:
					goScores[goID] = 1
				for dataset in goMatrix[goID]:
					goScores[goID] = min(goScores[goID], float(goMatrix[goID][dataset]))
			
		# define color function:
		def colorMaker(organismTag, context, colors, exception="#B8BEC2"): 
			if organismTag + "." + context in colors:
				return colors[organismTag + "." + context]
			else:
				return exception
		
		# define colors specifically:
		colors = {
			"ce.EE" : "#F3C520",
			"ce.LE" : "#F3D020",
			"ce.EM" : "#F3DB20",
			"ce.EX" : "#F3DB20",
			"ce.L1" : "#5BD4E7",
			"ce.L2" : "#5BD4E7",
			"ce.L3" : "#5BD4E7",
			"ce.L4" : "#5BD4E7",
			
			# Brasil:
			#"BP" : "#2BA01C",
			#"MF" : "#FFB700",
			#"CC" : "#1BA2E0"
			
			# Light-RGB:
			#"BP" : "#FF9987",
			#"MF" : "#87D785",
			#"CC" : "#85CAD7"
			
			# RGB:
			#"BP" : "#1A8191",
			#"MF" : "#40BBCE",
			#"CC" : "#2B7DDB"
			
			# Gray-scale:
			#"BP" : "#3C3C3C",
			#"MF" : "#B5B5B5",
			#"CC" : "#EEEEEE"
			
			# Black-n-White:
			#"BP" : "#FEFEFE",
			#"MF" : "#B5B5B5",
			#"CC" : "#1E1E1E"
			
			# Teals:
			#"BP" : "#1A8191",
			#"MF" : "#40BBCE",
			#"CC" : "#A7E8F2"
			
			# Wolfgang:
			"BP" : "#EDFFA3",
			"MF" : "#57CCD9",
			"CC" : "#2991E1"
			}
			
		# define blue-yellow color gradients:
		#gradients = {
		#	"EX" : ["#FFE600", "#FFB700", "#FF6F00"],
		#	"L1" : ["#8DECF3", "#1BA2E0", "#296DCC"],
		#	"L2" : ["#8DECF3", "#1BA2E0", "#296DCC"],
		#	"L3" : ["#8DECF3", "#1BA2E0", "#296DCC"],
		#	"L4" : ["#8DECF3", "#1BA2E0", "#296DCC"]
		#	}
		
		# define blue-red color gradients:
		#gradients = {
		#	"EX" : ["#F28E6F", "#F00000", "#B00E0E"],
		#	"L1" : ["#79AED1", "#477FBF", "#294380"],
		#	"L2" : ["#79AED1", "#477FBF", "#294380"],
		#	"L3" : ["#79AED1", "#477FBF", "#294380"],
		#	"L4" : ["#79AED1", "#477FBF", "#294380"]
		#	}
			
		# define blue-orange color gradients:
		gradients = {
			"EX" : ["#FF9500", "#FF421C", "#C40E0E"],
			"L1" : ["#4F88C9", "#2063B0", "#294380"],
			"L2" : ["#4F88C9", "#2063B0", "#294380"],
			"L3" : ["#4F88C9", "#2063B0", "#294380"],
			"L4" : ["#4F88C9", "#2063B0", "#294380"]
			}
			
		# generate sunburst matrix:
		print
		print "Generating sunburst diagram..."
		print "Combined GO terms:", len(goScores)
		sunburst = { "name":"sunburst", "color":"#B8BEC2", "children":list() }
		k = 0
		goIDs = list()
		goOrder = general.valuesort(goScores)
		for goID in goOrder:
			datasets = list()
			colorOntology = colors[goOntologies[goID]]
			#if goCounts[goID] >= 50 and goCounts[goID] <= 500:
			#if goCounts[goID] >= 5 and goCounts[goID] <= 250:
			if goCounts[goID] >= option.minCount and goCounts[goID] <= option.maxCount:
				k += 1
				inputSets = option.peaks.split(",")
				inputSets.reverse()
				for inputSet in inputSets:
					context = inputSet.split(".")[1].upper()
					colorContext = palette.Color(colorMaker(organismTag, context, colors))
					colorPalette = color.gradient("contextColors", gradients[context.upper()])
					
					if inputSet in goMaster[goID]:
						colorScore = min(goMaster[goID][inputSet].values())
						colorScore = -numpy.log10(colorScore)
						if colorScore > option.maxColor:
							colorScore = option.maxColor
						colorCode = matplotlib.colors.rgb2hex(colorPalette(float(colorScore)/option.maxColor))
						
					else:
						colorCode = palette.Color("#FFFFFF").hex
						#print "GO term not found in input set:", inputSet, goID
						#pdb.set_trace()
					
					if inputSet == inputSets[len(inputSets)-1]:
						label = str(goID.replace('"',''))
					else:
						label = ""
						
					if datasets == list():
						datasets = [{"name":context, "label":label, "color":colorCode}]
						#print "1:", context
					else:
						datasets = [{"name":context, "label":label, "color":colorCode, "children":datasets}]
						#print "2:", context
					#if colorValue.hex == "#FFFFFF":
					#	print inputSet, goID, colorScore, colorValue.hex
					#	pdb.set_trace()
						
				#print
				#print goID
				#print datasets
				#pdb.set_trace()	
				goIDs.append({"name":goID.replace('"', ''), "color":colorOntology, "children":datasets})
		sunburst["children"] = goIDs
		
		print "Exported GO terms:", k
		print "Writing out..."
		# export sunburst json file (to report path):
		f_output = open(sunburstpath + "mapgo_sunburst_" + option.target + ".json", "w")
		json.dump(sunburst, f_output)
		f_output.close()
	
		# export sunburst json file (to web directory):
		sunburstpath = "/Library/WebServer/Documents/metrn/data/test/"
		f_output = open(sunburstpath + "mapgo_sunburst_" + option.target + ".json", "w")
		json.dump(sunburst, f_output)
		f_output.close()
		print
	
	# compile mode:
	if option.mode == "compile":
	
		# specify master dictionary:
		masterDict = dict()
		
		# set P-value and configuration handles:
		hitPvalue_handle = "%.0e" % (float(option.hitPvalue))
		configuration_handle = "p" + str(option.minPvalue).split(".")[1] + "_hc" + str(option.hitCount) + "_hp" + hitPvalue_handle
			
		# process SOMs:
		for somName in option.peaks.split(","):
		
			# define report path:
			gopath = path_dict["go"] + somName + "/" + option.analysis + "/"
			
			# setup output key:
			s_outfile = "mapgo_complete_" + somName + "_" + configuration_handle + "_summary"
			m_outfile = "mapgo_complete_" + somName + "_" + configuration_handle + "_matrix"
			b_outfile = "mapgo_complete_" + somName + "_" + configuration_handle + "_matrix_bp"
			c_outfile = "mapgo_complete_" + somName + "_" + configuration_handle + "_matrix_cc"
			f_outfile = "mapgo_complete_" + somName + "_" + configuration_handle + "_matrix_mf"
			
			# generate input and output paths:
			resultspath = gopath + "results/"
			summarypath = gopath + "summary/"
			sunburstpath = gopath + "sunburst/"
			
			# load GO enrichment results:
			goDict = general.build2(summarypath + s_outfile, id_complex=["organism","id"], separator="-")
			
			# transfer GO results:
			for goKey in goDict:
				somNeuron, goID = goKey.split("-")
				if not goID in masterDict:
					masterDict[goID] = dict()
				if not somName in masterDict[goID]:
					masterDict[goID][somName] = dict()
				masterDict[goID][somName][somNeuron] = [goDict[goKey]["ontology"], goDict[goKey]["term"], goDict[goKey]["dataset.count"], goDict[goKey]["genome.count"], goDict[goKey]["adjusted.pvalue"]]
		
		# store counts per signature, and signatures per count:
		countDict, signDict, contextDict = dict(), dict(), dict()
		t = 0
		for goID in masterDict:
			n = 0
			c = len(masterDict[goID])
			for somName in masterDict[goID]:
				for somNeuron in masterDict[goID][somName]:
					n += 1
					t += 1
			if not n in countDict:
				countDict[n] = 0
			if not c in contextDict:
				contextDict[c] = 0
			countDict[n] += 1
			contextDict[c] += 1
			signDict[goID] = n
		counts = general.valuesort(signDict)
		counts.reverse()
			
		print
		print "Total GO Terms:", len(masterDict)
		print "Total Neurons:", t
		print "Specificity:", round(100*float(len(masterDict))/t, 2), "%"
		print
		print "Maximum Neurons per Signature:", max(countDict)
		print "Minimum Neurons per Signature:", min(countDict)
		print
		print "Neurons per GO Term:"
		for count in sorted(countDict.keys()):
			print count, ":", countDict[count], "GO terms"
		print
		print "Contexts per GO Term:"
		for count in sorted(contextDict.keys()):
			print count, ":", contextDict[count], "GO terms"
		print	
		
		# load neuron co-association signatures:
		neuronCompiler = neuronspath + "compile/" + option.name + "/mapneurons_compile_" + option.name + "_signs.txt"
		neuronSigns = general.build2(neuronCompiler, id_complex=["signature", "som.name"], separator=":")
		
		# store signatures per SOM, per neuron:
		neuronDict = dict()
		for neuronSign in neuronSigns:
			somSignature, somName = neuronSign.split(":")
			if not somName in neuronDict:
				neuronDict[somName] = dict()
			for somNeuron in neuronSigns[neuronSign]["neuron.ids"].split(","):
				neuronDict[somName][somNeuron] = somSignature
			somFactors = neuronSigns[neuronSign]["factor.ids"]
			
		# generate output file:
		compilepath = path_dict["go"] + "compile/" + option.name + "/"
		general.pathGenerator(compilepath)
		f_outfile = compilepath + "mapgo_compile_" + option.name + "_shared.txt"
		f_output = open(f_outfile, "w")
		print >>f_output, "\t".join(["go.id", "go.term", "go.ontology", "go.count", "som.count", "som.i", "som.j", "shared.factors", "shared.fraction", "shared.ids"])
			
		# examine correspondence between GO enrichments and SOM signatures:
		print
		somDict = dict()
		for goID in masterDict:
			somNames = sorted(masterDict[goID].keys())
			if len(somNames) > 1:
			
				print "Processing:", goID
			
				# collect target SOM:neurons:
				somTargets = list()
				for somName in somNames:
					for somNeuron in masterDict[goID][somName]:
						somTargets.append(somName + ":" + somNeuron)
				# Notes: These are the SOM:neurons with the specific GO term enrichment.
				
				# get info for the GO id:
				ontology, term, datasetCount, genomeCount, adjPvalue = masterDict[goID][somName][somNeuron]
				
				# gather signatures from SOM:neurons to compare:
				somSignatures = dict()
				for somTarget in somTargets:
					somName, somNeuron = somTarget.split(":")
					#print somTarget, somName, somNeuron
					if not somName in somSignatures:
						somSignatures[somName] = dict()
					somSignatures[somName][somNeuron] = neuronDict[somName][somNeuron]
				
				# select the least discrepant signatures between the two contexts (SOMs):
				somDistances, somSharing, somFraction = dict(), dict(), dict()
				for somNameA in somSignatures:
					for somNeuronA in somSignatures[somNameA]:
						for somNameB in somSignatures:
							for somNeuronB in somSignatures[somNameB]:
								
								if somNameA != somNameB:
									somCompareA, somCompareB = somNameA + ":" + somNeuronA, somNameB + ":" + somNeuronB 
									somShared = signatureSharing(somSignatures[somNameA][somNeuronA], somSignatures[somNameB][somNeuronB], somFactors)
									somUnion = signatureUnion(somSignatures[somNameA][somNeuronA], somSignatures[somNameB][somNeuronB], somFactors)
									somDistances[somCompareA + "-" + somCompareB] = len(somShared)
									somSharing[somCompareA + "-" + somCompareB] = somShared
									somFraction[somCompareA + "-" + somCompareB] = float(len(somShared))/len(somUnion)
				
				somMatches = general.valuesort(somDistances)
				somMatches.reverse()
				somMatch = somMatches[0]
				
				# export shared factor signature (score and IDs) for the GO term:
				somID1, somID2 = somMatch.split("-")
				output = [goID, term, ontology, genomeCount, len(somNames), somID1, somID2, somDistances[somMatch], somFraction[somMatch], ",".join(somSharing[somMatch])]
				print >>f_output, "\t".join(map(str, output))
				#print somDistances
				#pdb.set_trace()
				
		# close output file
		f_output.close()
		print
		#print neuronMatrix.keys()
		#pdb.set_trace()
		
	# compare mode:
	if option.mode == "compare":
	
		# specify master dictionary:
		masterDict = dict()
		
		# set P-value and configuration handles:
		hitPvalue_handle = "%.0e" % (float(option.hitPvalue))
		configuration_handle = "p" + str(option.minPvalue).split(".")[1] + "_hc" + str(option.hitCount) + "_hp" + hitPvalue_handle
		
		# pre-load peak set and SOM:
		peakSet, somName = option.peaks.split(",")
		inputDict = {
			"peakSet" : peakSet,
			"somName" : somName,
			}
		
		# load context, explicitly:
		mainContext = option.target
		codeContexts, targetContexts = metrn.options_dict["contexts.extended"][option.target]
		
		# process SOMs:
		print
		for inputType in inputDict:
			print "Loading", inputType, ":", inputDict[inputType]
			inputKey = inputDict[inputType]
			
			# define report path:
			gopath = path_dict["go"] + inputKey + "/" + option.analysis + "/"
			
			# setup output key:
			s_outfile = "mapgo_complete_" + inputKey + "_" + configuration_handle + "_summary"
			m_outfile = "mapgo_complete_" + inputKey + "_" + configuration_handle + "_matrix"
			b_outfile = "mapgo_complete_" + inputKey + "_" + configuration_handle + "_matrix_bp"
			c_outfile = "mapgo_complete_" + inputKey + "_" + configuration_handle + "_matrix_cc"
			f_outfile = "mapgo_complete_" + inputKey + "_" + configuration_handle + "_matrix_mf"
			
			# generate input and output paths:
			resultspath = gopath + "results/"
			summarypath = gopath + "summary/"
			sunburstpath = gopath + "sunburst/"
			
			# load GO enrichment results:
			goDict = general.build2(summarypath + s_outfile, id_complex=["factor","context","id"], separator=" ")
			
			# transfer GO results:
			for goKey in goDict:
				factor, context, goID = goKey.split(" ")
				
				# correct context names for SOMs:
				if context == "na":
					context = str(codeContexts)
				keyName = factor + " (" + context + ")"
				
				# sub-select matching context:
				if (context in targetContexts or inputType == "somName") and (option.include == "OFF" or goDict[goKey]["ontology"] in option.include.upper()):
					#if inputKey != "any.ex.som":
					#	print inputKey, factor, context, goID
					#	pdb.set_trace()
					if not goID in masterDict:
						masterDict[goID] = dict()
					if not inputKey in masterDict[goID]:
						masterDict[goID][inputKey] = dict()
					masterDict[goID][inputKey][keyName] = [goDict[goKey]["ontology"], goDict[goKey]["term"], goDict[goKey]["dataset.count"], goDict[goKey]["genome.count"], goDict[goKey]["adjusted.pvalue"]]
		
		# tally GO enrichment diversity:
		somCount, peakCount, dualCount = list(), list(), list()
		for goID in masterDict:
			if len(masterDict[goID]) == 2:
				dualCount.append(goID)
			elif masterDict[goID].keys() == [inputDict["somName"]]:
				somCount.append(goID)
			elif masterDict[goID].keys() == [inputDict["peakSet"]]:
				peakCount.append(goID)
		
		print
		print "GO terms found in SOM:", len(somCount)
		print "GO terms found in peak-data:", len(peakCount)
		print "GO terms found in both:", len(dualCount)
		print
		
		pdb.set_trace()
		
		#x = masterDict.keys()[0]
		#print x
		#print masterDict[x]
		#pdb.set_trace()
	
	
	# example finding mode:
	if option.mode == "example":
	
		# specify master dictionary:
		masterDict = dict()
		
		# set P-value and configuration handles:
		hitPvalue_handle = "%.0e" % (float(option.hitPvalue))
		configuration_handle = "p" + str(option.minPvalue).split(".")[1] + "_hc" + str(option.hitCount) + "_hp" + hitPvalue_handle
		
		# define report path:
		gopath = path_dict["go"] + option.peaks + "/" + option.analysis + "/"
			
		# setup output key:
		s_outfile = "mapgo_complete_" + option.peaks + "_" + configuration_handle + "_summary"
		m_outfile = "mapgo_complete_" + option.peaks + "_" + configuration_handle + "_matrix"
		b_outfile = "mapgo_complete_" + option.peaks + "_" + configuration_handle + "_matrix_bp"
		c_outfile = "mapgo_complete_" + option.peaks + "_" + configuration_handle + "_matrix_cc"
		f_outfile = "mapgo_complete_" + option.peaks + "_" + configuration_handle + "_matrix_mf"
			
		# generate input and output paths:
		resultspath = gopath + "results/"
		summarypath = gopath + "summary/"
		sunburstpath = gopath + "sunburst/"
			
		# load GO enrichment results:
		#goDict = general.build2(summarypath + s_outfile, id_complex=["factor","context","term"], separator=" ")
		
		# process inputs:
		print
		print "Scanning GO results for entrezID hits (" + hitPvalue_handle + ")"
		masterDict, scorerDict, entrezDict = dict(), dict(), dict()
		targets = option.target.split(",")
		for infile in os.listdir(resultspath):
			
			# get datasetID
			datasetID = infile.split("_enrichedGO")[0]
			datasetID = datasetID.split(option.peaks + "_")[1]
			
			# determine if factor is of interest:
			process = False
			for factor in option.include.split(","):
				if factor in datasetID:
					process = True
			
			if process:
	
				# load GO enrichment results:
				goDict = general.build2(resultspath + infile, id_complex=["go.id", "go.term", "EntrezID"], separator=";")
					
				# transfer GO results:
				for goKey in goDict:
					goID, goTerm, entrezID = goKey.split(";")
					
					# check whether this is a term of interest:
					process = False
					for target in targets:
						if target in goTerm:
							process = True
						
					# transfer results:
					if process and float(goDict[goKey]["BH.adjusted.p.value"]) < float(option.hitPvalue):
						if not goTerm in masterDict:
							masterDict[goTerm] = dict()
							scorerDict[goTerm] = dict()
							entrezDict[goTerm] = dict()
						if not entrezID in scorerDict[goTerm]:
							scorerDict[goTerm][entrezID] = dict()
						if not datasetID in masterDict[goTerm]:
							masterDict[goTerm][datasetID] = dict()
						masterDict[goTerm][datasetID][entrezID] = [goDict[goKey]["Ontology"], goDict[goKey]["go.term"], goDict[goKey]["count.InDataset"], goDict[goKey]["count.InGenome"], goDict[goKey]["BH.adjusted.p.value"]]
						scorerDict[goTerm][entrezID][datasetID] = float(goDict[goKey]["BH.adjusted.p.value"])
						entrezDict[goTerm][entrezID] = len(scorerDict[goTerm][entrezID])
						
		#print masterDict
		#print scorerDict
		#print entrezDict
		print "Printing candidate entrezID hits..."
		print
		for goTerm in sorted(scorerDict.keys()):
			entrezSort = general.valuesort(entrezDict[goTerm])
			entrezSort.reverse()
			entrezID = entrezSort[0]
			print goTerm, ":", entrezID, "(" + str(entrezDict[goTerm][entrezID]) + ")"
			for datasetID in scorerDict[goTerm][entrezID]:
				print datasetID, scorerDict[goTerm][entrezID][datasetID]
			print
		print
	
if __name__ == "__main__":
	main()
	
	print "Completed:", time.asctime(time.localtime())

#python mapGO.py --path ~/meTRN --peaks hs_selection_com_cx_hct --mode factor.context --minPvalue 0.01 --hitCount 1 --hitPvalue 0.01
#python mapGO.py --path ~/meTRN --peaks optimal_standard_final --mode factor.context --minPvalue 0.01 --hitCount 2 --hitPvalue 0.01
#python mapGO.py --path ~/meTRN --peaks optimal_standard_final --mode factor.context --minPvalue 0.01 --hitCount 3 --hitPvalue 0.01
#python mapGO.py --path ~/meTRN --peaks optimal_standard_final --mode factor.context --minPvalue 0.01 --hitCount 5 --hitPvalue 0.01


#python mapGO.py --path ~/meTRN --organism ce --mode sunburst --peaks any.ex.som,any.l1.som,any.l2.som,any.l3.som --analysis p5e-1 --target som.neurons --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05 --minCount 5 --maxCount 35 --maxColor 6