#!/usr/bin/env python
# SOM-related pre- and post-processing!

import sys
import time
import optparse
import general
import hyper
import numpy
import pickle
import pdb
import metrn
import modencode
import scipy
import random
import os

from scipy import stats
from scipy.stats.stats import pearsonr
from quantile import Quantile

print "Command:", " ".join(sys.argv)
print "Timestamp:", time.asctime(time.localtime())


""" define variables of internal use """
neuronMapping = {
	"any.ex.som"  : "mapneurons_matrix_context_any_ORC_FAC_<<PEAKS>>.txt",
	"any.l1.som"  : "mapneurons_matrix_context_any_ORC_FAC_<<PEAKS>>.txt",
	"any.l2.som"  : "mapneurons_matrix_context_any_ORC_FAC_<<PEAKS>>.txt",
	"any.l3.som"  : "mapneurons_matrix_context_any_ORC_FAC_<<PEAKS>>.txt",
	"any.l4.som"  : "mapneurons_matrix_context_any_ORC_FAC_<<PEAKS>>.txt",
	"all.elc.som" : "mapneurons_matrix_context_all_ELC_FAC_<<PEAKS>>.txt",
	"all.e1c.som" : "mapneurons_matrix_context_all_E1C_FAC_<<PEAKS>>.txt",
	"all.e2c.som" : "mapneurons_matrix_context_all_E2C_FAC_<<PEAKS>>.txt",
	"all.e3c.som" : "mapneurons_matrix_context_all_E3C_FAC_<<PEAKS>>.txt",
	"all.e4c.som" : "mapneurons_matrix_context_all_E4C_FAC_<<PEAKS>>.txt",
	"all.1v2.som" : "mapneurons_matrix_context_all_1V2_FAC_<<PEAKS>>.txt",
	"all.2v3.som" : "mapneurons_matrix_context_all_2V3_FAC_<<PEAKS>>.txt",
	"all.3v4.som" : "mapneurons_matrix_context_all_3V4_FAC_<<PEAKS>>.txt",
	"cell.ex.som" : "mapneurons_matrix_cells_all_xp1_FAC_<<PEAKS>>.txt"
	}


""" define functions of internal use """

def uppify(indict):
	output = dict()
	for key in indict:
		newkey = key.upper()
		output[newkey] = indict[key]
	return output


""" define a function to filter factors that are compliant, adjust and export signal matrix """
def exportMatrix(factors, regions, matrix, coverage, technique, outfile, filter="ON"):

	# adjust/scale signal matrix for 'signal' and 'rank' techniques:
	print
	print "Adjusting signal matrix..."
	if technique in ["signal", "normal", "rank"]:
		for factor in matrix:
			for region in matrix[factor]:
				if technique in ["signal", "normal"]:
					matrix[factor][region] = max(matrix[factor][region])
				elif technique == "rank":
					matrix[factor][region] = min(matrix[factor][region])
	
	# prepare output matrix values: 
	print "Preparing", technique, "matrix..."
	signal = dict()
	for factor in factors:
		values = matrix[factor].values()
		maximum = max(values)
		signal[factor] = dict()
		for region in regions:
			if region in matrix[factor]:
				if technique == "binary":
					value = 1
				elif technique in ["signal", "normal"]:
					value = float(matrix[factor][region])/maximum
				elif technique == "rank":
					value = float(maximum - matrix[factor][region] + 1)/maximum
			else:
				value = 0
			signal[factor][region] = value
	
	# sum signal per region for region normalization:
	if technique in ["normal"]:
		print "Normalizing region signals..."
		
		regionSums = dict()
		for factor in signal:
			for region in signal[factor]:
				if not region in regionSums:
					regionSums[region] = 0
				regionSums[region] += signal[factor][region]
		
		for factor in signal:
			for region in signal[factor]:
				if regionSums[region] == 0 and float(signal[factor][region]) == 0:
					signal[factor][region] = 0
				elif regionSums[region] != 0:
					signal[factor][region] = float(signal[factor][region])/regionSums[region]
				else:
					print "Error: Zero signal region with signal in factors!"
					pdb.set_trace()
				
	# export matrix data:
	print "Exporting matrix..."
	f_output = open(outfile, "w")
	print >>f_output, "\t".join([""] + factors)
	for region in regions:
		output = [region]
		for factor in factors:
			output.append(signal[factor][region])
		if filter == "OFF" or output.count(0) != len(output)-1:
			print >>f_output, "\t".join(map(str, output))
	f_output.close()
	
	print


""" define a function to filter factors that are compliant, adjust and export signal matrix """
def exportFaster(factors, regions, matrix, coverage, technique, outfile, filter="ON"):

	# check we're running in binary mode:
	if technique != "binary":
		raise "Error: This export mode can only be run for binary matrixes!"
	
	# invert the factors/region matrix:
	print "Inverting matrix for fast processing..."
	fastMatrix = dict()
	for factor in matrix:
		for region in matrix[factor]:
			if not region in fastMatrix:
				fastMatrix[region] = dict()
			fastMatrix[region][factor] = 1
	
	# prepare output matrix values: 
	print
	print "Exporting", technique, "matrix..."
	f_output = open(outfile, "w")
	print >>f_output, "\t".join([""] + factors)
	for region in regions:
		output = [region]
		for factor in factors:
			if factor in fastMatrix[region]:
				output.append(1)
			else:
				output.append(0)
		if filter == "OFF" or output.count(0) != len(output)-1:
			print >>f_output, "\t".join(map(str, output))
	f_output.close()
	
	print
	
	
""" define a function to load a binding/region matrix """
def loadMatrix(infile, coords="OFF", signatures="OFF", unique="ON"):
	factorMatrix, moduleMatrix, coordsMatrix, signatureMatrix = dict(), dict(), dict(), dict()
	indata = open(infile)
	factors = general.clean(indata.readline().strip().split("\t"))
	inline = indata.readline()
	codes = ",".join(factors)
	k = 0
	while inline:
		values = inline.strip().split("\t")
		module = values.pop(0)
		signature = ",".join(values)
		coords = module.split(":")
		coords.pop(0)
		coords = ":".join(coords)
		
		if not module in moduleMatrix:
			moduleMatrix[module] = dict()
		elif unique == "ON":
			raise "Error: Region already processed!"
		
		if coords == "ON":
			if not coords in coordsMatrix:
				coordsMatrix[coords] = list()
			coordsMatrix[coords].append(module)
		
		if signatures == "ON":
			if not signature in signatureMatrix:
				signatureMatrix[signature] = list()
			signatureMatrix[signature].append(module)
		
		index = 0
		for factor in factors:
			if not factor in factorMatrix:
				factorMatrix[factor] = dict()
			factorMatrix[factor][module] = values[index]
			moduleMatrix[module][factor] = values[index]
			index += 1
		inline = indata.readline()
		k += 1
		
	return factorMatrix, moduleMatrix, coordsMatrix, signatureMatrix


""" define a function to load a binding/region matrix (faster) """
def loadFaster(infile):
	factorMatrix, moduleMatrix, factorIndexes, moduleIndexes = dict(), dict(), dict(), dict()
	indata = open(infile)
	factors = general.clean(indata.readline().strip().split("\t"))
	index = 0
	for factor in factors:
		factorIndexes[factor] = index
		index += 1
	inline = indata.readline()
	index = 0
	while inline:
		values = inline.strip().split("\t")
		module = values.pop(0)
		for factor in factors:
			if not factor in factorMatrix:
				factorMatrix[factor] = list()
			factorMatrix[factor].append(module)
		moduleIndexes[module] = index
		moduleMatrix[module] = values
		inline = indata.readline()
		index += 1
	return factors, factorMatrix, moduleMatrix, factorIndexes, moduleIndexes


""" define a function to map regions to neurons """
def mapRegions(rawpath, regions=list()):
	neuronDict = dict()
	for neuron in os.listdir(rawpath):
		for inline in open(rawpath + neuron).readlines():
			neuron = neuron.split(".")[0]
			region = inline.strip().replace(" ", ":")
			if regions==list() or region in regions:
				if not neuron in neuronDict:
					neuronDict[neuron] = list()
				neuronDict[neuron].append(region)
	return neuronDict

	
def main():
	
	parser = optparse.OptionParser()
	parser.add_option("--path", action = "store", type = "string", dest = "path", help = "path from script to files")
	parser.add_option("--organism", action = "store", type = "string", dest = "organism", help = "Target organism for operations...", default="OFF")
	parser.add_option("--mode", action = "store", type = "string", dest = "mode", help = "type of operations to be performed: scan or filter")
	parser.add_option("--peaks", action = "store", type = "string", dest = "peaks", help = "basename for target peaks", default="OFF")
	parser.add_option("--infile", action = "store", type = "string", dest = "infile", help = "Input file for analysis", default="OFF")
	parser.add_option("--nuclear", action = "store", type = "string", dest = "nuclear", help = "Peaks are only nuclear?", default="ON")
	parser.add_option("--name", action = "store", type = "string", dest = "name", help = "Output file name-tag.", default="OFF")
	parser.add_option("--technique", action = "store", type = "string", dest = "technique", help = "What kind of matrix should I build? Binary, rank, or signal", default="binary")
	parser.add_option("--contexts", action = "store", type = "string", dest = "contexts", help = "What contexts of development should I track?", default="OFF")
	parser.add_option("--target", action = "store", type = "string", dest = "target", help = "What labels should be targeted?", default="dataset")
	parser.add_option("--coverage", action = "store", type = "string", dest = "coverage", help = "Filter factors for coverage of all target contexts?", default="all")
	parser.add_option("--sample", action = "store", type = "string", dest = "sample", help = "Sample comparable numbers of peaks per factor in comparisons?", default="OFF")
	parser.add_option("--filter", action = "store", type = "string", dest = "filter", help = "Cumulatively filter factors from matrix comparisons?", default="OFF")
	parser.add_option("--fraction", action = "store", type = "float", dest = "fraction", help = "Path to annotation file", default=0.75)
	parser.add_option("--remove", action = "store", type = "string", dest = "remove", help = "Individually remove factors from matrix comparisons?", default="OFF")
	parser.add_option("--expression", action = "store", type = "float", dest = "expression", help = "Fractional expression cutoff", default=0.1)
	parser.add_option("--minimum", action="store", type="float", dest="minimum", help="Minimum raw expression cutoff", default=2000)
	#parser.add_option("--selection", action = "store", type = "string", dest = "selection", help = "Report 'avg' or 'max' expression values?", default="avg")
	#parser.add_option("--limit", action="store", type="string", dest="limit", help="Distance limit for binding to TSS assignments...", default="OFF")
	parser.add_option("--neurons", action = "store", type = "string", dest = "neurons", help = "Name of SOM neurons/results...", default="OFF")
	parser.add_option("--regions", action = "store", type = "string", dest = "regions", help = "Regions to be examined...", default="OFF")
	parser.add_option("--include", action = "store", type = "string", dest = "include", help = "Forcefully include factors...", default="OFF")
	parser.add_option("--group", action = "store", type = "string", dest = "group", help = "Should cell groups be constructed from expression genes (expression.genes) or from those peak binding genes (peak.genes)?", default="peak.genes")
	parser.add_option("--source", action = "store", type = "string", dest = "source", help = "Path to annotation file", default="OFF")
	parser.add_option("--input", action = "store", type = "string", dest = "input", help = "Input annotation file", default="OFF")
	parser.add_option("--header", action = "store", type = "string", dest = "header", help = "Is there a header in input file?", default="OFF")
	parser.add_option("--id", action = "store", type = "string", dest = "id", help = "ID column/header", default="feature")
	parser.add_option("--tss", action = "store", type = "string", dest = "tss", help = "TSS reference file", default="OFF")
	parser.add_option("--query", action="store", type="string", dest="query", help="Query collections of cells whose enrichment will be searched in target cells", default="OFF")
	parser.add_option("--rename", action = "store", type = "string", dest = "rename", help = "Targets to rename. Comma-separated list of 'target:replacement' pairs to search and replace.", default="OFF")
	parser.add_option("--indexes", action = "store", type = "string", dest = "indexes", help = "values", default="feature:3,score:4")
	parser.add_option("--datatype", action = "store", type = "string", dest = "datatype", help = "values", default="float")
	parser.add_option("--overwrite", action="store", type="string", dest="overwrite", help="Overwrite outputs?", default="OFF")
	parser.add_option("--parameters", action = "store", type = "string", dest = "parameters", help = "Variable parameters...", default="")
	parser.add_option("--fast", action = "store", type = "string", dest = "fast", help = "Use fast loading?", default="ON")
	parser.add_option("--threads", action = "store", type = "int", dest = "threads", help = "Parallel processing threads", default=1)
	parser.add_option("--chunks", action = "store", type = "int", dest = "chunks", help = "", default=100)
	parser.add_option("--module", action = "store", type = "string", dest = "module", help = "", default="md1")
	parser.add_option("--qsub", action = "store", type = "string", dest = "qsub", help = "Qsub configuration header", default="OFF")
	parser.add_option("--server", action = "store", type = "string", dest = "server", help = "Are we on the server?", default="OFF")
	parser.add_option("--job", action = "store", type = "string", dest = "job", help = "Job name for cluster", default="OFF")
	parser.add_option("--copy", action = "store", type = "string", dest = "copy", help = "Copy simulated peaks to analysis folder?", default="OFF")
	parser.add_option("--tag", action = "store", type = "string", dest = "tag", help = "Add tag to TFBS?", default="")
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
	setuppath = path_dict["setup"]
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
	bindingpath = path_dict["binding"]
	neuronspath = path_dict["neurons"]
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
	
	# specify genome size file:
	if option.nuclear == "ON":
		chromosomes = metrn.chromosomes[organismTag]["nuclear"]
		genome_size_file = option.path + "/input/" + metrn.reference[organismTag]["nuclear_sizes"]
		genome_size_dict = general.build_config(genome_size_file, mode="single", separator="\t", spaceReplace=True)
	else:
		chromosomes = metrn.chromosomes[organismTag]["complete"]
		genome_size_file = option.path + "/input/" + metrn.reference[organismTag]["complete_sizes"]
		genome_size_dict = general.build_config(genome_size_file, mode="single", separator="\t", spaceReplace=True)
		
	# load gene ID dictionaries:
	id2name_dict, name2id_dict = modencode.idBuild(inpath + metrn.reference[organismTag]["gene_ids"], "Sequence Name (Gene)", "Gene Public Name", mode="label", header=True, idUpper=True, nameUpper=True)
	
	# define target/node type:
	if option.target == "organism":
		codeLabel = "ORG"
	elif option.target == "cellular":
		codeLabel = "CEL"
	else:
		codeLabel = metrn.labelGenerator(option.target, mode="handle")
	
	# prepare relevant paths:
	if option.mode != "neuron.compile":
		neuronspath = neuronspath + option.peaks + "/" + option.technique + "/"
		matrixpath = neuronspath + "matrix/"
		modulespath = neuronspath + "modules/"
		numberspath = neuronspath + "numbers/"
		resultspath = neuronspath + "results/"
		reportspath = neuronspath + "reports/"
		comparepath = neuronspath + "compare/"
		pearsonpath = neuronspath + "pearson/"
		general.pathGenerator(neuronspath)
		general.pathGenerator(matrixpath)
		general.pathGenerator(numberspath)
		general.pathGenerator(modulespath)
		general.pathGenerator(resultspath)
		general.pathGenerator(reportspath)
		general.pathGenerator(comparepath)
		general.pathGenerator(pearsonpath)
	else:
		compilepath = neuronspath + "compile/" + option.name + "/"
		general.pathGenerator(compilepath)
	
	# context matrix-construction mode:
	if option.mode == "matrix.context":
		
		print
		
		# load target contexts:
		codeContexts, targetContexts = metrn.options_dict["contexts.extended"][option.contexts]
		
		# define input files:
		completefile = peakspath + "mappeaks_" + option.peaks + "_complete.bed"
		compiledfile = peakspath + "mappeaks_" + option.peaks + "_compiled.bed"
		
		# prebuild peak dictionary:
		print "Preloading peak data..."
		peak_dict, ranks = dict(), list()
		inlines = open(completefile).readlines()
		for inline in inlines:
			chrm, start, end, peak, score, strand, signal, pvalue, qvalue, point = inline.strip().split("\t")[:10]
			dataset, rank = peak.split(":")
			dataset = dataset.replace("POL2", "AMA-1")
			rank = int(rank.split("P")[1])
			signal = float(signal)
			if option.technique == "rank":
				peak_dict[peak] = rank
			else:
				peak_dict[peak] = signal
			ranks.append(rank)
		
		# adjust peak ranks (in case you are working with HOT-filtered sets):
		if option.technique == "rank":
			minRank = min(ranks)
			for peak in peak_dict:
				peak_dict[peak] = peak_dict[peak] - minRank + 1
		
		# prebuild factors, regions, and matrix:
		print "Preloading factors and regions..."
		coverage, factors, regions, matrix, counts, places, labels = dict(), list(), list(), dict(), dict(), dict(), dict()
		contexts_processed = list()
		inlines = open(compiledfile).readlines()
		for inline in inlines:
			chrm, start, end, region, occupancy, strand, density, complexity, factorCount, contextCount, infos = inline.strip().split("\t")
			region = ":".join([chrm, start, end])
			peaks = infos.split(";")
			peaks.pop(0)
			for peak in peaks:
				dataset = peak.split(":")[0].replace("POL2", "AMA-1")
				organism, strain, factor, context, institute = metrn.labelComponents(dataset, target="subcomponents")
				process, context = metrn.resolve(context, option.contexts)
				if process:
					
					# record processed contexts:
					contexts_processed.append(context)
					
					# specialized factor and region labeling:
					if option.target == "organism":
						factor_label = factor
						region_label = organismTag.upper() + ":" + region
					elif option.target == "factor":
						factor_label = factor
						region_label = context + ":" + region
					elif option.target == "factor.context":
						factor_label = ".".join([factor, context])
						region_label = organismTag.upper() + ":" + region
					elif option.target == "factor.institute":
						factor_label = ".".join([factor, institute])
						region_label = organismTag.upper() + ":" + region
					elif option.target == "factor.context.institute":
						factor_label = ".".join([factor, context, institute])
						region_label = organismTag.upper() + ":" + region
					elif option.target == "dataset":
						factor_label = dataset
						region_label = organismTag.upper() + ":" + region
					
					# attach processed peak/region info:
					factors.append(factor_label)
					regions.append(region_label)
					if not factor_label in coverage:
						coverage[factor_label] = list()
					if not context in coverage[factor_label]:
						coverage[factor_label].append(context)
					if not factor_label in matrix:
						matrix[factor_label] = dict()
					if not region_label in matrix[factor_label]:
						matrix[factor_label][region_label] = list()
					matrix[factor_label][region_label].append(peak_dict[peak])
					
					# store counts per factor/context, factor-labels associations, and :
					if option.sample != "OFF" or option.filter != "OFF" or option.remove != "OFF":
						
						if not factor in counts:
							counts[factor] = dict()
							places[factor] = dict()
						if not context in counts[factor]:
							counts[factor][context] = 0
							places[factor][context] = list()
						counts[factor][context] += 1
						places[factor][context].append(region_label)
						labels[factor_label] = factor
					
		# prefilter complaint factors:
		print "Filtering compliant factors..."
		compliant = list()
		for factor in factors:
			if option.coverage == "any" and set(coverage[factor]).intersection(set(targetContexts)):
				compliant.append(factor)
			elif sorted(coverage[factor]) == sorted(targetContexts):
				compliant.append(factor)
		
		factors = sorted(list(set(compliant)))
		regions = sorted(list(set(regions)))
		print "Found", len(targetContexts), "target contexts..."
		print "Found", len(factors), "compliant factors..."
		
		# export observed matrix:
		if option.sample == "OFF" and option.filter == "OFF" and option.remove == "OFF":
		
			# define output file:
			f_outfile = matrixpath + "mapneurons_matrix_context_" + option.coverage + "_" + codeContexts + "_" + codeLabel + "_" + option.peaks + ".txt"
		
			# adjust signals and export matrix:
			exportMatrix(factors, regions, matrix, coverage, technique=option.technique, outfile=f_outfile)
		
		# export sub-sampled matrix(es):
		if option.sample != "OFF" and option.filter == "OFF" and option.remove == "OFF":
		
			# run through randomizations (iterations):
			print
			outfiles = list()
			for index in range(1, int(option.sample) + 1):
				
				print "Executing randomization:", index
			
				# subsample (without replacement) to match peak counts between factor/context (factor_labels) for the same factor:
				subset = dict()
				for factor_label in matrix:
					
					# determine factor and contexts:
					factor = labels[factor_label]
					contexts = counts[factor].keys()
					
					# determine limit (for sampling):
					limit = min(counts[factor].values())
					
					# subsample peaks for each context: 
					for context in contexts:
						for region_label in random.sample(places[factor][context], limit):
							if not factor_label in subset:
								subset[factor_label] = dict()
							subset[factor_label][region_label] = matrix[factor_label][region_label]
				
				# generate index tag:
				indexTag = general.indexTag(index, int(option.sample))
			
				# define output file:
				f_outfile = matrixpath + "mapneurons_sample_context_" + option.coverage + "_" + codeContexts + "_" + codeLabel + "_" + option.peaks + "." + indexTag
				outfiles.append(f_outfile)
				
				# adjust signals and export matrix:
				exportMatrix(factors, regions, subset, coverage, technique=option.technique, outfile=f_outfile)
		
			# define merged output file:
			m_outfile = matrixpath + "mapneurons_sample_context_" + option.coverage + "_" + codeContexts + "_" + codeLabel + "_" + option.peaks + ".sum"
			m_output = open(m_outfile, "w")
			
			# generate merged output file:
			print "Generating merged matrix..."
			start = True
			for outfile in outfiles:
				outdata = open(outfile)
				outline = outdata.readline().rstrip()
				if start:
					print >>m_output, outline
					start = False
				outline = outdata.readline().rstrip()
				while outline:
					print >>m_output, outline
					outline = outdata.readline().rstrip()
			
			# close merged output file:
			m_output.close()
		
		# export factor-filtered matrix(es):
		if option.sample == "OFF" and (option.filter != "OFF" or option.remove != "OFF"):
		
			# score factors for differences between contexts:
			fractional, difference = dict(), dict()
			for factor_label in sorted(matrix.keys()):
				
				# determine factor and contexts:
				factor = labels[factor_label]
				contexts = counts[factor].keys()
				
				# determine difference:
				if len(contexts) == 2:
					difference[factor] = abs(counts[factor][contexts[0]] - counts[factor][contexts[1]])
					fractional[factor] = float(abs(counts[factor][contexts[0]] - counts[factor][contexts[1]]))/max(counts[factor].values())
					masterContexts = sorted(list(contexts))
					
			# sort factors by the difference in peak counts:
			targets = general.valuesort(fractional)
			targets.reverse()
			
			# generate comparison report file:
			c_outfile = comparepath + "mapneurons_matrix_compare_" + option.coverage + "_" + codeContexts + "_" + codeLabel + "_" + option.peaks + ".txt"
			c_output = open(c_outfile, "w")
			print >>c_output, "\t".join(["rank", "factor", "index", "fraction", "difference"] + masterContexts)
			k = 1
			for target in targets:
				indexTag = general.indexTag(k, len(targets), minimum=3)
				output = [k, target, indexTag, fractional[target], difference[target]]
				for context in masterContexts:
					output.append(counts[target][context])
				print >>c_output, "\t".join(map(str, output))
				k += 1
			c_output.close()
			
			# determine iterations:
			if option.filter != "OFF":
				iterations = int(option.filter)
			if option.remove != "OFF":
				if option.remove != "ALL":
					iterations = int(option.remove)
				else:
					iterations = len(targets)
			
			# selectively remove factors from matrixes (iterations):
			print
			print "Targets for removal:", ", ".join(targets[:iterations])
			for index in range(1, iterations + 1):
				
				# generate index tag:
				indexTag = general.indexTag(index, iterations, minimum=3)
				
				# removal factors:
				if option.filter != "OFF":
					removal = list(targets[:index])
					print "Filtering factors:", index
					print "Excluding:", ", ".join(removal)
				if option.remove != "OFF":
					removal = [targets[index-1]]
					print "Removing factor:", index
					print "Excluding:", removal[0], difference[removal[0]], fractional[removal[0]]
				
				# generate new matrix (filtering removal factors):
				subset = dict()
				for factor_label in matrix:
					
					# determine factor and contexts:
					factor = labels[factor_label]
					
					# append factors (that are not targeted for removal):
					if not factor in removal:
						if not factor_label in subset:
							subset[factor_label] = dict()
							for region_label in matrix[factor_label]:
								subset[factor_label][region_label] = matrix[factor_label][region_label]
				
				# define output file:
				if option.filter != "OFF":
					f_outfile = matrixpath + "mapneurons_filter_context_" + option.coverage + "_" + codeContexts + "_" + codeLabel + "_" + option.peaks + "." + indexTag
				elif option.remove != "OFF":
					f_outfile = matrixpath + "mapneurons_remove_context_" + option.coverage + "_" + codeContexts + "_" + codeLabel + "_" + option.peaks + "." + indexTag
				
				# re-make factors list:
				subfactors = list(set(factors).difference(set(removal))) 
				
				# adjust signals and export matrix:
				exportMatrix(subfactors, regions, subset, coverage, technique=option.technique, outfile=f_outfile)
	
	
	# select matrix mode:
	elif option.mode == "matrix.select":
	
		# load target contexts:
		codeContexts, targetContexts = metrn.options_dict["contexts.extended"][option.contexts]
		
		# define matrix basename:
		f_observed = "mapneurons_matrix_context_" + option.coverage + "_" + codeContexts + "_" + codeLabel + "_" + option.peaks + ".txt"
		f_combined = "mapneurons_sample_context_" + option.coverage + "_" + codeContexts + "_" + codeLabel + "_" + option.peaks + ".sum"
		f_basename = "mapneurons_sample_context_" + option.coverage + "_" + codeContexts + "_" + codeLabel + "_" + option.peaks
		
		# remove previous selections:
		command = "rm -rf " + matrixpath + f_basename + ".max"
		os.system(command)
		
		print 
		# load observed signatures:
		factorMatrix, moduleMatrix, coordsMatrix, observedSignatures = loadMatrix(matrixpath + f_observed, coords="OFF", signatures="ON", unique="ON")
		print "Loaded observed signatures:", len(observedSignatures)
		
		# load combined signatures (from sampling analysis):
		factorMatrix, moduleMatrix, coordsMatrix, combinedSignatures = loadMatrix(matrixpath + f_combined, coords="OFF", signatures="ON", unique="OFF")
		print "Loaded combined signatures:", len(combinedSignatures)
		
		# find common signatures:
		commonSignatures = set(observedSignatures.keys()).intersection(set(combinedSignatures.keys()))
		print "Found common signatures:", len(commonSignatures)
		
		# define output file:
		f_outfile = pearsonpath + "mapneurons_pearson_context_" + option.coverage + "_" + codeContexts + "_" + codeLabel + "_" + option.peaks + ".txt"
		f_output = open(f_outfile, "w")
		print >>f_output, "\t".join(["index", "index.signatures", "observed.signatures", "observed.fraction", "observed.correlation", "combined.correlation", "observed.counts.avg", "orphaned.counts.avg", "observed.orphaned.ratio"])
		
		# load matrix files:
		f_infiles = os.listdir(matrixpath)
		
		# examine correlations:
		print
		print "Evaluating sampled signatures..."
		maxValue = 0
		for f_infile in f_infiles:
			if f_basename in f_infile and f_infile != f_combined:
				index = f_infile.replace(f_basename + "." , "")
				#print "Processing:", index
				factorMatrix, moduleMatrix, coordsMatrix, evaluateSignatures = loadMatrix(matrixpath + f_infile, coords="OFF", signatures="ON", unique="ON")
				
				observedVector, evaluateVector = list(), list()
				observedOverlap = set(observedSignatures.keys()).intersection(set(evaluateSignatures.keys()))
				observedFraction = float(len(observedOverlap))/len(evaluateSignatures)
				for signature in observedOverlap:
					observedVector.append(len(observedSignatures[signature]))
					evaluateVector.append(len(evaluateSignatures[signature]))
				observedCorrelation, observedPvalue = pearsonr(observedVector, evaluateVector)
				
				combinedVector, evaluateVector = list(), list()
				combinedOverlap = set(combinedSignatures.keys()).intersection(set(evaluateSignatures.keys()))
				combinedFraction = float(len(combinedOverlap))/len(evaluateSignatures)
				for signature in combinedOverlap:
					combinedVector.append(len(combinedSignatures[signature]))
					evaluateVector.append(len(evaluateSignatures[signature]))
				combinedCorrelation, combinedPvalue = pearsonr(combinedVector, evaluateVector)
				
				observedCounts = list()
				orphanedCounts = list()
				for signature in evaluateSignatures:
					if signature in observedSignatures:
						observedCounts.append(len(evaluateSignatures[signature]))
					else:
						orphanedCounts.append(len(evaluateSignatures[signature]))
				observedRatio = float(numpy.mean(observedCounts))/numpy.mean(orphanedCounts)
				print >>f_output, "\t".join(map(str, [index, len(evaluateSignatures), len(observedOverlap), observedFraction, observedCorrelation, combinedCorrelation, numpy.mean(observedCounts), numpy.mean(orphanedCounts), observedRatio]))
				if combinedCorrelation > maxValue:
					maxName, maxValue, maxSignatures = f_infile, float(combinedCorrelation), evaluateSignatures
		
		# export selected sample with highest correlation:
		print
		print "Selection:", maxName
		command = "cp " + matrixpath + maxName + " " + matrixpath + f_basename + ".max"
		os.system(command)
		print
	
	# cell matrix construction mode:
	elif option.mode == "matrix.cells":
	
		# define expression cutoff handle:
		codeExpression = "xp" + str(option.expression).split(".")[1]
		
		# load input expression path:
		expressionpath = cellspath + "expression/"
		
		# define input expression file:
		expressionfile = expressionpath + option.infile
		
		# define inclusion factors:
		if option.include != "OFF":
			inclusionGenes = ",".split(option.include)
		else:
			inclusionGenes = list()
			
		# load cellular expression matrix:
		expression_matrix = dict()
		hd = general.build_header_dict(expressionfile)
		tissueDict = dict()
		inlines = open(expressionfile).readlines()
		inlines.pop(0)
		for inline in inlines:
			initems = inline.strip().split("\t")
			gene, cell, cellExpression, fractionExpression, specificTissue, generalTissue, classTissue = initems[hd["gene"]], initems[hd["cell"]], float(initems[hd["cell.expression"]]), float(initems[hd["fraction.expression"]]), initems[hd["specific.tissue"]], initems[hd["general.tissue"]], initems[hd["class.tissue"]]
			if not cell in tissueDict:
				tissueDict[cell] = [classTissue, generalTissue, specificTissue]
			if not gene in expression_matrix:
				expression_matrix[gene] = dict()
			if (fractionExpression >= option.expression and cellExpression >= option.minimum)or gene in inclusionGenes:
				expression_matrix[gene][cell] = fractionExpression
		
		# define input files:
		completefile = peakspath + "mappeaks_" + option.peaks + "_complete.bed"
		compiledfile = peakspath + "mappeaks_" + option.peaks + "_compiled.bed"
		
		# define output file:
		f_outfile = matrixpath + "mapneurons_matrix_cells_" + option.name + "_" + codeExpression + "_" + codeLabel + "_" + option.peaks + ".txt"
		
		# prebuild peak dictionary:
		print
		print "Preloading peak data..."
		peak_dict, ranks = dict(), list()
		inlines = open(completefile).readlines()
		for inline in inlines:
			chrm, start, end, peak, score, strand, signal, pvalue, qvalue, point = inline.strip().split("\t")[:10]
			dataset, rank = peak.split(":")
			dataset = dataset.replace("POL2", "AMA-1")
			rank = int(rank.split("P")[1])
			signal = float(signal)
			if option.technique == "rank":
				peak_dict[peak] = rank
			else:
				peak_dict[peak] = signal
			ranks.append(rank)
		
		# adjust peak ranks (in case you are working with HOT-filtered sets):
		if option.technique == "rank":
			minRank = min(ranks)
			for peak in peak_dict:
				peak_dict[peak] = peak_dict[peak] - minRank + 1
				
		#if option.technique == "binary":
		#	f_output = open(outfile, "w")
		#	print >>f_output, "\t".join([""] + factors)
		
		# prebuild factors, regions, and matrix:
		print "Preloading factors and regions..."
		k, coverage, factors, regions, matrix = 0, dict(), list(), list(), dict()
		contexts_processed = list()
		inlines = open(compiledfile).readlines()
		for inline in inlines:
			genes, cells = list(), list()
			chrm, start, end, region, occupancy, strand, density, complexity, factorCount, contextCount, infos = inline.strip().split("\t")
			region = ":".join([chrm, start, end])
			altregion = ".".join([chrm, start, end])
			peaks = infos.split(";")
			peaks.pop(0)
			for peak in peaks:
				dataset = peak.split(":")[0].replace("POL2", "AMA-1")
				organism, strain, factor, context, institute = metrn.labelComponents(dataset, target="subcomponents")
				process, context = metrn.resolve(context, option.contexts)
				if factor in expression_matrix:
					genes.append(factor)
					for cell in expression_matrix[factor]:
						cells.append(cell)
			genes, cells = sorted(list(set(genes))), sorted(list(set(cells)))
			
			k += 1
			#print k, region, len(genes), len(cells)
			
			for cell in cells:
				
				# generate region label:
				if option.target == "organism":
					region_label = organismTag.upper() + ":" + region
				elif option.target == "factor":
					region_label = cell + ":" + region
				elif option.target == "factor.context":
					region_label = cell + ":" + region
				elif option.target == "factor.institute":
					region_label = cell + ":" + region
				elif option.target == "factor.context.institute":
					region_label = cell + ":" + region
				elif option.target == "dataset":
					region_label = cell + ":" + region
				elif option.target == "cellular":
					region_label = ":".join([classTissue, specificTissue, cell, altregion])
			
				"""
				# export for binaries!
				f_output = open(outfile, "w")
				print >>f_output, "\t".join([""] + factors)
				
				if option.technique == "binary":
					output = [region]
					
					for factor in factors:
						if region in matrix[factor]:
							value = 1
						else:
							value = 0
					output.append(value)
				
					if filter == "OFF" or output.count(0) != len(output)-1:
						print >>f_output, "\t".join(map(str, output))
				"""
				
				for peak in peaks:
					dataset = peak.split(":")[0].replace("POL2", "AMA-1")
					organism, strain, factor, context, institute = metrn.labelComponents(dataset, target="subcomponents")
					if factor in genes:
						if cell in expression_matrix[factor]:
							
							classTissue, generalTissue, specificTissue = tissueDict[cell]
							specificTissue = specificTissue.replace(" ", "_")
							
							# generate factor label:
							if option.target == "organism":
								factor_label = factor
							elif option.target == "factor":
								factor_label = factor
							elif option.target == "factor.context":
								factor_label = ".".join([factor, context])
							elif option.target == "factor.institute":
								factor_label = ".".join([factor, institute])
							elif option.target == "factor.context.institute":
								factor_label = ".".join([factor, context, institute])
							elif option.target == "dataset":
								factor_label = dataset
							elif option.target == "cellular":
								factor_label = factor
							
							factors.append(factor_label)
							regions.append(region_label)
							if not factor_label in coverage:
								coverage[factor_label] = list()
							if not cell in coverage[factor_label]:
								coverage[factor_label].append(cell)
							if not factor_label in matrix:
								matrix[factor_label] = dict()
							if not region_label in matrix[factor_label]:
								matrix[factor_label][region_label] = list()
							matrix[factor_label][region_label].append(peak_dict[peak])
		
		# simplify factor and region lists:
		factors = sorted(list(set(factors)))
		regions = sorted(list(set(regions)))
		
		# adjust signals and export matrix:
		if option.technique == "binary":
			exportFaster(factors, regions, matrix, coverage, technique=option.technique, outfile=f_outfile)
		else:
			exportMatrix(factors, regions, matrix, coverage, technique=option.technique, outfile=f_outfile)

	
	# expression group matrix construction mode:
	elif option.mode == "matrix.group":
	
		# define group handle:
		if option.group == "peak.genes":
			ghandle = "pk"
		elif option.group == "expression.genes":
			ghandle = "ex"
		
		# define inclusion factors:
		if option.include != "OFF":
			inclusionGenes = ",".split(option.include)
		else:
			inclusionGenes = list()
			
		# define group/expression cutoff handle:
		codeExpression = ghandle + str(option.expression).split(".")[1]
		
		# load input expression path:
		expressionpath = cellspath + "expression/"
		
		# define input expression file:
		expressionfile = expressionpath + option.infile
		
		# define input files:
		completefile = peakspath + "mappeaks_" + option.peaks + "_complete.bed"
		compiledfile = peakspath + "mappeaks_" + option.peaks + "_compiled.bed"
		
		# define output file:
		f_outfile = matrixpath + "mapneurons_matrix_group_" + option.name + "_" + codeExpression + "_" + codeLabel + "_" + option.peaks + ".txt"
		
		# prebuild peak dictionary:
		print
		print "Preloading peak data..."
		peak_dict, ranks = dict(), list()
		inlines = open(completefile).readlines()
		for inline in inlines:
			chrm, start, end, peak, score, strand, signal, pvalue, qvalue, point = inline.strip().split("\t")[:10]
			dataset, rank = peak.split(":")
			dataset = dataset.replace("POL2", "AMA-1")
			rank = int(rank.split("P")[1])
			signal = float(signal)
			if option.technique == "rank":
				peak_dict[peak] = rank
			else:
				peak_dict[peak] = signal
			ranks.append(rank)
		
		# adjust peak ranks (in case you are working with HOT-filtered sets):
		if option.technique == "rank":
			minRank = min(ranks)
			for peak in peak_dict:
				peak_dict[peak] = peak_dict[peak] - minRank + 1
				
		# collect factors for pre-filtering?
		if option.group == "peak.genes":
			print "Collecting factors targeted for grouping..."
			filterGenes = list()
			inlines = open(compiledfile).readlines()
			for inline in inlines:
				chrm, start, end, region, occupancy, strand, density, complexity, factorCount, contextCount, infos = inline.strip().split("\t")
				region = ":".join([chrm, start, end])
				peaks = infos.split(";")
				peaks.pop(0)
				for peak in peaks:
					dataset = peak.split(":")[0].replace("POL2", "AMA-1")
					organism, strain, factor, context, institute = metrn.labelComponents(dataset, target="subcomponents")
					if not factor in filterGenes:
						filterGenes.append(factor)
		
		# load cellular expression matrix:
		print "Loading cellular expression data..."
		expression_matrix = dict()
		hd = general.build_header_dict(expressionfile)
		inlines = open(expressionfile).readlines()
		inlines.pop(0)
		for inline in inlines:
			initems = inline.strip().split("\t")
			gene, cell, expression = initems[hd["gene"]], initems[hd["cell"]], float(initems[hd["fraction.expression"]])
			if option.group == "expression.genes" or gene in filterGenes:
				if expression >= option.expression:
					if not cell in expression_matrix or gene in inclusionGenes:
						expression_matrix[cell] = list()
					expression_matrix[cell].append(gene)
				
		# generate cellular expression groups: each is a group of genes expressed in a population of cells:
		print "Generating cellular expression groups..."
		group_matrix = dict()
		for cell in expression_matrix:
			group = ",".join(sorted(expression_matrix[cell]))
			if not group in group_matrix:
				group_matrix[group] = list()
			group_matrix[group].append(cell)
		
		# generate group codes:
		group2code_dict, i = dict(), 1
		for group in sorted(group_matrix.keys()):
			group2code_dict[group] = "G" + general.indexTag(i, len(group_matrix.keys()))
			i += 1
		print "Groups constructed:", len(group_matrix)
		
		# prebuild factors, regions, and matrix:
		print
		print "Preloading factors and regions..."
		k, coverage, factors, regions, matrix = 0, dict(), list(), list(), dict()
		contexts_processed = list()
		inlines = open(compiledfile).readlines()
		for inline in inlines:
			chrm, start, end, region, occupancy, strand, density, complexity, factorCount, contextCount, infos = inline.strip().split("\t")
			region = ":".join([chrm, start, end])
			altregion = ".".join([chrm, start, end])
			peaks = infos.split(";")
			peaks.pop(0)
			k += 1
			
			#print k, region, len(genes), len(cells)
			for group in sorted(group_matrix.keys()):
				code = group2code_dict[group]
				groupGenes = group.split(",")
				groupCells = group_matrix[group]
				for peak in peaks:
					dataset = peak.split(":")[0].replace("POL2", "AMA-1")
					organism, strain, factor, context, institute = metrn.labelComponents(dataset, target="subcomponents")
					if factor in groupGenes:
						
						if option.target == "organism":
							factor_label = factor
							region_label = organismTag.upper() + ":" + region
						elif option.target == "factor":
							factor_label = factor
							region_label = code + ":" + region
						elif option.target == "factor.context":
							factor_label = ".".join([factor, context])
							region_label = code + ":" + region
						elif option.target == "factor.institute":
							factor_label = ".".join([factor, institute])
							region_label = code + ":" + region
						elif option.target == "factor.context.institute":
							factor_label = ".".join([factor, context, institute])
							region_label = code + ":" + region
						elif option.target == "dataset":
							factor_label = dataset
							region_label = code + ":" + region
						elif option.target == "cellular":
							factor_label = factor
							region_label = ":".join([classTissue, specificTissue, cell, altregion])
							
						factors.append(factor_label)
						regions.append(region_label)
						if not factor_label in coverage:
							coverage[factor_label] = list()
						if not cell in coverage[factor_label]:
							coverage[factor_label].append(cell)
						if not factor_label in matrix:
							matrix[factor_label] = dict()
						if not region_label in matrix[factor_label]:
							matrix[factor_label][region_label] = list()
						matrix[factor_label][region_label].append(peak_dict[peak])
		
		# simplify factor and region lists:
		factors = sorted(list(set(factors)))
		regions = sorted(list(set(regions)))
		
		# adjust signals and export matrix:
		exportMatrix(factors, regions, matrix, coverage, technique=option.technique, outfile=f_outfile)
	
	
	# generate bed files from matrix files mode:
	elif option.mode == "matrix.regions":
	
		print
		for infile in os.listdir(matrixpath):
			
			print "Processing:", infile
			f_output = open(modulespath + infile.replace(".txt", ".bed"), "w")
			inlines = open(matrixpath + infile).readlines()
			header = inlines.pop(0).strip().split("\t")
			factors = ",".join(general.clean(header))
			incount = str(len(inlines))
			k = 1
			for inline in inlines:
				values = inline.strip().split("\t")
				module = values.pop(0)
				tag, chrm, start, stop = module.split(":")
				output = [chrm, start, stop, module, sum(map(float, values)), ".", tag, len(values), ",".join(values), factors]
				print >>f_output, "\t".join(map(str, output))
				k += 1
			f_output.close()
		print
	
	# annotation mode:
	elif option.mode == "matrix.reports":
	
		# Note: Used to require BEDOPS for function; now uses bedTools' "closestBed" function.
		#sort-bed unsortedData.bed > sortedData.bed
		#closest-features --closest repOrigins.bed TSS.bed
		
		# define treatment flags:
		#if option.selection == "avg":
		#	def expressionSelection(values):
		#		return numpy.mean(values)
		#elif option.selection == "max":
		#	def expressionSelection(values):
		#		return max(values)
		#else:
		#	raise "Error: Expression selection (option.selection) not recognized!"
		
		print
		print "Loading target contexts..."
		codeContexts, targetContexts = metrn.options_dict["contexts.extended"][option.contexts]
		if option.name != "OFF":
			codeContexts = option.name
		
		print "Loading binding modules..."
		factorMatrix, moduleMatrix, coordsMatrix, signatureMatrix = loadMatrix(matrixpath + option.infile.replace(".bed", ".txt"))
		modules = sorted(moduleMatrix.keys())
		
		# define input and output files:
		organism, selection, factors, stage, processing = option.peaks.split("_")
		bedfile = modulespath + option.infile
		tssfile = path_dict[option.source] + option.tss
		if option.parameters == "OFF":
			profile = path_dict[option.source] + option.tss
		else:
			upstream, dnstream = option.parameters.split(",")
			profile = path_dict[option.source] + option.tss.replace(".nh.bed", "_slopbed_up" + upstream + "_dn" + dnstream + ".nh.bed")
		
		# Note: The input files required for this analysis are as follows:
		# bedfile: input regions
		# tssfile: input TSS coordinates (for distance calculation)
		# profile: input promoter regions (for overlap identification)
		
		print "Finding overlapping TSSs..."
		ovpfile = numberspath + option.infile.replace(".bed", "_overlapTSS_" + stage + ".bed")
		command = "intersectBed -a " + bedfile + " -b " + profile  + " -wb > " + ovpfile
		os.system(command)
		
		print "Finding closest TSSs..."
		disfile = numberspath + option.infile.replace(".bed", "_closestTSS_" + stage + ".bed")
		command = "closestBed -d -a " + bedfile + " -b " + tssfile + " > " + disfile
		os.system(command)
		
		print "Loading overlapping TSSs..."
		tssMatrix, tssCount, tssMinus, hitCount = dict(), 0, 0, 0
		for inline in open(ovpfile).readlines():
			initems = inline.strip().split("\t")
			module = initems[3]
			if "Transcript:" in inline:
				if not module in tssMatrix:
					tssMatrix[module] = list()
				else:
					tssMinus += 1
				tssCount += 1
				for initem in initems:
					for tssName in initem.split(","):
						if "Transcript:" in tssName:
							tssName = tssName.replace("Transcript:", "")
							tssMatrix[module].append(tssName)
							#if len(initem.split(",")) > 1:
							#	print module, tssName
							#	pdb.set_trace()
							hitCount += 1
						
		print "Loading TSSs distances..."
		distMatrix = dict()
		for inline in open(disfile).readlines():
			initems = inline.strip().split("\t")
			module = initems[3]
			if "Transcript:" in inline:
				if not module in distMatrix:
					distMatrix[module] = dict()
				distance = int(initems[len(initems)-1])
				for initem in initems:
					for tssName in initem.split(","):
						if "Transcript:" in tssName:
							tssName = tssName.replace("Transcript:", "")
							distMatrix[module][tssName] = distance
						
		# load header:
		if option.header == "OFF" and option.indexes != "OFF":
			header_dict = dict()
			for details in option.indexes.split(","):
				column, index = details.split(":")
				header_dict[column] = int(index)
			header = False
		elif option.header == "ON":
			header_dict = "auto"
			header = True
		
		# load expression:
		print "Loading expression values (from overlapping TSSs)..."
		expressionDict = general.build2(path_dict[option.source] + option.input, id_column=option.id, header_dict=header_dict, header=header, datatype=option.datatype, skip=True)
		
		# check expression values per module:
		expressionMatrix = dict()
		print "Mapping expression to modules..."
		for module in tssMatrix:
			for tssName in tssMatrix[module]:
				tssParts = tssName.split("_")
				tssContext = tssParts[0]
				transcript = tssParts[len(tssParts)-1]
				if transcript in expressionDict:
					if not module in expressionMatrix:
						expressionMatrix[module] = dict()
					if not tssContext in expressionMatrix[module]:
						expressionMatrix[module][tssContext] = dict()
					expressionMatrix[module][tssContext][transcript] = float(expressionDict[transcript]["score"])
			
		print "Exporting expression values per module..."
		f_output = open(numberspath + option.infile.replace(".bed", "_expression_" + codeContexts + ".exp"), "w") # export information for modules with stage-matched expression measurements only
		c_output = open(numberspath + option.infile.replace(".bed", "_expression_" + codeContexts + ".com"), "w") # export information for modules with expression measurements only
		a_output = open(numberspath + option.infile.replace(".bed", "_expression_" + codeContexts + ".all"), "w") # export information for all modules
		print >>f_output, "\t".join(["module", "factors", "expression.avg", "expression.max", "expression.sum", "distance.med", "distance.top", "distance.low", "distance.std", "distance.min", "distance.max"])
		print >>c_output, "\t".join(["module", "factors", "expression.avg", "expression.max", "expression.sum", "distance.med", "distance.top", "distance.low", "distance.std", "distance.min", "distance.max"])
		print >>a_output, "\t".join(["module", "factors", "expression.avg", "expression.max", "expression.sum", "distance.med", "distance.top", "distance.low", "distance.std", "distance.min", "distance.max"])
		for module in modules:
			factorCount = sum(map(int, moduleMatrix[module].values()))
			#print >>c_output, "\t".join(map(str, [module, factorCount, value]))
			
			# calculate distance statistics:
			distances = distMatrix[module].values()
			distMean = numpy.mean(distances)
			distMedian = int(numpy.median(distances))
			distTop = Quantile(distances, 0.90)
			distLow = Quantile(distances, 0.10)
			distStDev = numpy.std(distances)
			distMin = int(min(distances))
			distMax = int(max(distances))
			
			# store expression values from complete (com) and stage-expressed (exp) TSSs only:
			comValues, expValues = list(), list()
			if module in expressionMatrix:
				for tssContext in expressionMatrix[module]:
					for value in expressionMatrix[module][tssContext].values():
						comValues.append(value)
						if tssContext in targetContexts:
							expValues.append(value)
				print >>c_output, "\t".join(map(str, [module, factorCount, numpy.mean(comValues), max(comValues), sum(comValues), distMedian, distTop, distLow, distStDev, distMin, distMax]))	
				if expValues != list():
					print >>f_output, "\t".join(map(str, [module, factorCount, numpy.mean(expValues), max(expValues), sum(expValues), distMedian, distTop, distLow, distStDev, distMin, distMax]))
				else:
					expValues = [0]
			else:
				comValues, expValues = [0], [0]
			print >>a_output, "\t".join(map(str, [module, factorCount, numpy.mean(comValues), max(comValues), sum(comValues), distMedian, distTop, distLow, distStDev, distMin, distMax]))
		
		# close output files:
		f_output.close()
		c_output.close()
		a_output.close()
		
		print
		print "Regulatory modules:", len(modules)
		print "Assigned TSS modules:", len(tssMatrix), "(" + str(round(float(len(tssMatrix))/len(modules), 5)) + ")"
		print "Expression modules:", len(expressionMatrix), "(" + str(round(float(len(expressionMatrix))/len(tssMatrix), 5)) + ")"
		print
		
	
	# report unique signatures mode:
	elif option.mode == "signatures":
	
		print
		print "Loading binding matrix..."
		
		k, regions, signatures = 0, list(), list()
		indata = open(matrixpath + option.infile)
		inline = indata.readline()
		inline = indata.readline()
		while inline:
			initems = inline.strip().split("\t")
			region = initems.pop(0)
			signature = ",".join(initems)
			regions.append(region)
			signatures.append(signature)
			inline = indata.readline()
			k += 1
		uniqueRegions = set(regions)
		uniqueSignatures = set(signatures)
		
		print "Loaded regions:", k
		print "Loaded signatures:", len(uniqueSignatures)
		print "Unique regions:", float(len(uniqueRegions))/len(uniqueSignatures)
		#print "Unique signatures:", float(len(signatureMatrix))/len(moduleMatrix)
		#print "Coords expansion (x):", float(len(moduleMatrix))/len(coordsMatrix)
		print
		#pdb.set_trace()
		
	
	# convert region to bed mode:
	elif option.mode == "convert.regions":
	
		# update results path, regions path and code path:
		resultspath = resultspath + option.neurons
		summarypath = resultspath + "/summary/"
		regionspath = resultspath + "/regions/"
		codespath = resultspath + "/codes/"
		rawpath = regionspath + "raw/"
		bedpath = regionspath + "bed/"
		general.pathGenerator(summarypath)
		general.pathGenerator(rawpath)
		general.pathGenerator(bedpath)
			
		# convert region files:
		print
		print "Converting regions to BED format..."
		x = 0
		rawfiles = general.clean(os.listdir(rawpath), ".DS_Store")
		for rawfile in rawfiles:
			name = rawfile.replace(".regions", "")
			f_output = open(bedpath + rawfile.replace(".regions", ".bed"), "w")
			rawlines = open(rawpath + rawfile).readlines()
			rawcount = str(len(rawlines))
			k = 1
			for rawline in rawlines:
				tag, chrm, start, stop = rawline.strip().split(" ")
				module = ":".join([tag, chrm, start, stop])
				print >>f_output, "\t".join([chrm, start, stop, module, "0", "+", tag, rawcount])
				k += 1
				x += 1
			f_output.close()
		print
	
	
	# generate neuron reports mode:
	elif option.mode == "neuron.summary":
	
		# load matrix file:
		matrixfile = neuronMapping[option.neurons].replace("<<PEAKS>>", option.peaks)
	
		# load matrix:
		print
		print "Loading matrix data..."
		if option.fast == "ON":
			factors, factorMatrix, moduleMatrix, factorIndexes, moduleIndexes = loadFaster(matrixpath + matrixfile)
		else:
			factorMatrix, moduleMatrix, coordsMatrix, signatureMatrix = loadMatrix(matrixpath + matrixfile, coords="OFF", signatures="OFF")
			factors = sorted(factorMatrix.keys())
		#print matrixfile
		#print sorted(factors)
		#pdb.set_trace()
			
		# factorMatrix > used only to list factors
		# moduleMatrix > is used critically
		# coordsMatrix > not used
		# signatureMatrix > not used
		
		# update results path, regions path and code path:
		resultspath = resultspath + option.neurons
		summarypath = resultspath + "/summary/"
		regionspath = resultspath + "/regions/"
		codespath = resultspath + "/codes/"
		rawpath = regionspath + "raw/"
		bedpath = regionspath + "bed/"
		general.pathGenerator(summarypath)
		general.pathGenerator(rawpath)
		general.pathGenerator(bedpath)
		
		# load region data from BED files:
		print "Loading regions and classes..."
		classLabels, classMatrix, sizeMatrix, signalMatrix, binaryMatrix = list(), dict(), dict(), dict(), dict()
		bedfiles = os.listdir(bedpath)
		for bedfile in bedfiles:
			neuron = int(bedfile.strip(".bed").strip("neuron"))
			regionCount, classDict = 0, dict()
			modules = list()
			for inline in open(bedpath + bedfile).readlines():
				chrm, start, stop, module, score, strand, classLabel, moduleCount = inline.strip().split("\t")
				modules.append(module)
				if not neuron in sizeMatrix:
					sizeMatrix[neuron] = int(moduleCount)
				if not classLabel in classDict:
					classLabels.append(classLabel)
					classDict[classLabel] = 0
				classDict[classLabel] += 1
				regionCount += 1
			if regionCount > 0:
				classMatrix[neuron] = classDict
				moduleCounts, signalSignature, binarySignature = dict(), list(), list()
				for module in modules:
					for factor in factors:
						if not factor in moduleCounts:
							moduleCounts[factor] = 0
						if option.fast == "ON":
							if not module in moduleMatrix:
								print module, factor, factorIndexes[factor]
								pdb.set_trace()
							moduleCounts[factor] += int(moduleMatrix[module][factorIndexes[factor]])
							
						else:
							moduleCounts[factor] += int(moduleMatrix[module][factor])
						#print factor, moduleCounts[factor]
						#pdb.set_trace()
				for factor in factors:
					factorSignal = float(moduleCounts[factor])/int(moduleCount)
					signalSignature.append(factorSignal)
					if factorSignal >= float(option.fraction):
						binarySignature.append(1)
					else:
						binarySignature.append(0)
				signalMatrix[neuron] = signalSignature
				binaryMatrix[neuron] = binarySignature
		classLabels = sorted(list(set(classLabels)))
		
		# define output file:
		f_output = open(summarypath + "mapneurons_summary.txt", "w")
		print >>f_output, "\t".join(["neuron", "modules", "classes", "factors", "max.class", "max.value", "max.fraction", "max.bias", "binary.signature", "signal.signature", "factor.ids", "class.ids"] + sorted(classLabels))
		
		print 
		print "Matrix factors:", len(factors)
		print "Classes:", ", ".join(classLabels), "(" + str(len(classLabels)) + ")"
		
		# process neurons:
		maxModules, maxFactors = 0, 0
		for neuron in sorted(classMatrix.keys()):
			classHits = general.valuesort(classMatrix[neuron])
			classHits.reverse()
			classMax = classHits[0]
			valueMax = classMatrix[neuron][classMax]
			valueSum = sum(classMatrix[neuron].values())
			fractionMax = float(valueMax)/valueSum
			
			if sum(binaryMatrix[neuron]) > maxFactors:
				maxFactors = sum(binaryMatrix[neuron])
			if sizeMatrix[neuron] > maxModules:
				maxModules = sizeMatrix[neuron]
			
			if len(classLabels) > 1:
				biasMin = float(1)/len(classLabels)
				biasMax = float(fractionMax-biasMin)/(1-biasMin)
			else:
				biasMax = 1
			output = ["neuron" + str(neuron), sizeMatrix[neuron], len(classLabels), sum(binaryMatrix[neuron]), classMax, valueMax, fractionMax, biasMax, ",".join(map(str, binaryMatrix[neuron])), ",".join(map(str, signalMatrix[neuron])), ",".join(map(str, factors)), ",".join(sorted(classMatrix[neuron].keys()))]
			for classLabel in sorted(classLabels):
				if classLabel in classMatrix[neuron]:
					classCount = classMatrix[neuron][classLabel]
				else:
					classCount = 0
				output.append(float(classCount)/valueSum)
			print >>f_output, "\t".join(map(str, output))
		f_output.close()
		print
		
		# define factor-count report file:
		f_output = open(summarypath + "mapneurons_factors.txt", "w")
		print >>f_output, "\t".join(["factors", "modules.avg", "modules.med", "modules.std", "modules.75", "modules.25", "fractions.avg", "fractions.med", "fractions.std", "fractions.75", "fractions.25", "bias.avg", "bias.med", "bias.std", "bias.75", "bias.25"])
		reportDict = general.build2(summarypath + "mapneurons_summary.txt", id_column="neuron")
		for factorCount in range(1, maxFactors+1):
			moduleCounts = list()
			fractionValues = list()
			biasValues = list()
			for neuron in reportDict:
				if int(reportDict[neuron]["factors"]) == factorCount:
					moduleCounts.append(int(reportDict[neuron]["modules"]))
					fractionValues.append(float(reportDict[neuron]["max.fraction"]))
					biasValues.append(float(reportDict[neuron]["max.bias"]))
			if moduleCounts == list():
				break
			else:
				output = [factorCount, 
							numpy.mean(moduleCounts), numpy.median(moduleCounts), numpy.std(moduleCounts), scipy.percentile(moduleCounts, 75), scipy.percentile(moduleCounts, 25), 
							numpy.mean(fractionValues), numpy.median(fractionValues), numpy.std(fractionValues), scipy.percentile(fractionValues, 75), scipy.percentile(fractionValues, 25), 
							numpy.mean(biasValues), numpy.median(biasValues), numpy.std(biasValues), scipy.percentile(biasValues, 75), scipy.percentile(biasValues, 25)]
				print >>f_output, "\t".join(map(str, output))
		f_output.close()

	
	# annotate numbers (expression, TSS distance, etc..) to neurons mode:
	elif option.mode == "neuron.reports":
	
		# Role: Associates values (columns) in BED files to neurons. 
		# Load: Loads setup (configuration) information indicating what to extract from where and what to call it. The target setup file format (tab-delimited):
		# 	name	format	source	infile	id	column	indexing
		# 	format: paired, bed, etc...
		# 	indexing: numeric, labels
		
		# load target setup dictionary and infiles:
		print
		print "Loading input setup..."
		setupDict = general.build2(setuppath + option.infile, id_column="name")
		
		# load numbers of interest:
		print "Loading input values..."
		valueDict = general.valueLoad(setupDict, pathDict=numberspath)
		
		# update results path, regions path and code path:
		resultspath = resultspath + option.neurons
		summarypath = resultspath + "/summary/"
		regionspath = resultspath + "/regions/"
		reportspath = resultspath + "/reports/" + option.infile.replace(".txt","") + "/"
		codespath = resultspath + "/codes/"
		rawpath = regionspath + "raw/"
		bedpath = regionspath + "bed/"
		general.pathGenerator(summarypath)
		general.pathGenerator(reportspath)
		general.pathGenerator(rawpath)
		general.pathGenerator(bedpath)

		# map regions to neurons:
		print "Loading regions to neurons..."
		neuronDict = mapRegions(rawpath, regions=valueDict.keys())
		
		# load neuron summary data:
		print "Loading neuron summaries..."
		summaryDict = general.build2(summarypath + "mapneurons_summary.txt", id_column="neuron")
		
		# start output file:
		f_output = open(reportspath + "mapneurons_signals.txt", "w")
		r_output = open(reportspath + "mapneurons_reports.txt", "w")
		e_output = open(reportspath + "mapneurons_extends.txt", "w")
		print >>f_output, "\t".join(["neuron", "factor", "signal"])
		
		# export expanded neuron data:
		print "Exporting neuron report..."
		inlines  = open(summarypath + "mapneurons_summary.txt").readlines()
		header = inlines.pop(0).strip().split("\t")
		headerDict = general.build_header_dict(header, mode="list")
		for name in sorted(setupDict.keys()):
			header.append(name + ".avg")
			header.append(name + ".med")
			header.append(name + ".std")
		print >>r_output, "\t".join(map(str, header))
		start = True
		for inline in inlines:
			output = inline.strip().split("\t")
			neuron = output[0]
			neuronCode = output[0].replace("neuron", "N")
			output[0] = neuronCode
			if start:
				factors = output[headerDict["factor.ids"]].split(",")
				header.extend(factors)
				print >>e_output, "\t".join(map(str, header))
			start = False
			if neuron in neuronDict:
				for name in sorted(setupDict.keys()):
					values = list()
					for region in neuronDict[neuron]:
						values.append(float(valueDict[region][name]))
					#if "EX:I:10185536:10185979" in neuronDict[neuron]:
					#	print "Found:", "EX:I:10185536:10185979"
					#	print len(neuronDict[neuron])
					#	print float(valueDict["EX:I:10185536:10185979"][name])
					#	pdb.set_trace()
					output.append(numpy.mean(values))
					output.append(numpy.median(values))
					output.append(numpy.std(values))
				print >>r_output, "\t".join(map(str, output))
				
				# extend to include factor signatures:
				signals = output[headerDict["signal.signature"]].split(",")
				output.extend(signals)
				print >>e_output, "\t".join(map(str, output))
				
				# export to signal matrix:
				for k in range(0, len(signals)):
					print >>f_output, "\t".join([neuronCode, factors[k], signals[k]])
		
		# close output file:
		f_output.close()
		r_output.close()
		e_output.close()
		
		# memory surveillance:
		#import resource
		#print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
		
		# Websites are (can be) pieces of art.
		print
		
	
	# mapping features to neurons mode:
	elif option.mode == "neuron.mapping":
	
		# Application: I don't exactly remember what this is used for!
	
		# update results path, regions path and code path:
		resultspath = resultspath + option.neurons
		summarypath = resultspath + "/summary/"
		regionspath = resultspath + "/regions/"
		overlappath = resultspath + "/overlap/" + option.name + "/"
		mappingpath = resultspath + "/mapping/" + option.name + "/"
		codespath = resultspath + "/codes/"
		rawpath = regionspath + "raw/"
		bedpath = regionspath + "bed/"
		general.pathGenerator(summarypath)
		general.pathGenerator(overlappath)
		general.pathGenerator(mappingpath)
		general.pathGenerator(rawpath)
		general.pathGenerator(bedpath)

		# load gene-expression similarities:
		fractionCorrelations = general.build2(cellspath + "overlap/" + option.parameters + "/mapcells_" + option.parameters + "_matrix_overlap", i="i", j="j", x="fraction.cor", mode="matrix")
		overlapsCorrelations = general.build2(cellspath + "overlap/" + option.parameters + "/mapcells_" + option.parameters + "_matrix_overlap", i="i", j="j", x="overlap.sum", mode="matrix")
		
		# load header:
		if option.header == "OFF":
			header = False
		elif option.header == "ON":
			header = True
		
		# load neurons per feature:
		index, matrix, quantx = 1, dict(), dict()
		for bedfile in os.listdir(bedpath):
			sys.stdout.write(str("\rProcessing BED file: {0}" + " (" + str(len(os.listdir(bedpath))) + ")").format(index))
			sys.stdout.flush()
			index += 1
					
			# determine overlap between features and neuron regions:
			queryfile = path_dict[option.source] + option.infile
			targetfile = bedpath + bedfile
			outputfile = overlappath + bedfile
			if not bedfile in os.listdir(overlappath) or option.overwrite == "ON":
				command = "intersectBed -a " + queryfile + " -b " + targetfile + " -wa >" + outputfile
				os.system(command)
				
			# store info:
			neuron = bedfile.strip(".bed")
			for inline in open(outputfile).readlines():
				feature = inline.strip().split("\t")[3]
				if not feature in matrix:
					matrix[feature] = list()
					quantx[feature] = dict()
				if not neuron in quantx[feature]:
					quantx[feature][neuron] = 0
				quantx[feature][neuron] += 1
				matrix[feature].append(neuron)
				matrix[feature] = list(set(matrix[feature]))
				#if quantx[feature][neuron] > 1:
				#	print "\n" + feature, neuron, quantx[feature][neuron]
				#	pdb.set_trace()
		
		# collect neurons:
		neurons = list()
		for i in matrix:
			neurons.extend(matrix[i])
		neurons = list(set(neurons))
				
		# load queries:
		p, f, n = 0, 0, 0
		queries = os.listdir(cellspath + "collect/" + option.query)
		for query in queries:
			if query in name2id_dict:
				f += 1
				if name2id_dict[query] in matrix:
					p += 1
			else:
				n += 1
		print
		print "Discovery:", p, f, n
		print
		
		# create overlaps:
		index = 0
		outputfile = mappingpath + "mapneurons_mapping_" + option.infile
		f_output = open(outputfile, "w")
		print >>f_output, "\t".join(["i", "j", "i.count", "j.count", "overlap", "total", "overlap.avg", "overlap.max", "overlap.sum", "overlap.pvalue", "adjusted.pvalue", "coassociation.cor", "coassociation.sig", "expression.cor", "expression.sum", "flag"])
		adjust = p*p
		for i in queries:
			for j in queries:
			
				process = False
				if i in name2id_dict and j in name2id_dict:
					iRNA = name2id_dict[i] 
					jRNA = name2id_dict[j]
					if iRNA in matrix and jRNA in matrix:
						process = True
				
				if process:
				
					sys.stdout.write(str("\rExporing matrix: {0}" + " (" + str(adjust) + ")").format(index))
					sys.stdout.flush()
					index += 1
			
					neuronFlag = 0
					ivector, jvector = list(), list()
					for neuron in neurons:
						#if neuron in quantx[iRNA] and neuron in quantx[jRNA]:
						#	ivector.append(quantx[iRNA][neuron])
						#	jvector.append(quantx[jRNA][neuron])
						#else:
						#	ivector.append(0)
						#	jvector.append(0)
						if neuron in quantx[iRNA] or neuron in quantx[jRNA]:
							neuronFlag += 1
						if neuron in quantx[iRNA]:
							ivector.append(quantx[iRNA][neuron])
						else:
							ivector.append(0)
						if neuron in quantx[jRNA]:
							jvector.append(quantx[jRNA][neuron])
						else:
							jvector.append(0)
					correlation, corPvalue = pearsonr(ivector, jvector)
					
					icounts = len(matrix[iRNA])
					jcounts = len(matrix[jRNA])
					overlap = len(set(matrix[iRNA]).intersection(set(matrix[jRNA])))
					total = len(neurons)
					
					setunion = len(set(matrix[iRNA]).union(set(matrix[jRNA])))
					ioverlap = float(overlap)/icounts
					joverlap = float(overlap)/jcounts
					overlap_avg = numpy.mean([ioverlap, joverlap])
					overlap_max = max([ioverlap, joverlap])
					overlap_sum = float(overlap)/setunion
				
					
					"""
					# Hypergeometric paramters:
					m = len(avalues) # number of white balls in urn
					n = len(universe) - len(avalues) # number of black balls in urn
					N = len(bvalues) # number of balls drawn from urn
					x = len(overlap) # number of white balls in drawn
					
					# If I pull out all balls with elephant tatoos (N), is the draw enriched in white balls?:
					pvalue = hyper.fishers(x, m+n, m, N, method="right")
					"""
					
					# calculate probability mass function (PMF):
					pvalue = hyper.fishers(overlap, total, icounts, jcounts, method="right")
					adjPvalue = hyper.limit(pvalue*adjust)
					
					# export overlap:
					output = [i, j, icounts, jcounts, overlap, total, overlap_avg, overlap_max, overlap_sum, pvalue, adjPvalue, correlation, corPvalue, fractionCorrelations[i][j], overlapsCorrelations[i][j], neuronFlag]
					print >>f_output, "\t".join(map(str, output))		
		
		# close output file:
		f_output.close()
		
		#print matrix.keys()
		#print len(matrix.keys())
		#pdb.set_trace()
	
	
	# collect features from different neurons mode:
	elif option.mode == "neuron.compile":
		
		# start master dictionary and factor list:
		masterDict = dict()
		factorList = list()
		
		# go through SOMs, preloading factors (from signatures):
		indexes = len(option.peaks.split(","))
		for index in range(0, indexes):
			
			# prepare relevant peakset and neurons:
			peakset = option.peaks.split(",")[index]
			somName = option.neurons.split(",")[index]
			
			# transfer factors:
			somNotes = path_dict["neurons"] + peakset + "/" + option.technique + "/results/" + somName + "/summary/mapneurons_summary.txt"
			notesDict = general.build2(somNotes, id_column="neuron")
			factorList.extend(notesDict[notesDict.keys()[0]]["factor.ids"].split(","))
		factorList = sorted(list(set(factorList)))
		
		# go through SOMs, preloading factors (from signatures):
		indexes = len(option.peaks.split(","))
		for index in range(0, indexes):
			
			# prepare relevant peakset and neurons:
			peakset = option.peaks.split(",")[index]
			somName = option.neurons.split(",")[index]
			
			# specify results:
			somNotes = path_dict["neurons"] + peakset + "/" + option.technique + "/results/" + somName + "/summary/mapneurons_summary.txt"
			
			# load results:
			notesDict = general.build2(somNotes, id_column="neuron")
			
			# transfer signatures:
			for neuron in notesDict:
				masterSignature = list()
				neuronSignature = notesDict[neuron]["signal.signature"]
				neuronFactorIDs = notesDict[neuron]["factor.ids"]
				neuronEncodings = dict()
				neuronElements = len(neuronSignature.split(","))
				for neuronElement in range(0, neuronElements):
					factor = neuronFactorIDs.split(",")[neuronElement]
					xscore = neuronSignature.split(",")[neuronElement]
					if float(xscore) >= option.fraction:
						neuronEncodings[factor] = "1"
					else:
						neuronEncodings[factor] = "0"
			
				for factor in factorList:
					if factor in neuronEncodings:
						masterSignature.append(neuronEncodings[factor])
					else:
						masterSignature.append("0")
				#print sum(map(int, notesDict[neuron][option.id].split(",")))
				#print sum(map(int, masterSignature))
				#pdb.set_trace()
				
				# rename neurons?
				for scheme in option.rename.split(","):
					target, replace = scheme.split(":")
					neuron = neuron.replace(target, replace)
				
				# store signature data:	
				masterSignature = ",".join(masterSignature)
				if not masterSignature in masterDict:
					masterDict[masterSignature] = dict()
				if not somName in masterDict[masterSignature]:
					masterDict[masterSignature][somName] = list()
				masterDict[masterSignature][somName].append(neuron)

		# store counts per signature, and signatures per count:
		countDict, signDict, contextDict = dict(), dict(), dict()
		t = 0
		for signature in masterDict:
			n = 0
			c = len(masterDict[signature])
			for somName in masterDict[signature]:
				for neuron in masterDict[signature][somName]:
					n += 1
					t += 1
			if not n in countDict:
				countDict[n] = 0
			if not c in contextDict:
				contextDict[c] = 0
			countDict[n] += 1
			contextDict[c] += 1
			signDict[signature] = n
		counts = general.valuesort(signDict)
		counts.reverse()
			
		print
		print "Total Factors:", len(factorList)
		print "Total Signatures:", len(masterDict)
		print "Total Neurons:", t
		print "Specificity:", round(100*float(len(masterDict))/t, 2), "%"
		print
		print "Maximum Neurons per Signature:", max(countDict)
		print "Minimum Neurons per Signature:", min(countDict)
		print
		print "Neurons per signature:"
		for count in sorted(countDict.keys()):
			print count, ":", countDict[count], "signatures"
		print
		print "Contexts per signature:"
		for count in sorted(contextDict.keys()):
			print count, ":", contextDict[count], "signatures"
		print
		
		# export signatures per neuron:
		f_outfile = compilepath + "mapneurons_compile_" + option.name + "_signs.txt"
		f_output = open(f_outfile, "w")
		print >>f_output, "\t".join(["signature", "som.name", "som.count", "neuron.count", "total.count", "neuron.ids", "factor.ids"])
		for signature in sorted(masterDict.keys()):
			for somName in sorted(masterDict[signature].keys()):
				output = [signature, somName, len(masterDict[signature]), len(masterDict[signature][somName]), signDict[signature], ",".join(sorted(masterDict[signature][somName])), ",".join(factorList)]
				print >>f_output, "\t".join(map(str, output))
		f_output.close()
		#pdb.set_trace()
		
	
if __name__ == "__main__":
	main()
	
	print "Completed:", time.asctime(time.localtime())

#python mapNeurons.py --path ~/ceTRN --mode collapse --peaks optimal_standard_totals_sx_filter

#python mapNeurons.py --path ~/ceTRN --mode matrix --peaks optimal_standard_totals_sx_filter --contexts shift.condensed --target factor
#python mapNeurons.py --path ~/ceTRN --mode matrix --peaks optimal_standard_totals_sx_filter --contexts basic.condensed --target factor
#python mapNeurons.py --path ~/ceTRN --mode matrix --peaks optimal_standard_totals_sx_filter --contexts embryonic --target factor

#python mapNeurons.py --path ~/ceTRN --mode matrix --peaks optimal_standard_totals_sx_filter --contexts basic.condensed --target factor.context
