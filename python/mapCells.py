#!/usr/bin/env python
# perform cellular-resolution expression analyses!

import sys
import time
import optparse
import general
import hyper
import numpy
import math
import pickle
import pdb
import metrn
import modencode
import itertools
import os
import re

import datetime
import calendar
#import simplejson as json
from scipy.stats.stats import pearsonr
from runner import *

from scipy import stats
from network import Network
from network import export

print "Command:", " ".join(sys.argv)
print "Timestamp:", time.asctime(time.localtime())


""" define functions of internal use """

""" define a function to recover cells in a time range """
def getTargetCells(inobject="", inpath="", mode="collection", timeRange=list()):

	# grab cells from collection:
	if mode == "collection":
		
		# load collection cells:
		cells = list()
		for gene in os.listdir(inpath):
			cells.extend(open(inpath + gene).read().split("\n"))
		cells = general.clean(sorted(list(set(cells))))
		print "Loading collection cells:", len(cells)
		
	# grab cells from time-points:
	elif mode == "time":
	
		# load time-point cells:
		cells = list()
		for timePoint in os.listdir(inpath):
			if int(timePoint) in timeRange:
				cells += general.clean(open(inpath + timePoint).read().split("\n"))
		cells = sorted(list(set(cells)))
		print "Loading time-point/range cells:", len(cells)
	
	# return collected cells:
	return cells
	

""" define a function to construct a cell-parent relationships, and pedigree cell list """
def expressionBuilder(expressionfile, path, cutoff, minimum, metric="fraction.expression"):
	
	# build header dict:
	hd = general.build_header_dict(path + expressionfile)
	
	# process input expression data:
	quantitation_matrix, expression_matrix, tracking_matrix, trackedCells = dict(), dict(), dict(), list()
	inlines = open(path + expressionfile).readlines()
	inlines.pop(0)
	for inline in inlines:
		initems = inline.strip().split("\t")
		cell, gene, rawSignal, metricSignal = initems[hd["cell.name"]], initems[hd["gene"]], initems[hd["cell.expression"]], initems[hd[metric]]
		trackedCells.append(cell)
		
		# store expression value
		if not gene in quantitation_matrix:
			quantitation_matrix[gene] = dict()
		quantitation_matrix[gene][cell] = float(metricSignal)
		
		# store tracked and expressing cells:
		if not gene in tracking_matrix:
			expression_matrix[gene] = list()
			tracking_matrix[gene] = list()
		tracking_matrix[gene].append(cell)
		if float(metricSignal) >= float(cutoff) and float(rawSignal) >= minimum:
			expression_matrix[gene].append(cell)

	trackedCells = list(set(trackedCells))
	return quantitation_matrix, expression_matrix, tracking_matrix, trackedCells


""" define a function to construct a cell-parent relationships, and pedigree cell list """
def relationshipBuilder(pedigreefile, path, trackedCells=list(), lineages="complete", mechanism="simple"):
	cell_dict, parent_dict = dict(), dict()
	inlines = open(path + pedigreefile).readlines()
	header = inlines.pop(0)
	for inline in inlines:
		cell, binCell, parent, binParent = inline.strip().split(",")[:4]
		tissues = inline.strip().split(",")[5]
		if not parent == "" and not cell == "":
			if mechanism == "simple" or lineages == "complete" or (lineages == "tracked" and parent in trackedCells and cell in trackedCells):
				if not parent in parent_dict:
					parent_dict[parent] = list()
				parent_dict[parent].append(cell)
				cell_dict[cell] = parent
	pedigreeCells = sorted(list(set(cell_dict.keys()).union(set(parent_dict.keys()))))
	return cell_dict, parent_dict, pedigreeCells


""" define a function to generate the underlying tree of a given parent """
def treeBuilder(parent_dict, cell_dict, highlights=list(), nodeColor="#FFFFFF", lineColor="#336699", textColor="#000000", highlightColor="#CC0000"):
	
	# set color rules:
	groups = { "unknown" : textColor, "highlight" : highlightColor }
	nodeColors = { "unknown" : nodeColor, "highlight" : highlightColor }
	lineColors = { "unknown" : lineColor, "highlight" : highlightColor }
	textColors = { "unknown" : textColor, "highlight" : highlightColor }
	
	# initialize tree:
	tree = {}
	for child in cell_dict:
		parent = cell_dict[child]
		
		# determine whether to 
		pkey, ckey = "unknown", "unknown"
		if parent in highlights:
			pkey = "highlight"
		if child in highlights:
			ckey = "highlight"
		
		# make an instance of a class for the parent if neccesary:
		if not tree.has_key(parent):
			tree[parent] = {'name':parent,'group':groups[pkey],'nodeColor':nodeColors[pkey],'lineColor':lineColors[pkey],'textColor':textColors[pkey],'children':[]}
			
		# make an instance of a class for the child if neccesary:
		if not tree.has_key(child):
			tree[child] = {'name':child,'group':groups[ckey],'nodeColor':nodeColors[ckey],'lineColor':lineColors[ckey],'textColor':textColors[ckey],'children':[]}
			
		# and child object to parent if necesary:
		if not tree[child] in tree[parent]['children']:
			tree[parent]['children'].append(tree[child])
			
	return tree
	

""" define a function to generate the list of cells that are parents to a given cell """
def ascendantsCollector(cell, parent_dict, cell_dict, ascendants=list(), sort=True):
	if not cell in ascendants:
		ascendants.append(cell)
	if cell in cell_dict:
		parent = cell_dict[cell]
		ascendants.append(parent)
		ascendants = ascendantsCollector(parent, parent_dict, cell_dict, ascendants, sort=sort)
	if sort:
		return sorted(list(set(ascendants)))
	else:
		return ascendants
		

""" define a function to generate the list of cells that are progeny to a given parent """
def descendantsCollector(parent, parent_dict, cell_dict, descendants=list(), sort=True):
	if not parent in descendants:
		descendants.append(parent)
	if parent in parent_dict:
		for cell in parent_dict[parent]:
			descendants.append(cell)
			descendants = descendantsCollector(cell, parent_dict, cell_dict, descendants, sort=sort)
	if sort:
		return sorted(list(set(descendants)))
	else:
		return descendants


""" define a function to generate the list of cells that are progeny to a given parent (using combinations function) """
def lineageGenerator(parent, parent_dict, cell_dict):
	descendants = descendantsCollector(parent, parent_dict, cell_dict, descendants=list())
	gList = list()
	for r in range(1, len(descendants)+1):
		for gCells in itertools.combinations(descendants, r):
			process = True
			gCells = list(gCells)
			for gCell in gCells:
				if gCell != parent:
					if not cell_dict[gCell] in gCells:
						process = False
			if process:
				gList.append(",".join(sorted(gCells)))
	return gList


""" define a function to generate the list of cells that are progeny to a given parent (using lineage growth) """
def lineageBuilder(parent, parent_dict, cell_dict, limit="OFF", descendants="ON"):
	mList = [parent]
	for mCells in mList:
		aCells, bCells, xCells, exit = str(mCells), str(mCells), str(mCells), False
		for mCell in mCells.split(","):
			if mCell in parent_dict and len(parent_dict[mCell]) == 2:
				aCell, bCell = parent_dict[mCell]
				
				if not aCell in aCells.split(","):
					aCells = ",".join(sorted(mCells.split(",") + [aCell]))
				if not bCell in bCells.split(","):
					bCells = ",".join(sorted(mCells.split(",") + [bCell]))
				if not aCell in xCells.split(",") and not bCell in xCells.split(","):
					xCells = ",".join(sorted(mCells.split(",") + [aCell, bCell]))
				
				if not aCells in mList:
					mList.append(aCells)
				if not bCells in mList:
					mList.append(bCells)
				if not xCells in mList:
					mList.append(xCells)
		
				if limit != "OFF" and len(mList) >= limit:
					if descendants == "ON":
						aCellx = sorted(list(set(mCells.split(",") + descendantsCollector(aCell, parent_dict, cell_dict, descendants=list()))))
						bCellx = sorted(list(set(mCells.split(",") + descendantsCollector(bCell, parent_dict, cell_dict, descendants=list()))))
						xCellx = sorted(list(set(aCellx).union(set(bCellx))))
						aCellx = ",".join(aCellx)
						bCellx = ",".join(bCellx)
						xCellx = ",".join(xCellx)
						if not aCellx in mList:
							mList.append(aCellx)
						if not bCellx in mList:
							mList.append(bCellx)
						if not xCellx in mList:
							mList.append(xCellx)
					exit = True
		if exit:
			break
	return sorted(mList)


""" define a function to generate lists of related-cells from a given set of of cells """
def lineageCollector(cells, parent_dict, cell_dict, siblings="ON"):
	collections, parent_tree, cell_tree = list(), dict(), dict()
	#ascendants = ascendantsCollector(descendant, parent_tree, cell_tree, ascendants=list())
	#descendants = descendantsCollector(parent, parent_dict, cell_dict, descendants=list())
	print len(cells), cells
	for cell in sorted(cells):
		found, relatives = False, [cell]
		if cell in cell_dict:
			relatives.append(cell_dict[cell])
		if cell in parent_dict:
			relatives.extend(parent_dict[cell])
		if siblings == "ON" and cell in cell_dict:
			relatives.extend(parent_dict[cell_dict[cell]])
		r, relatives = 0, list(set(relatives).intersection(set(cells)))
		print cell, relatives, "<-- relatives"
		updated = list()
		for collection in collections:
			if set(relatives).intersection(set(collection)):
				print collection, "<-- collection"
				collection.extend(relatives)
				collection = list(set(collection))
				print collection, "<-- updated"
				r += 1
				pdb.set_trace()
			updated.append(collection)
		if r == 0:
			updated.append(relatives)
		collections = updated
	return collections
			
	

""" define a function to calculate the number of possible subsets """
def combinationCalculator(n, R):
	combinations = 0
	for r in range(1,R):
		combinations += math.factorial(n)/(math.factorial(r)*math.factorial(n-r))
	return combinations


""" define a function to calculate the number of divisions between two cells """
def divisionCalculator(aCell, aParent, parent_dict, cell_dict):
	divisions = 0
	while aCell in cell_dict and aCell != aParent:
		if cell_dict[aCell] == aParent:
			divisions += 1
			break
		else:
			aCell = cell_dict[aCell]
			divisions += 1
	return divisions
	

""" define a function that calculates the lineage distance between two cells """
def lineageDistance(aCell, bCell, parent_dict, cell_dict):
	aParents = ascendantsCollector(aCell, parent_dict, cell_dict)
	bParents = ascendantsCollector(bCell, parent_dict, cell_dict)
	xParents = set(aParents).intersection(set(bParents))
	xDistances = dict()
	#print len(xParents), aCell, bCell, ":", ", ".join(xParents)
	for xParent in xParents:
		aDistance = divisionCalculator(aCell, xParent, parent_dict, cell_dict)
		bDistance = divisionCalculator(bCell, xParent, parent_dict, cell_dict)
		xDistances[xParent] = aDistance + bDistance
	xParents = general.valuesort(xDistances)
	distance, ancestor = xDistances[xParents[0]], xParents[0]
	return distance, ancestor

	
def main():
	
	parser = optparse.OptionParser()
	parser.add_option("--path", action="store", type="string", dest="path", help="Path from script to files")
	parser.add_option("--organism", action = "store", type = "string", dest = "organism", help = "Target organism for operations...", default="OFF")
	parser.add_option("--mode", action="store", type="string", dest="mode", help="Operation modes: import, map, or other...")
	parser.add_option("--peaks", action="store", type="string", dest="peaks", help="Peaks set to be used.", default="OFF")
	parser.add_option("--infile", action="store", type="string", dest="infile", help="Input file for abundance representation")
	parser.add_option("--nuclear", action = "store", type = "string", dest = "nuclear", help = "Peaks are only nuclear?", default="ON")
	parser.add_option("--expression", action="store", type="string", dest="expression", help="Input expression file for abundance representation", default="OFF")
	parser.add_option("--pedigree", action="store", type="string", dest="pedigree", help="Input pedigree file", default="OFF")
	parser.add_option("--mapping", action="store", type="string", dest="mapping", help="Input mapping file; associates tissue labels to more generic terms!", default="OFF")
	parser.add_option("--tissues", action="store", type="string", dest="tissues", help="Input tissues file", default="OFF")
	parser.add_option("--times", action="store", type="string", dest="times", help="Input cell times file", default="OFF")
	parser.add_option("--name", action="store", type="string", dest="name", help="Output file name", default="")
	parser.add_option("--nametag", action="store", type="string", dest="nametag", help="Output file name addition tag", default="")
	parser.add_option("--collection", action="store", type="string", dest="collection", help="Cell collection subset name", default="OFF")
	parser.add_option("--technique", action = "store", type = "string", dest = "technique", help = "What kind of matrix should I build? binary, fraction, or normal", default="binary")
	parser.add_option("--neurons", action="store", type="string", dest="neurons", help="Neurons to be used for 'collection' analysis...", default="OFF")
	parser.add_option("--factors", action="store", type="string", dest="factors", help="Infer factors (OFF) or load from file?", default="OFF")
	parser.add_option("--measure", action="store", type="string", dest="measure", help="Maximum (cells) or mean", default="avg.expression")
	parser.add_option("--fraction", action="store", type="float", dest="fraction", help="Fractional expression cutoff", default=0.1)
	parser.add_option("--minimum", action="store", type="float", dest="minimum", help="Minimum raw expression cutoff", default=2000)
	parser.add_option("--inherit", action="store", type="string", dest="inherit", help="Signal inheritance policy: 'max' or 'last' of ancestor expression signals...", default="last")
	parser.add_option("--overlap", action="store", type="float", dest="overlap", help="Cellular overlap cutoff", default=0.75)
	parser.add_option("--pvalue", action="store", type="float", dest="pvalue", help="Significance cutoff", default=0.01)
	parser.add_option("--header", action="store", type="string", dest="header", help="Is there a header?", default="OFF")
	parser.add_option("--format", action="store", type="string", dest="format", help="How should formatting be done?", default="bed")
	parser.add_option("--reference", action="store", type="string", dest="reference", help="Gene-coordinate reference file", default="in2shape_ce_wormbased_COM_gx.bed")
	parser.add_option("--up", action = "store", type = "int", dest = "up", help = "Upstream space", default=0)
	parser.add_option("--dn", action = "store", type = "int", dest = "dn", help = "Downstream space", default=0)
	parser.add_option("--method", action="store", type="string", dest="method", help="Should descendant cells or descendant lineages be examined?", default="lineages")
	parser.add_option("--cells", action="store", type="string", dest="cells", help="Reduce lineage cells to tracked cells (tracked) or use complete lineage cells (complete)?", default="tracked")
	parser.add_option("--lineages", action="store", type="string", dest="lineages", help="Reduce lineage tree to tracked cells (tracked) or use complete lineage tree (complete)?", default="tracked")
	parser.add_option("--descendants", action="store", type="string", dest="descendants", help="Apply descendants cutoff?", default="OFF")
	parser.add_option("--ascendants", action="store", type="string", dest="ascendants", help="Apply ascendants cutoff?", default="OFF")
	parser.add_option("--extend", action="store", type="string", dest="extend", help="Extend to include 0 signal expression values for cells not measured?", default="OFF")
	parser.add_option("--overwrite", action="store", type="string", dest="overwrite", help="Overwrite outputs?", default="OFF")
	parser.add_option("--parameters", action="store", type="string", dest="parameters", help="Optional parameters...", default="OFF")
	parser.add_option("--limit", action="store", type="string", dest="limit", help="Limit on lineage expansion? Numeric integer.", default="OFF")
	parser.add_option("--query", action="store", type="string", dest="query", help="Query collections of cells whose enrichment will be searched in target cells", default="OFF")
	parser.add_option("--source", action="store", type="string", dest="source", help="File source for inputs...", default="OFF")
	parser.add_option("--target", action="store", type="string", dest="target", help="Target collections of cells in which enrichment is searched for", default="OFF")
	parser.add_option("--domain", action="store", type="string", dest="domain", help="Domain of co-associations for hybrid-type analyses", default="OFF")
	parser.add_option("--A", action = "store", type = "string", dest = "a", help = "Paths to files of interest", default="OFF")
	parser.add_option("--B", action = "store", type = "string", dest = "b", help = "Files to be hybridized", default="OFF")
	parser.add_option("--indexes", action = "store", type = "string", dest = "indexes", help = "Indexes for matrix construction...", default="OFF")
	parser.add_option("--values", action = "store", type = "string", dest = "values", help = "Values for matrix construction...", default="OFF")
	parser.add_option("--contexts", action = "store", type = "string", dest = "contexts", help = "What contexts of development should I track?", default="OFF")
	parser.add_option("--exclude", action="store", type="string", dest="exclude", help="Are there items that should be excluded?", default="")
	parser.add_option("--start", action = "store", type = "int", dest = "start", help = "Start development time for cell search", default=1)
	parser.add_option("--stop", action = "store", type = "int", dest = "stop", help = "End development time for cell search", default=250)
	parser.add_option("--step", action = "store", type = "int", dest = "step", help = "Step size", default=1)
	parser.add_option("--total", action = "store", type = "int", dest = "total", help = "Total simulations (indexes) for 'master' operations ", default=1000)
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
	coassociationspath = path_dict["coassociations"]
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
	
	# update peaks path:
	peakspath = peakspath + option.peaks + "/"
	
	# define input/output folders:
	expressionpath = cellspath + "expression/"
	correctionpath = cellspath + "correction/"
	lineagepath = cellspath + "lineage/"
	bindingpath = cellspath + "peaks/"
	overlappath = cellspath + "overlap/"
	cellsetpath = cellspath + "cellset/"
	genesetpath = cellspath + "geneset/"
	reportspath = cellspath + "reports/"
	comparepath = cellspath + "compare/"
	matrixpath = cellspath + "matrix/"
	tissuespath = cellspath + "tissues/"
	distancepath = cellspath + "distance/"
	hybridpath = cellspath + "hybrid/"
	dynamicspath = cellspath + "dynamics/"
	cubismpath = cellspath + "cubism/"
	timepath = cellspath + "time/"
	cellnotationspath = cellspath + "annotations/"
	general.pathGenerator(expressionpath)
	general.pathGenerator(correctionpath)
	general.pathGenerator(lineagepath)
	general.pathGenerator(bindingpath)
	general.pathGenerator(overlappath)
	general.pathGenerator(cellsetpath)
	general.pathGenerator(genesetpath)
	general.pathGenerator(reportspath)
	general.pathGenerator(comparepath)
	general.pathGenerator(matrixpath)
	general.pathGenerator(tissuespath)
	general.pathGenerator(distancepath)
	general.pathGenerator(timepath)
	general.pathGenerator(hybridpath)
	general.pathGenerator(dynamicspath)
	general.pathGenerator(cubismpath)
	general.pathGenerator(cellnotationspath)
	
	# generate expression flag:
	if option.measure == "max.expression":
		expression_flag = "maxCel_"
	elif option.measure == "avg.expression":
		expression_flag = "avgExp_"
	
	# check that the index range is coherent:
	if option.stop > option.total:
		print
		print "Error: Range exceeded! Stop index is larger than total."
		print
		return
	
	
	# master mode:
	if "master" in option.mode:
	
		# capture master mode:
		master, mode = option.mode.split(":")
		
		# prepare for qsub:
		bash_path = str(option.path + "/data/cells/runs/").replace("//","/")
		bash_base = "_".join([mode, option.peaks, option.name]) + "-M"
		qsub_base = "_".join([mode, option.peaks, option.name])
		general.pathGenerator(bash_path)
		if option.qsub != "OFF":
			qsub_header = open(qsubpath + option.qsub).read()
			qsub = True
		else:
			qsub_header = ""
			qsub = False
		if option.job == "QSUB":
			qsub_header = qsub_header.replace("qsubRunner", "qsub-" + qsub_base)
		elif option.job != "OFF":
			qsub_header = qsub_header.replace("qsubRunner", "qsub-" + option.job)
			bash_base = option.job + "-M"
		
		# update server path:
		if option.qsub != "OFF":
			option.path = serverPath(option.path)
			
		# prepare slave modules:
		m, steps, modules, commands, sequences, chunks, start, complete = 1, 0, list(), list(), list(), option.chunks, option.start, False
		for index in range(option.start, option.stop+1, option.step):
			run = "rn" + general.indexTag(index, option.total)
			steps += 1
			
			# cellular peak generation mode:
			if mode == "cell.peaks":
			
				command = "python <<CODEPATH>>mapCells.py --path <<PATH>> --organism <<ORGANISM>> --mode <<MODE>> --peaks <<PEAKS>> --start <<START>> --stop <<STOP>> --total <<TOTAL>> --expression <<EXPRESSION>> --collection <<COLLECTION>> --times <<TIMES>> --fraction <<FRACTION>> --minimum <<MINIMUM>> --name <<NAME>> --qsub <<QSUB>> --server <<SERVER>> --module <<MODULE>>"
				command = command.replace("<<CODEPATH>>", option.path + "/python/")
				command = command.replace("<<PATH>>", option.path)
				command = command.replace("<<ORGANISM>>", option.organism)
				command = command.replace("<<MODE>>", mode)
				command = command.replace("<<PEAKS>>", option.peaks)
				command = command.replace("<<START>>", str(index))
				command = command.replace("<<STOP>>", str(index))
				command = command.replace("<<TOTAL>>", str(option.total))
				command = command.replace("<<EXPRESSION>>", option.expression)
				command = command.replace("<<COLLECTION>>", option.collection)
				command = command.replace("<<TIMES>>", option.times)
				command = command.replace("<<FRACTION>>", str(option.fraction))
				command = command.replace("<<MINIMUM>>", str(option.minimum))
				command = command.replace("<<NAME>>", option.name + general.indexTag(index, option.total))
				command = command.replace("<<QSUB>>", option.qsub)
				command = command.replace("<<SERVER>>", option.server)
				command = command.replace("<<MODULE>>", "md" + str(m))
			
			# cellular peak generation mode:
			if mode == "cell.annotation":
			
				command = "python <<CODEPATH>>mapCells.py --path <<PATH>> --organism <<ORGANISM>> --mode <<MODE>> --peaks <<PEAKS>> --start <<START>> --stop <<STOP>> --total <<TOTAL>> --infile <<INFILE>> --collection <<COLLECTION>> --times <<TIMES>> --name <<NAME>> --qsub <<QSUB>> --server <<SERVER>> --module <<MODULE>>"
				command = command.replace("<<CODEPATH>>", option.path + "/python/")
				command = command.replace("<<PATH>>", option.path)
				command = command.replace("<<ORGANISM>>", option.organism)
				command = command.replace("<<MODE>>", mode)
				command = command.replace("<<PEAKS>>", option.peaks)
				command = command.replace("<<START>>", str(index))
				command = command.replace("<<STOP>>", str(index))
				command = command.replace("<<TOTAL>>", str(option.total))
				command = command.replace("<<INFILE>>", option.infile)
				command = command.replace("<<COLLECTION>>", option.collection)
				command = command.replace("<<TIMES>>", option.times)
				command = command.replace("<<NAME>>", option.name + general.indexTag(index, option.total) + option.nametag)
				command = command.replace("<<QSUB>>", option.qsub)
				command = command.replace("<<SERVER>>", option.server)
				command = command.replace("<<MODULE>>", "md" + str(m))
				
			# cellular overlap mode:
			if mode == "cell.overlap":
			
				command = "python <<CODEPATH>>mapCells.py --path <<PATH>> --organism <<ORGANISM>> --mode <<MODE>> --peaks <<PEAKS>> --start <<START>> --stop <<STOP>> --total <<TOTAL>> --expression <<EXPRESSION>> --collection <<COLLECTION>> --times <<TIMES>> --fraction <<FRACTION>> --minimum <<MINIMUM>> --extend <<EXTEND>> --name <<NAME>> --qsub <<QSUB>> --server <<SERVER>> --module <<MODULE>>"
				command = command.replace("<<CODEPATH>>", option.path + "/python/")
				command = command.replace("<<PATH>>", option.path)
				command = command.replace("<<ORGANISM>>", option.organism)
				command = command.replace("<<MODE>>", mode)
				command = command.replace("<<PEAKS>>", option.peaks)
				command = command.replace("<<START>>", str(index))
				command = command.replace("<<STOP>>", str(index))
				command = command.replace("<<TOTAL>>", str(option.total))
				command = command.replace("<<EXPRESSION>>", option.expression)
				command = command.replace("<<COLLECTION>>", option.collection + general.indexTag(index, option.total) + option.nametag)
				command = command.replace("<<TIMES>>", option.times)
				command = command.replace("<<NAME>>", option.name)
				command = command.replace("<<FRACTION>>", str(option.fraction))
				command = command.replace("<<MINIMUM>>", str(option.minimum))
				command = command.replace("<<EXTEND>>", str(option.extend))
				command = command.replace("<<QSUB>>", option.qsub)
				command = command.replace("<<SERVER>>", option.server)
				command = command.replace("<<MODULE>>", "md" + str(m))
			
			# coassociations hybrid mode:
			if mode == "cell.hybrid":
				
				collection = option.collection + general.indexTag(index, option.total) + option.nametag
				command = "python <<CODEPATH>>mapCells.py --path <<PATH>> --organism <<ORGANISM>> --mode <<MODE>> --A <<A>> --B <<B>> --indexes <<INDEXES>> --values <<VALUES>> --contexts <<CONTEXTS>>"
				command = command.replace("<<CODEPATH>>", option.path + "/python/")
				command = command.replace("<<PATH>>", option.path)
				command = command.replace("<<ORGANISM>>", option.organism)
				command = command.replace("<<MODE>>", mode)
				command = command.replace("<<A>>", option.a)
				command = command.replace("<<B>>", collection + "/mapcells_" + collection + "_matrix_overlap")
				command = command.replace("<<INDEXES>>", option.indexes)
				command = command.replace("<<VALUES>>", option.values)
				command = command.replace("<<CONTEXTS>>", option.contexts)
				
			# is it time to export a chunk?
			if index-start+option.step == chunks:
				
				# update start, modules, commands, and module count (m):
				start = index + option.step
				commands.append(command)
				modules.append(commands)
				commands = list()
				complete = True
				m += 1
			
			# store whether the most recent index/command has been stored:
			else:
				complete = False
		
		# update if there are additional commands:
		if not complete:
			commands.append(command)
			modules.append(commands)
			m += 1
		
		# launch commands:
		print
		print "Launching comparisons:", len(modules)
		#for module in modules:
		#	for command in module:
		#		print command
		runCommands(modules, threads=option.threads, mode="module.run", run_mode="verbose", run_path=bash_path, run_base=bash_base, record=True, qsub_header=qsub_header, qsub=qsub)
		print "Analyses performed:", len(modules)
		print

	
	# filter cells :
	elif option.mode == "filter":
		
		# load cells to filter:
		filterCells = open(path_dict[option.source] + option.target).read().strip().split("\n")
		
		# generate output file:
		f_output = open(path_dict[option.source] + option.name, "w")
		
		# process input lines:
		f, k = 0, 0
		inlines = open(path_dict[option.source] + option.infile).readlines()
		for inline in inlines:
			process = True 
			items = inline.strip().split(",")
			for item in items:
				if item in filterCells:
					process = False
					f += 1
			if process:
				print >>f_output, inline.strip()
				k += 1	
		
		print
		print "Input lines:", len(inlines)
		print "Output lines:", k, "(" + str(f) + " filtered)"
		print
		
		# close output:
		f_output.close()
	
	
	# simplify cell annotations :
	elif option.mode == "simply":
		
		# generate output file:
		f_output = open(path_dict[option.source] + option.name, "w")
		
		# process input lines:
		f, k = 0, 0
		inlines = open(path_dict[option.source] + option.infile).read().strip().replace("\r","\n").split("\n")
		for inline in inlines:
			process = True
			if "cell_mapping" in option.infile:
				regExp, original, updated = inline.strip().split(",")
				if updated == "":
					annotation = str(original)
				else:
					annotation = str(updated)
				print >>f_output, ",".join([regExp,annotation])
				k += 1	
		
		print
		print "Input lines:", len(inlines)
		print "Output lines:", k, "(" + str(f) + " simplified)"
		print
		
		# close output:
		f_output.close()
		
	
	# robustness analysis mode:
	elif option.mode == "robust":
		
		import itertools
		
		print
		print "Loading input series data..."
		signalDict, replicateDict = dict(), dict()
		inlines = open(extraspath + option.infile).read().replace("\r","\n").split("\n")
		columnDict = dict()
		inline, index = inlines.pop(0), 0
		for column in inline.strip().split(","):
			columnDict[column] = index
			index += 1
		for inline in inlines:
			valueDict, initems = dict(), inline.strip().split(",")
			if initems != [""]:
				for column in columnDict:
					valueDict[column] = initems[columnDict[column]]
				gene, series, cell, value = valueDict["Gene"], valueDict["Series"], valueDict["Cell"], valueDict["Express"]
				if not gene in signalDict:
					signalDict[gene] = dict()
				if not cell in signalDict[gene]:
					signalDict[gene][cell] = dict()
				signalDict[gene][cell][series] = value
				if not gene in replicateDict:
					replicateDict[gene] = list()
				replicateDict[gene].append(series)
				replicateDict[gene] = sorted(list(set(replicateDict[gene])))
		
		# define output file:
		f_output = open(expressionpath + "mapcells_" + option.mode + "_" + option.infile.replace(".csv",".txt"), "w")
		s_output = open(expressionpath + "mapcells_" + option.mode + "_" + option.infile.replace(".csv",".sum"), "w")
		print >>f_output, "\t".join(["gene","series.count","i","j","cells","pearson.correlation","pearson.pvalue"])
		print >>s_output, "\t".join(["series.count", "gene.count"])
		
		print "Scoring replicate correlations .."
		countDict = dict()
		for gene in signalDict:
			if not len(replicateDict[gene]) in countDict:
				countDict[len(replicateDict[gene])] = list()
			countDict[len(replicateDict[gene])].append(gene)
			if len(replicateDict[gene]) > 1:
				#print gene, len(replicateDict[gene])
				for (i, j) in itertools.combinations(replicateDict[gene], 2):
					iValues, jValues = list(), list()
					for cell in signalDict[gene]:
						if i in signalDict[gene][cell] and j in signalDict[gene][cell]:
							iValues.append(float(signalDict[gene][cell][i]))
							jValues.append(float(signalDict[gene][cell][j]))
					correlation, corPvalue = pearsonr(iValues, jValues)
					output = [gene, len(replicateDict[gene]), i, j, len(iValues), correlation, corPvalue]
					print >>f_output, "\t".join(map(str, output))
				#pdb.set_trace()	
		
		for count in sorted(countDict.keys()):
			print >>s_output, "\t".join(map(str, [count, len(countDict[count])]))
				
		# close output file:
		f_output.close()
		s_output.close()
		print


	# fillin mode:
	elif option.mode == "fillin":
		
		print
		print "Loading annotation information..."
		annotationDict = general.build2(extraspath + option.infile, id_column="lineage", split=",")
		
		print "Checking parental annotation..."
		missingCells = list()
		for cell in annotationDict:
			parent = cell[:len(cell)-1]
			if not parent in annotationDict:
				if not parent in missingCells:
					missingCells.append(parent)
					print parent, cell
		print
		
	
	# import mode:
	elif option.mode == "import":
		
		# Cell annotations are cell-type and tissue-type (in the new Murray version):
		# specificDict:  cell > cell-type
		# generalDict: cell > tissue-type
		
		# construct tissue dictionary (if necessary):
		if option.tissues != "OFF":
			
			print
			print "Loading general and specific tissue information..."
			specificDict = general.build2(extraspath + option.tissues, i="lineage", x="cell", mode="values", split=",")
			specificTotal = specificDict.values()
			
			generalDict = general.build2(extraspath + option.tissues, i="lineage", x="tissue", mode="values", split=",")
			generalTotal = generalDict.values()
			
			print "Generating tissue classes..."
			classification = {
				"rectal" : "excretory",
				"na" : "other"
				}
			
			classDict, classTotal, classMissing = dict(), list(), 0
			for cell in generalDict:
				generalTissue = generalDict[cell]
				generalHits, classHits = list(), list()
				if generalTissue == "g":
					classTissue = "neuron/glial"
					generalHits.append(generalTissue)
					classHits.append(classTissue)
				else:
					for classTag in classification:
						if classTag in generalTissue:
							classTissue = classification[classTag]
							generalHits.append(generalTissue)
							classHits.append(classTissue)
				generalHits, classHits = list(set(generalHits)), list(set(classHits))
				
				#print generalTissue, ":", ", ".join(classHits)
				if len(classHits) > 1:
					classTissue = "mixed"
				elif len(classHits) == 1:
					classTissue = classHits[0]
				elif len(classHits) == 0:
					classTissue = generalTissue
					classMissing += 1
				classDict[cell] = classTissue
				classTotal.append(classTissue)
			classTotal = sorted(list(set(classTotal)))
			
			print
			print "Specific tissue terms:", len(set(specificDict.values()))
			print "General tissue terms:", len(set(generalDict.values()))
			generalCounts = dict()
			for cell in generalDict:
				generalTissue = generalDict[cell]
				if not generalTissue in generalCounts:
					generalCounts[generalTissue] = 0
				generalCounts[generalTissue] += 1
			generalTissues = general.valuesort(generalCounts)
			generalTissues.reverse()
			for generalTissue in generalTissues:
				print "\t" + generalTissue, ":", generalCounts[generalTissue]
				
			print
			print "Class tissue terms:", len(set(classDict.values()))
			classCounts = dict()
			for cell in classDict:
				classTissue = classDict[cell]
				if not classTissue in classCounts:
					classCounts[classTissue] = 0
				classCounts[classTissue] += 1
			classTissues = general.valuesort(classCounts)
			classTissues.reverse()
			for classTissue in classTissues:
				print "\t" + classTissue, ":", classCounts[classTissue]
			#pdb.set_trace()
			
		# prepare expression matrixes:
		series2cell_dict, gene2cell_dict, cell2gene_dict, gene2cell_list, allCells = dict(), dict(), dict(), dict(), list()
		
		# load expression data per series:
		print
		print "Loading cellular-expression data..."
		inlines = open(extraspath + option.infile).read().replace("\r","\n").split("\n")
		inheader = inlines.pop(0)
		for inline in inlines:
			if not inline == "":
				series, cell, gene, expression = inline.strip().split(",")
				gene = gene.upper()
				if not gene in option.exclude.split(","):
					if not cell in cell2gene_dict:
						cell2gene_dict[cell] = dict()
					if not gene in cell2gene_dict[cell]:
						cell2gene_dict[cell][gene] = dict()
					if not gene in gene2cell_dict:
						gene2cell_dict[gene] = dict()
						gene2cell_list[gene] = list()
					if not cell in gene2cell_dict[gene]:
						gene2cell_dict[gene][cell] = dict()
					if not series in series2cell_dict:
						series2cell_dict[series] = dict()
					gene2cell_dict[gene][cell][series] = float(expression)
					cell2gene_dict[cell][gene][series] = float(expression)
					series2cell_dict[series][cell] = float(expression)
					if not cell in gene2cell_list[gene]:
						gene2cell_list[gene].append(cell)
					if not cell in allCells:
						allCells.append(cell)
				
		# store cell-parent relationships:
		print "Loading cell-parent relationships..."
		cell_dict, parent_dict, pedigreeCells = relationshipBuilder(pedigreefile=option.pedigree, path=extraspath, mechanism="simple")
		
		# construct tissue dictionary (if necessary):
		if option.tissues != "OFF":
			
			print
			print "Expanding cell tissue information..."
			matchDict = { "specific":dict(), "general":dict(), "class":dict() }
			matchExpansion, matchTotal, matchMissing = list(), 0, 0
			for cell in pedigreeCells:
				
				if cell in generalDict and generalDict[cell] != "na":
					matchDict["specific"][cell] = specificDict[cell]
					matchDict["general"][cell] = generalDict[cell]
					matchDict["class"][cell] = classDict[cell]
				
				else:
					
					# find most closely-related, annotated cell (and use its associated tissue annotation):
					distanceDict = dict()
					queryDict, matchTissues = dict(), list(), 
					ancestorCells, descendantCells, matchCells, queryCells = list(), list(), list(), list()
					for queryCell in generalDict:
						relative = False
						if cell == queryCell[:len(cell)]:
							descendantCells.append(queryCell)
							relative = True
						if queryCell == cell[:len(queryCell)]:
							ancestorCells.append(queryCell)
							relative = True
						if relative:
							distance = abs(len(cell)-len(queryCell))
							if not distance in distanceDict:
								distanceDict[distance] = list()
							distanceDict[distance].append(queryCell)
					
					# determine which cells to obtain the annotations from:
					if descendantCells != list():
						queryCells = descendantCells
					else:
						queryCells = descendantCells + ancestorCells
					
					# find and weigh, most-related tissues:
					specificMatch, generalMatch, classMatch = dict(), dict(), dict()
					for distance in sorted(distanceDict.keys()):
						if distance != 0:
							for distanceCell in distanceDict[distance]:
								if distanceCell in queryCells:
									specificTissue = specificDict[distanceCell]
									generalTissue = generalDict[distanceCell]
									classTissue = classDict[distanceCell]
									if not specificTissue in specificMatch:
										specificMatch[specificTissue] = 0
									if not generalTissue in generalMatch:
										generalMatch[generalTissue] = 0
									if not classTissue in classMatch:
										classMatch[classTissue] = 0
									specificMatch[specificTissue] += float(1)/distance
									generalMatch[generalTissue] += float(1)/distance
									classMatch[classTissue] += float(1)/distance
						
					# Note: This section controls whether tissue annotations are obtained from
					# all related cells (parents and ancestors) or just subsets of these...
					
					""" define a function that returns the highest-likelihood tissue """
					def matchFunction(cell, matchDict, queryCells, verbose="OFF"):
						matchTissues = general.valuesort(matchDict)
						matchTissues.reverse()
						
						printFlag = False
						if len(matchTissues) > 1 and verbose == "ON":
							printFlag = True
							print cell, len(matchTissues), matchTissues, queryCells
							for matchTissue in matchTissues:
								print matchTissue, ":", matchDict[matchTissue]
						
						# Filter tissues associated with father/daughter cells:
						if len(matchTissues) > 1:
							matchTissues = general.clean(matchTissues, "death")
						if len(matchTissues) > 1:
							matchTissues = general.clean(matchTissues, "other")
							
						# Generate and store specific tissue label for cell:
						if len(matchTissues) == 0:
							matchTissue = "other"
						else:
							matchTissue = matchTissues[0]
						
						if printFlag and verbose == "ON":
							print ">", matchTissue
							print
					
						# return highest likelihood tissue match and ranked tissues:	
						return matchTissue, matchTissues
					
					# assign highest-scoring tissue types: 
					#specificDict[cell], specificTissues = matchFunction(cell, specificMatch, queryCells, verbose="OFF")
					#generalDict[cell], generalTissues = matchFunction(cell, generalMatch, queryCells, verbose="OFF")
					#classDict[cell], classTissues = matchFunction(cell, classMatch, queryCells, verbose="OFF")
					matchDict["specific"][cell], specificMatches = matchFunction(cell, specificMatch, queryCells, verbose="OFF")
					matchDict["general"][cell], generalMatches = matchFunction(cell, generalMatch, queryCells, verbose="OFF")
					matchDict["class"][cell], classMatches = matchFunction(cell, classMatch, queryCells, verbose="OFF")
					
					# update tissue counts:
					matchTotal += 1
					if matchDict["class"][cell] == "na":
						matchMissing += 1
					
					# Update/expand cell-tissue dictionary:
					matchTissue = matchDict["specific"][cell]
					if not matchTissue in matchExpansion:
						matchExpansion.append(matchTissue)
			
			# record counts for each type of tissue:
			specificCounts, generalCounts, classCounts = dict(), dict(), dict()
			for cell in specificDict:
				specificTissue = specificDict[cell]
				generalTissue = generalDict[cell]
				classTissue = classDict[cell]
				if not specificTissue in specificCounts:
					specificCounts[specificTissue] = 0
				specificCounts[specificTissue] += 1
				if not generalTissue in generalCounts:
					generalCounts[generalTissue] = 0
				generalCounts[generalTissue] += 1
				if not classTissue in classCounts:
					classCounts[classTissue] = 0
				classCounts[classTissue] += 1
				
			#print
			#print "Specific tissue terms:", len(set(specificDict.values()))
			#specificTissues = general.valuesort(specificCounts)
			#specificTissues.reverse()
			#for specificTissue in specificTissues:
			#	print "\t" + specificTissue, ":", specificCounts[specificTissue]
			
			print
			print "General tissue terms:", len(set(generalDict.values()))
			generalTissues = general.valuesort(generalCounts)
			generalTissues.reverse()
			for generalTissue in generalTissues:
				print "\t" + generalTissue, ":", generalCounts[generalTissue]
			
			print
			print "Class tissue terms:", len(set(classDict.values()))
			classTissues = general.valuesort(classCounts)
			classTissues.reverse()
			for classTissue in classTissues:
				print "\t" + classTissue, ":", classCounts[classTissue]
			
			print
			print "Tissue information expanded by:", len(matchExpansion)
			print "Tissue information expansion terms:", ", ".join(list(sorted(matchExpansion)))
			#pdb.set_trace()

		# calculate unique expression values for each gene/cell combination:
		print
		print "Generating per gene/cell expression values..."
		matrix, expression, expressing = dict(), dict(), dict()
		for gene in gene2cell_dict:
			for cell in gene2cell_list[gene]:
				values, maxSeries, maxValue = list(), "NA", 0
				for series in gene2cell_dict[gene][cell]:
					values.append(gene2cell_dict[gene][cell][series])
					if gene2cell_dict[gene][cell][series] >= maxValue:
						maxSeries, maxValue = series, gene2cell_dict[gene][cell][series]
				if not gene in matrix:
					matrix[gene] = dict()
					expression[gene] = dict()
				matrix[gene][cell] = [max(values), numpy.mean(values), numpy.median(values), numpy.std(values), len(gene2cell_dict[gene][cell]), ",".join(sorted(gene2cell_dict[gene][cell].keys())), maxSeries]
				if option.measure == "max.expression":
					expression[gene][cell] = max(values)
				elif option.measure == "avg.expression":
					expression[gene][cell] = numpy.mean(values)
		
		# calculate expression peaks...
		print "Generating per gene/cell expression statistics..."
		for gene in matrix:
			
			# find peak expression:
			peakCell, peakValue = "", 0
			for cell in matrix[gene]:
				maxValue, meanValue, medianValue, stdValue, seriesCount, seriesIDs, maxSeries = matrix[gene][cell]
				cellValue = expression[gene][cell]
				if cellValue > peakValue:
					peakCell, peakValue = cell, cellValue
					
			# calculate fractional expression, cell ranks, and add cells expressing the protein (above cutoff):
			cellRanks = general.valuesort(expression[gene])
			cellRanks.reverse()
			for cell in matrix[gene]:
				maxValue, meanValue, medianValue, stdValue, seriesCount, seriesIDs, maxSeries = matrix[gene][cell]
				cellValue = expression[gene][cell]
				fracValue = float(cellValue)/peakValue
				cellRank = cellRanks.index(cell) + 1
				if not gene in expressing:
					expressing[gene] = list()
				if fracValue >= option.fraction and cellValue >= option.minimum:
					expressing[gene].append(cell)
				matrix[gene][cell] = [cellValue, peakValue, fracValue, cellRank, maxValue, meanValue, medianValue, stdValue, seriesCount, seriesIDs, maxSeries]
		
		# define the ascendants cutoff:
		print
		print "Defining minimum ascendants across experiments..."
		cutAscendants = 0
		for gene in matrix:
			minAscendants, maxAscendants = 1000, 0
			for cell in matrix[gene]:
				ascendants = ascendantsCollector(cell, parent_dict, cell_dict, ascendants=list())
				if len(ascendants) < minAscendants:
					minAscendants = len(ascendants)
					minCell = cell
				if len(ascendants) > maxAscendants:
					maxAscendants = len(ascendants)
					maxCell = cell
			if minAscendants > cutAscendants:
				cutAscendants = minAscendants
		
		# define the set of cells tracked in target experiments:
		print "Defining cells focused: strict list of cells assayed in target experiments..."
		focusedCells = list()
		for gene in option.target.split(","):
			if focusedCells == list():
				focusedCells = gene2cell_list[gene]
			else:
				focusedCells = set(focusedCells).intersection(set(gene2cell_list[gene]))
		
		# define the set of cells tracked in all experiments:
		print "Defining cells tracked: strict list of cells assayed in all experiments..."
		trackedCells = list()
		for gene in gene2cell_dict:
			if trackedCells == list():
				trackedCells = gene2cell_list[gene]
			else:
				trackedCells = set(trackedCells).intersection(set(gene2cell_list[gene]))
		
		# define the set of ancestor or tracked cells:
		print "Defining cells started: parent-inclusive list of cells tracked in all experiments..."
		startedCells = list()
		for cell in pedigreeCells:
			ascendants = ascendantsCollector(cell, parent_dict, cell_dict, ascendants=list())
			if cell in trackedCells or len(ascendants) < int(option.ascendants):
				startedCells.append(cell)
				#if cell == "ABalaaaal":
				#	print cell, specificDict[cell], generalDict[cell], classDict[cell]
				#	pdb.set_trace()
		print "Ascendants cutoff:", cutAscendants
		
		# define output files:
		assayedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_assayed"
		startedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_started"
		trackedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_tracked"
		focusedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_focused"
		summaryfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_summary"
		tissuesfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_tissues"
		
		# define cellular expression header:
		expressionHeader = ["cell", "cell.name", "gene", "cell.expression", "peak.expression", "fraction.expression", "normal.expression", "rank", "max.expression", "avg.expression", "med.expression", "std.expression", "cells.expressing", "cells.count", "series.count", "time.series", "max.series", "specific.tissue", "general.tissue", "class.tissue", "match.tissue"]
		reportHeader = ["gene", "cells.expressing", "cells.assayed", "cells.tracked", "cell.expression", "peak.expression", "fraction.expression", "series.count", "time.series", "max.series"]
		tissueHeader = ["cell", "specific.tissue", "general.tissue", "class.tissue", "match.tissue"]
		
		# create output files:
		a_output = open(assayedfile, "w")
		s_output = open(startedfile, "w")
		t_output = open(trackedfile, "w")
		f_output = open(focusedfile, "w")
		r_output = open(summaryfile, "w")
		x_output = open(tissuesfile, "w")
		print >>a_output, "\t".join(expressionHeader)
		print >>s_output, "\t".join(expressionHeader)
		print >>t_output, "\t".join(expressionHeader)
		print >>f_output, "\t".join(expressionHeader)
		print >>r_output, "\t".join(reportHeader)
		print >>x_output, "\t".join(tissueHeader)
		
		# generate set-normalization values:
		maxAssayed, maxStarted, maxTracked, maxFocused = dict(), dict(), dict(), dict()
		for gene in sorted(matrix.keys()):
			cellsStarted = startedCells
			cellsAssayed = matrix[gene].keys()
			cellsTracked = trackedCells
			cellsFocused = focusedCells
			peakAssayed, peakStarted, peakTracked, peakFocused = 0, 0, 0, 0
			for cell in sorted(matrix[gene].keys()):
				cellValue, peakValue, fracValue, cellRank, maxValue, meanValue, medianValue, stdValue, seriesCount, seriesIDs, maxSeries = matrix[gene][cell]
				if cell in cellsAssayed and cellValue > peakAssayed:
					peakAssayed = cellValue
				if cell in cellsStarted and cellValue > peakStarted:
					peakStarted = cellValue
				if cell in cellsTracked and cellValue > peakTracked:
					peakTracked = cellValue
				if cell in cellsFocused and cellValue > peakFocused:
					peakFocused = cellValue
			maxAssayed[gene] = peakAssayed 
			maxStarted[gene] = peakStarted
			maxTracked[gene] = peakTracked
			maxFocused[gene] = peakFocused
		
		# export expression data:
		print "Exporting expression data..."
		for gene in sorted(matrix.keys()):
			cellsStarted = len(startedCells)
			cellsAssayed = len(matrix[gene].keys())
			cellsTracked = len(trackedCells)
			cellsFocused = len(focusedCells)
			cellsExpressingAssayed = len(set(expressing[gene]).intersection(set(matrix[gene].keys())))
			cellsExpressingTracked = len(set(expressing[gene]).intersection(set(trackedCells)))
			cellsExpressingFocused = len(set(expressing[gene]).intersection(set(focusedCells)))
			cellValues, fracValues = list(), list()
			for cell in sorted(matrix[gene].keys()):
				if option.tissues == "OFF" or not cell in specificDict:
					specificTissue = "*"
					generalTissue = "*"
					classTissue = "*"
					matchTissue = "*"
				else:
					specificTissue = specificDict[cell]
					generalTissue = generalDict[cell]
					classTissue = classDict[cell]
					matchTissue = matchDict["class"][cell]
				cellValue, peakValue, fracValue, cellRank, maxValue, meanValue, medianValue, stdValue, seriesCount, seriesIDs, maxSeries = matrix[gene][cell]
				print >>a_output, "\t".join(map(str, [cell, cell, gene, cellValue, peakValue, fracValue, float(cellValue)/maxAssayed[gene], cellRank, maxValue, meanValue, medianValue, stdValue, cellsExpressingAssayed, cellsAssayed, seriesCount, seriesIDs, maxSeries, specificTissue, generalTissue, classTissue, matchTissue]))
				if cell in startedCells:
					print >>s_output, "\t".join(map(str, [cell, cell, gene, cellValue, peakValue, fracValue, float(cellValue)/maxStarted[gene], cellRank, maxValue, meanValue, medianValue, stdValue, cellsExpressingTracked, cellsStarted, seriesCount, seriesIDs, maxSeries, specificTissue, generalTissue, classTissue, matchTissue]))
				if cell in trackedCells:
					print >>t_output, "\t".join(map(str, [cell, cell, gene, cellValue, peakValue, fracValue, float(cellValue)/maxTracked[gene], cellRank, maxValue, meanValue, medianValue, stdValue, cellsExpressingTracked, cellsTracked, seriesCount, seriesIDs, maxSeries, specificTissue, generalTissue, classTissue, matchTissue]))
				if cell in focusedCells:
					print >>f_output, "\t".join(map(str, [cell, cell, gene, cellValue, peakValue, fracValue, float(cellValue)/maxFocused[gene], cellRank, maxValue, meanValue, medianValue, stdValue, cellsExpressingFocused, cellsFocused, seriesCount, seriesIDs, maxSeries, specificTissue, generalTissue, classTissue, matchTissue]))
				if fracValue >= option.fraction and cellValue >= option.minimum:
					cellValues.append(cellValue)
					fracValues.append(fracValue)
			print >>r_output, "\t".join(map(str, [gene, cellsExpressingTracked, cellsAssayed, cellsTracked, numpy.mean(cellValues), peakValue, numpy.mean(fracValues), seriesCount, seriesIDs, maxSeries]))
			
		# export tissue annotations:
		print "Exporting tissue annotation data..."
		print "Annotated cells:", len(specificDict)
		for cell in sorted(specificDict.keys()):
			specificTissue = specificDict[cell]
			generalTissue = generalDict[cell]
			classTissue = classDict[cell]
			if cell in matchDict["class"]:
				matchTissue = matchDict["class"][cell]
			else:
				matchTissue = str(classTissue)
			print >>x_output, "\t".join([cell, specificTissue, generalTissue, classTissue, matchTissue])
	
		# close output:
		a_output.close()
		s_output.close()
		t_output.close()
		f_output.close()
		r_output.close()
		x_output.close()
		
		print
		print "Focused cells:", len(focusedCells)
		print "Tracked cells:", len(trackedCells)
		print "Started cells:", len(startedCells)
		print
		
	
	# inherit expression mode:
	elif option.mode == "inherit":
		
		# define input files:
		assayedinput = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_assayed"
		startedinput = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_started"
		trackedinput = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_tracked"
		focusedinput = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_focused"
		summaryinput = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_summary"
		tissuesinput = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_tissues"
		
		# define output files:
		assayedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_inassay"
		startedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_instart"
		trackedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_intrack"
		focusedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_infocus"
		inheritfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_inherit"
		inleafsfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_inleafs"
		maximalfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_maximal"
		mxleafsfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_mxleafs"
		
		# define cellular expression header:
		expressionHeader = ["cell", "cell.name", "gene", "cell.expression", "peak.expression", "fraction.expression", "normal.expression", "rank", "max.expression", "avg.expression", "med.expression", "std.expression", "cells.expressing", "cells.count", "series.count", "time.series", "max.series", "specific.tissue", "general.tissue", "class.tissue", "match.tissue"]
		reportHeader = ["gene", "cells.expressing", "cells.assayed", "cells.tracked", "cell.expression", "peak.expression", "fraction.expression", "series.count", "time.series", "max.series"]
		tissueHeader = ["cell", "specific.tissue", "general.tissue", "class.tissue", "match.tissue"]
		
		# create output files:
		a_output = open(assayedfile, "w")
		s_output = open(startedfile, "w")
		t_output = open(trackedfile, "w")
		f_output = open(focusedfile, "w")
		i_output = open(inheritfile, "w")
		l_output = open(inleafsfile, "w")
		m_output = open(maximalfile, "w")
		p_output = open(mxleafsfile, "w")
		print >>a_output, "\t".join(expressionHeader + ["inherited"])
		print >>s_output, "\t".join(expressionHeader + ["inherited"])
		print >>t_output, "\t".join(expressionHeader + ["inherited"])
		print >>f_output, "\t".join(expressionHeader + ["inherited"])
		print >>i_output, "\t".join(expressionHeader + ["inherited"])
		print >>l_output, "\t".join(expressionHeader + ["inherited"])
		print >>m_output, "\t".join(expressionHeader + ["inherited"])
		print >>p_output, "\t".join(expressionHeader + ["inherited"])
		
		# load terminal leaf cells:
		print
		print "Loading terminal cells..."
		inleafsCells = general.build2(extraspath + option.mapping, i="cell", x="cell.name", mode="values", skip=True)
		
		# store cell-parent relationships:
		print "Loading cell-parent relationships..."
		cell_dict, parent_dict, pedigreeCells = relationshipBuilder(pedigreefile=option.pedigree, path=extraspath, mechanism="simple")
		
		# loading tissue annotation data...
		print "Loading tissue annotation data..."
		tissuesAnnotation = general.build2(tissuesinput, id_column="cell", mode="table")
		
		# load expression data:
		print "Loading expression data..."
		assayedExpression = general.build2(assayedinput, id_complex=["gene","cell"], mode="table", separator=":")
		assayedMatrix = general.build2(assayedinput, i="gene", j="cell", x="cell.expression", mode="matrix")
		assayedCells = general.build2(assayedinput, i="cell", x="cell.name", mode="values", skip=True)
		startedCells = general.build2(startedinput, i="cell", x="cell.name", mode="values", skip=True)
		trackedCells = general.build2(trackedinput, i="cell", x="cell.name", mode="values", skip=True)
		focusedCells = general.build2(focusedinput, i="cell", x="cell.name", mode="values", skip=True)
		
		# define cellular space:
		print "Defining inheritance cells..."
		inheritCells = list()
		for inleafsCell in inleafsCells:
			inheritCells += ascendantsCollector(inleafsCell, parent_dict, cell_dict, ascendants=list())
		inheritCells = sorted(list(set(inheritCells)))
		
		# load header dictionary:
		hd = general.build_header_dict(assayedinput)
		header = general.valuesort(hd)
		
		# inherit peak expression from ancestors:
		print "Inheriting expression from ancestors..."
		inheritExpression, maximalExpression = dict(), dict()
		for gene in sorted(assayedMatrix.keys()):
			inheritExpression[gene] = dict()
			maximalExpression[gene] = dict()
			for inheritCell in inheritCells:
				ascendantCells, ascendantExpression = list(), dict()
				ascendants = ascendantsCollector(inheritCell, parent_dict, cell_dict, ascendants=list(), sort=False)
				
				#print inheritCell, ascendants
				if len(set(ascendants)) != len(ascendants):
					print "oh, oh: not a set!"
					pdb.set_trace()
				
				for ascendantCell in ascendants + [inheritCell]:
					if ascendantCell in assayedMatrix[gene]:
						ascendantExpression[ascendantCell] = float(assayedMatrix[gene][ascendantCell])
						ascendantCells.append(ascendantCell)
				
				if ascendantExpression != dict():
					
					# get inheritance cells for maximal expression and for last ancestor expression:
					maximalCells = general.valuesort(ascendantExpression)
					maximalCells.reverse()
					maximalCell = maximalCells[0]
					ascendantCell = ascendantCells[0]
					
					# store values for last ancestor expression:
					inheritExpression[gene][inheritCell] = dict(assayedExpression[gene + ":" + ascendantCell])
					inheritExpression[gene][inheritCell]["cell"] = str(inheritCell)
					inheritExpression[gene][inheritCell]["cell.name"] = str(inheritCell)
					inheritExpression[gene][inheritCell]["specific.tissue"] = tissuesAnnotation[inheritCell]["specific.tissue"]
					inheritExpression[gene][inheritCell]["general.tissue"] = tissuesAnnotation[inheritCell]["general.tissue"]
					inheritExpression[gene][inheritCell]["class.tissue"] = tissuesAnnotation[inheritCell]["class.tissue"]
					inheritExpression[gene][inheritCell]["match.tissue"] = tissuesAnnotation[inheritCell]["match.tissue"]
					inheritExpression[gene][inheritCell]["inherited"] = ascendantCell
					#if inheritCell != inheritExpression[gene][inheritCell]["cell"]:
					#	print cell, inheritExpression[gene][inheritCell]["cell"], 1
					#	pdb.set_trace()
					
					# store values for maximal ancestor expression:
					maximalExpression[gene][inheritCell] = dict(assayedExpression[gene + ":" + maximalCell])
					maximalExpression[gene][inheritCell]["cell"] = str(inheritCell)
					maximalExpression[gene][inheritCell]["cell.name"] = str(inheritCell)
					maximalExpression[gene][inheritCell]["specific.tissue"] = tissuesAnnotation[inheritCell]["specific.tissue"]
					maximalExpression[gene][inheritCell]["general.tissue"] = tissuesAnnotation[inheritCell]["general.tissue"]
					maximalExpression[gene][inheritCell]["class.tissue"] = tissuesAnnotation[inheritCell]["class.tissue"]
					maximalExpression[gene][inheritCell]["match.tissue"] = tissuesAnnotation[inheritCell]["match.tissue"]
					maximalExpression[gene][inheritCell]["inherited"] = ascendantCell
					
		
		# export inherited signals:
		print "Exporting inherited expression values..."
		for gene in sorted(inheritExpression):
			for cell in sorted(inheritExpression[gene].keys()):
				#if cell != inheritExpression[gene][cell]["cell"]:
				#	print cell, inheritExpression[gene][cell]["cell"], 2
				#	pdb.set_trace()
				output = list()
				for column in header + ["inherited"]:
					output.append(inheritExpression[gene][cell][column])
				if cell in assayedCells:
					print >>a_output, "\t".join(map(str, output))
				if cell in startedCells:
					print >>s_output, "\t".join(map(str, output))
				if cell in trackedCells:
					print >>t_output, "\t".join(map(str, output))
				if cell in focusedCells:
					print >>f_output, "\t".join(map(str, output))
				if cell in inheritCells:
					print >>i_output, "\t".join(map(str, output))
				if cell in inleafsCells:
					print >>l_output, "\t".join(map(str, output))
				#print "\t".join(map(str, output))
				#pdb.set_trace()
				
		# export inherited signals:
		print "Exporting maximal expression values..."
		for gene in sorted(maximalExpression):
			for cell in sorted(maximalExpression[gene].keys()):
				output = list()
				for column in header + ["inherited"]:
					output.append(maximalExpression[gene][cell][column])
				if cell in inheritCells:
					print >>m_output, "\t".join(map(str, output))
				if cell in inleafsCells:
					print >>p_output, "\t".join(map(str, output))
				#print "\t".join(map(str, output))
				#pdb.set_trace()
				
		print
		print "Total inherited cells:", len(inheritCells)
		print "Terminal (leaf) cells:",  len(inleafsCells)
		
		# close output files:
		a_output.close()
		s_output.close()
		t_output.close()
		f_output.close()
		i_output.close()
		l_output.close()
		m_output.close()
		p_output.close()
		print
		
		#k = inheritExpression.keys()[0]
		#print k
		#print inheritExpression[k][inleafsCell]
		#pdb.set_trace()
		
	
	# correct expression mode (detect outliers):
	elif option.mode == "correct":
		
		# load quantile functions
		from quantile import Quantile
		
		# define input files:
		startedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_started"
		trackedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_tracked"
		assayedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_assayed"
		summaryfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_summary"
		
		# load assayed expression data:
		print 
		print "Loading expression data..."
		expressionDict = general.build2(assayedfile, i="gene", j="cell", x="cell.expression", mode="matrix")
		
		# prepare to sort genes by quantile expression:
		print "Sorting genes by expression..."
		medianDict, quantDict = dict(), dict()
		for gene in expressionDict:
			values = map(float, expressionDict[gene].values())
			medianDict[gene] = numpy.median(values)
			quantDict[gene] = Quantile(values, 0.99)
		quantRanks = general.valuesort(quantDict)
		quantRanks.reverse()
		
		# store median rankings:
		rankDict = dict()
		medianRanks = general.valuesort(medianDict)
		medianRanks.reverse()
		k = 1
		for gene in medianRanks:
			rankDict[gene] = k
			k += 1
		
		# generate testing path:
		testingpath = correctionpath + "testing/"
		general.pathGenerator(testingpath)
		
		# Perform Gaussian Mixture Modeling (GMM):
		print "Performing GMM modeling..."
		gmmDict = dict()
		k = 1
		for gene in expressionDict:
			signals = map(int, map(float, expressionDict[gene].values()))
			signals = [1 if (x == 0) else x for x in signals]
			testingfile = testingpath + "mapCells-gmm_" + expression_flag + option.name + "_" + "temp"
			resultsfile = testingpath + "mapCells-gmm_" + expression_flag + option.name + "_" + gene
			#f_output = open(testingfile, "w")
			#print >>f_output, "\n".join(["signal"] + map(str, signals))
			#f_output.close()
			#command = " ".join(["Rscript", "~/meTRN/scripts/mapCells-gmm.r", testingfile, resultsfile, option.limit, option.parameters])
			#os.system(command)
			#Rscript ~/meTRN/scripts/mapCells-pilot.r ~/Desktop/data.test ~/Desktop/data.output 1000
			if "mapCells-gmm_" + expression_flag + option.name + "_" + gene in os.listdir(testingpath):
				gmmDict[gene] = open(resultsfile).readlines()[1].strip().split(" ")[2]
		#os.system("rm -rf " + testingfile)
		
		# export expression signals:
		rankingfile = correctionpath + "mapcells_" + expression_flag + option.name + "_correction_ranking" # rank information file
		percentfile = correctionpath + "mapcells_" + expression_flag + option.name + "_correction_percent" # gene-cell data, percentile-ranked genes
		mediansfile = correctionpath + "mapcells_" + expression_flag + option.name + "_correction_medians" # gene-cell data, median-ranked genes
		
		# define output headers:
		correctHeader = "\t".join(["index", "gene", "cell", "signal", "zscore", "nscore", "lscore", "rank", "median", "mean", "stdev", "alpha", "delta", "sigma", "gamma"])
		rankingHeader = "\t".join(["gene", "quantile.rank", "median.rank", "median", "mean", "stdev", "alpha", "delta", "sigma", "gamma"])
		
		# gather outputs:
		print "Generating expression thresholds..."
		r_output = open(rankingfile, "w")
		print >>r_output, rankingHeader
		outputDict = dict()
		k = 1
		for gene in quantRanks:
			signals = map(float, expressionDict[gene].values())
			maximal = max(signals)
			
			# calculate expression cutoffs:
			alpha = float(maximal)/10
			delta = float(quantDict[gene])/10
			sigma = float(quantDict[gene])/10
			
			# detect GMM expression cutoff:
			if gene in gmmDict:
				gamma = int(gmmDict[gene])
			else:
				gamma = int(option.limit)
		
			# threshold expression cutoffs:
			if alpha < int(option.limit):
				alpha = int(option.limit)
			if delta < int(option.limit):
				delta = int(option.limit)
			if gamma < int(option.limit):
				gamma = int(option.limit)
			
			# calculate general stats:
			median = numpy.median(signals)
			mean = numpy.mean(signals)
			stdev = numpy.std(signals)
			logMean = numpy.log10(mean)
			logStDev = numpy.log10(stdev)
			
			# store/export data:
			print >>r_output,  "\t".join(map(str, [gene, k, rankDict[gene], median, mean, stdev, alpha, delta, sigma, gamma]))
			if not gene in outputDict:
				outputDict[gene] = dict()
			for cell in sorted(expressionDict[gene].keys()):
				signal = float(expressionDict[gene][cell])
				if signal < 1:
					signal = 1
				zscore = float(signal-mean)/stdev
				nscore = float(signal)/maximal
				lscore = float(numpy.log10(signal) - logMean)/logStDev
				outputDict[gene][cell] = "\t".join(map(str, [k, gene, cell, signal, zscore, nscore, lscore, rankDict[gene], median, mean, stdev, alpha, delta, sigma, gamma]))
			k += 1
		r_output.close()
		
		# export expression signals, percentile-ranked genes:
		print "Exporting percentile-ranked expression signals..."
		f_output = open(percentfile, "w")
		print >>f_output, correctHeader
		for gene in quantRanks:
			for cell in sorted(outputDict[gene]):
				print >>f_output, outputDict[gene][cell]
		f_output.close()
	
		# export expression signals, median-ranked genes:
		print "Exporting median-ranked expression signals..."
		f_output = open(mediansfile, "w")
		print >>f_output, correctHeader
		for gene in medianRanks:
			for cell in sorted(outputDict[gene]):
				print >>f_output, outputDict[gene][cell]
		f_output.close()
		print
	
	
	# check status mode:
	elif option.mode == "check.status":
		
		# define input files:
		startedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_started"
		trackedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_tracked"
		assayedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_assayed"
		summaryfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_summary"
		
		# scan peak files:
		print
		print "Scanning peak files:"
		hd = general.build_header_dict(summaryfile)
		k, peaks, peak_files = 0, 0, os.listdir(peakspath)
		for inline in open(summaryfile).readlines()[1:]:
			gene, found = inline.strip().split("\t")[hd["gene"]], list()
			for peak_file in peak_files:
				dataset = peak_file.split("_peaks.bed")[0].replace("POL2", "AMA-1")
				if gene + "_" in dataset:
					found.append(dataset)
					peaks += general.countLines(peakspath + peak_file, header="OFF")
					
			if found != list():
				print gene, ":", ", ".join(sorted(found))
				k += 1
		print
		print "Found factors:", k
		print "Peaks called:", peaks
		print
	
		# scan expression files:
		print
		print "Scanning expression data:"
		caught = list()
		hd = general.build_header_dict(assayedfile)
		for inline in open(assayedfile).readlines()[1:]:
			initems = inline.strip().split("\t")
			gene, timeSeries = initems[hd["gene"]], initems[hd["time.series"]]
			for timeSerie in timeSeries.split(","):
				if not gene.lower() in timeSerie:
					if not gene in caught:
						print gene, timeSeries
					caught.append(gene)
		print
		print "Mismatched genes:", len(caught)
		print
	
	
	# lineage distance mode:
	elif option.mode == "cell.distance":
	
		# build cell-expression matrix:
		print
		print "Loading cellular expression..."
		quantitation_matrix, expression_matrix, tracking_matrix, trackedCells = expressionBuilder(expressionfile=option.expression, path=expressionpath, cutoff=0, minimum=0, metric="fraction.expression")
		
		# store cell-parent relationships:
		print "Loading cell-parent relationships..."
		cell_dict, parent_dict, pedigreeCells = relationshipBuilder(pedigreefile=option.pedigree, path=extraspath, mechanism="simple")
		
		print "Pedigree cells:", len(pedigreeCells)
		print "Tracked cells:", len(trackedCells)
		print
		
		# define output files:
		signalsmatrixfile = str(option.expression + "_distance_signals")
		lineagematrixfile = str(option.expression + "_distance_lineage")
		combinematrixfile = str(option.expression + "_distance_combine")
		
		# build cell-cell expression correlation matrix
		if not signalsmatrixfile in os.listdir(distancepath) or option.overwrite == "ON":
			
			print "Calculating expression correlation matrix..."
			correlation_matrix, index = dict(), 1
			f_output = open(distancepath + signalsmatrixfile, "w")
			print >>f_output, "\t".join(["i", "j", "correlation", "correlation.pvalue", "correlation.adjusted.pvalue"])
			for aCell in sorted(trackedCells):
				print index, aCell
				for bCell in sorted(trackedCells):
					aValues, bValues = list(), list()
					for gene in sorted(quantitation_matrix.keys()):
						aValues.append(quantitation_matrix[gene][aCell])
						bValues.append(quantitation_matrix[gene][bCell])
					correlation, corPvalue = pearsonr(aValues, bValues)
					adjCorPvalue = corPvalue*len(trackedCells)*len(trackedCells)
					if adjCorPvalue > 1:
						adjCorPvalue = 1
					if not aCell in correlation_matrix:
						correlation_matrix[aCell] = dict()
					correlation_matrix[aCell][bCell] = [correlation, corPvalue, adjCorPvalue]
					print >>f_output, "\t".join(map(str, [aCell, bCell] + correlation_matrix[aCell][bCell]))
				index += 1
			f_output.close()
			print
			
		else:
			print "Loading expression correlation matrix..."
			correlation_matrix = general.build2(distancepath + signalsmatrixfile, i="i", j="j", x=["correlation","correlation.pvalue","correlation.adjusted.pvalue"], datatype="float", mode="matrix", header_dict="auto")
		
		# build lineage distance matrix:
		if not lineagematrixfile in os.listdir(distancepath) or option.overwrite == "ON":
			
			print "Calculating lineage distance matrix..."
			lineage_matrix, index = dict(), 1
			f_output = open(distancepath + lineagematrixfile, "w")
			print >>f_output, "\t".join(["i", "j", "distance", "parent"])
			for aCell in sorted(trackedCells):
				print index, aCell 
				for bCell in sorted(trackedCells):
					distance, ancestor = lineageDistance(aCell, bCell, parent_dict, cell_dict)
					if not aCell in lineage_matrix:
						lineage_matrix[aCell] = dict()
					lineage_matrix[aCell][bCell] = [distance, ancestor]
					print >>f_output, "\t".join(map(str, [aCell, bCell] + lineage_matrix[aCell][bCell]))
				index += 1
			f_output.close()
			print
			
		else:
			print "Loading lineage distance matrix..."
			lineage_matrix = general.build2(distancepath + lineagematrixfile, i="i", j="j", x=["distance","parent"], datatype="list", mode="matrix", header_dict="auto", listtypes=["int", "str"])
		
		#print correlation_matrix["ABal"]["ABal"]
		#print lineage_matrix["ABal"]["ABal"]
		#pdb.set_trace()
		
		# build expression distance matrix (as a function of fraction expression):
		print "Generating combined distance matrix (at fraction range):"
		f_output = open(distancepath + combinematrixfile, "w")
		print >>f_output, "\t".join(["i", "j", "minimal", "fraction", "distance", "parent", "expression.correlation", "expression.correlation.pvalue", "expression.correlation.adjusted.pvalue", "i.genes", "j.genes", "overlap", "total", "overlap.max", "overlap.sum", "pvalue", "adjusted.pvalue", "flag"])
		fraction_matrix, genes = dict(), sorted(tracking_matrix.keys())
		for minimal in [1500, 1750, 2000]:
			for fraction in general.drange(0.10, 0.50, 0.10):
				print "...", minimal, fraction
				fraction_matrix[fraction] = dict()
				
				# find genes expressed per cell (using fraction cutoff):
				cellular_matrix = dict()
				fraction_quantitation_matrix, fraction_expression_matrix, fraction_tracking_matrix, fraction_trackedCells = expressionBuilder(expressionfile=option.expression, path=expressionpath, cutoff=fraction, minimum=minimal, metric="fraction.expression")
				for gene in fraction_expression_matrix:
					for cell in fraction_expression_matrix[gene]:
						if not cell in cellular_matrix:
							cellular_matrix[cell] = list()
							cellular_matrix[cell].append(gene)
							
				# find multiple hypothesis adjustment factor:
				adjust = 0
				for aCell in sorted(fraction_trackedCells):
					for bCell in sorted(fraction_trackedCells):
						if aCell in cellular_matrix and bCell in cellular_matrix:
							adjust += 1
							
				# find gene expression overlap between cells:
				overlap_matrix = dict()
				universe = len(quantitation_matrix.keys())
				for aCell in sorted(trackedCells):
					for bCell in sorted(trackedCells):
						if aCell in cellular_matrix and bCell in cellular_matrix:
							aGenes = cellular_matrix[aCell]
							bGenes = cellular_matrix[bCell]
							union = set(aGenes).union(set(bGenes))
							overlap = set(aGenes).intersection(set(bGenes))
							
							maxOverlap = float(len(overlap))/min(len(aGenes), len(bGenes))
							sumOverlap = float(len(overlap))/len(union)
							
							# Hypergeometric paramters:
							m = len(aGenes) # number of white balls in urn
							n = universe - len(bGenes) # number of black balls in urn
							N = len(bGenes) # number of balls drawn from urn
							x = len(overlap) # number of white balls in drawn
							
							# If I pull out all balls with elephant tatoos (N), is the draw enriched in white balls?:
							pvalue = hyper.fishers(x, m+n, m, N, method="right")
							adjPvalue = hyper.limit(pvalue*adjust)
							
							# Store overlap and significance:
							if not aCell in overlap_matrix:
								overlap_matrix[aCell] = dict()
							overlap_matrix[aCell][bCell] = [len(aGenes), len(bGenes), len(overlap), universe, maxOverlap, sumOverlap, pvalue, adjPvalue]
			
				# generate combined distance output line:
				for aCell in sorted(trackedCells):
					for bCell in sorted(trackedCells):
					
						# load lineage distances:
						distance, ancestor = lineage_matrix[aCell][bCell]
						
						# load correlation distances:
						correlation, corPvalue, adjCorPvalue = correlation_matrix[aCell][bCell] 
						
						# load expresssion distances:
						if aCell in cellular_matrix and bCell in cellular_matrix:
							aGenes, bGenes, overlap, universe, maxOverlap, sumOverlap, pvalue, adjPvalue = overlap_matrix[aCell][bCell]
							madeFlag = "both.observed"
						elif aCell in cellular_matrix:
							aGenes, bGenes, overlap, universe, maxOverlap, sumOverlap, pvalue, adjPvalue = len(cellular_matrix[aCell]), 0, 0, len(trackedCells), 0, 0, 1, 1
							madeFlag = "only.observed"
						elif bCell in cellular_matrix:
							aGenes, bGenes, overlap, universe, maxOverlap, sumOverlap, pvalue, adjPvalue = 0, len(cellular_matrix[bCell]), 0, len(trackedCells), 0, 0, 1, 1
							madeFlag = "only.observed"
						else:
							aGenes, bGenes, overlap, universe, maxOverlap, sumOverlap, pvalue, adjPvalue = 0, 0, 0, len(trackedCells), 0, 0, 1, 1
							madeFlag = "none.observed"
						
						# export data:
						print >>f_output, "\t".join(map(str, [aCell, bCell, minimal, fraction, distance, ancestor, correlation, corPvalue, adjCorPvalue, aGenes, bGenes, overlap, universe, maxOverlap, sumOverlap, pvalue, adjPvalue, madeFlag]))
		
		# close output file:
		f_output.close()
		print
	
	
	# cell time mode:
	elif option.mode == "cell.times":
	
		# define input expression files:
		assayedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_assayed"
		startedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_started"
		trackedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_tracked"
		focusedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_focused"
		inheritfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_inherit"
		maximalfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_maximal"
		
		# load cell times:
		print
		print "Loading cellular times..."
		time_matrix = dict()
		inlines = open(extraspath + option.times).readlines()
		for inline in inlines:
			cell, start, stop = inline.strip().split(",")
			for time in range(int(start), int(stop)+1):
				if not time in time_matrix:
					time_matrix[time] = list()
				time_matrix[time].append(cell)
		
		# export cell times:
		populationDict = {
			"assayed" : assayedfile,
			"started" : startedfile,
			"tracked" : trackedfile,
			"focused" : focusedfile,
			"inherit" : inheritfile,
			"maximal" : maximalfile
			}
		print "Exporting cells per time point..."
		for population in populationDict:
			populationCells = general.build2(populationDict[population], id_column="cell", skip=True, mute=True).keys()
			for time in sorted(time_matrix.keys()):
				general.pathGenerator(timepath + population + "/cells/")
				f_output = open(timepath + population + "/cells/" + str(time), "w")
				timedCells = sorted(set(time_matrix[time]).intersection(set(populationCells)))
				if len(timedCells) > 0:
					print >>f_output, "\n".join(timedCells)
				f_output.close()
		
		# generate reports:
		print "Generating reports..."
		for population in populationDict:
			general.pathGenerator(timepath + population + "/report/")
			f_output = open(timepath + population + "/report/mapcells_" + population + "_time_report.txt", "w")
			print >>f_output, "\t".join(["time", "cell.count", "cell.percent", "cell.ids"])
			for time in sorted(time_matrix.keys()):
				general.pathGenerator(timepath + population + "/report/")
				timedCount = general.countLines(timepath + population + "/cells/" + str(time))
				timedPercent = round(100*float(timedCount)/len(time_matrix[time]), 2)
				timedCells = open(timepath + population + "/cells/" + str(time)).read().split("\n")
				print >>f_output, "\t".join([str(time), str(timedCount), str(timedPercent), ",".join(timedCells).rstrip(",")])
			f_output.close()
		print
			
	# cubism graph mode:
	elif option.mode == "cell.cubism":
	
		# define input expression files:
		assayedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_assayed"
		startedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_started"
		trackedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_tracked"
		focusedfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_focused"
		inheritfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_inherit"
		maximalfile = expressionpath + "mapcells_" + expression_flag + option.name + "_expression_maximal"
		
		# define cell populations:
		populationDict = {
			"assayed" : assayedfile,
			"started" : startedfile,
			"tracked" : trackedfile,
			"focused" : focusedfile,
			"inherit" : inheritfile,
			"maximal" : maximalfile
			}
		
		# parse reports:
		print
		print "Exporting per gene, per timepoint expression cells:"
		for population in populationDict:
			
			print "Processing:", population
			
			# define output paths:
			factorpath = cubismpath + population + "/factor/"
			matrixpath = cubismpath + population + "/matrix/"
			general.pathGenerator(factorpath)
			general.pathGenerator(matrixpath)
			
			# build cell-expression matrix:
			quantitation_matrix, expression_matrix, tracking_matrix, trackedCells = expressionBuilder(expressionfile=populationDict[population], path="", cutoff=option.fraction, minimum=option.minimum, metric="fraction.expression")
			
			# load timepoint data:
			timeDict = general.build2(timepath + population + "/report/mapcells_" + population + "_time_report.txt", id_column="time")
			
			# load calendar months:
			monthDict = dict((k,v) for k,v in enumerate(calendar.month_abbr))
			
			# process genes x timepoints:
			m_output = open(matrixpath + "mapcells_cubism_matrix.txt", "w")
			print >>m_output, "\t".join(["gene", "time", "cells", "gene.cells", "time.cells"])
			for gene in sorted(expression_matrix.keys()):
				
				timeStamp = 1001856000000
				timeAdded = 100000000
				factorLines = list()
				for time in sorted(map(int, timeDict.keys())):
					geneCells = expression_matrix[gene]
					timeCells = timeDict[str(time)]["cell.ids"].split(",")
					dateCells = len(set(geneCells).intersection(set(timeCells)))
					date = datetime.datetime.fromtimestamp(timeStamp / 1e3)
					day, month, year = date.day, date.month, str(date.year)[2:]
					date = "-".join(map(str, [day, monthDict[month], year])) 
					factorLines.append(",".join(map(str, [date, dateCells, dateCells, dateCells, len(trackedCells)])))
					print >>m_output, "\t".join(map(str, [gene, time, dateCells, len(geneCells), len(timeCells)]))
					timeStamp += timeAdded
				factorLines.reverse()
				
				f_output = open(factorpath + gene + ".csv", "w")
				print >>f_output, "Date,Open,High,Low,Close,Volume"
				for factorLine in factorLines:
					print >>f_output, factorLine.strip()
				f_output.close()
			
			# close output matrix:
			m_output.close()
			
		print

	# cell annotation mode:
	elif option.mode == "cell.annotation":
	
		# load target cells from time-points:
		print
		if option.times != "OFF":
		
			# define output file:
			f_output = open(cellnotationspath + "mapcells_" + option.name + "_" + option.infile, "w")
		
			# load time-point cells:
			cells = getTargetCells(inpath=timepath + option.times + "/cells/", mode="time", timeRange=range(option.start, option.stop + 1, option.step))
		
		# load target cells from collection:
		elif option.collection != "OFF":
			
			# define output file:
			f_output = open(cellnotationspath + "mapcells_" + option.collection + "_" + option.infile, "w")
		
			# load collection cells:
			cells = getTargetCells(inpath=cellsetpath + option.collection + "/", mode="collection")
		
		# export features per cell:
		print "Exporting features per cell..."
		k = 0
		inlines = open(annotationspath + option.infile).readlines()
		if option.header == "ON":
			inlines.pop(0)
		for cell in cells:
			for inline in inlines:
				if option.format == "bed":
					print >>f_output, cell + ":" + inline.strip()
					k += 1
		f_output.close()
	
		print "Features scaled from", len(inlines), "to", k, ": " + str(round(float(k)/len(inlines), 0)) + "x"
		print
		
	
	# build matrix mode:
	elif option.mode == "cell.matrix":
		
		# update overlappath:
		matrixpath = matrixpath + option.collection + "/"
		general.pathGenerator(matrixpath)

		# define input files:
		infile = expressionpath + option.expression
		
		# define output files:
		matrixfile = matrixpath + str(option.expression + "_" + option.name + "_matrix")
		
		# load header dictionary:
		hd = general.build_header_dict(infile)
		
		# build cellular expression matrix:
		matrix, cells, genes, tissueDict = dict(), list(), list(), dict()
		inlines = open(infile).readlines()
		inlines.pop(0)
		for inline in inlines:
			initems = inline.strip().split("\t")
			cell, gene, cellExpression, fractionExpression, normalExpression, specificTissue, generalTissue, classTissue = initems[hd["cell"]], initems[hd["gene"]], initems[hd["cell.expression"]], initems[hd["fraction.expression"]], initems[hd["normal.expression"]], initems[hd["specific.tissue"]], initems[hd["general.tissue"]], initems[hd["class.tissue"]]
			
			# extract expression value (using specified technique):
			if option.technique == "binary":
				if float(fractionExpression) >= option.fraction and float(cellExpression) >= option.minimum:
					value = 1
				else:
					value = 0
			elif option.technique == "signal":
				value = float(cellExpression)
			elif option.technique == "fraction":
				value = float(fractionExpression)
			elif option.technique == "normal":
				value = float(normalExpression)
			
			# store cells, genes, and values:	
			if not cell in cells:
				cells.append(cell)
			if not gene in genes:
				genes.append(gene)
			if not cell in tissueDict:
				tissueDict[cell] = [classTissue, generalTissue, specificTissue]
			if not cell in matrix:
				matrix[cell] = dict()
			matrix[cell][gene] = value
		
		# export the cellular expression matrix!
		f_output = open(matrixfile, "w")
		cells, genes = sorted(cells), sorted(genes)
		print >>f_output, "\t".join([""] + genes)
		for cell in cells:
			values = list()
			for gene in genes:
				if gene in matrix[cell]:
					values.append(matrix[cell][gene])
				else:
					values.append(0)
			valueCount = len(values) - values.count(0)
			classTissue, generalTissue, specificTissue = tissueDict[cell]
			specificTissue = specificTissue.replace(" ", "_")
			label = ":".join([classTissue, generalTissue, specificTissue, cell])
			print >>f_output, "\t".join([label] + map(str, values))
		f_output.close()
		
	
	# build in silico binding peaks mode:
	elif option.mode == "cell.peaks":
	
		# define the target contexts:
		if option.contexts != "OFF":
			shandle, target_context_dict = metrn.options_dict["contexts.condensed"][option.contexts]
			target_contexts = list()
			for target in target_context_dict:
				target_contexts.append(target_context_dict[target])
			target_contexts = sorted(list(set(target_contexts)))
			
		# generate output paths:
		insilicopath = bindingpath + option.name + "/"
		general.pathGenerator(insilicopath)
		
		# load header dictionary:
		hd = general.build_header_dict(expressionpath + option.expression)
		
		# load expression matrix:
		print
		print "Loading expression matrix..."
		matrix = dict()
		inlines = open(expressionpath + option.expression).readlines()
		inlines.pop(0)
		for inline in inlines:
			initems = inline.strip().split("\t")
			gene, cell, cellExpression, fractionExpression = initems[hd["gene"]], initems[hd["cell"]], float(initems[hd["cell.expression"]]), float(initems[hd["fraction.expression"]])
			if not gene in matrix:
				matrix[gene] = dict()
			if fractionExpression >= option.fraction and float(cellExpression) >= option.minimum:
				matrix[gene][cell] = fractionExpression
		
		# load target cells:
		if option.times != "OFF":
		
			# load time-point cells:
			timedCells = getTargetCells(inpath=timepath + option.times + "/cells/", mode="time", timeRange=range(option.start, option.stop + 1, option.step))
		
		# scan peak files:
		print
		print "Generating cell-resolution peaks..."
		k, peak_files, insilico_files = 0, os.listdir(peakspath), list()
		for peak_file in peak_files:
			dataset = peak_file.split("_peaks.bed")[0].replace("POL2", "AMA-1")
			organism, strain, factor, context, institute, method = metrn.labelComponents(dataset, target="components")
			if factor in matrix:
				if option.contexts == "OFF" or context in target_contexts:
					print "Processing:", dataset
					insilico_file = peak_file.replace("POL2", "AMA-1")
					f_output = open(insilicopath + insilico_file, "w")
					for cell in sorted(matrix[factor].keys()):
						if option.times == "OFF" or cell in timedCells:
							for inline in open(peakspath + peak_file).readlines():
								print >>f_output, cell + ":" + inline.strip()
					f_output.close()
					insilico_files.append(insilico_file)
				
		# define output peak files:
		unsortedfile = bindingpath + "mapcells_silico_" + option.name + "_unsorted.bed"
		completefile = bindingpath + "mapcells_silico_" + option.name + "_complete.bed"
		compiledfile = bindingpath + "mapcells_silico_" + option.name + "_compiled.bed"
		
		# generate compilation files:
		if not "mapcells_silico_" + option.peaks + "_complete.bed" in os.listdir(bindingpath) or option.overwrite == "ON":
		
			# gather peak files and compiled into a single file:
			print
			print "Gathering peaks into single file..."
			joint = " " + insilicopath
			command = "cat " + insilicopath + joint.join(insilico_files) + " > " + unsortedfile
			os.system(command)
		
			print "Sorting peaks in single file..."
			command = "sortBed -i " + unsortedfile + " > " + completefile
			os.system(command)
		
			# merge peaks into single file:
			print "Collapsing peaks in sorted file..."
			command = "mergeBed -nms -i " + completefile + " > " + compiledfile
			os.system(command)
			
			# remove unsorted file:
			command = "rm -rf " + unsortedfile
			os.system(command)
		
		print
	
	
	# gene and cell collection reporting mode (gene expressed per cell):
	elif option.mode == "reports":
		
		taskDict = {
			"gene" : [cellsetpath, "cells"],
			"cell" : [genesetpath, "genes"]
			}
		
		print
		for task in taskDict:
			inputpath, column = taskDict[task]
			for collection in os.listdir(inputpath):
				if collection in option.collection.split(",") or option.collection == "OFF":
					print "Processing:", task, collection
					f_output = open(reportspath + "mapcell_report_" + task + "_" + collection, "w")
					print >>f_output, "\t".join([task, "count", column])
					for item in sorted(os.listdir(inputpath + collection)):
						contents = open(inputpath + collection + "/" + item).read().strip().split("\n")
						contents = general.clean(contents)
						print >>f_output, "\t".join(map(str, [item, len(contents), ",".join(sorted(contents))]))
					f_output.close()	
			print
			
		print
	
	
	# cell collection mode (cells expressed per gene):
	elif option.mode == "cell.collection":
		
		# establish descendants cutoff:
		if option.descendants == "OFF":
			descendants_cutoff = 1000000
			descendants_handle = "XX"
		else:
			descendants_cutoff = int(option.descendants)
			descendants_handle = option.descendants
		
		# establish ascendants cutoff:
		if option.ascendants == "OFF":
			ascendants_cutoff = 0
			ascendants_handle = "XX"
		else:
			ascendants_cutoff = int(option.ascendants)
			ascendants_handle = option.ascendants
		
		# establish limit cutoff:
		if option.limit == "OFF":
			limit_cutoff = "OFF"
			limit_handle = "XX"
		else:
			limit_cutoff = int(option.limit)
			limit_handle = option.limit
		
		# define output folder:
		cellsetpath = cellsetpath + option.collection + "/"
		general.pathGenerator(cellsetpath)
		
		# export expressing-cells for each gene:
		if option.expression != "OFF":
	
			# build cell-expression matrix:
			quantitation_matrix, expression_matrix, tracking_matrix, trackedCells = expressionBuilder(expressionfile=option.expression, path=expressionpath, cutoff=option.fraction, minimum=option.minimum, metric="fraction.expression")
			
			# export cells per gene:
			for gene in sorted(expression_matrix.keys()):
				f_output = open(cellsetpath + gene, "w")
				for cell in sorted(expression_matrix[gene]):
					print >>f_output, cell
				f_output.close()
		
		# export cells for SOM neurons:
		if option.neurons != "OFF":
			
			# update path to neurons:
			neuronspath = neuronspath + option.peaks + "/"
			
			# define input path:
			sumpath = neuronspath + option.technique + "/results/" + option.neurons + "/summary/" 
			sumfile = "mapneurons_summary.txt"
			
			# build header dict:
			hd = general.build_header_dict(sumpath + sumfile)
		
			# build SOM-cell matrix:
			collection_matrix, trackedCells = dict(), list()
			inlines = open(sumpath + sumfile).readlines()
			inlines.pop(0)
			for inline in inlines:
				initems = inline.rstrip("\n").split("\t")
				neuron, cells = initems[hd["neuron"]], initems[hd["class.ids"]]
				collection_matrix[neuron] = general.clean(cells.split(","), "")
				trackedCells.extend(cells.split(","))
			trackedCells = general.clean(list(set(trackedCells)), "")
			
			# export cells per gene:
			for neuron in sorted(collection_matrix.keys()):
				f_output = open(cellsetpath + neuron, "w")
				for cell in sorted(collection_matrix[neuron]):
					print >>f_output, cell
				f_output.close()
		
	
	# gene collection mode (gene expressed per cell):
	elif option.mode == "gene.collection":
		
		# establish descendants cutoff:
		if option.descendants == "OFF":
			descendants_cutoff = 1000000
			descendants_handle = "XX"
		else:
			descendants_cutoff = int(option.descendants)
			descendants_handle = option.descendants
		
		# establish ascendants cutoff:
		if option.ascendants == "OFF":
			ascendants_cutoff = 0
			ascendants_handle = "XX"
		else:
			ascendants_cutoff = int(option.ascendants)
			ascendants_handle = option.ascendants
		
		# establish limit cutoff:
		if option.limit == "OFF":
			limit_cutoff = "OFF"
			limit_handle = "XX"
		else:
			limit_cutoff = int(option.limit)
			limit_handle = option.limit
		
		# define output folder:
		genesetpath = genesetpath + option.collection + "/"
		general.pathGenerator(genesetpath)
		
		# export expressing-cells for each gene:
		if option.expression != "OFF":
	
			# build cell-expression matrix:
			quantitation_matrix, expression_matrix, tracking_matrix, trackedCells = expressionBuilder(expressionfile=option.expression, path=expressionpath, cutoff=option.fraction, minimum=option.minimum, metric="fraction.expression")
			
			# gather cells:
			cells = list()
			for gene in sorted(quantitation_matrix.keys()):
				for cell in sorted(quantitation_matrix[gene]):
					cells.append(cell)
			cells = sorted(list(set(cells)))
			
			# invert matrix:
			inverted_matrix = dict()
			for gene in sorted(expression_matrix.keys()):
				for cell in sorted(expression_matrix[gene]):
					if not cell in inverted_matrix:
						inverted_matrix[cell] = list()
					inverted_matrix[cell].append(gene)
			
			# export cells per gene:
			for cell in cells:
				f_output = open(genesetpath + cell, "w")
				if cell in inverted_matrix:
					for gene in sorted(inverted_matrix[cell]):
						print >>f_output, gene
				f_output.close()
	
	
	# cell transfer mode:
	elif option.mode == "cell.transfer":
		
		# define time-range to examine:
		timeRange=range(option.start, option.stop + 1, option.step)
		
		# generate new collections:
		for timePoint in timeRange:
			
			# define output path:
			outpath = cellsetpath + option.name + general.indexTag(timePoint, option.total) + option.nametag + "/"
			general.pathGenerator(outpath)
			
			# load timePoint cells:
			timedCells = getTargetCells(inpath=timepath + option.times + "/cells/", mode="time", timeRange=[timePoint])
			
			# parse per gene signatures in collection: 
			for gene in os.listdir(cellsetpath + option.collection):
				
				# load expression cells:
				expressionCells = open(cellsetpath + option.collection + "/" + gene).read().split("\n")
				
				# export timed, expressionCells:
				f_output = open(outpath + gene, "w")
				print >>f_output, "\n".join(sorted(list(set(timedCells).intersection(set(expressionCells)))))
				f_output.close()
		
	
	# mapping overlap mode:
	elif option.mode == "cell.overlap":

		# update overlappath:
		overlappath = overlappath + option.collection + "/"
		general.pathGenerator(overlappath)

		# build cell-expression matrix:
		print
		print "Loading cellular expression..."
		signal_matrix, expression_matrix, tracking_matrix, trackedCells = expressionBuilder(expressionfile=option.expression, path=expressionpath, cutoff=option.fraction, minimum=option.minimum, metric="cell.expression")
		fraction_matrix, expression_matrix, tracking_matrix, trackedCells = expressionBuilder(expressionfile=option.expression, path=expressionpath, cutoff=option.fraction, minimum=option.minimum, metric="fraction.expression")
		normal_matrix, expression_matrix, tracking_matrix, trackedCells = expressionBuilder(expressionfile=option.expression, path=expressionpath, cutoff=option.fraction, minimum=option.minimum, metric="fraction.expression")
		rank_matrix, expression_matrix, tracking_matrix, trackedCells = expressionBuilder(expressionfile=option.expression, path=expressionpath, cutoff=option.fraction, minimum=option.minimum, metric="rank")
		
		# load collection cells:
		#print "Loading target cells..."
		targetCells = getTargetCells(inpath=cellsetpath + option.collection + "/", mode="collection")
		
		# create output file:
		o_output = open(overlappath + "mapcells_" + option.collection + "_matrix_overlap", "w")
		o_header = ["i", "j", "i.cells", "j.cells", "overlap.cells", "total.cells", "overlap.max", "overlap.sum", "overlap.avg", "expression.cor", "fraction.cor", "normal.cor", "rank.cor", "pvalue", "pvalue.adj", "score", "i.only.ids", "j.only.ids", "overlap.ids"]
		print >>o_output, "\t".join(o_header)
		
		# find maximum rank, if necessary:
		if option.extend == "ON":
			maxRank = 0 
			for gene in rank_matrix:
				for targetCell in rank_matrix[gene]:
					if int(rank_matrix[gene][targetCell]) > maxRank:
						maxRank = int(rank_matrix[gene][targetCell])
						
		# load gene-expressing cells data:
		print 
		print "Build expression matrix per gene..."
		genes = os.listdir(cellsetpath + option.collection)
		matrix = dict()
		for gene in genes:
			matrix[gene] = dict()
			matrix[gene]["cells"], matrix[gene]["signals"], matrix[gene]["fractions"], matrix[gene]["normals"], matrix[gene]["ranks"] = list(), list(), list(), list(), list()
			expressionCells = open(cellsetpath + option.collection + "/" + gene).read().split("\n")
			for targetCell in targetCells:
				if targetCell in expressionCells:
					matrix[gene]["cells"].append(targetCell)
				if targetCell in signal_matrix[gene]:
					matrix[gene]["signals"].append(signal_matrix[gene][targetCell])
					matrix[gene]["fractions"].append(fraction_matrix[gene][targetCell])
					matrix[gene]["normals"].append(normal_matrix[gene][targetCell])
					matrix[gene]["ranks"].append(rank_matrix[gene][targetCell])
				elif option.extend == "ON":
					matrix[gene]["signals"].append(0)
					matrix[gene]["fractions"].append(0)
					matrix[gene]["normals"].append(0)
					matrix[gene]["ranks"].append(maxRank)
		
		# print a matrix of cell expression overlap between genes:
		print "Exporting cellular overlap matrix..."
		adjust = len(matrix)*len(matrix)
		universe = len(targetCells)
		for geneX in sorted(matrix.keys()):
			cellsX = matrix[geneX]["cells"]
			signalsX, fractionsX, normalsX, ranksX = numpy.array(matrix[geneX]["signals"]), numpy.array(matrix[geneX]["fractions"]), numpy.array(matrix[geneX]["normals"]), numpy.array(matrix[geneX]["ranks"])		
			
			for geneY in sorted(matrix.keys()):
				cellsY = matrix[geneY]["cells"]
				signalsY, fractionsY, normalsY, ranksY = numpy.array(matrix[geneY]["signals"]), numpy.array(matrix[geneY]["fractions"]), numpy.array(matrix[geneY]["normals"]), numpy.array(matrix[geneY]["ranks"])
				
				signalCor = numpy.corrcoef(signalsX, signalsY)[0][1]
				fractionCor = numpy.corrcoef(fractionsX, fractionsY)[0][1]
				normalCor = numpy.corrcoef(normalsX, normalsY)[0][1]
				rankCor = numpy.corrcoef(ranksX, ranksY)[0][1]
				
				cellsXo = sorted(set(cellsX).difference(set(cellsY))) # X-only cells
				cellsYo = sorted(set(cellsY).difference(set(cellsX))) # Y-only cells
				cellsO = sorted(set(cellsX).intersection(set(cellsY))) # overlap
				cellsU = sorted(set(cellsX).union(set(cellsY))) # union
				cellsT = targetCells
				
				# Hypergeometric paramters:
				m = len(cellsX) # number of white balls in urn
				n = universe - len(cellsX) # number of black balls in urn
				N = len(cellsY) # number of balls drawn from urn
				x = len(cellsO) # number of white balls in drawn
					
				# If I pull out all balls with elephant tatoos (N), is the draw enriched in white balls?:
				pvalue = hyper.fishers(x, m+n, m, N, method="right")
				adjPvalue = hyper.limit(pvalue*adjust)
				score = hyper.directional(x, m+n, m, N, adjust=adjust)
				
				output = [geneX, geneY]
				output.append(len(cellsX))
				output.append(len(cellsY))
				output.append(len(cellsO))
				output.append(len(cellsT))
				
				if len(cellsO) > 0:
					output.append(float(len(cellsO))/min(len(cellsX), len(cellsY)))
					output.append(float(len(cellsO))/len(cellsU))
					output.append(float(len(cellsO))/numpy.mean([len(cellsX), len(cellsY)]))
				else:
					output.append(0)
					output.append(0)
					output.append(0)
				
				output.append(signalCor)
				output.append(fractionCor)
				output.append(normalCor)
				output.append(rankCor)
				
				output.append(pvalue)
				output.append(adjPvalue)
				output.append(score)
				
				output.append(";".join(cellsXo))
				output.append(";".join(cellsYo))
				output.append(";".join(cellsO))
				
				if len(output) != len(o_header):
					print len(o_header), len(output)
					print output
					print
					pdb.set_trace()
				
				if " " in "\t".join(map(str, output)):
					print output
					pdb.set_trace()
				
				print >>o_output, "\t".join(map(str, output))
		
		# close output:
		o_output.close()
		print
	
	# hybrid (datatypes) matrix mode:
	elif option.mode == "cell.hybrid":
		
		# get comparison properties:
		peaks, domain = option.a.split("/")[:2]
		collection = option.b.split("/")[0]
		
		# load target contexts:
		codeContexts, targetContexts = metrn.options_dict["contexts.extended"][option.contexts]
		
		# make comparison output folders:
		hybridpath = hybridpath + collection + "/" + peaks + "/" + domain + "/" + codeContexts + "/"
		general.pathGenerator(hybridpath)
		
		# define input files:
		ainfile = str(coassociationspath + option.a).replace("//","/")
		binfile = str(cellspath + "overlap/" + option.b).replace("//","/")
		
		# load input headers:
		aheader = general.build_header_dict(ainfile)
		bheader = general.build_header_dict(binfile)
		
		# load co-association results:
		print
		print "Loading co-associations..."
		cobindingFrames = general.build2(ainfile, id_complex=["i", "j"], separator=":")
		
		# load cellular expression overlap:
		print "Loading co-expression..."
		coexpressionFrames = general.build2(binfile, id_complex=["i", "j"], separator=":")
		coexpressionMatrix = general.build2(binfile, i="i", j="j", x="overlap.sum", mode="matrix")
		
		# characterize input file basenames:
		abasename = option.a.split("/")[len(option.a.split("/"))-1].replace(".txt","").replace(".bed","")
		bbasename = option.b.split("/")[len(option.b.split("/"))-1].replace(".txt","").replace(".bed","")
		
		# define output file:
		f_outfile = hybridpath + "mapcells_hybrid_" + collection + "-" + peaks + "-" + domain + "_combined.txt"
		f_output = open(f_outfile, "w")
		
		# generate output header:
		header = ["i", "j", "label"]
		acolumns = list()
		for acolumn in general.valuesort(aheader):
			if not acolumn in header:
				acolumns.append(acolumn)
		bcolumns = list()
		for bcolumn in general.valuesort(bheader):
			if not bcolumn in header and not bcolumn in ["i.only.ids", "j.only.ids", "overlap.ids"]:
				bcolumns.append(bcolumn)
		print >>f_output, "\t".join(header + acolumns + bcolumns)
		
		# filter-match results:
		print "Merging co-binding and co-expression..."
		ifactor, jfactor = option.indexes.split(",")
		icontext, jcontext = option.values.split(",")
		comparisons = list()
		for cobindingComparison in sorted(cobindingFrames.keys()):
			iFactor, jFactor = cobindingFrames[cobindingComparison][ifactor], cobindingFrames[cobindingComparison][jfactor]
			iContext, jContext = cobindingFrames[cobindingComparison][icontext], cobindingFrames[cobindingComparison][jcontext]
			if iContext in targetContexts and jContext in targetContexts:
				if iFactor in coexpressionMatrix and jFactor in coexpressionMatrix:
					coexpressionComparison = iFactor + ":" + jFactor
					label = ":".join(sorted([iFactor, jFactor]))
					if not coexpressionComparison in comparisons:
						output = [iFactor, jFactor, label]
						for acolumn in acolumns:
							output.append(cobindingFrames[cobindingComparison][acolumn])
						for bcolumn in bcolumns:
							output.append(coexpressionFrames[coexpressionComparison][bcolumn])
						print >>f_output, "\t".join(map(str, output))
						comparisons.append(coexpressionComparison)
						
						# NOTE: this filtering for unique comparisons ensures that only one 
						# of the RNA PolII datasets gets used.
					
		# close output file:
		f_output.close()
		
		print "Merged comparisons:", len(comparisons)
		print


	# dynamics (hybrid) matrix mode:
	elif option.mode == "cell.dynamics":
	
		# are working with hybrid co-binding and co-expression data?
		if option.peaks != "OFF" and option.domain != "OFF":
			hybridMode = "ON"
		else:
			hybridMode = "OFF"
		
		# make comparison output folders:
		if hybridMode == "ON":
			dynamicspath = dynamicspath + option.name + "/" + option.peaks + "/" + option.domain + "/"
			general.pathGenerator(dynamicspath)
			f_outfile = dynamicspath + "mapcells_hybrid_" + option.name + "-" + option.peaks + "-" + option.domain + "_dynamics.txt"
			f_output = open(f_outfile, "w")
			u_outfile = dynamicspath + "mapcells_hybrid_" + option.name + "-" + option.peaks + "-" + option.domain + "_uniqueID.txt"
			u_output = open(u_outfile, "w")
		else:
			dynamicspath = dynamicspath + option.name + "/overlap/"
			general.pathGenerator(dynamicspath)
			f_outfile = dynamicspath + "mapcells_direct_" + option.name + "_dynamics.txt"
			f_output = open(f_outfile, "w")
			u_outfile = dynamicspath + "mapcells_direct_" + option.name + "_uniqueID.txt"
			u_output = open(u_outfile, "w")
			
		# define output file:
		k = 0
		
		# load target contexts:
		codeContexts, targetContexts = metrn.options_dict["contexts.extended"][option.contexts]
		
		# load overlap data from collections:
		print
		print "Transfer dynamic co-binding and co-expression analysis..."
		for index in range(option.start, option.stop+1, option.step):
			collection = option.collection + general.indexTag(index, option.total) + option.nametag
			labels = list()
			
			# process hybrid co-binding and co-expression data:
			if hybridMode == "ON" and collection in os.listdir(hybridpath):
				if option.peaks in os.listdir(hybridpath + collection):
					if option.domain in os.listdir(hybridpath + collection + "/" + option.peaks):
					
						infile = hybridpath + collection + "/" + option.peaks + "/" + option.domain + "/" + codeContexts + "/mapcells_hybrid_" + collection + "-" + option.peaks + "-" + option.domain + "_combined.txt"
						inlines = open(infile).readlines()
						header = inlines.pop(0)
						if k == 0:
							print >>f_output, "timepoint" + "\t" + header.strip()
							print >>u_output, "timepoint" + "\t" + header.strip()
							k += 1
						for inline in inlines:
							label = inline.strip().split("\t")[2]
							print >>f_output, str(index) + "\t" + inline.strip()
							if not label in labels:
								print >>u_output, str(index) + "\t" + inline.strip()	
								labels.append(label)
								
			# process direct co-expression data:
			if hybridMode == "OFF" and collection in os.listdir(overlappath):
				infile = overlappath + collection + "/mapcells_" + collection + "_matrix_overlap"
				inlines = open(infile).readlines()
				header = inlines.pop(0)
				if k == 0:
					headerItems = header.strip().split("\t")[:15]
					print >>f_output, "\t".join(["timepoint", "label"] + headerItems)
					print >>u_output, "\t".join(["timepoint", "label"] + headerItems)
					k += 1
				for inline in inlines:
					initems = inline.strip().split("\t")[:15]
					label = ":".join(sorted([initems[0], initems[1]]))
					print >>f_output, "\t".join([str(index), label] + initems)
					if not label in labels:
						print >>u_output, "\t".join([str(index), label] + initems)	
						labels.append(label)
						
		f_output.close()
		u_output.close()
		print
	
	
	# hypergeometric tissue-testing mode:
	elif option.mode == "test.tissues":
	
		# load tissue annotation matrixes:
		print
		print "Loading tissue annotations..."
		specificTissues = general.build2(expressionpath + option.infile, i="specific.tissue", j="cell", mode="matrix", counter=True)
		generalTissues = general.build2(expressionpath + option.infile, i="general.tissue", j="cell", mode="matrix", counter=True)
		classTissues = general.build2(expressionpath + option.infile, i="class.tissue", j="cell", mode="matrix", counter=True)
		totalCells = general.build2(expressionpath + option.expression, i="cell", x="specific.tissue", mode="values", skip=True)
		totalCells = sorted(totalCells.keys())
		
		# define a function for testing:
		def tissueTesting(queryCells, tissueMatrix, totalCells, adjust=1, match=True):
			if match:
				queryCells = set(queryCells).intersection(set(totalCells))
			tissueOverlap = dict()
			for tissue in sorted(tissueMatrix.keys()):
				tissueCells = sorted(tissueMatrix[tissue].keys())
				if match:
					tissueCells = set(tissueCells).intersection(set(totalCells))
				overlapCells = set(queryCells).intersection(set(tissueCells))
				
				m = len(queryCells)
				n = len(totalCells) - len(queryCells)
				U = len(totalCells)
				N = len(tissueCells)
				x = len(overlapCells)
				
				unionized = len(set(queryCells).union(set(tissueCells)))
				maximized = min(len(queryCells), len(tissueCells))
					
				# determine overlap fractions:
				if maximized > 0:
					maxOverlap = float(x)/maximized
				else:
					maxOverlap = 0
				if unionized > 0:
					sumOverlap = float(x)/unionized
				else:
					sumOverlap = 0
					
				# calculate probability mass function (PMF):
				pvalue = hyper.fishers(x, U, m, N, adjust=1, method="right")
				adjPvalue = hyper.limit(pvalue*adjust)
				
				# calculate enrichment/depletion score:
				score = hyper.directional(x, U, m, N, adjust=adjust)
				
				# store overlap scores:
				tissueOverlap[tissue] = [len(queryCells), len(tissueCells), len(overlapCells), len(totalCells), maxOverlap, sumOverlap, pvalue, adjPvalue, score]
			
			# return overlap scores:
			return tissueOverlap
		
		# load genes:
		genes = sorted(os.listdir(cellsetpath + option.collection))
		
		# determine Bonferroni correction factors
		adjustSpecific = len(genes)*len(specificTissues)
		adjustGeneral = len(genes)*len(generalTissues)
		adjustClass = len(genes)*len(classTissues)
		
		#print adjustSpecific
		#print adjustGeneral
		#print adjustClass
		#pdb.set_trace()
		
		# load cellular expression patterns per gene:
		print "Loading per gene expression matrix..."
		specificMatrix, generalMatrix, classMatrix = dict(), dict(), dict()
		for gene in genes:
			cells = open(cellsetpath + option.collection + "/" + gene).read().split("\n")
			specificMatrix[gene] = tissueTesting(cells, specificTissues, totalCells, adjust=adjustSpecific)
			generalMatrix[gene] = tissueTesting(cells, generalTissues, totalCells, adjust=adjustGeneral)
			classMatrix[gene] = tissueTesting(cells, classTissues, totalCells, adjust=adjustClass)
			
		
		# load cellular expression patterns per gene:
		print "Exporting overlap scores..."
		s_output = open(tissuespath + "mapcells_" + option.collection + "_matrix_specific.txt", "w")
		g_output = open(tissuespath + "mapcells_" + option.collection + "_matrix_general.txt" , "w")
		c_output = open(tissuespath + "mapcells_" + option.collection + "_matrix_class.txt"	  , "w")
		print >>s_output, "\t".join(["gene", "tissue", "gene.cells", "tissue.cells", "overlap.cells", "total.cells", "overlap.max", "overlap.sum", "pvalue", "pvalue.adj", "score"])
		print >>g_output, "\t".join(["gene", "tissue", "gene.cells", "tissue.cells", "overlap.cells", "total.cells", "overlap.max", "overlap.sum", "pvalue", "pvalue.adj", "score"])
		print >>c_output, "\t".join(["gene", "tissue", "gene.cells", "tissue.cells", "overlap.cells", "total.cells", "overlap.max", "overlap.sum", "pvalue", "pvalue.adj", "score"])
		for gene in sorted(specificMatrix.keys()):
			for tissue in sorted(specificMatrix[gene].keys()):
				print >>s_output, "\t".join(map(str, [gene, tissue] + specificMatrix[gene][tissue]))
			for tissue in sorted(generalMatrix[gene].keys()):
				print >>g_output, "\t".join(map(str, [gene, tissue] + generalMatrix[gene][tissue]))
			for tissue in sorted(classMatrix[gene].keys()):
				print >>c_output, "\t".join(map(str, [gene, tissue] + classMatrix[gene][tissue]))
		
		# close outputs:
		s_output.close()
		g_output.close()
		c_output.close()
		print
		
	# lineage construction/generation mode:
	elif option.mode == "build.lineages":
		
		import time
		
		# establish descendants cutoff:
		if option.descendants == "OFF":
			descendants_cutoff = 1000000
			descendants_handle = "XX"
		else:
			descendants_cutoff = int(option.descendants)
			descendants_handle = option.descendants
		
		# establish ascendants cutoff:
		if option.ascendants == "OFF":
			ascendants_cutoff = 0
			ascendants_handle = "XX"
		else:
			ascendants_cutoff = int(option.ascendants)
			ascendants_handle = option.ascendants
		
		# establish limit cutoff:
		if option.limit == "OFF":
			limit_cutoff = "OFF"
			limit_handle = "XX"
		else:
			limit_cutoff = int(option.limit)
			limit_handle = option.limit
		
		# define output paths:
		logpath = lineagepath + option.name + "/" + option.method + "/lineage." + option.lineages + "/ascendants." + ascendants_handle + "/descendants." + descendants_handle + "/limit." + limit_handle + "/log/"
		buildpath = lineagepath + option.name + "/" + option.method + "/lineage." + option.lineages + "/ascendants." + ascendants_handle + "/descendants." + descendants_handle + "/limit." + limit_handle + "/build/"
		general.pathGenerator(logpath)
		general.pathGenerator(buildpath)
		
		
		# prepare log file:
		l_output = open(logpath + "mapcells_build_" + option.cells + ".log", "w")
		
		# clear output folder contents:
		command = "rm -rf " + buildpath + "*"
		os.system(command)
		
		# build cell-expression matrix:
		print
		print "Loading cellular expression..."
		print >>l_output, "Loading cellular expression..."
		quantitation_matrix, expression_matrix, tracking_matrix, trackedCells = expressionBuilder(expressionfile=option.expression, path=expressionpath, cutoff=option.fraction, minimum=option.minimum, metric="fraction.expression")
		
		# store cell-parent relationships:
		print "Loading cell-parent relationships..."
		print >>l_output, "Loading cell-parent relationships..."
		cell_dict, parent_dict, pedigreeCells = relationshipBuilder(pedigreefile=option.pedigree, path=extraspath, trackedCells=trackedCells, lineages=option.lineages)
		
		print "Pedigree cells:", len(pedigreeCells)
		print "Tracked cells:", len(trackedCells)
		print >>l_output, "Pedigree cells:", len(pedigreeCells)
		print >>l_output, "Tracked cells:", len(trackedCells)
		
		# generate lineages for enrichment:
		print
		print "Generating lineages..."
		print >>l_output, ""
		print >>l_output, "Generating lineages..."
		i, j, maxDN, minUP = 0, 0, 0, 10000
		for parent in pedigreeCells:
			i += 1
			
			# define descendant cells:
			descendants = descendantsCollector(parent, parent_dict, cell_dict, descendants=list())
			
			# define ascendants cells:
			ascendants = ascendantsCollector(parent, parent_dict, cell_dict, ascendants=list())
			
			# calculate combinations possible:
			combinations = combinationCalculator(len(descendants), len(descendants))
			
			# apply descendants cutoff:
			if len(descendants) <= descendants_cutoff and len(ascendants) >= ascendants_cutoff:
				j += 1
				
				print parent, len(ascendants), len(descendants), time.asctime(time.localtime())
				print >>l_output, parent, len(ascendants), len(descendants), time.asctime(time.localtime())
				
				# record max and min cutoffs:
				if len(ascendants) < minUP:
					minUP = len(ascendants)
				if len(descendants) > maxDN:
					maxDN = len(descendants)
				
				# define lineage cells:
				if option.method == "descender":
					subtrees = [",".join(descendants)]
				elif option.method == "generator":
					subtrees = lineageGenerator(parent, parent_dict, cell_dict)
				elif option.method == "builder":
					subtrees = lineageBuilder(parent, parent_dict, cell_dict, limit=limit_cutoff)
				elif option.method == "collector":
					subtrees = lineageCollector(expression_matrix[gene], parent_dict, cell_dict)
					print subtrees
					pdb.set_trace() # not implemented yet
				
				# export lineage cells:
				f_output = open(buildpath + parent, "w")
				index = 1
				for subtree in subtrees:
					print >>f_output, "\t".join([parent, parent + "." + str(index), subtree])
					index += 1
				f_output.close()
		
		print
		print "Pedigree nodes lineaged:", i
		print "Pedigree nodes examined:", j, "(" + str(round(100*float(j)/i, 2)) + "%)"
		print "Maximum descendants:", maxDN
		print "Minimum ascendants:", minUP
		print
		
		print >>l_output, ""
		print >>l_output, "Pedigree nodes lineaged:", i
		print >>l_output, "Pedigree nodes examined:", j, "(" + str(round(100*float(j)/i, 2)) + "%)"
		print >>l_output, "Maximum descendants:", maxDN
		print >>l_output, "Minimum ascendants:", minUP
		
		# close output files:
		l_output.close()
		
		#pdb.set_trace()
		
					
	# hypergeometric lineage-testing mode:
	elif option.mode == "test.lineages":
		
		# establish descendants cutoff:
		if option.descendants == "OFF":
			descendants_cutoff = 1000000
			descendants_handle = "XX"
		else:
			descendants_cutoff = int(option.descendants)
			descendants_handle = option.descendants
		
		# establish ascendants cutoff:
		if option.ascendants == "OFF":
			ascendants_cutoff = 0
			ascendants_handle = "XX"
		else:
			ascendants_cutoff = int(option.ascendants)
			ascendants_handle = option.ascendants
		
		# establish limit cutoff:
		if option.limit == "OFF":
			limit_cutoff = "OFF"
			limit_handle = "XX"
		else:
			limit_cutoff = int(option.limit)
			limit_handle = option.limit
		
		# define output paths:
		logpath = lineagepath + option.name + "/" + option.method + "/lineage." + option.lineages + "/ascendants." + ascendants_handle + "/descendants." + descendants_handle + "/limit." + limit_handle + "/log/"
		buildpath = lineagepath + option.name + "/" + option.method + "/lineage." + option.lineages + "/ascendants." + ascendants_handle + "/descendants." + descendants_handle + "/limit." + limit_handle + "/build/"
		hyperpath = lineagepath + option.name + "/" + option.method + "/lineage." + option.lineages + "/ascendants." + ascendants_handle + "/descendants." + descendants_handle + "/limit." + limit_handle + "/hyper/"
		cellsetpath = cellsetpath + option.collection + "/"
		general.pathGenerator(logpath)
		general.pathGenerator(buildpath)
		general.pathGenerator(hyperpath)
		#general.pathGenerator(cellsetpath)
		
		# prepare log file:
		l_output = open(logpath + "mapcells_hyper_" + option.collection + "_" + option.cells + ".log", "w")
		
		# build cell-expression matrix:
		print
		print "Loading cellular expression..."
		print >>l_output, "Loading cellular expression..."
		quantitation_matrix, expression_matrix, tracking_matrix, trackedCells = expressionBuilder(expressionfile=option.expression, path=expressionpath, cutoff=option.fraction, minimum=option.minimum, metric="fraction.expression")
		
		# load cell-parent relationships:
		print "Loading cell-parent relationships..."
		print >>l_output, "Loading cell-parent relationships..."
		cell_dict, parent_dict, pedigreeCells = relationshipBuilder(pedigreefile=option.pedigree, path=extraspath, trackedCells=trackedCells, lineages=option.lineages)
		
		print "Pedigree cells:", len(pedigreeCells)
		print "Tracked cells:", len(trackedCells)
		print >>l_output, "Pedigree cells:", len(pedigreeCells)
		print >>l_output, "Tracked cells:", len(trackedCells)
		
		# prepare for scanning...
		i, j, k = 0, 0, 0
		nodes = general.clean(os.listdir(buildpath), '.DS_Store')
		overlap_dict, pvalue_dict, score_dict = dict(), dict(), dict()
		
		# prepare output file:
		f_output = open(hyperpath + "mapcells_hyper_" + option.collection + "_" + option.cells + ".txt", "w")
		header = ["gene", "node", "lineage", "experiment.cells", "lineage.cells", "overlap.sum", "overlap.max", "overlap.count", "total.count", "lineage.count", "expressed.count", "pvalue", "pvalue.adj", "score", "cells"]
		print >>f_output, "\t".join(map(str, header))
		
		# load target cells:
		print
		print "Loading target cells..."
		print >>l_output, ""
		print >>l_output, "Loading target cells..."
		collection_matrix = dict()
		for collection in os.listdir(cellsetpath):
			collectionCells = general.clean(open(cellsetpath + collection).read().split("\n"), "")
			collection_matrix[collection] = collectionCells
			#print collection, collectionCells
		
		# define multiple-hypothesis correction factor:
		lineageTotal = 0
		for node in nodes:
			lineageTotal += general.countLines(buildpath + node)
		adjust = len(collection_matrix)*lineageTotal
		
		# check background cell population:
		if option.cells == "tracked":
			pedigreeCells = list(trackedCells)
		
		# scan cells for enrichment:
		print "Scanning cells for lineage enrichments..."
		print >>l_output, ""
		print >>l_output, "Scanning cells for lineage enrichments..."
		collectionsEnriched = list()
		for collection in collection_matrix:
			collectionCells = collection_matrix[collection]
			
			# filter (complete) pedigree cells to reduce to tracked cells?
			if option.cells == "tracked" and collection in tracking_matrix:
				completeCells = set(tracking_matrix[collection]).intersection(set(pedigreeCells))
			else:
				completeCells = pedigreeCells
			# Note: These are cells for which we have expression measurements for gene ('collection')...
			# Note: It is not necessary to filter expression cells because these are by definition a subset of the tracked cells.
					
			# scan lineages for enrichment:
			nodesEnriched, linesEnriched = list(), list()
			for node in nodes:
			
				# load node-specific lineages:
				lineageLines = open(buildpath + node).readlines()
				for lineageLine in lineageLines:
					lineageNode, lineageName, lineageCells = lineageLine.strip().split("\t")
					lineageCells = lineageCells.split(",")
					lineageCount = len(lineageCells)
					
					# filter lineage cells to reduce to tracked cells?
					if option.cells == "tracked" and collection in tracking_matrix:
						lineageCells = set(tracking_matrix[collection]).intersection(set(lineageCells))
						
					#print collection, node, len(tracking_matrix[collection]), len(lineageCells), ",".join(lineageCells)
					#pdb.set_trace()
				
					# test enrichment in lineage:
					i += 1
					completed = len(completeCells)
					descended = len(lineageCells)
					collected = len(collectionCells)
					overlaped = len(set(collectionCells).intersection(set(lineageCells)))
					unionized = len(set(collectionCells).union(set(lineageCells)))
					maximized = min(descended, collected)
					
					# determine overlaps:
					if maximized > 0:
						maxOverlap = float(overlaped)/maximized
					else:
						maxOverlap = 0
					if unionized > 0:
						sumOverlap = float(overlaped)/unionized
					else:
						sumOverlap = 0
					
					# check overlap:
					if maxOverlap >= float(option.overlap):
						j += 1
						
						# calculate probability mass function (PMF):
						pvalue = hyper.fishers(overlaped, completed, descended, collected, adjust=1, method="right")
						adjPvalue = hyper.limit(pvalue*adjust)
					
						# calculate enrichment/depletion score:
						score = hyper.directional(overlaped, completed, descended, collected, adjust=adjust)
						
						# should we store this result?
						if adjPvalue < float(option.pvalue):
							k += 1
							if not collection in overlap_dict:
								overlap_dict[collection], pvalue_dict[collection], score_dict[collection] = dict(), dict(), dict()
							overlap_dict[collection][node] = maxOverlap
							pvalue_dict[collection][node] = adjPvalue
							score_dict[collection][node] = score
			
							output = [collection, node, lineageName, len(collectionCells), lineageCount, sumOverlap, maxOverlap, overlaped, completed, descended, collected, pvalue, adjPvalue, score, ','.join(lineageCells)]
							print >>f_output, "\t".join(map(str, output))
							
							if not collection in collectionsEnriched:
								collectionsEnriched.append(collection)
							if not node in nodesEnriched:
								nodesEnriched.append(node)
							linesEnriched.append(lineageName)
			
			print collection, i, j, k, len(collectionsEnriched), len(nodesEnriched), len(linesEnriched)
			print >>l_output, collection, i, j, k, len(collectionsEnriched), len(nodesEnriched), len(linesEnriched)
			
		
		print
		print "Lineages examined:", i
		print "Lineages overlapped:", j
		print "Lineages significant:", k, "(" + str(round(100*float(k)/i, 2)) + "%)"
		print
		
		print >>l_output, ""
		print >>l_output, "Lineages examined:", i
		print >>l_output, "Lineages overlapped:", j
		print >>l_output, "Lineages significant:", k, "(" + str(round(100*float(k)/i, 2)) + "%)"
		
		# close output file
		f_output.close()
		l_output.close()
		
		#pdb.set_trace()

	
	# hypergeometric testing between sets of cells mode:
	elif option.mode == "test.comparison":
		
		# define output paths:
		querypath = cellsetpath + option.query + "/"
		targetpath = cellsetpath + option.target + "/"
		hyperpath = comparepath + option.name + "/hyper/"
		logpath = comparepath + option.name + "/log/"
		general.pathGenerator(hyperpath)
		general.pathGenerator(logpath)
		
		# prepare log file:
		l_output = open(logpath + "mapcells_comparison_" + option.name + "_" + option.cells + ".log", "w")
		
		# build cell-expression matrix:
		print
		print "Loading cellular expression..."
		print >>l_output, "Loading cellular expression..."
		quantitation_matrix, expression_matrix, tracking_matrix, trackedCells = expressionBuilder(expressionfile=option.expression, path=expressionpath, cutoff=option.fraction, minimum=option.minimum, metric="fraction.expression")
		
		# store cell-parent relationships:
		print "Loading cell-parent relationships..."
		print >>l_output, "Loading cell-parent relationships..."
		cell_dict, parent_dict, pedigreeCells = relationshipBuilder(pedigreefile=option.pedigree, path=extraspath, trackedCells=trackedCells, lineages=option.cells)
		# Note that here the lineage-filtering uses the indicated cells option!
		
		print "Pedigree cells:", len(pedigreeCells)
		print "Tracked cells:", len(trackedCells)
		print >>l_output, "Pedigree cells:", len(pedigreeCells)
		print >>l_output, "Tracked cells:", len(trackedCells)
		
		# prepare for scanning...
		overlap_dict, pvalue_dict = dict(), dict()
		i, j, k = 0, 0, 0
		
		# prepare output file:
		f_output = open(hyperpath + "mapcells_test_" + option.name + "_" + option.cells + "_comparison.txt", "w")
		header = ["query", "target", "lineage", "query.cells", "target.cells", "overlap.sum", "overlap.max", "overlap.count", "total.count", "query.count", "target.count", "pvalue", "pvalue.adj", "score", "cells"]
		print >>f_output, "\t".join(map(str, header))
		
		# load query cells:
		print
		print "Loading query cells..."
		print >>l_output, ""
		print >>l_output, "Loading query cells..."
		query_matrix = dict()
		for query in os.listdir(querypath):
			queryCells = general.clean(open(querypath + query).read().split("\n"), "")
			query_matrix[query] = queryCells
			#print query, queryCells
		
		# load target cells:
		print
		print "Loading target cells..."
		print >>l_output, ""
		print >>l_output, "Loading target cells..."
		target_matrix = dict()
		for target in os.listdir(targetpath):
			targetCells = general.clean(open(targetpath + target).read().split("\n"), "")
			target_matrix[target] = targetCells
			#print target, targetCells
		
		# define multiple-hypothesis correction factor:
		adjust = len(query_matrix)*len(target_matrix)
		
		# check background cell population:
		if option.cells == "tracked":
			pedigreeCells = list(trackedCells)
		
		# scan query cells for enrichment:
		print "Scanning target cells for query cells enrichment..."
		print >>l_output, ""
		print >>l_output, "Scanning target cells for query cells enrichment..."
		queriesEnriched = list()
		for query in sorted(query_matrix.keys()):
			queryCells = list(set(query_matrix[query]))
			
			# filter query cells to reduce to tracked cells?
			if option.cells == "tracked":
				queryCells = set(queryCells).intersection(set(pedigreeCells))
					
			# scan target cells for enrichment:
			targetsEnriched, linesEnriched = list(), list()
			for target in sorted(target_matrix.keys()):
				targetCells = list(set(target_matrix[target]))
				
				# filter target cells to reduce to tracked cells?
				if option.cells == "tracked":
					targetCells = set(targetCells).intersection(set(pedigreeCells))
						
				#print query, target, len(queryCells), len(targetCells), ",".join(targetCells)
				#pdb.set_trace()
				
				# test enrichment in lineage:
				i += 1
				completed = len(pedigreeCells)
				descended = len(targetCells)
				collected = len(queryCells)
				overlaped = len(set(queryCells).intersection(set(targetCells)))
				unionized = len(set(queryCells).union(set(targetCells)))
				maximized = min(descended, collected)
					
				# determine overlaps:
				if maximized > 0:
					maxOverlap = float(overlaped)/maximized
				else:
					maxOverlap = 0
				if unionized > 0:
					sumOverlap = float(overlaped)/unionized
				else:
					sumOverlap = 0
				
				# check overlap:
				if maxOverlap >= float(option.overlap):
					j += 1
					
					# calculate probability mass function (PMF):
					pvalue = hyper.fishers(overlaped, completed, descended, collected, adjust=1, method="right")
					adjPvalue = hyper.limit(pvalue*adjust)
				
					# calculate enrichment/depletion score:
					score = hyper.directional(overlaped, completed, descended, collected, adjust=adjust)
				
					if adjPvalue < float(option.pvalue):
						k += 1
						if not query in overlap_dict:
							overlap_dict[query], pvalue_dict[query] = dict(), dict()
						overlap_dict[query][target] = maxOverlap
						pvalue_dict[query][target] = pvalue
			
						output = [query, target, option.target, len(queryCells), len(targetCells), sumOverlap, maxOverlap, overlaped, completed, descended, collected, pvalue, adjPvalue, score, ','.join(targetCells)]
						print >>f_output, "\t".join(map(str, output))
						
						if not query in queriesEnriched:
							queriesEnriched.append(query)
						if not target in targetsEnriched:
							targetsEnriched.append(target)
			
			print query, i, j, k, len(queriesEnriched), len(targetsEnriched), len(linesEnriched)
			print >>l_output, query, i, j, k, len(queriesEnriched), len(targetsEnriched), len(linesEnriched)
			
		
		print
		print "Lineages examined:", i
		print "Lineages overlapped:", j
		print "Lineages significant:", k, "(" + str(round(100*float(k)/i, 2)) + "%)"
		print
		
		print >>l_output, ""
		print >>l_output, "Lineages examined:", i
		print >>l_output, "Lineages overlapped:", j
		print >>l_output, "Lineages significant:", k, "(" + str(round(100*float(k)/i, 2)) + "%)"
		
		# close output file
		f_output.close()
		l_output.close()

		
	# filter testing results to neurons where the region is contained mode:
	elif option.mode == "test.regions":
	
		# update path to neurons:
		neuronspath = neuronspath + option.peaks + "/"
			
		# define input/output paths:
		bedpath = neuronspath + option.technique + "/results/" + option.neurons + "/regions/bed/" 
		querypath = cellsetpath + option.query + "/"
		targetpath = cellsetpath + option.target + "/"
		hyperpath = comparepath + option.name + "/hyper/"
		logpath = comparepath + option.name + "/log/"
		general.pathGenerator(hyperpath)
		general.pathGenerator(logpath)
		
		# load region coordinates per neuron:
		print
		print "Loading regions per neuron matrix..."
		neuron_matrix = dict()
		for bedfile in general.clean(os.listdir(bedpath), ".DS_Store"):
			neuron = bedfile.replace(".bed", "")
			neuron_matrix[neuron] = dict()
			for bedline in open(bedpath + bedfile).readlines():
				chrm, start, stop, region = bedline.strip("\n").split("\t")[:4]
				neuron_matrix[neuron][region] = [chrm, int(start), int(stop)]
		
		# load gene coordinates:
		print "Loading gene/feature coordinates..."
		coord_dict = dict()
		ad = general.build_header_dict(annotationspath + option.reference)
		inlines = open(annotationspath + option.reference).readlines()
		inlines.pop(0)
		for inline in inlines:
			initems = inline.strip("\n").split("\t")
			chrm, start, stop, feature, strand, name = initems[ad["chrm"]], initems[ad["start"]], initems[ad["end"]], initems[ad["feature"]], initems[ad["strand"]], initems[ad["name"]]
			if strand == "+":
				coord_dict[name] = [chrm, int(start)-option.up, int(start)+option.dn]
				coord_dict[feature] = [chrm, int(start)-option.up, int(start)+option.dn]
			elif strand == "-":
				coord_dict[name] = [chrm, int(stop)-option.dn, int(stop)+option.up]
				coord_dict[feature] = [chrm, int(stop)-option.dn, int(stop)+option.up]
		
		# prepare output file:
		f_output = open(hyperpath + "mapcells_test_" + option.name + "_" + option.cells + "_regions.txt", "w")
		
		# define hypergeometric results file:
		hyperfile = hyperpath + "mapcells_test_" + option.name + "_" + option.cells + "_comparison.txt"
		
		# build header dict:
		hd = general.build_header_dict(hyperfile)
	
		# scan hypergeometric results file for cases of overlap:
		print "Scanning hypergeometric results..."
		i, j, k = 0, 0, 0
		inlines = open(hyperfile).readlines()
		print >>f_output, inlines.pop(0).strip("\n")
		queriesMissed, queriesFound, targetsFound = list(), list(), list()
		for inline in inlines:
			initems = inline.strip("\n").split("\t")
			query, target, pvalue = initems[hd["query"]], initems[hd["target"]], initems[hd["pvalue.adj"]]
			if query in coord_dict:
				i += 1
				qchrm, qstart, qstop = coord_dict[query]
				hits = False
				for region in neuron_matrix[target]:
					j += 1
					rchrm, rstart, rstop = neuron_matrix[target][region]
					if qchrm == rchrm:
						if qstart <= rstart and qstop >= rstop:
							k += 1
							hits = True
				if hits:
					print >>f_output, inline.strip("\n")
					queriesFound.append(query)
					targetsFound.append(target)
			else:
				queriesMissed.append(query)
				queriesMissed = sorted(list(set(queriesMissed)))
				
		# close output file
		f_output.close()
		
		queriesFound = sorted(list(set(queriesFound)))
		targetsFound = sorted(list(set(targetsFound)))
		
		#pdb.set_trace()
		print 
		print "Queries found in neurons:", len(queriesFound)
		print "Neurons found in queries:", len(targetsFound)
		print "Searches performed:", i
		print "Searches performed (x Regions):", j
		print "Searches with hits (x Regions):", k
		print "Queries with coordinates and found:", ", ".join(queriesFound)
		print "Queries missed (no coordinates):", len(queriesMissed)
		print "\n".join(queriesMissed)
		print
		

	# false discovery rate mode:
	elif option.mode == "test.fdr":
	
		# update path to neurons:
		neuronspath = neuronspath + option.peaks + "/"
			
		# define input/output paths:
		bedpath = neuronspath + option.technique + "/results/" + option.neurons + "/regions/bed/" 
		querypath = cellsetpath + option.query + "/"
		targetpath = cellsetpath + option.target + "/"
		hyperpath = comparepath + option.name + "/hyper/"
		logpath = comparepath + option.name + "/log/"
		general.pathGenerator(hyperpath)
		general.pathGenerator(logpath)
		
		# load region coordinates per neuron:
		print
		print "Loading regions per neuron matrix..."
		neuron_matrix = dict()
		for bedfile in general.clean(os.listdir(bedpath), ".DS_Store"):
			neuron = bedfile.replace(".bed", "")
			neuron_matrix[neuron] = dict()
			for bedline in open(bedpath + bedfile).readlines():
				chrm, start, stop, region = bedline.strip("\n").split("\t")[:4]
				neuron_matrix[neuron][region] = [chrm, int(start), int(stop)]
		
		# load gene coordinates:
		print "Loading gene/feature coordinates..."
		coord_dict = dict()
		ad = general.build_header_dict(annotationspath + option.reference)
		inlines = open(annotationspath + option.reference).readlines()
		inlines.pop(0)
		for inline in inlines:
			initems = inline.strip("\n").split("\t")
			chrm, start, stop, feature, strand, name = initems[ad["chrm"]], initems[ad["start"]], initems[ad["end"]], initems[ad["feature"]], initems[ad["strand"]], initems[ad["name"]]
			if strand == "+":
				coord_dict[name] = [chrm, int(start)-option.up, int(start)+option.dn]
				coord_dict[feature] = [chrm, int(start)-option.up, int(start)+option.dn]
			elif strand == "-":
				coord_dict[name] = [chrm, int(stop)-option.dn, int(stop)+option.up]
				coord_dict[feature] = [chrm, int(stop)-option.dn, int(stop)+option.up]
		
		# prepare output file:
		f_output = open(hyperpath + "mapcells_test_" + option.name + "_" + option.cells + "_fdr.txt", "w")
		
		# define hypergeometric results file:
		hyperfile = hyperpath + "mapcells_test_" + option.name + "_" + option.cells + "_comparison.txt"
		
		# build header dict:
		hd = general.build_header_dict(hyperfile)
	
		# load positive hypergeometric results:
		print "Loading hypergeometric results (hits)..."
		inlines = open(hyperfile).readlines()
		inlines.pop(0)
		hyper_matrix, hyperTargets = dict(), list()
		for inline in inlines:
			initems = inline.strip("\n").split("\t")
			query, target, pvalue = initems[hd["query"]], initems[hd["target"]], initems[hd["pvalue.adj"]]
			if not query in hyper_matrix:
				hyper_matrix[query] = dict()
			hyper_matrix[query][target] = float(pvalue)
			if not target in hyperTargets:
				hyperTargets.append(target)
			
		# select the best matching neuron for each query:
		match_matrix = dict()
		print "Scanning hypergeometric results per query..."
		i, j, k = 0, 0, 0
		positiveRate, negativeRate, matchTargets = list(), list(), list()
		for query in hyper_matrix:
			if query in coord_dict:
				i += 1
				qchrm, qstart, qstop = coord_dict[query]
				for target in neuron_matrix:
					j += 1
					hits = 0
					for region in neuron_matrix[target]:
						rchrm, rstart, rstop = neuron_matrix[target][region]
						if qchrm == rchrm:
							if qstart <= rstart and qstop >= rstop:
								hits += 1
					if hits != 0:
						if not query in match_matrix:
							match_matrix[query] = dict()
						match_matrix[query][target] = float(hits)/len(neuron_matrix[target])
					if not target in matchTargets:
						matchTargets.append(target)
						
		#print hyper_matrix.keys()
		#print match_matrix.keys()
		#print query
		#print target
		#print match_matrix[query][target]
		#pdb.set_trace()
		
		# Test A
		"""
		print
		print "Testing positive and negative hits..."
		positiveRate, negativeRate, unknownRate = list(), list(), list()
		for query in match_matrix:
			hits = general.valuesort(match_matrix[query])
			hits.reverse()
			target = hits[0]
			if query in hyper_matrix and target in hyper_matrix[query]:
				positiveRate.append(query + ":" + target)
				print "+", query, target, match_matrix[query][target]
			else:
				print "-", query, target, match_matrix[query][target]
				negativeRate.append(query + ":" + target)
				if query in hyper_matrix:
					unknownRate.append(query + ":" + target)
		
		print "True Positive Rate:", len(positiveRate), 100*float(len(positiveRate))/(len(positiveRate)+len(negativeRate))
		print "False Positive Rate:", len(negativeRate), 100*float(len(negativeRate))/(len(positiveRate)+len(negativeRate))
		print "False Unknown Rate:", len(unknownRate), 100*float(len(unknownRate))/(len(positiveRate)+len(unknownRate))
		print
		"""
		
		# Test B
		"""
		print
		print "Testing positive and negative hits..."
		positiveRate, negativeRate, unknownRate = list(), list(), list()
		for query in hyper_matrix:
			hits = 0
			for target in general.valuesort(hyper_matrix[query]):
				if query in match_matrix and target in match_matrix[query]:
					hits += 1
			if hits != 0:
				positiveRate.append(query + ":" + target)
			else:
				negativeRate.append(query + ":" + target)
				if query in match_matrix:
					unknownRate.append(query + ":" + target)
		
		print "True Positive Rate:", len(positiveRate), 100*float(len(positiveRate))/(len(positiveRate)+len(negativeRate))
		print "False Positive Rate:", len(negativeRate), 100*float(len(negativeRate))/(len(positiveRate)+len(negativeRate))
		print "False Unknown Rate:", len(unknownRate), 100*float(len(unknownRate))/(len(positiveRate)+len(unknownRate))
		print
		"""
		
		# Test C
		print
		print "Testing positive and negative hits..."
		positiveRate, negativeRate, unknownRate = list(), list(), list()
		for query in match_matrix:
			hits = 0
			for target in general.valuesort(match_matrix[query]):
				if query in hyper_matrix and target in hyper_matrix[query]:
					hits += 1
			if hits != 0:
				positiveRate.append(query + ":" + target)
			else:
				negativeRate.append(query + ":" + target)
				if query in hyper_matrix:
					unknownRate.append(query + ":" + target)
		
		print "Genes enriched in SOM neurons:", len(hyper_matrix)
		print "Genes with promoter in SOM neurons:", len(match_matrix)
		print "Neurons enriched in gene expression:", len(hyperTargets)
		print "Neurons with gene promoter matches:", len(matchTargets)
		print "True Positive Rate:", len(positiveRate), 100*float(len(positiveRate))/(len(positiveRate)+len(negativeRate))
		print "False Positive Rate:", len(negativeRate), 100*float(len(negativeRate))/(len(positiveRate)+len(negativeRate))
		print "False Unknown Rate (not enriched in any neuron):", len(unknownRate), 100*float(len(unknownRate))/(len(positiveRate)+len(unknownRate))
		print
		
		# scan each gene for cellular overlap in neurons where the promoter is found:
		"""
		print "Scanning positive and negative hits..."
		i, j, k = 0, 0, 0
		positiveRate, negativeRate = list(), list()
		for query in hyper_matrix:
			if query in coord_dict:
				i += 1
				qchrm, qstart, qstop = coord_dict[query]
				for target in neuron_matrix:
					j += 1
					hits = 0
					for region in neuron_matrix[target]:
						rchrm, rstart, rstop = neuron_matrix[target][region]
						if qchrm == rchrm:
							if qstart <= rstart and qstop >= rstop:
								hits += 1
					if hits != 0:
						if target in hyper_matrix[query]:
							positiveRate.append(query + ":" + target)
						else:
							negativeRate.append(query + ":" + target)
		"""
		
		#print >>f_output, inlines.pop(0).strip("\n")
		
		# close output file
		f_output.close()
		
		
		#print 
		#print "Queries found in neurons:", len(queriesFound)
		#print "Neurons found in queries:", len(targetsFound)
		#print "Searches performed:", i
		#print "Searches performed (x Regions):", j
		#print "Searches with hits (x Regions):", k
		#print "Queries with coordinates and found:", ", ".join(queriesFound)
		#print "Queries missed (no coordinates):", len(queriesMissed)
		#print

	
	# annotate tissue composition in neurons:
	elif option.mode == "test.composition":
	
		# update path to neurons:
		neuronspath = neuronspath + option.peaks + "/"
			
		# define input/output paths:
		bedpath = neuronspath + option.technique + "/results/" + option.neurons + "/regions/bed/" 
		codespath = neuronspath + option.technique + "/results/" + option.neurons + "/codes/" 
		summarypath = neuronspath + option.technique + "/results/" + option.neurons + "/summary/" 
		querypath = cellsetpath + option.query + "/"
		targetpath = cellsetpath + option.target + "/"
		compositionpath = comparepath + option.name + "/composition/"
		hyperpath = comparepath + option.name + "/hyper/"
		logpath = comparepath + option.name + "/log/"
		general.pathGenerator(compositionpath)
		general.pathGenerator(hyperpath)
		general.pathGenerator(logpath)
		
		# load codes:
		inlines = open(codespath + option.neurons + ".codes").readlines()
		codes = inlines.pop(0).strip().split("\t")
		codeDict = dict()
		for inline in inlines:
			initems = inline.strip().split("\t")
			neuron = initems.pop(0)
			codeDict["neuron" + neuron] = initems
		
		# load cellular expression data:
		print
		print "Loading cellular annotation..."
		annotationDict = general.build2(expressionpath + option.expression, id_column="cell", value_columns=["specific.tissue", "general.tissue", "class.tissue", "match.tissue"], skip=True, verbose=False)
		
		# load tissue annotation matrixes:
		print "Loading tissue annotations..."
		#specificCounts = general.build2(expressionpath + option.infile, i="specific.tissue" , mode="values", skip=True, counter=True)
		#generalCounts = general.build2(expressionpath + option.infile, i="general.tissue", mode="values", skip=True, counter=True)
		#classCounts = general.build2(expressionpath + option.infile, i="class.tissue", mode="values", skip=True, counter=True)
		totalCells = general.build2(expressionpath + option.expression, i="cell", x="specific.tissue", mode="values", skip=True)
		totalCells = sorted(totalCells.keys())
		
		# gather tissue labels
		specificTissues, generalTissues, classTissues, matchTissues = list(), list(), list(), list()
		for cell in annotationDict:
			if not annotationDict[cell]["specific.tissue"] in specificTissues:
				specificTissues.append(annotationDict[cell]["specific.tissue"])
			if not annotationDict[cell]["general.tissue"] in generalTissues:
				generalTissues.append(annotationDict[cell]["general.tissue"])
			if not annotationDict[cell]["class.tissue"] in classTissues:
				classTissues.append(annotationDict[cell]["class.tissue"])
			if not annotationDict[cell]["match.tissue"] in matchTissues:
				matchTissues.append(annotationDict[cell]["match.tissue"])
			
		# load cells identified in each neuron:
		print
		print "Loading cell identities per neuron..."
		neuronDict = general.build2(summarypath + "mapneurons_summary.txt", id_column="neuron", value_columns=["class.ids"])
		
		# load cells counted in each neuron:
		print "Loading cell counts per neuron..."
		countMatrix, binaryMatrix = dict(), dict()
		for neuron in os.listdir(bedpath):
			inlines = open(bedpath + neuron).readlines()
			neuron = neuron.replace(".bed", "")
			countMatrix[neuron] = dict()
			binaryMatrix[neuron] = dict()
			for inline in inlines:
				chrm, start, end, feature, score, strand, cell, regions = inline.strip().split("\t")
				if not cell in countMatrix[neuron]:
					countMatrix[neuron][cell] = 0
					binaryMatrix[neuron][cell] = 1
				countMatrix[neuron][cell] += 1
		
		# generate tissue class scores:
		cellList = list()
		cellMatrix, specificMatrix, generalMatrix, classMatrix, matchMatrix = dict(), dict(), dict(), dict(), dict()
		for neuron in neuronDict:
			if not neuron in cellMatrix:
				cellMatrix[neuron] = dict()
				specificMatrix[neuron] = dict()
				generalMatrix[neuron] = dict()
				classMatrix[neuron] = dict()
				matchMatrix[neuron] = dict()
			for cell in neuronDict[neuron]["class.ids"].split(","):
				specificTissue, generalTissue, classTissue, matchTissue = annotationDict[cell]["specific.tissue"], annotationDict[cell]["general.tissue"], annotationDict[cell]["class.tissue"], annotationDict[cell]["match.tissue"]
				
				if not cell in cellMatrix[neuron]:
					cellMatrix[neuron][cell] = 0
				if not specificTissue in specificMatrix[neuron]:
					specificMatrix[neuron][specificTissue] = 0
				if not generalTissue in generalMatrix[neuron]:
					generalMatrix[neuron][generalTissue] = 0
				if not classTissue in classMatrix[neuron]:
					classMatrix[neuron][classTissue] = 0
				if not matchTissue in matchMatrix[neuron]:
					matchMatrix[neuron][matchTissue] = 0
				
				cellList.append(cell)
				cellMatrix[neuron][cell] += binaryMatrix[neuron][cell]
				specificMatrix[neuron][specificTissue] += binaryMatrix[neuron][cell]
				generalMatrix[neuron][generalTissue] += binaryMatrix[neuron][cell]
				classMatrix[neuron][classTissue] += binaryMatrix[neuron][cell]
				matchMatrix[neuron][matchTissue] += binaryMatrix[neuron][cell]
		cellList = sorted(list(set(cellList)))
		
		# Note: The above dictionaries record how many of the cell (ids)
		# in a given neuron have correspond to a given tissue.
		
		# prepare class tallies for normalization:
		specificTallies, generalTallies, classTallies, matchTallies = dict(), dict(), dict(), dict()
		for cell in cellList:
			if not annotationDict[cell]["specific.tissue"] in specificTallies:
				specificTallies[annotationDict[cell]["specific.tissue"]] = 0
			if not annotationDict[cell]["general.tissue"] in generalTallies:
				generalTallies[annotationDict[cell]["general.tissue"]] = 0
			if not annotationDict[cell]["class.tissue"] in classTallies:
				classTallies[annotationDict[cell]["class.tissue"]] = 0
			if not annotationDict[cell]["match.tissue"] in matchTallies:
				matchTallies[annotationDict[cell]["match.tissue"]] = 0
			specificTallies[annotationDict[cell]["specific.tissue"]] += 1
			generalTallies[annotationDict[cell]["general.tissue"]] += 1
			classTallies[annotationDict[cell]["class.tissue"]] += 1
			matchTallies[annotationDict[cell]["match.tissue"]] += 1
		
		# Note: The above tallies record the number of cells (observed, 
		# in neurons) that correspond to each tissue.***
		
		# prepare output files:
		f_output = open(compositionpath + "mapcells_composition_codes.txt", "w")
		c_output = open(compositionpath + "mapcells_composition_cellular.txt", "w")
		s_output = open(compositionpath + "mapcells_composition_specific.txt", "w")
		g_output = open(compositionpath + "mapcells_composition_general.txt", "w")
		l_output = open(compositionpath + "mapcells_composition_class.txt", "w")
		m_output = open(compositionpath + "mapcells_composition_match.txt", "w")
		
		# print out headers:
		print >>f_output, "\t".join(["neuron", "id", "fraction.ids"])
		print >>c_output, "\t".join(["neuron", "id", "id.found", "id.cells", "fraction.ids", "fraction.sum", "fraction.max", "fraction.nrm", "pvalue", "pvalue.adj", "score"])
		print >>s_output, "\t".join(["neuron", "id", "id.found", "id.cells", "fraction.ids", "fraction.sum", "fraction.max", "fraction.nrm", "pvalue", "pvalue.adj", "score"])
		print >>g_output, "\t".join(["neuron", "id", "id.found", "id.cells", "fraction.ids", "fraction.sum", "fraction.max", "fraction.nrm", "pvalue", "pvalue.adj", "score"])
		print >>l_output, "\t".join(["neuron", "id", "id.found", "id.cells", "fraction.ids", "fraction.sum", "fraction.max", "fraction.nrm", "pvalue", "pvalue.adj", "score"])
		print >>m_output, "\t".join(["neuron", "id", "id.found", "id.cells", "fraction.ids", "fraction.sum", "fraction.max", "fraction.nrm", "pvalue", "pvalue.adj", "score"])
		
		# Note: We will now output the following information:
		# id.found : is ID found in neuron?
		# id.cells : number of cells (diversity) that match ID.
		# fraction.ids: fraction of ID diversity in neuron.
		# fraction.sum: fraction of cellular diversity in neuron that matches ID.
		# fraction.rat: fraction of cellular diversity in neuron that matches ID, normalized by the representation of the ID.
		# fraction.max: fraction of cellular diversity in neuron as normalized by the ID with the highest cellular diversity in neuron. 
		# fraction.nrm: fraction of cellular diversity in neuron as normalized by the total number of cells with said ID.
		
		# determine missed tissues:
		print
		specificMissed, generalMissed, classMissed, matchMissed = set(specificTissues).difference(set(specificTallies.keys())), set(generalTissues).difference(set(generalTallies.keys())), set(classTissues).difference(set(classTallies.keys())), set(matchTissues).difference(set(matchTallies.keys()))
		print "Specific tissues not found:", str(len(specificMissed)) + " (" + str(len(specificTissues)) + ") ; " + ",".join(sorted(specificMissed))
		print "General tissues not found:", str(len(generalMissed)) + " (" + str(len(generalTissues)) + ") ; " + ",".join(sorted(generalMissed))
		print "Class tissues not found:", str(len(classMissed)) + " (" + str(len(classTissues)) + ") ; " + ",".join(sorted(classMissed))
		print "Match tissues not found:", str(len(matchMissed)) + " (" + str(len(matchTissues)) + ") ; " + ",".join(sorted(matchMissed))
		print
		
		# export the fractions:
		print "Exporting representation per neuron..."
		for neuron in sorted(neuronDict.keys()):
			
			if neuron in codeDict:
			
				# export factor signals:			
				index = 0
				for code in codes:
					print >>f_output, "\t".join(map(str, [neuron, code, codeDict[neuron][index]]))
					index += 1
				
				# export cell counts:
				for cell in cellList:
					adjust = len(neuronDict.keys())*len(cellList)
					types = len(cellMatrix[neuron].keys())
					total = sum(cellMatrix[neuron].values())
					maxxx = max(cellMatrix[neuron].values())
					if cell in cellMatrix[neuron]:
						count = float(cellMatrix[neuron][cell])
						index = 1
					else:
						count = 0
						index = 0
					print >>c_output, "\t".join(map(str, [neuron, cell, index, count, float(index)/types, float(count)/total, float(count)/maxxx, 1, 1, 1, 0]))
				
				# export specific tissue enrichment:
				for specificTissue in sorted(specificTallies.keys()):
					types = len(specificMatrix[neuron].keys())
					total = sum(specificMatrix[neuron].values())
					maxxx = max(specificMatrix[neuron].values())
					tally = specificTallies[specificTissue]
					if specificTissue in specificMatrix[neuron]:
						count = float(specificMatrix[neuron][specificTissue])
						index = 1
					else:
						count = 0
						index = 0
					adjust = len(neuronDict.keys())*len(specificTallies.keys())
					universe = sum(specificTallies.values())
					pvalue = hyper.fishers(count, universe, total, tally, adjust=1, method="right")
					adjPvalue = hyper.limit(pvalue*adjust)
					score = hyper.directional(count, universe, total, tally, adjust=adjust)
					print >>s_output, "\t".join(map(str, [neuron, specificTissue, index, count, float(index)/types, float(count)/total, float(count)/maxxx, float(count)/tally, pvalue, adjPvalue, score]))
				
				# export general tissue enrichment:
				for generalTissue in sorted(generalTallies.keys()):
					types = len(generalMatrix[neuron].keys())
					total = sum(generalMatrix[neuron].values())
					maxxx = max(generalMatrix[neuron].values())
					tally = generalTallies[generalTissue]
					if generalTissue in generalMatrix[neuron]:
						count = float(generalMatrix[neuron][generalTissue])
						index = 1
					else:
						count = 0
						index = 0
					adjust = len(neuronDict.keys())*len(generalTallies.keys())
					universe = sum(generalTallies.values())
					pvalue = hyper.fishers(count, universe, total, tally, adjust=1, method="right")
					adjPvalue = hyper.limit(pvalue*adjust)
					score = hyper.directional(count, universe, total, tally, adjust=adjust)
					print >>g_output, "\t".join(map(str, [neuron, generalTissue, index, count, float(index)/types, float(count)/total, float(count)/maxxx, float(count)/tally, pvalue, adjPvalue, score]))
					
				# export class tissue enrichment:
				for classTissue in sorted(classTallies.keys()):
					types = len(classMatrix[neuron].keys())
					total = sum(classMatrix[neuron].values())
					maxxx = max(classMatrix[neuron].values())
					tally = classTallies[classTissue]
					if classTissue in classMatrix[neuron]:
						count = float(classMatrix[neuron][classTissue])
						index = 1
					else:
						count = 0
						index = 0
					adjust = len(neuronDict.keys())*len(classTallies.keys())
					universe = sum(classTallies.values())
					pvalue = hyper.fishers(count, universe, total, tally, adjust=1, method="right")
					adjPvalue = hyper.limit(pvalue*adjust)
					score = hyper.directional(count, universe, total, tally, adjust=adjust)
					print >>l_output, "\t".join(map(str, [neuron, classTissue, index, count, float(index)/types, float(count)/total, float(count)/maxxx, float(count)/tally, pvalue, adjPvalue, score]))
		
				# export match tissue enrichment:
				for matchTissue in sorted(matchTallies.keys()):
					types = len(matchMatrix[neuron].keys())
					total = sum(matchMatrix[neuron].values())
					maxxx = max(matchMatrix[neuron].values())
					tally = matchTallies[matchTissue]
					if matchTissue in matchMatrix[neuron]:
						count = float(matchMatrix[neuron][matchTissue])
						index = 1
					else:
						count = 0
						index = 0
					adjust = len(neuronDict.keys())*len(matchTallies.keys())
					universe = sum(matchTallies.values())
					pvalue = hyper.fishers(count, universe, total, tally, adjust=1, method="right")
					adjPvalue = hyper.limit(pvalue*adjust)
					score = hyper.directional(count, universe, total, tally, adjust=adjust)
					print >>m_output, "\t".join(map(str, [neuron, matchTissue, index, count, float(index)/types, float(count)/total, float(count)/maxxx, float(count)/tally, pvalue, adjPvalue, score]))
		
		# close outputs:
		f_output.close()
		c_output.close()
		s_output.close()
		g_output.close()
		l_output.close()
		m_output.close()
		
		print
		print "Combining cell and factor (mix) information.."
		
		# load input factor information:
		factorDict = general.build2(compositionpath + "mapcells_composition_codes.txt", i="neuron", j="id", x="fraction.ids", mode="matrix")
		
		# define input cell/tissue files:
		infiles = ["mapcells_composition_cellular.txt", "mapcells_composition_specific.txt", "mapcells_composition_general.txt", "mapcells_composition_class.txt", "mapcells_composition_match.txt"]
		for infile in infiles:
			print "Processing:", infile
			
			# initiate neuron data extraction:
			f_output = open(compositionpath + infile.replace(".txt", ".mix"), "w")
			inheader = open(compositionpath + infile).readline().strip().split("\t")
			inlines = open(compositionpath + infile).readlines()
			print >>f_output, inlines.pop(0)
			
			# append factor information to neuron data:
			processed = list()
			for inline in inlines:
				neuron, label = inline.strip().split("\t")[:2]
				if not neuron in processed:
					processed.append(neuron)
					for factor in factorDict[neuron]:
						output = list()
						for column in inheader:
							if column == "neuron":
								output.append(neuron)
							elif column == "id":
								output.append(factor)
							elif column in ["pvalue", "pvalue.adj"]:
								output.append("1")
							else:
								output.append(factorDict[neuron][factor])
						print >>f_output, "\t".join(output)
				print >>f_output, inline.strip()
				
			# close outputs:
			f_output.close()
		
		print
		
	
	# examine co-association correspondence between genes:
	elif option.mode == "test.similarity":
	
		# update path to neurons:
		neuronspath = neuronspath + option.peaks + "/"
			
		# define input/output paths:
		bedpath = neuronspath + option.technique + "/results/" + option.neurons + "/regions/bed/" 
		querypath = cellsetpath + option.query + "/"
		targetpath = cellsetpath + option.target + "/"
		hyperpath = comparepath + option.name + "/hyper/"
		logpath = comparepath + option.name + "/log/"
		general.pathGenerator(hyperpath)
		general.pathGenerator(logpath)
		
		# load query cells:
		print
		print "Loading query cells..."
		query_matrix = dict()
		for query in os.listdir(querypath):
			queryCells = general.clean(open(querypath + query).read().split("\n"), "")
			query_matrix[query] = queryCells
			#print query, queryCells
		
		print "Generating merged region file..."
		#queryfile = hyperpath  + "query.bed"
		#regionsfile = hyperpath + "regions.bed"
		#overlapfile = hyperpath + "overlap.bed"
			
		joint = " " + bedpath
		command = "cat " + bedpath + joint.join(os.listdir(bedpath)) + " > " + regionsfile
		os.system(command)
			
		# load gene coordinates:
		print "Loading gene/feature coordinates..."
		coord_dict = dict()
		ad = general.build_header_dict(annotationspath + option.reference)
		inlines = open(annotationspath + option.reference).readlines()
		inlines.pop(0)
		for inline in inlines:
			initems = inline.strip("\n").split("\t")
			chrm, start, stop, feature, strand, name = initems[ad["#chrm"]], initems[ad["start"]], initems[ad["end"]], initems[ad["feature"]], initems[ad["strand"]], initems[ad["name"]]
			if strand == "+":
				start, end = int(start)-option.up, int(start)+option.dn
			elif strand == "-":
				start, end = int(stop)-option.dn, int(stop)+option.up
			
			for query in query_matrix:
				if query == feature or query == name:
					f_output = open(queryfile, "w")
					print >>f_output, "\t".join(map(str, [chrm, start, end, feature, 0, strand]))
					f_output.close()
					overlaps = list()
					command = "intersectBed -u -a " + regionsfile + " -b " + queryfile + " > " + overlapfile
					os.system(command)
					for inline in open(overlapfile).readlines():
						overlaps.append(inline.strip())
					print query, len(overlaps)
					if len(overlaps) > 0:
						pdb.set_trace()
					break
	
	# tree building mode:
	elif option.mode == "tree.build":
		
		# build cell-expression matrix:
		print
		print "Loading cellular expression..."
		quantitation_matrix, expression_matrix, tracking_matrix, trackedCells = expressionBuilder(expressionfile=option.expression, path=expressionpath, cutoff=option.fraction, minimum=option.minimum, metric="fraction.expression")
		
		# store cell-parent relationships:
		print "Loading cell-parent relationships..."
		cell_dict, parent_dict, pedigreeCells = relationshipBuilder(pedigreefile=option.pedigree, path=extraspath, mechanism="simple")
		
		print "Pedigree cells:", len(pedigreeCells)
		print "Tracked cells:", len(trackedCells)
		
		# trim tree:
		cell_tree, parent_tree = dict(), dict()
		for parent in parent_dict:
			for cell in parent_dict[parent]:
				ascendants = ascendantsCollector(cell, parent_dict, cell_dict, ascendants=list())
				process = False
				if option.lineages == "complete":
					process = True
				elif parent in trackedCells and cell in trackedCells:
					process = True
				elif option.ascendants != "OFF" and len(ascendants) < int(option.ascendants):
					process = True
				
				if process:
					if not parent in parent_tree:
						parent_tree[parent] = list()
					parent_tree[parent].append(cell)
					cell_tree[cell] = parent
		
		tree = treeBuilder(parent_tree, cell_tree)
		#print sorted(tree.keys())
		#print tree["P0"]
		#pdb.set_trace()
		f_output = open(cellspath + "mapcells_tree_" + option.name + ".json", "w")
		json.dump(tree["P0"], f_output)
		f_output.close()
	
	
	# tree coloring mode:
	elif option.mode == "tree.color":
		
		# build cell-expression matrix:
		print
		print "Loading cellular expression..."
		quantitation_matrix, expression_matrix, tracking_matrix, trackedCells = expressionBuilder(expressionfile=option.expression, path=expressionpath, cutoff=option.fraction, minimum=option.minimum, metric="fraction.expression")
		
		# store cell-parent relationships:
		print "Loading cell-parent relationships..."
		cell_dict, parent_dict, pedigreeCells = relationshipBuilder(pedigreefile=option.pedigree, path=extraspath, mechanism="simple")
		
		print "Pedigree cells:", len(pedigreeCells)
		print "Tracked cells:", len(trackedCells)
		
		# trim tree:
		cell_tree, parent_tree = dict(), dict()
		for parent in parent_dict:
			for cell in parent_dict[parent]:
				ascendants = ascendantsCollector(cell, parent_dict, cell_dict, ascendants=list())
				process = False
				if option.lineages == "complete":
					process = True
				elif parent in trackedCells and cell in trackedCells:
					process = True
				elif option.ascendants != "OFF" and len(ascendants) < int(option.ascendants):
					process = True
				
				if process:
					if not parent in parent_tree:
						parent_tree[parent] = list()
					parent_tree[parent].append(cell)
					cell_tree[cell] = parent
		
		# build header dict:
		hd = general.build_header_dict(option.infile)
		
		# load input lines:
		pvalue_matrix, cells_matrix = dict(), dict()
		inlines = open(option.infile).readlines()
		inlines.pop(0)
		for inline in inlines:
			initems = inline.strip("\n").split("\t")
			query, target, pvalue, cells = initems[hd["query"]], initems[hd["target"]], initems[hd["pvalue"]], initems[hd["cells"]]
			if not query in pvalue_matrix:
				pvalue_matrix[query] = dict()
				cells_matrix[query] = dict()
			pvalue_matrix[query][target] = float(pvalue)
			cells_matrix[query][target] = cells.split(",")
		
		# scan inputs, selecting the targets of highest enrichment and generating color tree for each:
		k = 0
		print
		print "Scanning queries..."
		for query in cells_matrix:
			target = general.valuesort(pvalue_matrix[query])[0]
			cells = cells_matrix[query][target]
			print query, target, pvalue_matrix[query][target], len(cells)
			tree = treeBuilder(parent_tree, cell_tree, highlights=cells)
			#print sorted(tree.keys())
			#print tree["P0"]
			#pdb.set_trace()
			f_output = open(cellspath + "mapcells_tree_" + option.name + "_" + query + "-" + target + ".json", "w")
			json.dump(tree["P0"], f_output)
			f_output.close()
			k += 1
		
		print
		print "Queries processed:", k
		print 
					
if __name__ == "__main__":
	main()
	
	print "Completed:", time.asctime(time.localtime())

#python mapCells.py --path ~/meTRN --mode import --infile murray_2012_supplemental_dataset_1_per_gene.txt --name murray # Retired!

#python mapCells.py --path ~/meTRN --mode import --infile waterston_avgExpression.csv --name waterston --measure max.expression
#python mapCells.py --path ~/meTRN --mode import --infile waterston_avgExpression.csv --name waterston --measure avg.expression

#python mapCells.py --path ~/meTRN --mode check.status --peaks optimal_standard_factor_sx_rawraw --name waterston --measure avg.expression
#python mapCells.py --path ~/meTRN --mode check.status --peaks optimal_standard_factor_ex_rawraw --name waterston --measure avg.expression

#python mapCells.py --path ~/meTRN/ --mode build.lineages --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_tracked --name waterston.tracked --method builder --lineages tracked --descendants OFF --ascendants OFF --limit 10000
#python mapCells.py --path ~/meTRN/ --mode build.lineages --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_tracked --name waterston.tracked --method builder --lineages complete --descendants OFF --ascendants OFF --limit 10000

#python mapCells.py --path ~/meTRN/ --mode test.lineages --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_tracked --name waterston.tracked --method builder --lineages tracked --descendants OFF --ascendants OFF --limit 10000
#python mapCells.py --path ~/meTRN/ --mode test.lineages --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_assayed --name waterston.assayed --method builder --lineages tracked --descendants OFF --ascendants OFF --limit 10000


#python mapCells.py --path ~/meTRN --organism ce --mode robust --infile waterston_avgExpression.csv