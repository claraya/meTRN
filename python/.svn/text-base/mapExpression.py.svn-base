#!/usr/bin/env python
# convert expression data into expressed and not-expressed genes

import sys
import time
import optparse
import general
import numpy
import pickle
import pdb
import metrn
import modencode
import multiprocessing
import os

from runner import *

print "Command:", " ".join(sys.argv)
print "Timestamp:", time.asctime(time.localtime())


""" define a function that converts a local user path to SCG3 or GS server paths """
def serverPath(inpath, server="ON"):
	if server="ON":
		return inpath.replace("/Users/claraya/", "/srv/gs1/projects/snyder/claraya/")
	elif server="GS":
		return inpath.replace("/Users/claraya/", "/net/fields/vol1/home/araya/")

				
def main():
	
	parser = optparse.OptionParser()
	parser.add_option("--path", action = "store", type = "string", dest = "path", help = "path from script to files")
	parser.add_option("--organism", action = "store", type = "string", dest = "organism", help = "Target organism for operations...", default="OFF")
	parser.add_option("--mode", action = "store", type = "string", dest = "mode", help = "launch")
	parser.add_option("--infile", action = "store", type = "string", dest = "infile", help = "input expression file", default="OFF")
	parser.add_option("--folder", action = "store", type = "string", dest = "folder", help = "folder location of input files", default="OFF")
	parser.add_option("--nuclear", action = "store", type = "string", dest = "nuclear", help = "Peaks are only nuclear?", default="ON")
	parser.add_option("--genes", action = "store", type = "string", dest = "genes", help = "reference gene file", default="OFF")
	parser.add_option("--transcripts", action = "store", type = "string", dest = "transcripts", help = "reference transcript file")
	parser.add_option("--coord", action = "store", type = "string", dest = "coord", help = "reference coordinates to export: RNA or TSS", default="RNA")
	parser.add_option("--cutoff", action = "store", type = "float", dest = "cutoff", help = "expression cutoff", default=0.05)
	parser.add_option("--name", action = "store", type = "string", dest = "name", help = "name for output file", default="OFF")
	parser.add_option("--threads", action = "store", type = "string", dest = "threads", help = "multiprocessing threads", default="1")
	parser.add_option("--qsub", action = "store", type = "string", dest = "qsub", help = "qsub configuration header", default="OFF")
	parser.add_option("--server", action = "store", type = "string", dest = "server", help = "are we on the server?", default="OFF")
	parser.add_option("--job", action = "store", type = "string", dest = "job", help = "job name for cluster", default="OFF")
	(option, args) = parser.parse_args()
	
	# import paths:
	if option.server == "OFF":
		path_dict = modencode.configBuild(option.path + "/input/" + "configure_path.txt")
	elif option.server == "ON":
		path_dict = modencode.configBuild(option.path + "/input/" + "configure_server.txt")
	elif option.server == "GS":
		path_dict = modencode.configBuild(option.path + "/input/" + "configure_nexus.txt")
	
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
	elif option.organism == "mm" or option.organism == "m.musculus":
		organismTag = "mm"
	elif option.organism == "ce" or option.organism == "c.elegans":
		organismTag = "ce"
	elif option.organism == "dm" or option.organism == "d.melanogaster":
		organismTag = "dm"
	
	# specify genome size file:
	if option.organism != "OFF":
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
	
	# define expression cutoff handle:
	expression = "%.0e" % (float(option.cutoff))
	
	# parse expression data:
	if option.mode == "parse":
		
		# determine path of input files:
		loadpath = path_dict[option.folder]
		
		# specify expression header dict:
		HD = {
			"chrm" : 0,
			"start" : 1,
			"stop" : 2,
			"feature" : 3,
			"dcpm" : 4,
			"strand" : 5,
			"id" : 6,
			"expressed" : 7,
			"normalized" : 8,
			"overlap.gene" : 9,
			"status" : 10
			}
		
		# load transcript expression data:
		measured_dict = general.build2(loadpath + option.infile, id_column="feature", header_dict=HD)
		
		# load reference gene data:
		gene_dict = general.build2(annotationspath + option.genes, id_column="feature", skip=True, mute=True)
		
		# load reference transcript data:
		transcript_dict = general.build2(annotationspath + option.transcripts, id_column="feature", skip=True, mute=True)
		
		# load gene expression matrix:
		k, m = 0, 0
		expression_dict, conversion_dict, genes = dict(), dict(), list()
		for transcript in measured_dict:
			gene, dcpm, status = measured_dict[transcript]["overlap.gene"], float(measured_dict[transcript]["dcpm"]), measured_dict[transcript]["status"]
			genes.append(gene)
			if status == "Confirmed" and float(dcpm) > float(option.cutoff):
				k += 1
				if not gene in expression_dict:
					expression_dict[gene] = 0
					conversion_dict[gene] = transcript
					m += 1
				if float(dcpm) > expression_dict[gene]:
					expression_dict[gene] = float(dcpm)
					conversion_dict[gene] = transcript
		genes = set(genes)
		
		print
		print "Input transcripts:", len(measured_dict.keys())
		print "Input genes:", len(genes)
		print "Reference genes:", len(id2name_dict.keys())
		print "Expressed genes:", len(expression_dict.keys())
		print "Expressed genes matched (in name dict):", len(set(id2name_dict.keys()).intersection(set(expression_dict.keys())))
		print "Expressed genes matched (in gene dict):", len(set(gene_dict.keys()).intersection(set(expression_dict.keys())))
		print "Expressed genes matched (in RNAs dict):", len(set(transcript_dict.keys()).intersection(set(expression_dict.keys())))
		print
		
		# define output files:
		e_output = open(annotationspath + option.name.replace("RENAME", "expressed"), "w")
		r_output = open(annotationspath + option.name.replace("RENAME", "repressed"), "w")
		f_output = open(annotationspath + option.name.replace("RENAME", "reference"), "w")
			
		# export expressed genes:
		for gene in expression_dict:
			transcript = conversion_dict[gene]
			chrm, start, stop, score, strand = measured_dict[transcript]["chrm"], measured_dict[transcript]["start"], measured_dict[transcript]["stop"], measured_dict[transcript]["dcpm"], measured_dict[transcript]["strand"]
			if option.coord == "TSS":
				if strand == "+":
					start, stop = start, int(start) + 1
				else:
					start, stop = int(stop) - 1, stop
			print >>e_output, "\t".join(map(str, [chrm, start, stop, gene, score, strand, transcript]))
			print >>f_output, "\t".join(map(str, [chrm, start, stop, gene, score, strand, transcript]))
		
		# export repressed genes:
		for gene in sorted(gene_dict.keys()):
			if not gene in expression_dict:
				chrm, start, stop, score, strand = gene_dict[gene]["chrm"], gene_dict[gene]["start"], gene_dict[gene]["end"], 0, gene_dict[gene]["strand"]
				
				# select largest transcript:
				matching_dict = dict()
				for transcript in transcript_dict:
					if gene in transcript:
						tstart, tstop = int(transcript_dict[transcript]["start"]), int(transcript_dict[transcript]["end"])
						matching_dict[transcript] = tstop-tstart
				transcripts = general.valuesort(matching_dict)
				transcripts.reverse()
				transcript = transcripts[0]			
				#print gene, transcript
				#pdb.set_trace()
				
				# export data:
				chrm, start, stop, score, strand = transcript_dict[transcript]["chrm"], transcript_dict[transcript]["start"], transcript_dict[transcript]["end"], 0, transcript_dict[transcript]["strand"]
				if option.coord == "TSS":
					if strand == "+":
						start, stop = start, int(start) + 1
					else:
						start, stop = int(stop) - 1, stop
				print >>r_output, "\t".join(map(str, [chrm, start, stop, gene, score, strand, transcript]))
				print >>f_output, "\t".join(map(str, [chrm, start, stop, gene, score, strand, transcript]))
		
		# close output files:
		e_output.close()
		r_output.close()
		f_output.close()
		

if __name__ == "__main__":
	main()
	
	print "Completed:", time.asctime(time.localtime())

#python mapExpression.py --path ~/ceTRN --mode parse --infile in2shape_expression_ee.bed --folder annotations --cutoff 0 --genes in2shape_wbGene.bed --transcripts in2shape_wbTrans.bed --coord TSS