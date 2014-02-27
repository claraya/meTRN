#!/usr/bin/env python
# prepare configuration files for Circos!

import sys
import time
import optparse
import general
import numpy
import pickle
import pdb
import metrn
import fasta
import modencode
import os

print "Command:", " ".join(sys.argv)
print "Timestamp:", time.asctime(time.localtime())

""" define a function that converts a local user path to SCG3 or GS server paths """
def serverPath(inpath, server="ON"):
	if server == "ON":
		return inpath.replace("/Users/claraya/", "/srv/gs1/projects/snyder/claraya/")
	elif server == "GS":
		return inpath.replace("/Users/claraya/", "/net/fields/vol1/home/araya/")
				
def main():
	
	parser = optparse.OptionParser()
	parser.add_option("--path", action = "store", type = "string", dest = "path", help = "Path from script to files")
	parser.add_option("--organism", action = "store", type = "string", dest = "organism", help = "Model organism targeted for analysis")
	parser.add_option("--infile", action = "store", type = "string", dest = "infile", help = "Input file (with path)")
	parser.add_option("--mode", action = "store", type = "string", dest = "mode", help = "Type of operations to be performed: karyotype, import, extend...")
	parser.add_option("--name", action = "store", type = "string", dest = "name", help = "Name of Circos visualization", default="OFF")
	parser.add_option("--track", action = "store", type = "string", dest = "track", help = "Name of Circos data track", default="OFF")
	parser.add_option("--nuclear", action = "store", type = "string", dest = "nuclear", help = "Peaks are only nuclear?", default="ON")
	parser.add_option("--color", action = "store", type = "string", dest = "color", help = "Color for the karyotype track", default="red")
	parser.add_option("--scale", action = "store", type = "string", dest = "scale", help = "Scale for chromosome labeling", default="1000000")
	parser.add_option("--min", action = "store", type = "string", dest = "min", help = "Min value on track", default="1")
	parser.add_option("--max", action = "store", type = "string", dest = "max", help = "Max value on track", default="50")
	parser.add_option("--hi", action = "store", type = "string", dest = "hi", help = "Hi-point of track", default="0.99")
	parser.add_option("--lo", action = "store", type = "string", dest = "lo", help = "Lo-point of track", default="0.50")
	parser.add_option("--fillcolor", action = "store", type = "string", dest = "fillcolor", help = "Fill color for the karyotype track", default="optred")
	parser.add_option("--fillunder", action = "store", type = "string", dest = "fillunder", help = "Fill under the line?", default="yes")
	parser.add_option("--parameters", action = "store", type = "string", dest = "parameters", help = "Variable parameters...", default="")
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
	circospath = path_dict["circos"]
	
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
		
	# define organism parameters:
	if option.organism == "h.sapiens" or option.organism == "human" or option.organism == "hs":
		contextTag = "cells"
		idColumns = ["name", "code", "hgcn","ensembl"]
		idComplexList = list()
	elif option.organism == "m.musculus" or option.organism == "mouse" or option.organism == "mm":
		contextTag = "cells"
		idColumns = ["name", "code", "hgcn","ensembl"]
		idComplexList = list()
	elif option.organism == "c.elegans" or option.organism == "worm" or option.organism == "ce":
		contextTag = "stage"
		idColumns = ["name", "code", "wormbase","ensembl"]
		idComplexList = list()
	elif option.organism == "d.melanogaster" or option.organism == "fly" or option.organism == "dm":
		contextTag = "stage"
		idColumns = ["name", "code", "flybase","ensembl"]
		idComplexList = list()
	
	# prepare karyotype file:
	if option.mode == "karyotype":
	
		# generate Circos path:
		general.pathGenerator(circospath)
		
		# setup output file:
		f_output = open(circospath + "mapcircos_" + organismTag +"_karyotype.txt", "w")
		index = 1
		
		print
		print "Generating karyotype file..."
		for chrm in chromosomes:
			#chr - hs1 1 0 249250621 chr1
			print >>f_output, " ".join(["chr", "-", organismTag + str(index), chrm, "0", str(genome_size_dict[chrm]), "chr" + str(index)])
			index += 1
		
		print "Karyotype: ", len(chromosomes), "chromosomes"
		print 
				
		# close output files:
		f_output.close()
	
	# import bed features into circos format:
	if option.mode == "import":
	
		print
		print "Loading chromosomes..."
		index, chrm_dict = 1, dict()
		for chrm in chromosomes:
			chrm_dict[chrm.upper()] = str(index)
			index += 1
		
		# determine output name:
		if option.name == "OFF":
			filename = option.infile.split("/")
			filename = filename[len(filename)-1]
			option.name = filename.replace(".bed","").replace(".txt","")
		
		# remove previous outputs:
		#command = "rm -rf " + circospath + option.name
		#os.system(command)
			
		# prepare output path:
		analysispath = circospath + option.name + "/"
		datapath = circospath + option.name + "/data/"
		etcpath = circospath + option.name + "/etc/"
		general.pathGenerator(analysispath)
		general.pathGenerator(datapath)
		general.pathGenerator(etcpath)
		
		# copy circos templates:
		command = "cp -r " + circospath + "templates/* " + circospath + option.name
		os.system(command)
			
		command = "cp " + circospath + "mapcircos_" + organismTag +"_karyotype.txt " + datapath + "mapcircos_karyotype.txt"
		os.system(command)
		
		# replace color line:
		command = "cp " + circospath + "templates/etc/colors-" + organismTag + ".conf " + etcpath + "colors.conf"
		os.system(command)
		
		# mess with configuration file:
		configuration = open(circospath + "templates/circos.conf").read()
		configuration = configuration.replace("chromosomes_units = 1000000", "chromosomes_units = " + option.scale)
		configuration = configuration.replace("min = 1", "min = " + option.min)
		configuration = configuration.replace("max = 50", "max = " + option.max)
		configuration = configuration.replace("r0 = 0.50r", "r0 = " + option.hi + "r")
		configuration = configuration.replace("r1 = 0.99r", "r1 = " + option.lo + "r")
		configuration = configuration.replace("color = red", "color = " + option.color)
		configuration = configuration.replace("fill_color = optred", "fill_color = " + option.fillcolor)
		configuration = configuration.replace("fill_under = yes", "fill_under = " + option.fillunder)
		c_output = open(analysispath + "circos.conf", "w")
		print >>c_output, configuration
		c_output.close
		
		# setup output file:
		f_output = open(datapath + "mapcircos_data.txt", "w")
		
		# import the features:
		print "Creating Circos data file..."
		inlines = open(option.infile).readlines()
		inlines.pop(0)
		for inline in inlines:
			chrm, start, stop, feature, score, strand = inline.strip().split("\t")[:6]
			if chrm in chromosomes:
				chrm = organismTag + chrm_dict[chrm.upper()]
				print >>f_output, " ".join([chrm, start, stop, score])
		
		# close output files:
		f_output.close()
		print

	# extend bed features into circos format:
	if option.mode == "extend":
	
		# determine track name:
		if option.track == "OFF":
			trackname = option.infile.split("/")
			trackname = trackname[len(trackname)-1]
			option.track = trackname.replace(".bed","").replace(".txt","")
		option.track = "mapcircos_track_" + option.track + ".txt"
		
		print
		print "Track:", option.track
		
		# load chromosomes:
		index, chrm_dict = 1, dict()
		for chrm in chromosomes:
			chrm_dict[chrm.upper()] = str(index)
			index += 1
		
		# load output path:
		analysispath = circospath + option.name + "/"
		datapath = circospath + option.name + "/data/"
		etcpath = circospath + option.name + "/etc/"
		
		# load configuration segment:
		configuration = open(circospath + "templates/circos.conf").read()
		configuration = configuration.split("<plots>")[1].split("</plots>")[0]
		configuration = configuration.replace("min = 1", "min = " + option.min)
		configuration = configuration.replace("max = 50", "max = " + option.max)
		configuration = configuration.replace("r0 = 0.50r", "r0 = " + option.hi + "r")
		configuration = configuration.replace("r1 = 0.99r", "r1 = " + option.lo + "r")
		configuration = configuration.replace("color = red", "color = " + option.color)
		configuration = configuration.replace("fill_color = optred", "fill_color = " + option.fillcolor)
		configuration = configuration.replace("fill_under = yes", "fill_under = " + option.fillunder)
		configuration = configuration.replace("file = data/mapcircos_data.txt", "file = data/" + option.track)
		
		# update configuration file:
		print "Updating configuration..."
		updated = open(analysispath + "circos.conf").read()
		updated = updated.replace("</plots>", configuration + "</plots>")
		c_output = open(analysispath + "circos.conf", "w")
		print >>c_output, updated
		c_output.close
		
		# setup output file:
		f_output = open(datapath + option.track, "w")
		
		# import the features:
		print "Creating track data file..."
		inlines = open(option.infile).readlines()
		inlines.pop(0)
		for inline in inlines:
			chrm, start, stop, feature, score, strand = inline.strip().split("\t")[:6]
			if chrm in chromosomes:
				chrm = organismTag + chrm_dict[chrm.upper()]
				print >>f_output, " ".join([chrm, start, stop, score])
		
		# close output files:
		f_output.close()
		print
		
	
	# launch circos visualization mode:
	if option.mode == "launch":
	
		# load output path:
		analysispath = circospath + option.name + "/"
		
		# launch circos:
		print
		print "Generating:", option.name
		os.chdir(analysispath)
		command = "circos -conf circos.conf"
		os.system(command)
		print
		
		
if __name__ == "__main__":
	main()
	
	print "Completed:", time.asctime(time.localtime())

#python mapCircos.py --path ~/meTRN --mode karyotype --organism ce
#python mapCircos.py --path ~/meTRN --mode import --organism ce --infile ~/meTRN/data/hot/regions/maphot_ce_selection_reg_cx_occP05_any.bed --scale 100000 --hi 0.90 --lo 0.98 --color dgrey --fillcolor dgrey
#python mapCircos.py --path ~/meTRN --mode extend --organism ce --name maphot_ce_selection_reg_cx_occP05_any --infile ~/meTRN/data/hot/regions/maphot_ce_selection_reg_ex_occP05_hot.bed --track ex --hi 0.80 --lo 0.88 --color optblue --fillcolor optblue
#python mapCircos.py --path ~/meTRN --mode extend --organism ce --name maphot_ce_selection_reg_cx_occP05_any --infile ~/meTRN/data/hot/regions/maphot_ce_selection_reg_l1_occP05_hot.bed --track l1 --hi 0.70 --lo 0.78 --color optgreen --fillcolor optgreen
#python mapCircos.py --path ~/meTRN --mode extend --organism ce --name maphot_ce_selection_reg_cx_occP05_any --infile ~/meTRN/data/hot/regions/maphot_ce_selection_reg_l2_occP05_hot.bed --track l2 --hi 0.60 --lo 0.68 --color optyellow --fillcolor optyellow
#python mapCircos.py --path ~/meTRN --mode extend --organism ce --name maphot_ce_selection_reg_cx_occP05_any --infile ~/meTRN/data/hot/regions/maphot_ce_selection_reg_l3_occP05_hot.bed --track l3 --hi 0.50 --lo 0.58 --color optorange --fillcolor optorange
#python mapCircos.py --path ~/meTRN --mode extend --organism ce --name maphot_ce_selection_reg_cx_occP05_any --infile ~/meTRN/data/hot/regions/maphot_ce_selection_reg_l4_occP05_hot.bed --track l4 --hi 0.40 --lo 0.48 --color optred --fillcolor optred
#circos -conf ./circos.conf