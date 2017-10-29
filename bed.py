#!/usr/bin/env python

import sys
import os
import csv

import interval

class CommentedFileReader:
	"""
	Helper class for file reading.
	Skips lines starting with '#'

	tsv_file = csv.reader(CommentedFileReader("inputfile.txt"),
		delimiter='\t')
	for row in tsv_file:
		print row[2] # prints column 3 of each line
	"""
    
	def __init__(self, f, commentstring="#"):
		self.f = open(f, 'rU')
		self.commentstring = commentstring
		return

	def next(self):
		line = self.f.next()
		while line.startswith(self.commentstring):
			line = self.f.next()
		return line

	def __iter__(self):
		return self

csv.register_dialect('bed', delimiter = '\t', \
	quoting = csv.QUOTE_NONE, skipinitialspace = True)

csv.register_dialect('bedgraph', delimiter = '\t', \
	quoting = csv.QUOTE_NONE, skipinitialspace = True)

class BEDReader(csv.DictReader):
	"""
	Read BED files into a DictReader.
	See BEDReader.FIELDS for field names

	Example:
		bed = BEDReader("file.bed")
		for line in bed:
			# print the chromStart
			print(line['chromStart'])
	"""
	FIELDS = ('chrom', 'chromStart', 'chromEnd', \
		'name', 'score', 'strand', 'thickStart', 'thickEnd', \
		'itemRgb', 'blockCount', 'blockSizes', 'blockStarts')
	
	def __init__(self, filename):
		csv.DictReader.__init__(self, CommentedFileReader(filename), dialect='bed', \
			fieldnames=self.FIELDS)


class BedGraphReader(csv.DictReader):
	"""
	Read BedGraph files into a DictReader.
	See BedGraphReader.FIELDS for field names

	Example:
		bedgraph = BedGraphReader("file.bed")
		for line in bedgraph:
			# print the chromStart
			print(line['chromStart'])
	"""
	FIELDS = ('chrom', 'chromStart', 'chromEnd', 'signal')

	def __init__(self, filename):
		csv.DictReader.__init__(self, CommentedFileReader(filename), dialect='bedgraph', \
			fieldnames=self.FIELDS)
		return




class Bed():

	def __init__(self):
		self.features = [] #of type interval.Interval
		return

	def readFromFile(self):
		#This applies to bed files to keep in memory
		for line in self.bedFile:
			#print line["chromStart"], line["chromEnd"], len(line)
			currInterval = interval.Interval()
			fields = line.rstrip().split("\t")
			currInterval.chr = fields[0]
			currInterval.start = int(fields[1])
			currInterval.end = int(fields[2])
			self.features.append(currInterval)
		return

	def sortByChromosomeAndStartAndEnd(self):
		self.features.sort(key=lambda Interval:Interval.end)
		self.features.sort(key=lambda Interval:Interval.start)
		self.features.sort(key=lambda Interval:Interval.chr)
		return

	def readLineBed(self):
		for line in self.bedFile:
			currInterval = interval.Interval()
			currInterval.extractFeatureBed(line)
			#print currInterval.start
			break

		return currInterval

class Peak():

	def __init__(self, filename):
		self.features = [] #of type interval.Interval
		if filename != None:
			try:
				self.bedFile = open(filename, "r")
			except:
				sys.stderr.write("Not able to open peak file " + filename + "\n")
				sys.exit()
		return

	def readFromFile(self):
		for line in self.bedFile:
			if "#" in line:
				continue
			fields = line.rstrip().split("\t")
			currInterval = interval.Interval()
			currInterval.chr = fields[0]
			currInterval.start = int(fields[1])
			currInterval.end = int(fields[2])
			self.features.append(currInterval)
		return

	def sortByChromosomeAndStartAndEnd(self):
		self.features.sort(key=lambda Interval:Interval.end)
		self.features.sort(key=lambda Interval:Interval.start)
		self.features.sort(key=lambda Interval:Interval.chr)
		return

class BedGraph():

	def __init__(self, filename):
		self.features = [] #of type Interval
		try:
			self.bedgraphFile = BedGraphReader(filename)
		except:
			sys.stderr.write("Not able to open bedgraph file " + filename + "\n")
			sys.exit()
		return

	def readLineBG(self):

		currInterval = interval.Interval()
		for line in self.bedgraphFile:
			
			currInterval.extractFeatureBedGraph(line)
			#print currInterval.start
			break

		return currInterval