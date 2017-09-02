#!/usr/bin/env python
#Requirements:
# python2.7
# numpy, matplotlib, scipy, pybedtools, metaseq

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
import random
from pylab import xticks
from scipy import interpolate
from matplotlib.backends.backend_pdf import PdfPages
from collections import OrderedDict

import pybedtools as pbt
import metaseq
from scipy.stats import norm
import seaborn as sns

import bed
import interval

class shapeCharacteristics():
	def __init__(self):
		self.h1 = -1
		self.h2 = -1
		self.h3 = -1
		self.w1 = -1
		self.w2 = -1
		self.w = -1
		self.ratio1 = -1
		self.ratio2 = -1
		self.dis = -1
		return

	def parseCharacteristics(self, line):
		fields = line.rstrip().split("\t")
		self.h1 = float(fields[0])
		self.h2 = float(fields[1])
		self.h3 = float(fields[2])
		self.w1 = float(fields[3])
		self.w2 = float(fields[4])
		self.w = self.w1 + self.w2
		self.ratio1 = float(fields[5])
		self.ratio2 = float(fields[6])
		self.dis = float(fields[7])
		self.ratioW = self.w1/self.w
		return

def extractFeature(line):
	currInterval = interval.Interval()
	fields = line.strip().split("\t")
	currInterval.chr = fields[0]
	currInterval.start = int(fields[1])
	currInterval.end = int(fields[2])
	if len(fields) >= 6:
		#currInterval.signal = float(fields[3])
		currInterval.strand = fields[3]
	return currInterval


def calculateNumberNegatives(nonrandomFile, histone, opPrefix):
	"""
	The calculateNumberNegatives function calculates the number of nonrandom regions with signal. 
	Args:
		nonrandomFile: nonrandom file (of datatype: string)
		histone: is the histone Signal file (of datatype:string)
	Returns:
		numberNegativeIntervals is the number of possible negative intervals possible
	"""

	os.system("bigWigToBedGraph " + histone + " " + opPrefix + "_signal.bedgraph")
	os.system("subtractBed -a " + opPrefix + "_signal.bedgraph -b " + nonrandomFile + " > " + opPrefix + "_signal2.bedgraph")
	histoneFile =  bed.BedGraph(opPrefix + "_signal2.bedgraph")
	numberNegativeIntervals = 0
	idx = 0
	sigFileEnds = False
	while not(sigFileEnds):
		currInterval = histoneFile.readLineBG()
		if currInterval.chr == None:
			sigFileEnds = True
		numCurrIntervals = (currInterval.end - currInterval.start)/25
		if (currInterval.end - currInterval.start) % 25 != 0:
			numCurrIntervals += 1
		numberNegativeIntervals += numCurrIntervals
	os.system("rm " + opPrefix + "_signal.bedgraph")

	return numberNegativeIntervals

def smoothSignal(currSignal):
	numBins = np.size(currSignal)
	smoothedSignal = np.zeros(numBins)
	smoothedSignal[0] = currSignal[0]
	smoothedSignal[1] = currSignal[1]
	for idx in range(2, numBins-2):
		smoothedSignal[idx]=  np.mean(currSignal[idx-2:idx+3])
	smoothedSignal[numBins-2] = currSignal[numBins-2]
	smoothedSignal[numBins-1] = currSignal[numBins-2]
	return smoothedSignal

def calculateDifferential(signal):
	diffSig = np.zeros(np.size(signal))
	diffSig[0] = 0
	for idx in range(1, np.size(signal)):
		diffSig[idx] = signal[idx] - signal[idx - 1]
	#print signal, diffSig
	return diffSig

def findMinima(diffSig):
	middleIdx = np.size(diffSig)/2
	for idx in range(middleIdx, 3, -1):
		if diffSig[idx-3] < 0 and diffSig[idx-2] < 0 and diffSig[idx-1] < 0 \
			and diffSig[idx] >= 0 and diffSig[idx + 1] >= 0 and diffSig[idx + 2] >= 0:
			#print "minima", intersectingIntervals[idx].start, intersectingIntervals[idx].end
			return idx

	for idx in range(middleIdx, np.size(diffSig) - 2):
		if diffSig[idx-3] < 0 and diffSig[idx-2] < 0 and diffSig[idx-1] < 0 \
			and diffSig[idx] >= 0 and diffSig[idx + 1] >= 0 and diffSig[idx + 2] >= 0:
			#print "minima", intersectingIntervals[idx].start, intersectingIntervals[idx].end
			#print middleIdx, idx
			return idx
	return -1

def findSurroundingMaxima(diffSig, minimaIdx):
	if minimaIdx == 0:
		return [-1, -1]

	for idx in range(minimaIdx, 2, -1):
		if diffSig[idx-3] > 0 and diffSig[idx-2] > 0 and diffSig[idx-1] > 0 \
			and diffSig[idx] <= 0 and diffSig[idx + 1] <= 0 and diffSig[idx + 2] <= 0:
			#print "maxima1", intersectingIntervals[idx].start, intersectingIntervals[idx].end
			maximaIdx1 = idx
			break	
	if idx == 3:
		return[-1, -1]

	for idx in range(minimaIdx, np.size(diffSig) - 2):
		if diffSig[idx-3] > 0 and diffSig[idx-2] > 0 and diffSig[idx-1] > 0 \
			and diffSig[idx] <= 0 and diffSig[idx + 1] <= 0 and diffSig[idx + 2] <= 0:
			#print "maxima2", intersectingIntervals[idx+1].start, intersectingIntervals[idx+1].end
			maximaIdx2 = idx + 1
			break
	if idx == np.size(diffSig) - 3:
		return [-1, -1]
	else:
		return [maximaIdx1, maximaIdx2]

def calculatePeakCharacteristics(filteredPositives, smoothedSignal, currFeature, maximaIndices, minimaIdx, metaMaxima1, metaMaxima2, metaMinima, metaIntersectingIntervals, doublePeakRegions, op, op2):
	binWidth = 25
	h1 = smoothedSignal[maximaIndices[0]]
	h2 = smoothedSignal[minimaIdx]
	h3 = smoothedSignal[maximaIndices[1]]
	coord1 = currFeature.start + binWidth * maximaIndices[0] + binWidth/2 
	coord2 = currFeature.start + binWidth * minimaIdx + binWidth/2 
	coord3 = currFeature.start + binWidth * maximaIndices[1] + binWidth/2 
	w1 = coord2 - coord1
	w2 = coord3 - coord2
	maximum = max([h1, h3])
	maximum2 = min([h1, h3])
	deltaH1 = maximum - maximum2
	deltaH2 = maximum - h2
	dis = coord2 - (currFeature.start + currFeature.end)/2
	
	if dis <= 500 and dis >= -500:
		op.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(h1, h2, h3, w1, w2, deltaH1, deltaH2, dis))
		metaIntersectingIntervals.append(smoothedSignal)
		metaMaxima1.append(maximaIndices[0])
		metaMaxima2.append(maximaIndices[1])
		metaMinima.append(minimaIdx)
		currDP = interval.Interval()
		currDP.chr = currFeature.chrom
		currDP.start = currFeature.start + binWidth * maximaIndices[0]
		currDP.end = currFeature.start + (binWidth + 1)* maximaIndices[1]
		op2.write(currDP.chr + "\t" + str(currDP.start) + "\t" + str(currDP.end) + "\n")
		doublePeakRegions.append(currDP)
		filteredPositives.append(currFeature)
		
	return 
				

def getDoublePeakRegions(MPRApeaks, signal, metaMaxima1, metaMaxima2, metaMinima, op, op2, metaIntersectingIntervals, doublePeakRegions, otherMarkFile=""):
	"""
	This function finds the double peak regions around MPRA peaks. 
	Args:
		MPRApeaks: MPRApeaks (datatype:list of Intervals)
		signal: is the histone Signal file (of datatype:bigWig)
		metaMaxima1: is the list with indices of all maxima in double peak regions
		metaMaxima2: is the list with all indices of 2nd maxima in double peak regions
		metaMinima: is the list with all indices of minima in double peak regions
		op: is file in which double peak characteristics are written
		op2: is bed file in which double peak intervals are written
		metaIntersectingIntervals: is list of signals within double peak regions
	Returns:
		None
	"""
	width = 2000
	binWidth = 25
	numberMinima = 0
	numberDoublePeaks = 0
	filteredPositives = []
	strands = []
	for currFeature in MPRApeaks:
		middle = (currFeature.start + currFeature.end)/2
		currFeature.start = middle - width/2
		currFeature.end = middle + width/2
		numBins = width/binWidth
		currSignal = signal.array([currFeature], bins=numBins)
		smoothedSignal = smoothSignal(currSignal[0])
		diffSig = calculateDifferential(smoothedSignal)
		minimaIdx = findMinima(diffSig)
		if minimaIdx != -1:
			numberMinima += 1
			maximaIndices = findSurroundingMaxima(diffSig, minimaIdx)
			if -1 not in maximaIndices:
				calculatePeakCharacteristics(filteredPositives, smoothedSignal, currFeature, maximaIndices, minimaIdx, metaMaxima1, metaMaxima2, metaMinima, metaIntersectingIntervals, doublePeakRegions, op, op2)
				numberDoublePeaks += 1

	print numberDoublePeaks, numberMinima
	return filteredPositives

def readShapeCharacteristics(allShapeCharacteristics, inp):
	ipFile = open(inp, "r")
	for line in ipFile:
		if "Maxima1" in line:
			continue
		currShapeChar = shapeCharacteristics()
		currShapeChar.parseCharacteristics(line.rstrip())
		allShapeCharacteristics.append(currShapeChar)
	ipFile.close()
	return

def plotCharacteristics(allShapeCharacteristics, pp):
	h1 = [currShapeChar.h1 for currShapeChar in allShapeCharacteristics]
	h2 = [currShapeChar.h2 for currShapeChar in allShapeCharacteristics]
	h3 = [currShapeChar.h3 for currShapeChar in allShapeCharacteristics]
	w1 = [currShapeChar.w1 for currShapeChar in allShapeCharacteristics]
	w2 = [currShapeChar.w2 for currShapeChar in allShapeCharacteristics]
	w = [currShapeChar.w for currShapeChar in allShapeCharacteristics]
	deltaH1 = [currShapeChar.ratio1 for currShapeChar in allShapeCharacteristics]
	deltaH2 = [currShapeChar.ratio2 for currShapeChar in allShapeCharacteristics]
	dis = [currShapeChar.dis for currShapeChar in allShapeCharacteristics]
	ratioW = [currShapeChar.ratioW for currShapeChar in allShapeCharacteristics]
	plt.figure()
	sns.distplot(h1, bins=25);
	plt.xlabel('height of maxima1', fontsize=14)
	plt.ylabel('# of doublepeaks', fontsize=14)
	plt.savefig(pp, facecolor='w', edgecolor='w', format='pdf')
	#plt.show()
	plt.figure()
	sns.distplot(h2, bins=25)
	plt.xlabel('height of minima', fontsize=14)
	plt.ylabel('# of doublepeaks', fontsize=14)
	plt.savefig(pp, facecolor='w', edgecolor='w', format='pdf')
	#plt.show()
	plt.figure()
	sns.distplot(h3, bins=25)
	plt.xlabel('height of maxima2', fontsize=14)
	plt.ylabel('# of doublepeaks', fontsize=14)
	plt.savefig(pp, facecolor='w', edgecolor='w', format='pdf')
	#plt.show()
	plt.figure()
	sns.distplot(w1, bins=25)
	plt.xlabel('width1 (bp)', fontsize=14)
	plt.ylabel('# of doublepeaks', fontsize=14)
	plt.savefig(pp, facecolor='w', edgecolor='w', format='pdf')
	#plt.show()
	plt.figure()
	sns.distplot(w2, bins=25)
	plt.xlabel('width2 (bp)', fontsize=14)
	plt.ylabel('# of doublepeaks', fontsize=14)
	plt.savefig(pp, facecolor='w', edgecolor='w', format='pdf')
	#plt.show()
	plt.figure()
	plt.hist(w, bins=25, alpha=0.5)
	plt.xlabel('total width (bp)', fontsize=14)
	plt.ylabel('# of doublepeaks', fontsize=14)
	plt.savefig(pp, facecolor='w', edgecolor='w', format='pdf')
	#plt.show()
	plt.figure()
	sns.distplot(deltaH1, bins=25)
	plt.xlabel('difference in height between maxima', fontsize=14)
	plt.ylabel('# of doublepeaks', fontsize=14)
	plt.grid(True)
	plt.savefig(pp, facecolor='w', edgecolor='w', format='pdf')
	#plt.show()
	plt.figure()
	sns.distplot(deltaH2, bins=25)
	plt.xlabel('difference in height between maxima and minima', fontsize=14)
	plt.ylabel('# of doublepeaks', fontsize=14)
	plt.savefig(pp, facecolor='w', edgecolor='w', format='pdf')
	#plt.show()
	plt.figure()
	sns.distplot(dis, bins=25)
	plt.xlabel('difference between minima and STARR-seq peak', fontsize=14)
	plt.ylabel('# of doublepeaks', fontsize=14)
	plt.savefig(pp, facecolor='w', edgecolor='w', format='pdf')
	#plt.show()
	plt.figure()
	plt.hist(ratioW, bins=25, alpha=0.5)
	plt.xlabel('w1/w', fontsize=14)
	plt.ylabel('# of doublepeaks', fontsize=14)
	plt.savefig(pp, facecolor='w', edgecolor='w', format='pdf')
	#plt.show()
	return

def getCurrentProfile(lowerBound, upperBound, currProfile, maxima1, maxima2, bins):
	x = np.arange(float(lowerBound), float(upperBound) + 1.0, 1.0)
	f = interpolate.splrep(x, np.array(currProfile[lowerBound:upperBound+1]), s=0)
	x1new = np.linspace(float(lowerBound), float(maxima1), num=20, endpoint=True)
	x2new = np.linspace(float(maxima1), float(maxima2), num=bins, endpoint=True)
	x3new = np.linspace(float(maxima2), float(upperBound), num=20, endpoint=True)

	y1new = list(interpolate.splev(x1new, f, der=0))
	y2new = list(interpolate.splev(x2new, f, der=0))
	y3new = list(interpolate.splev(x3new, f, der=0))
	ynew = y1new[0:-1] + y2new + y3new[1:]

	return ynew

def calculateMetaProfile(filteredPositives, metaIntersectingIntervals, bins, metaMaxima1, metaMaxima2, metaMinima, pp, op, op2, dependent=None):
	smoothedMetaProfile = []
	asymMetaProfile = []
	for idx in range(0, bins + 10):
		smoothedMetaProfile.append(0.0)
		asymMetaProfile.append(0.0)

	for idx in range(0, len(metaMaxima1)):
		if dependent == None:
			currIntervals = list(metaIntersectingIntervals[idx])
		else:
			currIntervals = list(dependent[idx])
		
		if metaMaxima1[idx] < 5:
			for i in range(0, 5-metaMaxima1[idx]):
				currIntervals.insert(0, currIntervals[0])
				lowerBound = 0
				maxima1 = 5
				maxima2 = metaMaxima2[idx] + (5 - metaMaxima1[idx])
		else:
			lowerBound = metaMaxima1[idx] - 5
			maxima1 = metaMaxima1[idx]
			maxima2 = metaMaxima2[idx]

		if metaMaxima2[idx] + 5 >= len(currIntervals):
			for i in range(0, metaMaxima2[idx] + 5 - len(currIntervals)):
				currIntervals.append(currIntervals[-1])
			upperBound = len(currIntervals) - 1
		else:
			upperBound = metaMaxima2[idx] + 5
		
		ynew = getCurrentProfile(lowerBound, upperBound, currIntervals, maxima1, maxima2, bins)
		for currIdx in range(10, bins + 20):
			smoothedMetaProfile[currIdx - 10] += ynew[currIdx]

		#if filteredPositives[idx].strand == "-":
		if metaIntersectingIntervals[idx][metaMaxima1[idx]] > metaIntersectingIntervals[idx][metaMaxima2[idx]]:
			currSignals = currIntervals[::-1]
			lowerBound2 = len(currIntervals) - upperBound
			upperBound2 = len(currIntervals) - lowerBound
			if lowerBound == 0:
				upperBound2 -= 1
			maximaIdx2 = len(currIntervals) - maxima1
			maximaIdx1 = len(currIntervals) - maxima2
			
			ynew = getCurrentProfile(lowerBound2, upperBound2, currSignals, maximaIdx1, maximaIdx2, bins)
			for currIdx in range(10, bins + 20):
				asymMetaProfile[currIdx - 10] += ynew[currIdx]
		else:
			for currIdx in range(10, bins + 20):
				asymMetaProfile[currIdx - 10] += ynew[currIdx]

	for currIdx in range(0, bins + 10):
		smoothedMetaProfile[currIdx] /= len(metaIntersectingIntervals)
		asymMetaProfile[currIdx] /= len(metaIntersectingIntervals)

	plt.figure()
	plt.plot(smoothedMetaProfile,'o')
	plt.plot(asymMetaProfile, 'x')
	plt.savefig(pp, format='pdf')
	for currProfilePoint in smoothedMetaProfile:
		op.write(str(currProfilePoint) + "\n")
	for currProfilePoint in asymMetaProfile:
		op2.write(str(currProfilePoint) + "\n")

	return smoothedMetaProfile

def getCorrespondingIntervals(signal, featureList):
	width = 2000
	binWidth = 25
	numBins = width/binWidth
	intervals = []
	for currFeature in featureList:
		currSignal = signal.array([currFeature], bins=numBins)
		intervals.append(currSignal[0])
	return intervals


def calculateDependentProfile(otherMarkFile, filteredPositives, metaIntersectingIntervals, bins, metaMaxima1, metaMaxima2, metaMinima, pp, opPrefix):
	markFiles = OrderedDict()
	corrIntervals = OrderedDict()
	ip = open(otherMarkFile, "r")
	for line in ip:
		fields = line.strip().split("\t")
		markFiles[fields[0]] = metaseq.genomic_signal(fields[1], "bigWig")
		corrIntervals[fields[0]] = getCorrespondingIntervals(markFiles[fields[0]], filteredPositives)
	ip.close()

	for currMark in markFiles:
		op = open(opPrefix + "_" + currMark + "_metaProfile.dat", "w")
		op2 = open(opPrefix + "_" + currMark + "_asymProfile.dat", "w")
		calculateMetaProfile(filteredPositives, metaIntersectingIntervals, bins, metaMaxima1, metaMaxima2, metaMinima, pp, op, op2, dependent=corrIntervals[currMark])
		op.close()
		op2.close()

	return



def main(histoneFile, MPRApeakFile, nonrandomFile, opPrefix, pp, otherMarkFile=None, plotChar=True):
	"""
	The main function of the program. 
	Calculates the pattern of troughs between STARR-seq peaks.
	Args:
		histoneFile: bigWigfile of histone signal
		MPRApeaks: contains peaks from massively parallel reporter assay for regulatory regions
		nonrandomFile: contains regions on which negatives should not intersect
		opPrefix: output prefix for all output files
	Returns:
		statistics for output
	"""

	#Checking input file
	try:
		MPRApeaks = pbt.BedTool(MPRApeakFile)
	except:
		sys.stderr.write("ERROR: Cannot open MPRA peak file " + MPRApeakFile + "\n")
		sys.exit()
	try:
		histoneSignal = metaseq.genomic_signal(histoneFile, "bigWig")
	except:
		sys.stderr.write("ERROR: Cannot open histone signal file " + histoneFile + "\n")
		sys.exit()
	try:
		nrFile = open(nonrandomFile, "r")
	except:
		sys.stderr.write("ERROR: " + nonrandomFile + " does not open\n")
		sys.exit()

	#Checking output files
	try:
		op = open(opPrefix + "_doublePeakStats.txt", "w")
		op2 = open(opPrefix + "_doublePeak.bed", "w")
	except:
		sys.stderr.write("ERROR: Cannot create output files\n")
		sys.exit()

	#Read nonrandom intervals
	nonrandom = bed.Bed()
	for line in nrFile:
		currInterval = extractFeature(line.rstrip())
		nonrandom.features.append(currInterval)
		del currInterval
	nonrandom.sortByChromosomeAndStartAndEnd()

	#Getting total number of 25 bp intervals that can be negative
	numNegativesTotal = calculateNumberNegatives(nonrandomFile, histoneFile, opPrefix)
	print "Number of negatives are", numNegativesTotal
	del histoneFile
	
	smoothingWindow = 2
	numberDoublePeaks = 0
	numberMinima = 0
	metaIntersectingIntervals = []
	metaMinima = []
	metaMaxima1 = []
	metaMaxima2 = []
	op.write("Maxima1\tMinima\tMaxima2\tdistanceMaxima1\tdistanceMaxima2\tratioMaxima\tratioMaximaMinima\tdistancePeakMinima\n")
	doublePeakRegions = []
	filteredPositives = getDoublePeakRegions(MPRApeaks, histoneSignal, metaMaxima1, metaMaxima2, metaMinima, op, op2, metaIntersectingIntervals, doublePeakRegions, otherMarkFile)
	op.close()
	op2.close()

	#Getting shape characteristics
	allShapeCharacteristics = []
	readShapeCharacteristics(allShapeCharacteristics, opPrefix + "_doublePeakStats.txt")
	if plotChar:
		plotCharacteristics(allShapeCharacteristics, pp)
	length1 = []
	length2 = []
	for idx in range(0, len(metaMaxima1)):
		length1.append(metaMinima[idx] - metaMaxima1[idx])
		length2.append(metaMaxima2[idx] - metaMinima[idx])
		if length1[-1] <= 0 or length2[-1] <= 0:
			print metaMaxima1[idx], metaMinima[idx], metaMaxima2[idx]
	
	#Calculating metaprofile
	try:
		op = open(opPrefix + "_metaProfile.dat", "w")
		op2 = open(opPrefix + "_asymProfile.dat", "w")
	except:
		sys.stderr.write("ERROR: Could not open " + opPrefix + "_metaProfile.dat\n")
		sys.exit()
	bins = 2 * max([max(length1), max(length2)])
	smoothedMetaProfile = calculateMetaProfile(filteredPositives, metaIntersectingIntervals, bins, metaMaxima1, metaMaxima2, metaMinima, pp, op, op2)
	op.close()
	op2.close()

	if otherMarkFile != None:
		calculateDependentProfile(otherMarkFile, filteredPositives, metaIntersectingIntervals, bins, metaMaxima1, metaMaxima2, metaMinima, pp, opPrefix)

	return

if __name__ == "__main__":

	if not(4 < len(sys.argv) < 7):
		sys.stderr.write("Usage: " + sys.argv[0] + " <histoneFile.bigWig> <MPRApeaks.bed> <nonrandomFile.bed> <opPrefix> [<otherMarks>]\n")
		sys.stderr.write("where:\n")
		sys.stderr.write("	<histoneFile.bigWig> is the histone signal file (bigWig format)\n")
		sys.stderr.write("	<MPRApeaks.bed> is the  peak file for MPRA (bed4 format with 4th column containing -log(p-value))\n")
		sys.stderr.write("	<nonrandomFile.bed> is the nonrandom file containing regions that negatives should not intersect with\n")
		sys.stderr.write("	<opPrefix> is the output prefix\n")
		sys.stderr.write("  <otherMarks> is a tab delimited file with information about dependent epigenetic marks\n")
		sys.exit()
	
	pp = PdfPages(sys.argv[4] + '_plots.pdf')
	if len(sys.argv) == 5:
		main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], pp)
	elif len(sys.argv) == 6:
		main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], pp, otherMarkFile=sys.argv[5])
	pp.close()
