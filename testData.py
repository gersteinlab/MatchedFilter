import sys
import os
from collections import OrderedDict

import metaseq
import pybedtools as pbt
import numpy as np
from scipy import interpolate

import scanMatchedFilter

def getTestRegions(bigWigFiles, bigWigList, currRegions):
	#Getting signals for test regions
	signalList = OrderedDict()
	for currSignal in bigWigFiles:
		signalList[currSignal] = scanMatchedFilter.getSignals(bigWigFiles, bigWigList, currRegions, currSignal)		
	 	print len(signalList[currSignal])
	 	del bigWigList[currSignal] 
		bigWigList[currSignal] =  metaseq.genomic_signal(bigWigFiles[currSignal], "bigWig")
		print len(signalList)
	
	return signalList

def getPossiblePairings(signalList, currRegions, metaprofiles, meanDist, stdDist, op):
	#Find pairs of maxima with appropriate distance constraints and their scores for the ML model
	binSize = 25
	minWidth = 350
	maxWidth = 1100
	y = OrderedDict()
	f = OrderedDict()
	Zscore = OrderedDict()
	numPoints = len(metaprofiles["H3K27ac"])
	profile = OrderedDict()
	op.write("#Region")
	for currSignal in signalList:
		profile[currSignal] = np.array(metaprofiles[currSignal][::-1])
		op.write("\t" + currSignal)
	op.write("\n")

	for idx in range(0, len(signalList["H3K27ac"])):
		if idx % 1000 ==0:
			print idx
		region = currRegions[idx]
		#Smoothing signal
		smoothedSignal = scanMatchedFilter.smoothSignals(signalList["H3K27ac"][idx], win=10)
		#Get maxima
		maximaIndices = scanMatchedFilter.getMaxima(smoothedSignal)
		if len(maximaIndices) < 2:
			continue
		#Get possible pairings of maxima for calculating MF scores
		pairings = scanMatchedFilter.findPossiblePairings(maximaIndices, binSize, minWidth, maxWidth)
		
		#Interpolating for curr region in different signals
		numBins = len(signalList["H3K27ac"][idx])
		x = np.linspace(region.start - 487.5, region.end + 487.5, num=numBins)
		
		for currSignal in signalList:
			y[currSignal] = signalList[currSignal][idx]
			#print len(x), len(y[currSignal])
			f[currSignal] = interpolate.interp1d(x, y[currSignal])
			#print len(x), len(y[currSignal])
		
		#Calculating MF score in different signals
		for currPairing in pairings:
			currWidth = (currPairing[1] - currPairing[0])*25
			finalWidth = 10.0/(numPoints - 10) * currWidth + currWidth
			start = region.start - 487.5 + (currPairing[1] + currPairing[0])*12.5 - float(finalWidth)/2
			end = start + finalWidth
			x1 = np.linspace(start, end, num=numPoints)	
			op.write(region.chrom + ":" + str(np.floor(start)) + "-" + str(np.ceil(end)))
			for currSignal in signalList:
				y1 = f[currSignal](x1)
				y1 = y1.T
				score = np.dot(profile[currSignal], y1)
				Zscore[currSignal] = (score - meanDist[currSignal][currWidth])/stdDist[currSignal][currWidth]
				op.write("\t" + str(Zscore[currSignal]))
			op.write("\n")
	return

def scoreTestData(opPrefix, SVM_model):
	#Score test data
	ip = open(opPrefix + "_testScores.dat" ,"r")
	X = []
	regions = []
	op1 = open(opPrefix + "_SVMpredScores.dat" ,"w")
	op2 = open(opPrefix + "_SVMpredictions.dat" ,"w")
	for line in ip:
		if "#" in line:
			continue
		fields = line.strip().split("\t")
		regions.append(fields[0])
		X.append([])
		for idx in range(1, len(fields)):
			X[-1].append(float(fields[idx]))
		if len(X) == 1000:
			y_pred = SVM_model.predict_proba(X)[:, 1]	
			for idx in range(0, len(X)):
				op1.write(regions[idx] + "\t" + str(y_pred[idx]) + "\n")
				if y_pred[idx] >= 0.85:
					op2.write(regions[idx] + "\t" + str(y_pred[idx]) + "\n")
			X = []
			regions = []
	if len(X) != 0:
		y_pred = SVM_model.predict_proba(X)[:, 1]	
		for idx in range(0, len(X)):
			op1.write(regions[idx] + "\t" + str(y_pred[idx]) + "\n")
			if y_pred[idx] >= 0.85:
				op2.write(regions[idx] + "\t" + str(y_pred[idx]) + "\n")
	
	op1.close()
	op2.close()
	return