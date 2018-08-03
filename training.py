import sys
import os

import numpy as np
from sklearn import svm
from scipy.optimize import curve_fit

def getScores(trainingScoreFile):
	#This function reads the scores for a given training file
	features = []
	try:
		ip = open(trainingScoreFile, "r")
	except:
		sys.stderr.write("Issue: Cannot open training file " + trainingScoreFile + "\n")
		sys.exit()
	
	header = ip.readline()
	fields = header.strip().split("\t")
	## Assuming the fields are chrom start end bestStart(MF) bestEnd(MF) followed by list of features
	features = fields[5:]
	featureScores = []
	for line in ip:
		fields = line.strip().split("\t")
		featureScores.append([])
		for idx in range(5, len(fields)):
			featureScores[-1].append(float(fields[idx]))

	return features, featureScores

def chooseRelevantColumns(scores, features, inputFiles):
	#This function chooses features based on input
	idx = 0
	relFeatures = []
	for currFeature in features:
		if currFeature in inputFiles:
			relFeatures.append(idx)
		idx += 1
	for currFeature in inputFiles:
		if currFeature not in features:
			sys.stderr.write("The feature " + currFeature + " does not have training data\n")
			sys.exit()

	#Now choosing just the relevant scores
	relScores = []
	for idx1 in range(0, len(scores)):
		relScores.append([])
		for idx2 in range(0, len(scores[0])):
			if idx2 in relFeatures:
				relScores[-1].append(scores[idx1][idx2])

	return relScores

def gauss_function(x, a, x0, sigma):
	#Gaussian function
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def calculateZscores(positiveScores, negativeScores, trainingNegatives):
	#Get Z-scores for training with respect to negative base Gaussian curve
	Zpos = [[] for i in range(len(positiveScores))]
	Zneg = [[] for i in range(len(negativeScores))]
	for idx in range(0, len(positiveScores[0])):
		subset = []
		idx2 = 0
		for currX in trainingNegatives:
			if currX[idx] == 0:
				continue
			subset.append(currX[idx])

		y, binEdges = np.histogram(subset, bins=50)
		x = np.zeros(50)
		y = y/float(len(subset))
		for i in range(0, len(binEdges) - 1):
			x[i] = (binEdges[i] + binEdges[i + 1])/2
			x = np.array(x)
		n = len(x)

		median = np.median(subset) #Slightly right-skewed. So median might be better starting point for fit
		sigma = np.std(subset)
		nopt,ncov = curve_fit(gauss_function, x, y, p0=[max(y), median, sigma])

		idx2 = 0
		for currScore in positiveScores:
			Zpos[idx2].append((currScore[idx] - nopt[1])/nopt[2])
			idx2 += 1

		idx2 = 0
		for currScore in negativeScores:
			Zneg[idx2].append((currScore[idx] - nopt[1])/nopt[2])
			idx2 += 1
		
	return Zpos, Zneg

def performSVM(trainingScores, trainingResults):
	print "->SVM training"

	clf = svm.SVC(kernel='linear', probability=True)
	clf.fit(trainingScores, trainingResults)
	#print clf.coef_
	print "<-SVM training"

	return clf
	
