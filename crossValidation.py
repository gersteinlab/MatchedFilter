#!/usr/bin/env python2.7

import sys
import os

import matplotlib
matplotlib.use('Agg')
import numpy as np
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import math
from matplotlib.backends.backend_pdf import PdfPages
import pylab as P
import metaseq
import matplotlib.pyplot as plt
import pybedtools as pbt
from scipy import interpolate
from collections import OrderedDict
import random
from sklearn import metrics
from matplotlib.font_manager import FontProperties
import seaborn as sns
import six
from matplotlib import colors
from sklearn import linear_model
from sklearn.ensemble import RandomForestClassifier 
from sklearn import svm
from sklearn.naive_bayes import GaussianNB


import bed
import interval
import matchedFilter
import createMetaPattern
import scorePredictions


def createSetsAndPattern(histoneFile, positiveFile, negativeFile, opPrefix, n, pp, otherMarkFile):
	
	positives = pbt.BedTool(positiveFile)
	negatives = pbt.BedTool(negativeFile)
	numNeg = len(negatives)

	p = 1.0/n
	for idx in range(0, n):
		op = open(opPrefix + "_positives" + str(idx) + ".bed", "w")
		op2 = open(opPrefix + "_testPos" + str(idx) + ".bed", "w")
		numPos = 0
		for currFeature in positives:
			if random.random() >= p:
				op.write(currFeature.chrom + "\t" + str(currFeature.start) + "\t" + str(currFeature.end) + "\n")
			else:
				op2.write(currFeature.chrom + "\t" + str(currFeature.start) + "\t" + str(currFeature.end) + "\n")
				numPos += 1
		op2.close()
		op.close()

		chosenNumNeg = min(numNeg, 10*numPos)
		pNeg = float(chosenNumNeg)/numNeg
		op = open(opPrefix + "_negatives" + str(idx) + ".bed", "w")
		op2 = open(opPrefix + "_testNeg" + str(idx) + ".bed", "w")
		for currFeature in negatives:
			if random.random() <= pNeg:
				op.write(currFeature.chrom + "\t" + str(currFeature.start) + "\t" + str(currFeature.end) + "\n")
			else:
				op2.write(currFeature.chrom + "\t" + str(currFeature.start) + "\t" + str(currFeature.end) + "\n")
		op2.close()
		op.close()
		if otherMarkFile == None:
			createMetaPattern.main(histoneFile, opPrefix + "_positives" + str(idx) + ".bed", opPrefix + "_negatives" + str(idx) + ".bed", opPrefix + str(idx), pp, otherMarkFile=None, plotChar=False)
		else:
			createMetaPattern.main(histoneFile, opPrefix + "_positives" + str(idx) + ".bed", opPrefix + "_negatives" + str(idx) + ".bed", opPrefix + str(idx), pp, otherMarkFile=otherMarkFile, plotChar=False)
		
	return

def scorePositivesAndNegatives(histoneFile, positiveFile, negativeFile, opPrefix, idx, pp, otherMarkSignal):
	
	dirName = os.path.dirname(os.path.realpath(__file__))
	if otherMarkSignal == None:
		metaProfileFile = opPrefix + str(idx) + "_metaProfile.dat"
		asymFile = opPrefix + str(idx) + "_asymProfile.dat"
		positives = opPrefix + "_positives" + str(idx) + ".bed"
		negatives = opPrefix + "_negatives" + str(idx) + ".bed"
		#metaProfileFile = "/Users/anurag/Projects/Enhancer/STARRseq/BG3/H3K27ac/BG3_H3K27ac_metaProfile.dat"
		#asymFile = "/Users/anurag/Projects/Enhancer/STARRseq/BG3/H3K27ac/BG3_H3K27ac_asymProfile.dat"
		#positives = "/Users/anurag/Projects/Enhancer/STARRseq/BG3/BG3_STARRseq_H3K27ac_relPeaks.bed"
		#negatives = "/Users/anurag/Projects/Enhancer/STARRseq/BG3/H3K27ac/BG3_H3K27ac_profile_negatives0.bed"
		os.system("python2.7 " + dirName + "/scorePredictions.py " + histoneFile + " " + opPrefix + "_positives" + str(idx) + ".bed " + metaProfileFile + " " + opPrefix + "_testPos" + str(idx) + " 0 1500")
		os.system("python2.7 " + dirName + "/scorePredictions.py " + histoneFile + " " + positives + " " + metaProfileFile + " " + opPrefix + "_positives" + str(idx) + " 0 1500")
		os.system("python2.7 " + dirName + "/scorePredictions.py " + histoneFile + " " + negatives + " " + metaProfileFile + " " + opPrefix + "_negatives" + str(idx) + " 0 0")
		os.system("python2.7 " + dirName + "/scorePredictions.py " + histoneFile + " " + opPrefix + "_negatives" + str(idx) + ".bed " + metaProfileFile + " " + opPrefix + "_testNeg" + str(idx) + " 0 0")
		os.system("python2.7 " + dirName + "/scorePredictions.py " + histoneFile + " " + opPrefix + "_positives" + str(idx) + ".bed " + asymFile + " " + opPrefix + "_testPos_asym" + str(idx) + " 1 1500")
		os.system("python2.7 " + dirName + "/scorePredictions.py " + histoneFile + " " + positives + " " + asymFile + " " + opPrefix + "_positives_asym" + str(idx) + " 1 1500")
		os.system("python2.7 " + dirName + "/scorePredictions.py " + histoneFile + " " + opPrefix + "_negatives" + str(idx) + ".bed " + asymFile + " " + opPrefix + "_testNeg_asym" + str(idx) + " 1 0")
		os.system("python2.7 " + dirName + "/scorePredictions.py " + histoneFile + " " + negatives + " " +  asymFile + " " + opPrefix + "_negatives_asym" + str(idx) + " 1 0")
	else:
		op = open(opPrefix + "_otherMarks_sym", "w")
		op2 = open(opPrefix + "_otherMarks_asym", "w")
		for currMark in otherMarkSignal:
			metaProfileFile = opPrefix + str(idx) + "_" + currMark + "_metaProfile.dat"
			asymFile = opPrefix + str(idx) + "_" + currMark + "_asymProfile.dat"
			#metaProfileFile = "/Users/anurag/Projects/Enhancer/STARRseq/S2/H3K27ac/S2_H3K27ac_" + currMark + "_metaProfile.dat"
			#asymFile = "/Users/anurag/Projects/Enhancer/STARRseq/S2/H3K27ac/S2_H3K27ac_" + currMark + "_asymProfile.dat"
			op.write(currMark + "\t" + otherMarkSignal[currMark] + "\t" + metaProfileFile + "\n")
			op2.write(currMark + "\t" + otherMarkSignal[currMark] + "\t" + asymFile + "\n")
		op.close()
		op2.close()
		metaProfileFile = opPrefix + str(idx) + "_metaProfile.dat"
		asymFile = opPrefix + str(idx) + "_asymProfile.dat"
		positives = opPrefix + "_positives" + str(idx) + ".bed"
		negatives = opPrefix + "_negatives" + str(idx) + ".bed"
		#metaProfileFile = "/Users/anurag/Projects/Enhancer/STARRseq/BG3/H3K27ac/BG3_H3K27ac_metaProfile.dat"
		#asymFile = "/Users/anurag/Projects/Enhancer/STARRseq/BG3/H3K27ac/BG3_H3K27ac_asymProfile.dat"
		#positives = "/Users/anurag/Projects/Enhancer/STARRseq/BG3/BG3_STARRseq_H3K27ac_relPeaks.bed"
		#negatives = "/Users/anurag/Projects/Enhancer/STARRseq/BG3/H3K27ac/BG3_H3K27ac_profile_negatives0.bed"
		print "python2.7 " + dirName + "/scorePredictions.py " + histoneFile + " " + opPrefix + "_testPos" + str(idx) + ".bed " +  metaProfileFile + " " +  opPrefix + "_testPos" + str(idx) + " 0 1500 " + opPrefix + "_otherMarks_sym"
		print "python2.7 " + dirName + "/scorePredictions.py " + histoneFile + " " + positives + " " +  metaProfileFile + " " +  opPrefix + "_positives" + str(idx) + " 0 1500 " + opPrefix + "_otherMarks_sym"
		os.system("python2.7 " + dirName + "/scorePredictions.py " + histoneFile + " " + opPrefix + "_testPos" + str(idx) + ".bed " +  metaProfileFile + " " +  opPrefix + "_testPos" + str(idx) + " 0 1500 " + opPrefix + "_otherMarks_sym")
		os.system("python2.7 " + dirName + "/scorePredictions.py " + histoneFile + " " + positives + " " +  metaProfileFile + " " +  opPrefix + "_positives" + str(idx) + " 0 1500 " + opPrefix + "_otherMarks_sym")
		os.system("python2.7 " + dirName + "/scorePredictions.py " + histoneFile + " " + opPrefix + "_testNeg" + str(idx) + ".bed " +  metaProfileFile + " " +  opPrefix + "_testNeg" + str(idx) + " 0 0 " + opPrefix + "_otherMarks_sym")
		os.system("python2.7 " + dirName + "/scorePredictions.py " + histoneFile + " " + negatives + " " +  metaProfileFile + " " +  opPrefix + "_negatives" + str(idx) + " 0 0 " + opPrefix + "_otherMarks_sym")
		os.system("python2.7 " + dirName + "/scorePredictions.py " + histoneFile + " " + opPrefix + "_testNeg" + str(idx) + ".bed " +  asymFile + " " +  opPrefix + "_testNeg_asym" + str(idx) + " 1 1500 " + opPrefix + "_otherMarks_asym" )
		os.system("python2.7 " + dirName + "/scorePredictions.py " + histoneFile + " " + positives + " " + asymFile + " " +  opPrefix + "_positives_asym" + str(idx) + " 1 1500 " + opPrefix + "_otherMarks_asym" )
		os.system("python2.7 " + dirName + "/scorePredictions.py " + histoneFile + " " + opPrefix + "_testPos" + str(idx) + ".bed " +  asymFile + " " +  opPrefix + "_testPos_asym" + str(idx) + " 1 0 " + opPrefix + "_otherMarks_asym")
		os.system("python2.7 " + dirName + "/scorePredictions.py " + histoneFile + " " + negatives + " " + asymFile + " " +  opPrefix + "_negatives_asym" + str(idx) + " 1 0 " + opPrefix + "_otherMarks_asym")


	return

def scoresCurrFile(filename, currMark):
	ip = open(filename, "r")
	header = ip.readline()
	fields = header.strip().split("\t")
	masterField = fields.index(currMark)
	currScores = []
	for line in ip:
		fields = line.strip().split("\t")
		currScores.append(float(fields[masterField]))
	ip.close()
	return currScores

def sortByChromosomeAndStartAndEnd(features):
	features.sort(key=lambda Interval:Interval.end)
	features.sort(key=lambda Interval:Interval.start)
	features.sort(key=lambda Interval:Interval.chr)
	return

def parseBedFile(peakFile, filetype="bed"):
	ip = open(peakFile, "r")
	intervalList = []
	for line in ip:
		currInterval = interval.Interval()
		fields = line.strip().split("\t")
		#print fields
		currInterval.chr = fields[0]
		currInterval.start = int(fields[1])
		currInterval.end = int(fields[2])
		if filetype=="bed":
			currInterval.prob = float(fields[3])
		intervalList.append(currInterval)
	sortByChromosomeAndStartAndEnd(intervalList)
	ip.close()
	return intervalList

def calculateRankList(intervalList):

	#print "->Calculating Rank list"
	region = OrderedDict()
	allKeys = []
	currPrediction = OrderedDict()
	for currInterval in intervalList:
		keyStr = currInterval.chr + ":" + str(currInterval.start) + "-" + str(currInterval.end)
		currPrediction[keyStr] = currInterval.prob
		allKeys.append(keyStr)


	probs = [value for (key, value) in currPrediction.iteritems()]
	probs.sort(reverse=True)
	idx = 0
	for currInterval in intervalList:
		currKeyStr = currInterval.chr + ":" + str(currInterval.start) + "-" + str(currInterval.end)
		currInterval.rank = probs.index(currPrediction[currKeyStr])
		idx += 1
	#print "<-Calculating Rank list"
	
	return

def findIntersectingRegions(intFile, peakList, filetype="bed"):
	ip = open(intFile, "r")
	prevRegion = ""
	regions = []
	currInterval = None
	for line in ip:
		fields = line.split("\t")
		currRegion = fields[0] + ":" + fields[1] + "-" + fields[2]
		if currRegion == prevRegion:
			if filetype == "narrowPeak":
				prob = float(fields[12])
			elif filetype == "peak":
				prob = float(fields[11])
			elif filetype == "bed":
				prob = float(fields[7])
			if prob > currInterval.prob:
				currInterval.chr = fields[4]
				currInterval.start = int(fields[5])
				currInterval.end = int(fields[6])
				currInterval.prob = prob
		else:
			if currInterval != None:
				regions.append(currInterval)
				currInterval = None
			currInterval = interval.Interval()
			currInterval.chr = fields[4]
			currInterval.start = int(fields[5])
			currInterval.end = int(fields[6])
			if filetype == "narrowPeak":
				currInterval.prob = float(fields[12])
			elif filetype == "peak":
				currInterval.prob = float(fields[11])
			elif filetype == "bed":
				currInterval.prob = float(fields[7])
			currInterval.result = int(fields[3])
			prevRegion = currRegion
			#print currRegion, len(regions)
	if currInterval != None:
		regions.append(currInterval)
	#print len(regions)
	sortByChromosomeAndStartAndEnd(regions)

	idx = 0
	for currInterval in regions:
		while idx < len(peakList) and peakList[idx].chr < currInterval.chr:
			idx += 1
		while idx < len(peakList) and peakList[idx].chr == currInterval.chr and peakList[idx].start < currInterval.start:
			idx += 1

		if peakList[idx].chr == currInterval.chr and peakList[idx].start == currInterval.start and peakList[idx].end == currInterval.end:
			currInterval.rank = peakList[idx].rank
			#print idx
		else:
			print "Issue", idx, peakList[idx].chr, peakList[idx].start, peakList[idx].end, currInterval.chr, currInterval.start, currInterval.end

	return regions

def statisticsIntersectingRegions(filename):
	ip = open(filename, "r")
	numPos = 0
	numNeg = 0
	for line in ip:
		fields = line.strip().split("\t")
		if int(fields[3]) == 1:
			numPos += 1
		elif int(fields[3]) == 0:
			numNeg += 1
		else:
			print "Issue : " + fields[3]
	return numPos, numNeg


def calculatePeakAUC(positiveFile, negativeFile, peakFile, opPrefix, currMark):
	allPeaks = parseBedFile(peakFile, filetype="bed")
	calculateRankList(allPeaks)

	os.system("awk '{print $1 \"\t\" $2 \"\t\" $3 \"\t1\"}' " + positiveFile + " |grep -v \# > temp1.bed")
	os.system("awk '{print $1 \"\t\" $2 \"\t\" $3 \"\t0\"}' " + negativeFile + " |grep -v \# | cat temp1.bed - > annotation.bed")

	if "DHS" in currMark:
		os.system("intersectBed -a annotation.bed -b " + peakFile + " -wb -wa -f 0.10 | sortBed -i - > " + opPrefix + "_intersect.bed")
		os.system("intersectBed -a annotation.bed -b " + peakFile + " -wa -v -f 0.10 | sortBed -i -  > " + opPrefix + "_notIntersect.bed")
	else:
		os.system("intersectBed -a annotation.bed -b " + peakFile + " -wb -wa -f 0.50 | sortBed -i - > " + opPrefix + "_intersect.bed")
		os.system("intersectBed -a annotation.bed -b " + peakFile + " -wa -v -f 0.50 | sortBed -i -  > " + opPrefix + "_notIntersect.bed")
	intersectingPeaks = findIntersectingRegions(opPrefix + "_intersect.bed", allPeaks, filetype="bed")
	donotIntersectPos, donotIntersectNeg = statisticsIntersectingRegions(opPrefix + "_notIntersect.bed")
	rankingPeaks = [currInterval.rank for currInterval in intersectingPeaks]
	result = [currInterval.result for currInterval in intersectingPeaks]
	nonintersectingResult = donotIntersectPos * [1] + donotIntersectNeg * [0]
	nonintersectingRank = len(nonintersectingResult) * [100000000]
	fpr, tpr, auc, precision, recall, average_precision = calculateMetrics(1.0/np.array(rankingPeaks + nonintersectingRank), result + nonintersectingResult, "Peaks")
	print currMark + " " + str(auc) + " " + str(average_precision)
	return fpr, tpr, auc, precision, recall, average_precision

def calculateMetrics(scores, results, currMark):
	numPts = len(results)
	
	fpr, tpr, _ = metrics.roc_curve(np.array(results).reshape(numPts,1), np.array(scores).reshape(numPts,1))
	auc = metrics.auc(fpr,tpr)
	precision, recall, thresholds = metrics.precision_recall_curve(np.array(results).reshape(numPts,1), np.array(scores).reshape(numPts,1))	
	average_precision = metrics.average_precision_score(np.array(results).reshape(numPts,1), np.array(scores).reshape(numPts,1))


	return fpr, tpr, auc, precision, recall, average_precision


def performLinearRegression(trainingScores, trainingResults, testScores):
	
	y_pred = OrderedDict()
	y_pred2 = OrderedDict()
	for currMark in trainingScores:
		if "Master" in currMark:
			continue
		if "Asym" in currMark:
			continue
		X = []
		for idx in range(0, len(trainingScores["MasterAsym"])):
			X.append([trainingScores["MasterAsym"][idx], trainingScores[currMark][idx]])

		regr = linear_model.Ridge(fit_intercept=True,  copy_X=True, normalize=False, n_jobs=1)
		regr.fit(X, np.array(trainingResults))

		X_test = []
		y_pred2["Master + " + currMark] = []
		for idx in range(0, len(testScores["MasterAsym"])):
			X_test.append([testScores["MasterAsym"][idx], testScores[currMark][idx]])
			y_pred2["Master + " + currMark].append(testScores["MasterAsym"][idx] + testScores[currMark][idx])
		y_pred["Master + " + currMark] = list(regr.predict(X_test))

	return y_pred2

def performRandomForest(trainingScores, trainingResults, testScores):

	X = []
	for idx in range(0, len(trainingScores["MasterAsym"])):
		X.append([])

	for currMark in trainingScores:
		#if "Master" in currMark:
		#	continue
		if "Asym" in currMark:
			continue
		#print currMark, trainingScores[currMark][0]
		print currMark, 
		for idx in range(0, len(trainingScores[currMark])):
			X[idx].append(trainingScores[currMark][idx])

	X_test = []
	for idx in range(0, len(testScores["MasterAsym"])):
		X_test.append([])

	for currMark in trainingScores:
		if "Asym" in currMark:
			continue
		#print currMark, testScores[currMark][0]
		for idx in range(0, len(testScores[currMark])):
			X_test[idx].append(testScores[currMark][idx])

	forest = RandomForestClassifier(n_estimators=1000, criterion="gini")
	forest = forest.fit(X, np.array(trainingResults))
	y_pred = forest.predict_proba(X_test)[:, 1]
	#print y_pred
	print forest.feature_importances_
	return y_pred, forest.feature_importances_

def performSVM(trainingScores, trainingResults, testScores):
	print "->SVM"
	X = []
	for currMark in trainingScores:
		pass
	for idx in range(0, len(trainingScores[currMark])):
		X.append([])

	for currMark in trainingScores:
		if "Asym" in currMark:
			continue
		print currMark, 
		for idx in range(0, len(trainingScores[currMark])):
			X[idx].append(trainingScores[currMark][idx])

	X_test = []
	for idx in range(0, len(testScores[currMark])):
		X_test.append([])

	for currMark in trainingScores:
		if "Asym" in currMark:
			continue
		for idx in range(0, len(testScores[currMark])):
			X_test[idx].append(testScores[currMark][idx])
	clf = svm.SVC(kernel='linear', probability=True)
	clf.fit(X, np.array(trainingResults))
	print clf.coef_
	y_pred = clf.predict_proba(X_test)[:, 1]
	print "<-SVM"

	return y_pred, clf.coef_[0]

def performLR(trainingScores, trainingResults, testScores):
	print "->LR"
	X = []
	for currMark in trainingScores:
		pass
	for idx in range(0, len(trainingScores[currMark])):
		X.append([])

	for currMark in trainingScores:
		if "Asym" in currMark:
			continue
		print currMark, 
		for idx in range(0, len(trainingScores[currMark])):
			X[idx].append(trainingScores[currMark][idx])

	X_test = []
	for idx in range(0, len(testScores[currMark])):
		X_test.append([])

	for currMark in trainingScores:
		if "Asym" in currMark:
			continue
		for idx in range(0, len(testScores[currMark])):
			X_test[idx].append(testScores[currMark][idx])
	regr = linear_model.LinearRegression(fit_intercept=True,  copy_X=True, normalize=False, n_jobs=1)
	regr.fit(X, np.array(trainingResults))
	y_pred = list(regr.predict(X_test))
	print regr.coef_
	print "<-LR"
	return y_pred, regr.coef_

def performNB(trainingScores, trainingResults, testScores):
	print "->Gaussian NB"
	X = []
	for currMark in trainingScores:
		pass
	for idx in range(0, len(trainingScores[currMark])):
		X.append([])

	for currMark in trainingScores:
		if "Asym" in currMark:
			continue
		print currMark, 
		for idx in range(0, len(trainingScores[currMark])):
			X[idx].append(trainingScores[currMark][idx])

	X_test = []
	for idx in range(0, len(testScores[currMark])):
		X_test.append([])

	for currMark in trainingScores:
		if "Asym" in currMark:
			continue
		for idx in range(0, len(testScores[currMark])):
			X_test[idx].append(testScores[currMark][idx])
	gnb = GaussianNB()
	gnb.fit(X, np.array(trainingResults))
	y_pred = gnb.predict_proba(X_test)[:, 1]
	print "->Gaussian NB"
	return y_pred

def gauss_function(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def renormalizeScores(positiveScores, negativeScores, title):
	minimumScore = min(positiveScores + negativeScores)
	maximumScore = max(positiveScores + negativeScores)
	renPos = []
	for currScore in positiveScores:
		renPos.append((currScore-minimumScore)/(maximumScore-minimumScore))
	renNeg = []
	for currScore in negativeScores:
		renNeg.append((currScore-minimumScore)/(maximumScore-minimumScore))
	plt.figure()
	plt.hist(renPos, bins=50)
	plt.xlabel("score")
	plt.hist(renNeg,bins=50)
	plt.title(title)
	plt.savefig(pp, format="pdf")
	return renPos, renNeg

def calculateZscores(positiveScores, negativeScores, trainingNegatives, pp, title):
	subset = []
	for currX in trainingNegatives:
		if currX == 0:
			continue
		subset.append(currX)

	y, binEdges = np.histogram(trainingNegatives, bins=50)
	x = np.zeros(50)
	y = y/float(len(trainingNegatives))
	for i in range(0, len(binEdges) - 1):
		x[i] = (binEdges[i] + binEdges[i + 1])/2
		x = np.array(x)
	#plt.show()
	n = len(x)

	median = np.median(trainingNegatives) #Slightly right-skewed. So median might be better starting point for fit
	sigma = np.std(trainingNegatives)

	#nopt,ncov = curve_fit(gauss_function, x, y, p0=[1, median, sigma])

	Zpos = []
	for currScore in positiveScores:
		Zpos.append((currScore - median)/sigma)

	Zneg = []
	for currScore in negativeScores:
		Zneg.append((currScore - median)/sigma)
	plt.figure()
	plt.hist(Zpos, bins=50, weights=np.array([1.0/len(Zpos)] * len(Zpos)))
	plt.xlabel("score")
	plt.hist(Zneg,bins=50, weights=np.array([1.0/len(Zneg)] * len(Zneg)))
	plt.title(title)
	plt.savefig(pp, format="pdf")
	return Zpos, Zneg


def readScores(opPrefix, n, positiveFile, negativeFile, peakFile, otherMarkSignal, otherMarkPeak):
	scores = []
	results = []
	trainingScores = []
	trainingResults = []
	totalAuc, totalAupr = OrderedDict(), OrderedDict()
	fpr = OrderedDict()
	tpr = OrderedDict()
	prec = OrderedDict()
	rec = OrderedDict()
	colorIdx = OrderedDict()
	totalAuc["Master"], fpr["Master"], tpr["Master"] = 0, [], []
	totalAuc["MasterAsym"], fpr["MasterAsym"], tpr["MasterAsym"] = 0, [], []
	totalAuc["MasterPeak"], fpr["MasterPeak"], tpr["MasterPeak"] = 0, [], []
	totalAuc["RandomForest"], fpr["RandomForest"], tpr["RandomForest"] = 0, [], []
	totalAupr["Master"], prec["Master"], rec["Master"] = 0, [], []
	totalAupr["MasterAsym"], prec["MasterAsym"], rec["MasterAsym"] = 0, [], []
	totalAupr["MasterPeak"], prec["MasterPeak"], rec["MasterPeak"] = 0, [], []
	totalAupr["RandomForest"], prec["RandomForest"], rec["RandomForest"] = 0, [], []
	colorIdx["Master"] = 0
	colorIdx["MasterAsym"] = 0
	colorIdx["MasterPeak"] = 0
	colorIdx["RandomForest"] = 0
	colorIdx["SVM"] = 1
	totalAuc["SVM"], fpr["SVM"], tpr["SVM"] = 0, [], []
	totalAupr["SVM"], prec["SVM"], rec["SVM"] = 0, [], []
	colorIdx["LR"] = 2
	totalAuc["LR"], fpr["LR"], tpr["LR"] = 0, [], []
	totalAupr["LR"], prec["LR"], rec["LR"] = 0, [], []
	colorIdx["NB"] = 3
	totalAuc["NB"], fpr["NB"], tpr["NB"] = 0, [], []
	totalAupr["NB"], prec["NB"], rec["NB"] = 0, [], []
	idx = 1
	peakMarks = {"MasterPeak": peakFile}
	featureNames = ["Master"]
	if otherMarkSignal != None:
		for currMark in otherMarkSignal:
			print currMark
			featureNames.append(currMark)
			totalAuc[currMark], fpr[currMark], tpr[currMark] = 0, [], []
			totalAuc[currMark + "Asym"], fpr[currMark + "Asym"], tpr[currMark + "Asym"] = 0, [], []
			totalAuc[currMark + "Peak"], fpr[currMark + "Peak"], tpr[currMark + "Peak"] = 0, [], []
			totalAupr[currMark], prec[currMark], rec[currMark] = 0, [], []
			totalAupr[currMark + "Asym"], prec[currMark + "Asym"], rec[currMark + "Asym"] = 0, [], []
			totalAupr[currMark + "Peak"], prec[currMark + "Peak"], rec[currMark + "Peak"] = 0, [], []
			colorIdx[currMark] = idx
			colorIdx[currMark + "Asym"] = idx
			colorIdx[currMark + "Peak"] = idx
			peakMarks[currMark + "Peak"] = otherMarkPeak[currMark]
			totalAuc["Master + " + currMark], fpr["Master + " + currMark], tpr["Master + " + currMark] = 0, [], []
			totalAupr["Master + " + currMark], prec["Master + " + currMark], rec["Master + " + currMark] = 0, [], []
			colorIdx["Master + " + currMark] = idx
			idx += 1

	colors_ = list(six.iteritems(colors.cnames))
	for name, rgb in six.iteritems(colors.ColorConverter.colors):
		hex_ = colors.rgb2hex(rgb)
		colors_.append((name, hex_))
	hex_ = [color[1] for color in colors_]
	rgb = [colors.hex2color(color) for color in hex_]
	hsv = [colors.rgb_to_hsv(color) for color in rgb]
	hue = [color[0] for color in hsv]
	sat = [color[1] for color in hsv]
	val = [color[2] for color in hsv]
	ind = np.lexsort((val, sat, hue))
	sorted_colors = [colors_[i] for i in ind]
	usedColors = []
	n1 = len(sorted_colors)
	current_palette = sns.color_palette(n_colors=len(otherMarkSignal) + 1)
	sns.palplot(current_palette)
	for i in range(0, len(otherMarkSignal) + 1):
		usedColors.append(sorted_colors[int(i*float(n1)/(len(otherMarkSignal) + 1))])
	c = current_palette

	for idx in range(0, n):
		#print idx, n
		scores.append(OrderedDict())
		trainingScores.append(OrderedDict())
		trainingResults.append([])
		results.append([])
		trainingScores1 = scoresCurrFile(opPrefix + "_positives" + str(idx) + "_MFscores.bed", "Master")
		trainingResults1 = [1] * len(trainingScores1) 
		trainingScores3 = scoresCurrFile(opPrefix + "_positives_asym" + str(idx) + "_MFscores.bed", "Master")
		trainingScores2 = scoresCurrFile(opPrefix + "_negatives" + str(idx) + "_MFscores.bed", "Master")
		trainingResults2 = [0] * len(trainingScores2) 
		trainingScores4 = scoresCurrFile(opPrefix + "_negatives_asym" + str(idx) + "_MFscores.bed", "Master")
		#renTrainScores1, renTrainScores2 = renormalizeScores(trainingScores1, trainingScores2, "Master")
		renTrainScores1, renTrainScores2 = calculateZscores(trainingScores1, trainingScores2, trainingScores2, pp, "trainMaster")
		renTrainScores3, renTrainScores4 = calculateZscores(trainingScores3, trainingScores4, trainingScores4, pp, "trainMasterAsym")
		trainingScores[idx]["Master"] = renTrainScores1 + renTrainScores2
		trainingScores[idx]["MasterAsym"] = trainingScores3 + trainingScores4
		trainingResults[idx] = trainingResults1 + trainingResults2

		scores1 = scoresCurrFile(opPrefix + "_testPos" + str(idx) + "_MFscores.bed", "Master")
		results1 = [1] * len(scores1) 
		scores3 = scoresCurrFile(opPrefix + "_testPos_asym" + str(idx) + "_MFscores.bed", "Master")
		scores2 = scoresCurrFile(opPrefix + "_testNeg" + str(idx) + "_MFscores.bed", "Master")
		results2 = [0] * len(scores2) 
		scores4 = scoresCurrFile(opPrefix + "_testNeg_asym" + str(idx) + "_MFscores.bed", "Master")
		#renTestScores1, renTestScores2 = renormalizeScores(scores1, scores2, "Master")
		renTestScores1, renTestScores2 = calculateZscores(scores1, scores2, scores2, pp, "testMaster")
		renTestScores3, renTestScores4 = calculateZscores(scores3, scores4, scores4, pp, "testMasterAsym")
		scores[idx]["Master"] = renTestScores1 + renTestScores2
		scores[idx]["MasterAsym"] = scores3 + scores4
		results[idx] = results1 + results2

		if otherMarkSignal != None:
			for currMark in otherMarkSignal:
				#print opPrefix + "_testPos" + str(idx) + "_MFscores.bed"
				scores1 = scoresCurrFile(opPrefix + "_testPos" + str(idx) + "_MFscores.bed", currMark)
				scores3 = scoresCurrFile(opPrefix + "_testPos_asym" + str(idx) + "_MFscores.bed", currMark)
				scores2 = scoresCurrFile(opPrefix + "_testNeg" + str(idx) + "_MFscores.bed", currMark)
				scores4 = scoresCurrFile(opPrefix + "_testNeg_asym" + str(idx) + "_MFscores.bed", currMark)		
				trainingScores1 = scoresCurrFile(opPrefix + "_positives" + str(idx) + "_MFscores.bed", currMark)
				trainingScores3 = scoresCurrFile(opPrefix + "_positives_asym" + str(idx) + "_MFscores.bed", currMark)
				trainingScores2 = scoresCurrFile(opPrefix + "_negatives" + str(idx) + "_MFscores.bed", currMark)
				trainingScores4 = scoresCurrFile(opPrefix + "_negatives_asym" + str(idx) + "_MFscores.bed", currMark)
				#renTrainScores1, renTrainScores2 = renormalizeScores(trainingScores1, trainingScores2, currMark)
				renTrainScores1, renTrainScores2 = calculateZscores(trainingScores1, trainingScores2, trainingScores2, pp, "train" + currMark)
				renTrainScores3, renTrainScores4 = calculateZscores(trainingScores3, trainingScores4, trainingScores4, pp, "train" + currMark + "Asym")
				#trainingScores[idx][currMark] = renTrainScores1 + renTrainScores2
				#trainingScores[idx][currMark + "Asym"] = renTrainScores3 + renTrainScores4
				trainingScores[idx][currMark] = renTrainScores1 + renTrainScores2
				trainingScores[idx][currMark + "Asym"] = trainingScores3 + trainingScores4
				renTestScores1, renTestScores2 = calculateZscores(scores1, scores2, scores2, pp, "test" + currMark)
				renTestScores3, renTestScores4 = calculateZscores(scores3, scores4, scores4, pp, "test" + currMark + "Asym")
				scores[idx][currMark] = renTestScores1 + renTestScores2
				scores[idx][currMark + "Asym"] = scores3 + scores4

		for currMark in scores[idx]:

			currFpr, currTpr, currRoc_auc, currPrec, currRec, curr_aupr = calculateMetrics(scores[idx][currMark], results[idx], currMark)
			currFpr = list(currFpr)
			currTpr = list(currTpr)
			currPrec = list(currPrec)
			currRec = list(currRec)
			if idx == 0:
				for element in currFpr:
					fpr[currMark].append(element)
				for element in currTpr:
					tpr[currMark].append(element)
				for element in currPrec:
					prec[currMark].append(element)
				for element in currRec:
					rec[currMark].append(element)
			totalAuc[currMark] += currRoc_auc
			totalAupr[currMark] += curr_aupr


	for idx in range(0, n):
		for currMark in peakMarks:
			currFpr, currTpr, currRoc_auc, currPrec, currRec, curr_aupr = calculatePeakAUC(opPrefix + "_testPos" + str(idx) + "_MFscores.bed", opPrefix + "_testNeg" + str(idx) + "_MFscores.bed", peakMarks[currMark], opPrefix, currMark)
			currFpr = list(currFpr)
			currTpr = list(currTpr)
			currPrec = list(currPrec)
			currRec = list(currRec)
			if idx == 0:
				for element in currFpr:
					fpr[currMark].append(element)
				for element in currTpr:
					tpr[currMark].append(element)
				for element in currPrec:
					prec[currMark].append(element)
				for element in currRec:
					rec[currMark].append(element)
			totalAuc[currMark] += currRoc_auc
			totalAupr[currMark] += curr_aupr

	if otherMarkSignal != None:
		scores_pred = []
		featureImportance = []
		SVMcoeff = []
		LRcoeff = []
		for idx in range(0, n):
			scores_pred.append(OrderedDict())
			
			#scores_pred.append(performLinearRegression(scores[idx], results[idx], scores[idx]))
			scores_pred[idx]["RandomForest"], currFeatureImportance = performRandomForest(trainingScores[idx], trainingResults[idx], scores[idx])
			featureImportance.append(currFeatureImportance)
			scores_pred[idx]["SVM"], currFeatureImportance = performSVM(trainingScores[idx], trainingResults[idx], scores[idx])
			SVMcoeff.append(currFeatureImportance)
			scores_pred[idx]["LR"], currFeatureImportance = performLR(trainingScores[idx], trainingResults[idx], scores[idx])
			LRcoeff.append(currFeatureImportance)
			scores_pred[idx]["NB"] = performNB(trainingScores[idx], trainingResults[idx], scores[idx])
			
			for currMark in scores_pred[idx]:
				currFpr, currTpr, currRoc_auc, currPrec, currRec, curr_aupr = calculateMetrics(scores_pred[idx][currMark], results[idx], currMark)
				currFpr = list(currFpr)
				currTpr = list(currTpr)
				currPrec = list(currPrec)
				currRec = list(currRec)
				if idx == 0:
					for element in currFpr:
						fpr[currMark].append(element)
					for element in currTpr:
						tpr[currMark].append(element)
					for element in currPrec:
						prec[currMark].append(element)
					for element in currRec:
						rec[currMark].append(element)
				totalAuc[currMark] += currRoc_auc
				totalAupr[currMark] += curr_aupr
		featureImportanceMean = np.mean(featureImportance, axis=0)
		featureImportanceStd = np.std(featureImportance, axis = 0)
		SVMcoeffMean = np.mean(SVMcoeff, axis=0)
		SVMcoeffStd = np.std(SVMcoeff, axis = 0)
		LRcoeffMean = np.mean(LRcoeff, axis=0)
		LRcoeffStd = np.std(LRcoeff, axis = 0)
		plt.figure()
		fig, ax = plt.subplots()
		ind = np.arange(len(featureNames))
		width = 0.5
		rects1 = ax.bar(ind, featureImportanceMean, width, color='g', yerr=featureImportanceStd)
		ax.set_xlabel('Feature')
		ax.set_ylabel('Feature Importance')
		ax.set_title('Importance of feature in Random Forest Model')
		ax.set_xticks(ind + width/2)
		ax.set_xticklabels(tuple(featureNames))
		plt.savefig(pp, format="pdf")
		fig, ax = plt.subplots()
		ind = np.arange(len(featureNames))
		width = 0.5
		rects1 = ax.bar(ind, SVMcoeffMean, width, color='g', yerr=SVMcoeffStd)
		ax.set_xlabel('Feature')
		ax.set_ylabel('Feature Importance')
		ax.set_title('Coefficients in SVM Model')
		ax.set_xticks(ind + width/2)
		ax.set_xticklabels(tuple(featureNames))
		plt.savefig(pp, format="pdf")
		fig, ax = plt.subplots()
		ind = np.arange(len(featureNames))
		width = 0.5
		rects1 = ax.bar(ind, LRcoeffMean, width, color='g', yerr=LRcoeffStd)
		ax.set_xlabel('Feature')
		ax.set_ylabel('Feature Importance')
		ax.set_title('Coefficients in LR Model')
		ax.set_xticks(ind + width/2)
		ax.set_xticklabels(tuple(featureNames))
		plt.savefig(pp, format="pdf")

	plt.figure()
	ax = plt.subplot(111)
	for currMark in totalAuc:
		totalAuc[currMark] /= n
		if "Master" in currMark and "Asym" not in currMark and "Peak" not in currMark and "+" not in currMark:
			#ax.plot(fpr[currMark], tpr[currMark], "--", linewidth=1, c=c[colorIdx[currMark]], label=currMark+' (area = %0.2f)' % totalAuc[currMark])
			pass
		elif "Peak" in currMark:
			#ax.plot(fpr[currMark], tpr[currMark], "-.", linewidth=1, c=c[colorIdx[currMark]], label=currMark+' Peak (area = %0.2f)' % totalAuc[currMark])
			pass
		elif "Asym" not in currMark and "Peak" not in currMark and "+" not in currMark and currMark not in ["RandomForest", "SVM", "LR", "NB"]:
			#ax.plot(fpr[currMark], tpr[currMark], "--", linewidth=1, c=c[colorIdx[currMark]], label=currMark+' (area = %0.2f)' % totalAuc[currMark])
			pass
		elif "Asym" not in currMark and "+" not in currMark and currMark not in ["RandomForest", "SVM", "LR", "NB"]:
			pass
			#ax.plot(fpr[currMark], tpr[currMark], "-.", linewidth=1, c=c[colorIdx[currMark]], label=currMark+' (area = %0.2f)' % totalAuc[currMark])
		elif currMark in ["RandomForest", "SVM", "LR", "NB"]:# or "Random" in currMark:
			#pass
			ax.plot(fpr[currMark], tpr[currMark], ":", linewidth=1, c=c[colorIdx[currMark]], label=currMark+' (area = %0.2f)' % totalAuc[currMark])
		else:
			pass
			#ax.plot(fpr[currMark], tpr[currMark], "--", linewidth=1, c=c[colorIdx[currMark]], label=currMark+' (area = %0.2f)' % totalAuc[currMark])
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.75, box.height * 0.75])
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=7)
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.title("ROC Plot")
	#plt.legend(loc="lower right")
	plt.savefig(pp, format="pdf")

	plt.figure()
	ax = plt.subplot(111)
	for currMark in totalAupr:
		totalAupr[currMark] /= n
		if "Master" in currMark and "Asym" not in currMark and "Peak" not in currMark and "+" not in currMark:
			#ax.plot(rec[currMark], prec[currMark], "--", linewidth=1, c=c[colorIdx[currMark]], label=currMark+' (area = %0.2f)' % totalAupr[currMark])
			pass
		elif "Peak" in currMark:
			pass
			#ax.plot(rec[currMark], prec[currMark], "-.", linewidth=1, c=c[colorIdx[currMark]], label=currMark+' Peak (area = %0.2f)' % totalAupr[currMark])
		elif "Asym" not in currMark and "Peak" not in currMark and "+" not in currMark and currMark not in ["RandomForest", "SVM", "LR", "NB"]:
			pass
			#ax.plot(rec[currMark], prec[currMark], "--", linewidth=1, c=c[colorIdx[currMark]], label=currMark+' (area = %0.2f)' % totalAupr[currMark])
		elif "Asym" not in currMark and "+" not in currMark and currMark not in ["RandomForest", "SVM", "LR", "NB"]:
			pass
			#ax.plot(rec[currMark], prec[currMark], "-.", linewidth=1, c=c[colorIdx[currMark]], label=currMark+' (area = %0.2f)' % totalAupr[currMark])
		elif currMark in ["RandomForest", "SVM", "LR", "NB"]: #or "Random" in currMark:
			#pass
			ax.plot(rec[currMark], prec[currMark], ":", linewidth=2, c=c[colorIdx[currMark]], label=currMark+' (area = %0.2f)' % totalAupr[currMark])
		else:
			pass
			#ax.plot(rec[currMark], prec[currMark], "--", linewidth=1, c=c[colorIdx[currMark]], label=currMark+' (area = %0.2f)' % totalAupr[currMark])

		#plt.plot(rec[currMark], prec[currMark], linewidth=2, label=currMark+' (area = %0.2f)' % totalAupr[currMark])
	plt.xlabel('Recall')
	plt.ylabel('Precision')
	plt.title("PR Plot")
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.0])
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.75, box.height * 0.75])
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=7)
	#plt.legend(loc="upper right")
	plt.savefig(pp, format="pdf")
	#sys.exit()
	#print "Final AUC (symmetric) = ", totalAuc["Master"]
	#print "Final AUC (asymmetric) = ", totalAuc["MasterAsym"]

	return 


def main(histoneFile, positiveFile, negativeFile, peakFile, opPrefix, n, pp, otherMarkFile=None):

	#createSetsAndPattern(histoneFile, positiveFile, negativeFile, opPrefix, n, pp, otherMarkFile)
	n = 1
	
	if otherMarkFile != None:
		otherMarkSignal = OrderedDict()
		otherMarkPeaks = OrderedDict()
		ip = open(otherMarkFile)
		for line in ip:
			fields = line.strip().split("\t")
			otherMarkSignal[fields[0]] = fields[1]
			otherMarkPeaks[fields[0]] = fields[2]
		ip.close()

	for idx in range(0, n):
	  	if otherMarkFile == None:
	  		scorePositivesAndNegatives(histoneFile, positiveFile, negativeFile, opPrefix, idx, pp, otherMarkFile)
	  	else:
	  		scorePositivesAndNegatives(histoneFile, positiveFile, negativeFile, opPrefix, idx, pp, otherMarkSignal)
	
	if otherMarkFile != None:
		readScores(opPrefix, n, positiveFile, negativeFile, peakFile, otherMarkSignal=otherMarkSignal, otherMarkPeak=otherMarkPeaks)
	else:
		readScores(opPrefix, n, positiveFile, negativeFile, peakFile, otherMarkSignal=None, otherMarkPeak=None)

	return


if __name__ == "__main__":
	if not(5 < len(sys.argv) < 9):
		sys.stderr.write("Usage: " + sys.argv[0] + " <histoneFile.bigWig> <positives.bed> <negatives.bed> <peaks.bed> <opPrefix> [<n> [<otherFiles>]] \n")
		sys.stderr.write("where:\n")
		sys.stderr.write("	<histoneFile.bigWig> is the histone signal file (bigWig format)\n")
		sys.stderr.write("	<positives.bed> is the bed file with positives\n")
		sys.stderr.write("	<negatives.bed> is the bed file with negatives\n")
		sys.stderr.write("	<peaks.bed> is the bed file with peaks for master mark\n")
		sys.stderr.write("	<opPrefix> is the output prefix\n")
		sys.stderr.write("	<n> is the number of blocks to divide the positives into (Default: n = 10)\n")
		sys.stderr.write("	<otherMarkFile> is the file with other marks\n")
		sys.exit()

	if len(sys.argv) >= 7:
		n = int(sys.argv[6])
	else:
		n = 10

	pp = PdfPages(sys.argv[5] + "_figures.pdf")
	
	if len(sys.argv) == 6 or len(sys.argv) == 7:
		main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], n, pp)
	elif len(sys.argv) == 8:
		main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], n, pp, otherMarkFile=sys.argv[7])
	
	pp.close()