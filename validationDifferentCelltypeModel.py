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
import crossValidation


def createSets(positiveFile, negativeFile, opPrefix, n):
	positives = pbt.BedTool(positiveFile)
	negatives = pbt.BedTool(negativeFile)
	numNeg = len(negatives)
	n = 1
	for idx in range(0, n):
		op = open(opPrefix + "_positives" + str(idx) + ".bed", "w")
		numPos = 0
		for currFeature in positives:
			
			op.write(currFeature.chrom + "\t" + str(currFeature.start) + "\t" + str(currFeature.end) + "\n")
			numPos += 1
		op.close()

		op = open(opPrefix + "_negatives" + str(idx) + ".bed", "w")
		for currFeature in negatives:
			op.write(currFeature.chrom + "\t" + str(currFeature.start) + "\t" + str(currFeature.end) + "\n")
		op.close()

	return

def scorePositivesAndNegatives(opPrefix, idx, profileFiles, allMarkSignal):	
	dirName = os.path.dirname(os.path.realpath(__file__))
	positives = opPrefix + "_positives" + str(idx) + ".bed"
	negatives = opPrefix + "_negatives" + str(idx) + ".bed"
	
	op = open(opPrefix + "_marks", "w")
	markIdx = 0
	for currMark in profileFiles:
		if markIdx == 0:
			masterProfile = profileFiles[currMark]
			masterSignal = allMarkSignal[currMark]
			markIdx +=1 
			continue
		metaProfileFile = profileFiles[currMark]
		signalFile = allMarkSignal[currMark]
		op.write(currMark + "\t" + signalFile + "\t" + metaProfileFile + "\n")
	op.close()
	print("python2.7 " + dirName + "/scorePredictions.py " + masterSignal + " " + positives + " " + masterProfile + " " + opPrefix + "_testPos" + str(idx) + " 0 1500 " + opPrefix + "_marks")
	os.system("python2.7 " + dirName + "/scorePredictions.py " + masterSignal + " " + positives + " " + masterProfile + " " + opPrefix + "_testPos" + str(idx) + " 0 1500 " + opPrefix + "_marks")
	os.system("python2.7 " + dirName + "/scorePredictions.py " + masterSignal + " " + negatives + " " + masterProfile + " " + opPrefix + "_testNeg" + str(idx) + " 0 1500 " + opPrefix + "_marks")
	return

def scoreTrainingSamples(trainingPosFile, trainingNegFile, opPrefix, profileFiles, trainingSignals):
	dirName = os.path.dirname(os.path.realpath(__file__))
	markIdx = 0
	op = open(opPrefix + "_marks", "w")
	for currMark in profileFiles:
		if markIdx == 0:
			masterProfile = profileFiles[currMark]
			masterSignal = trainingSignals[currMark]
			markIdx +=1 
			continue
		metaProfileFile = profileFiles[currMark]
		signalFile = trainingSignals[currMark]
		op.write(currMark + "\t" + signalFile + "\t" + metaProfileFile + "\n")
	op.close()

	print("python2.7 " + dirName + "/scorePredictions.py " + masterSignal + " " + trainingPosFile + " " + masterProfile + " " + opPrefix + "_trainPos" + " 0 1500 "  + opPrefix + "_marks")
	print("python2.7 " + dirName + "/scorePredictions.py " + masterSignal + " " + trainingNegFile + " " + masterProfile + " " + opPrefix + "_trainNeg" + " 0 0 "  + opPrefix + "_marks")

	os.system("python2.7 " + dirName + "/scorePredictions.py " + masterSignal + " " + trainingPosFile + " " + masterProfile + " " + opPrefix + "_trainPos" + " 0 1500 "  + opPrefix + "_marks")
	os.system("python2.7 " + dirName + "/scorePredictions.py " + masterSignal + " " + trainingNegFile + " " + masterProfile + " " + opPrefix + "_trainNeg" + " 0 0 "  + opPrefix + "_marks")
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
		os.system("intersectBed -a annotation.bed -b " + peakFile + " -wb -wa -f 0.05 | sortBed -i - > " + opPrefix + "_intersect.bed")
		os.system("intersectBed -a annotation.bed -b " + peakFile + " -wa -v -f 0.05 | sortBed -i -  > " + opPrefix + "_notIntersect.bed")
	else:
		os.system("intersectBed -a annotation.bed -b " + peakFile + " -wb -wa -f 0.25 | sortBed -i - > " + opPrefix + "_intersect.bed")
		os.system("intersectBed -a annotation.bed -b " + peakFile + " -wa -v -f 0.25 | sortBed -i -  > " + opPrefix + "_notIntersect.bed")
	intersectingPeaks = findIntersectingRegions(opPrefix + "_intersect.bed", allPeaks, filetype="bed")
	donotIntersectPos, donotIntersectNeg = statisticsIntersectingRegions(opPrefix + "_notIntersect.bed")
	rankingPeaks = [currInterval.rank for currInterval in intersectingPeaks]
	result = [currInterval.result for currInterval in intersectingPeaks]
	nonintersectingResult = donotIntersectPos * [1] + donotIntersectNeg * [0]
	nonintersectingRank = len(nonintersectingResult) * [100000000]
	fpr, tpr, auc, precision, recall, average_precision = calculateMetrics(1.0/np.array(rankingPeaks + nonintersectingRank), result + nonintersectingResult, "Peaks")
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

		regr = linear_model.LinearRegression(fit_intercept=True,  copy_X=True, normalize=False, n_jobs=1)
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
	for currMark in trainingScores:
		pass
	for idx in range(0, len(trainingScores[currMark])):
		X.append([])

	for currMark in trainingScores:
		for idx in range(0, len(trainingScores[currMark])):
			X[idx].append(trainingScores[currMark][idx])

	X_test = []
	for idx in range(0, len(testScores[currMark])):
		X_test.append([])

	for currMark in trainingScores:
		for idx in range(0, len(testScores[currMark])):
			X_test[idx].append(testScores[currMark][idx])

	forest = RandomForestClassifier(n_estimators=1000, criterion="gini")
	forest = forest.fit(X, np.array(trainingResults))
	y_pred = forest.predict_proba(X_test)[:, 1]
	#print forest.feature_importances_
	return y_pred, forest.feature_importances_

def performSVM(trainingScores, trainingResults, testScores):

	X = []
	for currMark in trainingScores:
		pass
	for idx in range(0, len(trainingScores[currMark])):
		X.append([])

	for currMark in trainingScores:
		print currMark, 
		for idx in range(0, len(trainingScores[currMark])):
			X[idx].append(trainingScores[currMark][idx])

	X_test = []
	for idx in range(0, len(testScores[currMark])):
		X_test.append([])

	for currMark in trainingScores:
		for idx in range(0, len(testScores[currMark])):
			X_test[idx].append(testScores[currMark][idx])

	clf = svm.SVC(kernel='linear', probability=True)
	clf.fit(X, np.array(trainingResults))
	y_pred = clf.predict_proba(X_test)[:, 1]
	print clf.coef_

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
	plt.hist(renPos, bins=50, weights=np.array([1.0/len(renPos)] * len(renPos)))
	plt.xlabel("score")
	plt.hist(renNeg,bins=50, weights=np.array([1.0/len(renNeg)] * len(renNeg)))
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

def readScores(opPrefix, n, trainingPosFile, trainingNegFile, allMarkPeaks, allMarkSignal):
	testScores = []
	testResults = []
	trainingScores = []
	trainingResults = []
	totalAuc, totalAupr = OrderedDict(), OrderedDict()
	fpr, tpr = OrderedDict(), OrderedDict()
	prec, rec = OrderedDict(), OrderedDict()
	colorIdx = OrderedDict()
	idx = 0
	colorIdx["RandomForest"] = 0
	totalAuc["RandomForest"], fpr["RandomForest"], tpr["RandomForest"] = 0, [], []
	totalAupr["RandomForest"], prec["RandomForest"], rec["RandomForest"] = 0, [], []
	colorIdx["SVM"] = 1
	totalAuc["SVM"], fpr["SVM"], tpr["SVM"] = 0, [], []
	totalAupr["SVM"], prec["SVM"], rec["SVM"] = 0, [], []
	colorIdx["LR"] = 2
	totalAuc["LR"], fpr["LR"], tpr["LR"] = 0, [], []
	totalAupr["LR"], prec["LR"], rec["LR"] = 0, [], []
	colorIdx["NB"] = 3
	totalAuc["NB"], fpr["NB"], tpr["NB"] = 0, [], []
	totalAupr["NB"], prec["NB"], rec["NB"] = 0, [], []
	featureNames = []
	for currMark in allMarkSignal:
		featureNames.append(currMark)
		totalAuc[currMark], fpr[currMark], tpr[currMark] = 0, [], []
		totalAupr[currMark], prec[currMark], rec[currMark] = 0, [], []
		colorIdx[currMark] = idx
		totalAuc[currMark + "Peaks"], fpr[currMark + "Peaks"], tpr[currMark + "Peaks"] = 0, [], []
		totalAupr[currMark + "Peaks"], prec[currMark + "Peaks"], rec[currMark + "Peaks"] = 0, [], []
		colorIdx[currMark + "Peaks"] = idx
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
	current_palette = sns.color_palette(n_colors=len(allMarkSignal))
	sns.palplot(current_palette)
	for i in range(0, len(allMarkSignal)):
		usedColors.append(sorted_colors[int(i*float(n1)/(len(allMarkSignal)))])
	c = current_palette
	
	for idx in range(0, n):
		testScores.append(OrderedDict())
		trainingScores.append(OrderedDict())
		trainingResults.append([])
		testResults.append([])
		markIdx = 0
		for currMark in allMarkSignal:
			if markIdx == 0:
				markIdx += 1
				markName = "Master"
			else:
				markIdx += 1
				markName = currMark
			testScores1 = scoresCurrFile(opPrefix + "_testPos" + str(idx) + "_MFscores.bed", markName)
			testScores2 = scoresCurrFile(opPrefix + "_testNeg" + str(idx) + "_MFscores.bed", markName)
			trainScores1 = scoresCurrFile(opPrefix + "_trainPos_MFscores.bed", markName)
			trainScores2 = scoresCurrFile(opPrefix + "_trainNeg_MFscores.bed", markName)
			renTrainScores1, renTrainScores2 = calculateZscores(trainScores1, trainScores2, trainScores2, pp, "train_" + markName)
			renTestScores1, renTestScores2 = calculateZscores(testScores1, testScores2, testScores2, pp, "test_" + markName)
			#renormalizeScores(trainScores1, trainScores2, currMark + "Train")
			#renormalizeScores(testScores1, testScores2, currMark + "Test")
			testScores[idx][currMark] = renTestScores1 + renTestScores2
			trainingScores[idx][currMark] = renTrainScores1 + renTrainScores2
		testResults1 = [1] * len(testScores1)
		testResults2 = [0] * len(testScores2)
		trainResults1 = [1] * len(trainScores1)
		trainResults2 = [0] * len(trainScores2)
		testResults[idx] = testResults1 + testResults2
		trainingResults[idx] = trainResults1 + trainResults2

		for currMark in testScores[idx]:
			currFpr, currTpr, currRoc_auc, currPrec, currRec, curr_aupr = calculateMetrics(testScores[idx][currMark], testResults[idx], currMark)
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

	# for idx in range(0, n):
	#  	for currMark in allMarkPeaks:
	#  		currFpr, currTpr, currRoc_auc, currPrec, currRec, curr_aupr = calculatePeakAUC(opPrefix + "_testPos" + str(idx) + "_MFscores.bed", opPrefix + "_testNeg" + str(idx) + "_MFscores.bed", allMarkPeaks[currMark], opPrefix, currMark)
	#  		currFpr = list(currFpr)
	#  		currTpr = list(currTpr)
	#  		currPrec = list(currPrec)
	#  		currRec = list(currRec)
	#  		if idx == 0:
	#  			for element in currFpr:
	#  				fpr[currMark].append(element)
	#  			for element in currTpr:
	#  				tpr[currMark].append(element)
	#  			for element in currPrec:
	#  				prec[currMark].append(element)
	#  			for element in currRec:
	#  				rec[currMark].append(element)
	#  		totalAuc[currMark + "Peaks"] += currRoc_auc
	#  		totalAupr[currMark + "Peaks"] += curr_aupr

	scores_pred = []
	featureImportance = []
	SVMcoeff = []
	LRcoeff = []
	for idx in range(0, n):
		scores_pred.append(OrderedDict())
		scores_pred[idx]["RandomForest"], currFeatureImportance = performRandomForest(trainingScores[idx], trainingResults[idx], testScores[idx])
		featureImportance.append(currFeatureImportance)
		scores_pred[idx]["SVM"], currFeatureImportance = performSVM(trainingScores[idx], trainingResults[idx], testScores[idx])
		SVMcoeff.append(currFeatureImportance)
		scores_pred[idx]["LR"], currFeatureImportance = performLR(trainingScores[idx], trainingResults[idx], testScores[idx])
		LRcoeff.append(currFeatureImportance)
		scores_pred[idx]["NB"] = performNB(trainingScores[idx], trainingResults[idx], testScores[idx])
			
		#print len(testScores[idx]), np.size(scores_pred[idx]["SVM"]), np.size(testResults[idx])
		for currMark in scores_pred[idx]:
			currFpr, currTpr, currRoc_auc, currPrec, currRec, curr_aupr = calculateMetrics(scores_pred[idx][currMark], testResults[idx], currMark)
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
			print currRoc_auc, curr_aupr
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
		plt.figure()
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
		if "Asym" not in currMark and "Peak" not in currMark and "+" not in currMark and currMark not in ["RandomForest", "SVM", "LR", "NB"]:
			print currMark
			ax.plot(fpr[currMark], tpr[currMark], linewidth=2, c=c[colorIdx[currMark]], label=currMark+' (area = %0.2f)' % totalAuc[currMark])
		elif "Asym" not in currMark and "+" not in currMark and currMark not in ["RandomForest", "SVM", "LR", "NB"]:
			pass
			#ax.plot(fpr[currMark], tpr[currMark], "-.", linewidth=2, c=c[colorIdx[currMark]], label=currMark+' (area = %0.2f)' % totalAuc[currMark])
		elif currMark in ["RandomForest", "SVM", "LR", "NB"]:# or "Random" in currMark:
			ax.plot(fpr[currMark], tpr[currMark], ":", linewidth=1, c=c[colorIdx[currMark]], label=currMark+' (area = %0.2f)' % totalAuc[currMark])
		else:
			pass
			#ax.plot(fpr[currMark], tpr[currMark], "--", linewidth=1, c=c[colorIdx[currMark]], label=currMark+' (area = %0.2f)' % totalAuc[currMark])
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.75, box.height * 0.75])
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=7)
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.title("ROC Plot")
	#plt.legend(loc="lower right")
	plt.savefig(pp, format="pdf")

	plt.figure()
	ax = plt.subplot(111)
	for currMark in totalAupr:
		totalAupr[currMark] /= n
		if "Asym" not in currMark and "Peak" not in currMark and "+" not in currMark and currMark not in ["RandomForest", "SVM", "LR", "NB"]:
			ax.plot(rec[currMark], prec[currMark], linewidth=2, c=c[colorIdx[currMark]], label=currMark+' (area = %0.2f)' % totalAupr[currMark])
		elif "Asym" not in currMark and "+" not in currMark and currMark not in ["RandomForest", "SVM", "LR", "NB"]:
			pass
			#ax.plot(rec[currMark], prec[currMark], "-.", linewidth=1, c=c[colorIdx[currMark]], label=currMark+' (area = %0.2f)' % totalAupr[currMark])
		elif currMark in ["RandomForest", "SVM", "LR", "NB"]: #or "Random" in currMark:
			#pass
			print currMark
			ax.plot(rec[currMark], prec[currMark], ":", linewidth=2, c=c[colorIdx[currMark]], label=currMark+' (area = %0.2f)' % totalAupr[currMark])
		else:
			pass
			#ax.plot(rec[currMark], prec[currMark], "--", linewidth=1, c=c[colorIdx[currMark]], label=currMark+' (area = %0.2f)' % totalAupr[currMark])

		#plt.plot(rec[currMark], prec[currMark], linewidth=2, label=currMark+' (area = %0.2f)' % totalAupr[currMark])
	plt.xlabel('Recall')
	plt.ylabel('Precision')
	plt.title("PR Plot")
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.75, box.height * 0.75])
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=7)
	#plt.legend(loc="upper right")
	plt.savefig(pp, format="pdf")
	#sys.exit()
	#print "Final AUC (symmetric) = ", totalAuc["Master"]
	#print "Final AUC (asymmetric) = ", totalAuc["MasterAsym"]

	return


def main(positiveFile, negativeFile, trainingPosFile, trainingNegFile, opPrefix, n, pp, profileFile, allMarkFile):

	print "-> Create Sets"
	createSets(positiveFile, negativeFile, opPrefix, n)
	n = 1
	print "<- Create Sets"
	
	allMarkSignal = OrderedDict()
	allMarkPeaks = OrderedDict()
	ip = open(allMarkFile, "r")
	for line in ip:
		fields = line.strip().split("\t")
		#print fields[1], fields[2]
		allMarkSignal[fields[0]] = fields[1]
		allMarkPeaks[fields[0]] = fields[2]
	ip.close()

	ip = open(profileFile, "r")
	profileFiles = OrderedDict()
	trainingSignals = OrderedDict()
	for line in ip:
		fields = line.strip().split("\t")
		profileFiles[fields[0]] = fields[1]
		trainingSignals[fields[0]] = fields[2]
	ip.close()

	for idx in range(0, n):
	 	scorePositivesAndNegatives(opPrefix, idx, profileFiles, allMarkSignal)
	scoreTrainingSamples(trainingPosFile, trainingNegFile, opPrefix, profileFiles, trainingSignals)

	readScores(opPrefix, n, trainingPosFile, trainingNegFile, allMarkPeaks, allMarkSignal)

	return


if __name__ == "__main__":
	if not( len(sys.argv) == 8):
		sys.stderr.write("Usage: " + sys.argv[0] + " <positives.bed> <negatives.bed> <trainingPos.bed> <trainingNeg.bed> <opPrefix> <trainingProfiles> <allMarks>\n")
		sys.stderr.write("where:\n")
		sys.stderr.write("	<positives.bed> is the bed file with positives\n")
		sys.stderr.write("	<negatives.bed> is the bed file with negatives\n")
		sys.stderr.write("	<trainingPos.bed> is the bed file with training positives\n")
		sys.stderr.write("	<trainingNeg.bed> is the bed file with training negatives\n")
		sys.stderr.write("	<opPrefix> is the output prefix\n")
		sys.stderr.write("	<profileFile> is the file with training profiles\n")
		sys.stderr.write("	<allMarks> is the file with all marks\n")
		sys.exit()


	pp = PdfPages(sys.argv[5] + "_figures.pdf")
	
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], n, pp, sys.argv[6], sys.argv[7])

	pp.close()