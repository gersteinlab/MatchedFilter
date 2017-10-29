import sys
import os

formatWig  = ""
wigChr = ""
wigPosn = -1
wigSpan = -1
wigStep = -1

class Interval():
	def __init__(self):
		self.chr = None
		self.start = -1
		self.end = -1
		self.signal = 0.0
		return

	def extractFeatureBed(self, line):
		if line["chrom"] == None:
			currInterval.chr = None
			return
		self.chr = line["chrom"]
		self.start = int(line["chromStart"])
		self.end = int(line["chromEnd"])
		#print line
		if line["name"] != None:
			self.name = line["name"]
		if line["score"] != None:
			self.signal = float(line["score"])
		return

	def extractFeatureBedGraph(self, line):
		if line["chrom"] == None:
			currInterval.chr = None
			return
		self.chr = line["chrom"]
		self.start = int(line["chromStart"])
		self.end = int(line["chromEnd"])
		if line["signal"] != None:
			self.signal = float(line["signal"])
		return

	def extractFeatureWig(self, ipWig):
		global formatWig, wigChr, wigPosn, wigStep, wigSpan
		
		line = ipWig.readline().rstrip()
		if not line:
			self.chr = None
			return
		if "fixedStep" in line:
			formatWig = "fixedStep"
			wigChr = line.split()[1].split("chrom=")[1]
			wigPosn = int(line.split()[2].split("start=")[1]) - 1 #Consistent with bed format
			wigStep = int(line.split()[3].split("step=")[1])
			if "span" in line:
				wigSpan = int(line.split()[4].split("span=")[1])
			else:
				wigSpan = 1
			line = ipWig.readline().rstrip()
		elif "variableStep" in line:
			formatWig = "variableStep"
			wigChr = line.split()[1].split("chrom=")[1]
			if "span" in line:
				wigSpan = int(line.split()[2].split("span=")[1])
			else:
				wigSpan = 1
			line = ipWig.readline().rstrip()

		if formatWig == "fixedStep":
			self.chr = wigChr
			self.start = wigPosn
			self.end = self.start + wigSpan
			wigPosn += wigStep
			self.signal = float(line)
		elif formatWig == "variableStep":
			self.chr = wigChr
			fields = line.split("\t")
			self.start = int(fields[0]) - 1
			self.signal = float(fields[1])
			self.end = self.start + wigSpan

		return

	def copyInterval(self, oldInterval):
		self.chr = oldInterval.chr
		self.start = oldInterval.start
		self.end = oldInterval.end
		self.signal = oldInterval.signal
		return

	def areTheyContinuous(self, oldInterval, midWidth):
		if self.chr == oldInterval.chr and self.start <= oldInterval.end + minWidth:
			return True
		else:
			return False

def makeIntervalsSame(currInterval, newInterval):
	currInterval.chr = newInterval.chr
	currInterval.start = newInterval.start
	currInterval.end = newInterval.end
	currInterval.signal = newInterval.signal
	return