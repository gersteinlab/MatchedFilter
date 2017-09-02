

There are 4 scripts within this tool:

1) createMetaPattern.py creates the metaprofile given peaks from a massively parallel reporter assay for promoter and enhancer activity. To run this script:

python createMetaPattern.py <histoneFile.bigWig> <MPRApeaks.bed> <nonrandomFile.bed> <opPrefix> [<otherMarks>]

where:
	<histoneFile.bigWig> is the histone signal over the genome for the primary mark (H3K27ac in journal article) in bigWig format (https://genome.ucsc.edu/goldenpath/help/bigWig.html). The signal is in log-foldchange between histone and control ChIP-seq datasets.

	<MPRApeaks.bed> is the peaks from a massively parallel assay in bed format (chromosome, start, and end positions are most important and are necessary).

	<nonrandomFile.bed> is a bed file containing all nonrandom regions of the genome in bed3 format (that is regions in the genome with non-zero ChIP-seq signal and not overlapping with a MPRA peak).

	<opPrefix> is the prefix for all output files.

	<otherMarks> is an optional file that contains any other histone marks. The format of this file is tab delimited 2-columns with experimental dataset name (for example, H3K4me1) and 2nd column containing file name for the experimental dataset (in bigWig format - log foldchange signal enrichment over control).


2) crossValidation.py performs the k-fold cross validation of Matched Filter for different chromatin marks. To run this script:

python crossValidation.py <histoneFile.bigWig> <positives.bed> <negatives.bed> <peaks.bed> <opPrefix> [<n> [<otherFiles>]]

where:

	<histoneFile.bigWig> is the histone signal over the genome for the primary mark in bigWig format (https://genome.ucsc.edu/goldenpath/help/bigWig.html). The signal is in log-foldchange between histone and control ChIP-seq datasets.

	<positives.bed> is the list of positives (for example, MPRA peaks in bed format)

	<negatives.bed> is the list of negatives (chosen from nonrandom regions of the genome)

	<peaks.bed> is the bed file with peaks for primary chromatin mark (histone or DHS).

	<opPrefix> is the prefix for all output files.

	<n> is the value for k in k-fold cross validation. Default value for k=10.

	<otherFiles> contains the information about secondary chromatin marks to be tested with matched filter. The format of this file is tab delimited 3-columns with experimental dataset name (for example, H3K4me1), the 2nd column containing file name for the experimental dataset (in bigWig format - log foldchange signal enrichment over control), and the 3rd column contains the name of the file for the corresponding chromatin peaks (to compare AUC and AUPR for peaks against MF score).

3) validationDifferentCelltypeModel.py performs cross validation of models learned in training cell-type (all training and metaprofiles learned in this cell-type) to predict regulatory elements in a test cell-type

python validationDifferentCelltypeModel.py <positives.bed> <negatives.bed> <trainingPos.bed> <trainingNeg.bed> <opPrefix> <trainingProfiles> <allMarks>

where:

	<positives.bed> labeled positives in test cell-type - for example, from MPRA assay.

	<negatives.bed> labeled negatives in test cell-type.

	<trainingPositives.bed> labeled positives in training cell-type. These files are provided in training data directory.

	<trainingNegatives.bed> labeled negatives in training cell-type. These files are provided in training data directory.

	<opPrefix> is the prefix for all output files.

	<trainingProfiles> is a file with the matched filter learned from training data. This is a tab-delimited 2 column file with the first column containing experimental dataset name (for example, H3K4me1) and the 2nd column containing file name with metaprofile (provided in subdirectory metaprofiles).

	<allMarks> is a file with the histone datasets in bigWig format for the test cell-type. The training profiles will be used to give best matched filter score over positives and negatives in test cell-type using these datasets. The datasets should match training positives and negatives datasets. 

4) scanMatchedFilter.py scans the whole genome with matched filter and applies the SVM model using training data provided with the code to predict enhancers and promoters in a cell-type specific manner.

python scanMatchedFilter.py <fileList> <metaProfileList> <chrNameList> <peakFileList> <opPrefix>

where:

	<fileList> is a file with the list of chromatin signals in the format (2 column tab delimited with experimental dataset name in column 1 and filename in column 2). The chromatin signals are in log foldchange signal signal enrichment over control.

	<metaProfileList> is the list with training profiles. This is a tab-delimited 2 column file with the first column containing experimental dataset name (for example, H3K4me1) and the 2nd column containing file name with metaprofile (provided in subdirectory metaprofiles).

	<chrNameList> is the list with chromosome names. The first column contains chromosome name and 2nd column contains length of chromosome.
	
	<peakFileList> is the file with chromatin peaks. These regions are removed during background model fitting.
	
	<positiveScores> is the file containing the scores for all training positives, provided in the training data directory.
	
	<negativeScores> is the file containing the scores for all training negatives, provided in the training data directory.

	<opPrefix> is the prefix for all output files.

	The final output file test_SVMpredScores.dat contains the SVM scores and predictions.
	
