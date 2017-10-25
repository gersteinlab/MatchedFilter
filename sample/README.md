The directory contains sample datasets for testing the tool. 

The narrow peak files are under peakFiles directory in bed format. The signal files are log-scaled bigwig files under signalFile directory. All peak files and signal files are downloaded from ENCODE portal.

Our script generates predictions based on these input files, which are summarized in Input_peakFile and Input_signalFile.
To run the script, use:

    python scanMatchedFilter.py Input_signalFile metaprofileFile sample.chrom.sizes Input_peakFile S2_H3K27ac_positivesAll_MFscores.bed S2_H3K27ac_negatives_MFscores.bed samplePrefix

The sample output file is provided for reference
