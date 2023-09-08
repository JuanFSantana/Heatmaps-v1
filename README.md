# heatmap.py #
This script will create heatmaps from bigWig files.K-means clustering can also be applied for in-depth data analysis. 

# File requirements #
The input regions file should be a six-column, tab-delimited bed file containing chromosome, start and end positions, and the strand information for each region. The regions can be of any length as long as it is an even number and the center is a feature under study (e.g. transcription start site). 
 
| chr6 | 142946246 | 142946446 | Gene_A | 255 | - |
|:----:|:---------:|:---------:|:------:|:---:|:-:|

BigWig files of data to be abalayzed. If stranded data, separate forward and reverse bigWigs are needed.


# Behavior #
It can make heatmaps for stranded data (e.g. PRO-Seq) or non-stranded data (e.g. ChIP). It requires two bigWig files (e.g. control and experimental) creating heatmaps for each of them and a log2 fold change color heatmap. 
The program is capable of applying a k-means clustering algorithm on a chosen interval from the log2 fold change heatmaps.

# Dependencies #
### Python libraries ###
Pandas: https://pypi.org/project/pandas/

Numpy: https://pypi.org/project/numpy/

Matplotlib: https://matplotlib.org/stable/users/installing/index.html

pyBigWig: https://github.com/deeptools/pyBigWig

scikit_learn: https://scikit-learn.org/stable/

# Example command usag #
```
python heatmap.py maxtss-plus-minus-500bp.bed \
                  - numerator POLII-CHIP-DEPLETION.bw \
                  - denominator POLII-CHIP-DMSO.bw \
                  -chip \
                  -x 3 \
                  -id Depleted DMSO \
                  -o /home/usr/

python heatmap.py maxtss-plus-minus-500bp.bed \
                  - numerator POLII-PRO-SEQ-DEPLETION-forward.bw POLII-PRO-SEQ-DEPLETION-reverse.bw\
                  - denominator POLII-PRO-SEQ-DMSO-forward.bw POLII-PRO-SEQ-DMSO-reverse.bw\
                  -x 3 \
                  -id Depleted DMSO \
                  -o /home/usr/
		  -k 10
 	          -r 400 800
```
# Parameter description #
```
regions: <str> Bed file of genomic regions of chosen length with the format described above

-numerator: <str> If stranded data, forward and reverse (in that order) bigwigs for condition to be used as numerator in log2FC. If non-stranded data, only one file is required.

-denominator: <str> If stranded data, forward and reverse (in that order) bigwigs for condition to be used as denominator in log2FC. If non-stranded data, only one file is required.

-mb: <int> Sets the chosen value as black, default is the largest number in the matrix

-mc: <int> Sets the chosen value as most red/blue in the log2 fold change heatmap, default is the largest absolute fold change

-n: If the argument is invoked, the average total number of reads for all regions will be calculated between the numerator and denominator datasets and the reads per base will be normalized to this value 

-chip: If bigwigs are ChIP-seq data, argument must be invoked

-k: <int> Number of K-means clusters to be used for clustering the data

-r: <int> <int> Location in the regions file to be used for k-means clustering. The first argument is the start position and the second is the end position. For example, in regions of length 1000 bp, if you want to cluster the middle 500 bp, the arguments would be -r 250 750. Clustering will be done based on the data from log2 fold change.

-antiSense: If the argument is invoked, the counts from the anti-sense strand will be plotted instead of the sense strand. This is only applicable for transcriptional data.

-y: <int> (value greater than or equal to 1) Horizontal lines/bp for each fragment length, default is 1

-x: <float> or <int> (value less than or equal to 1) Vertical lines/bp for each genomic interval displayed, for example, -x 1 is one vertical line/bp; -x 0.1 is one vertical line/averaged 10 bp, default is 1

-g: <float> Gamma correction factor, default is 1 but 0.5 provides an image which is more interpretable by the human eye. For more information: https://en.wikipedia.org/wiki/Gamma_correction

-o: <str> Ouput directory

-id: <str> Image output names for numerator and denominator (in that order). The log2FC heatmaps will be named by combining the numerator and denominator file names
```

###### Example ouput applying k-means clustering analysis (k=10) with command: -numerator Control-H3K4me3-DFF-ChIP-centers.bw -denominator Treatment-H3K4me3-DFF-ChIP-centers.bw -k 10 -r 400 900 -x 6 

&nbsp; Control-H3K4me3 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Treatment-H3K4me3 &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; log2FC
<p float="left">
    <img src="https://github.com/JuanFSantana/DNA-and-RNA-seq-analysis-essentials/blob/main/Heatmaps/images/Cluster-10-chip-H3K4me3-Max-10.0-X-3.0-Y-1.0.png?raw=true" width="250" height="400" />
    <img src="https://github.com/JuanFSantana/DNA-and-RNA-seq-analysis-essentials/blob/main/Heatmaps/images/Cluster-10-chip-2h-Max-10.0-X-3.0-Y-1.0.png?raw=true" width="250" height="400" />
    <img src="https://github.com/JuanFSantana/DNA-and-RNA-seq-analysis-essentials/blob/main/Heatmaps/images/Cluster-10-chip-H3K4me3-vs-2h-Max-1.67-X-3.0-Y-1.0.png?raw=true" width="250" height="400" />
</p>
