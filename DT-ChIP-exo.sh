#!/usr/bin/env bash
#
################################################
# PROGRAM:
# DT-ChIP-exo.sh
#
# DESCRIPTION:
# I wrote a very basic script to generate the plots from the Rossi et al. (2021) ChIP-exo data mapped to 
# annotated TSS and TES positions from Steinmetz et al. (2013). Additionally, the data is plotted over
# two sets of gene quartile lists. One is a proxy of expression, ranked by Rpb3 occupancy around the TSS
# (performed by myself), and the other is a gene list ranked by length. 
#
# AUTHOR:
# Andrew J. Tonsager
#
# START DATE:
# January 7, 2024
#
# DEPENDENCIES:
# 	Requires the installation of the following software: 
#		deeptools
#		bedGraphToBigWig
#	
#
#
# REQUIRES:
#    INPUT: _filtered.bam file:    			Download the entire folder of files from yeastepigenome.org for a factor
#											of your choice. The entire folder should be unzipped. This script will
#											generate a new folder in a parent directory named after the factor of your
#											choice containing /images/ and /data/ folders which will contain the 
#											data and figures produced from the pipeline.
#
#    GENOME SIZES: sacCer3.chrom.sizes: 	File obtainable from UCSC Genome Browser Downloads which contains the sizes
#											of each yeast chromosome, which is required for producing a bigWig from the
#											bed file produced during the pipeline. This should be in "input_files". 
#											For the download, see: 
#											https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/
#
#
#    GENOME ANNOTATION: .gtf or .txt files: This is a genome annotation file of gene features. For this pipeline, I 
#											wrote it to input a unranked annotation and ranked gene lists for expression
#											and length. If desired, you could input other gene lists of interest. This
#											should also be in "input_files".
#
# USAGE:
# 	1)	conda create -n DT-ChIP-exo deeptools ucsc-bedgraphtobigwig (only need to do this once)
# 	2) 	source activate DT-ChIP-exo 
#	3) 	bash DT-ChIP-exo.sh <factor name>
#
# OUTPUT:
#	This pipeline will produce a new folder titled <factor name> containing /images/ and /data/
#
# KNOWN BUGS:
#	Currently I'm having issues with having bedGraphToBigWig in a different folder. As a result, I installed it to my environment.
#	It would be handy to have it just added to your path rather than in the conda environment.
#	Also, my plotProfile doesn't work with the samplelabel... the quotes are used to name the variable in the command line but it
#	no longer will read the variable.
#
# THINGS TO IMPROVE:
#	Currently, I've written this pipeline to pull the .bam file from the current folder. Navigating to /Rossi-2021/<number>_YEP/
#	folder for each analysis might be a bit tedious in the future. Additionally, if I could write the code to be more high throughput
#	to analyze multiple factors at the same time, that might be advantageous. It would also be great to include the option to average
#	replicates if needed. For now, I've been analyzing single replicates.
#
################################################

####### MODIFY THIS SECTION #############

#Factor name. This pulls the name from the command line
name=$1

#This is where the annotation files live:
annotation="../../input_files/"

#This is the name of the full gene list file:
genes="Sc.cerevisiae.annotation_Steinmetz_2013.gtf"

#This is the name of your length quartiles file:
length="Length_quartiles.txt"

#This is the name of your expression quartiles file:
expression="Rpb3_quartiles.txt"

#This is the output_directory:
DATE=`date +%Y-%m-%d`
outputdir="../../"$name"/"$DATE"_output/"

#This is your desired file type for images (options: png, pdf, svg, eps, plotly)
image="png"


########## DONE MODIFYING ###############





########## BEGIN CODE ###############
#
echo -e ">>> INITIATING DT-ChIP-exo.sh with command:\n\t$0 $@"

# Make output directories
echo -e ">>> MAKING output directory..."
echo -e "\tmkdir $outputdir"
mkdir -p $outputdir
mkdir -p ""$outputdir"data/"
mkdir -p ""$outputdir"images/"

# Use samtools to index the .bam file in the folder you're in, then bamCoverage to convert to a .bed file.
echo -e ">>> INDEXING bam and creating .bed file..."
samtools index *.bam
bamCoverage -b *_filtered.bam -o $outputdir"data"/$name"_filtered".bed -of bedgraph
sed -i -e 's/chr1/chrI/g' -e 's/chr2/chrII/g' -e 's/chr3/chrIII/g' -e 's/chr4/chrIV/g' -e 's/chr5/chrV/g' -e 's/chr6/chrVI/g' -e's/chr7/chrVII/g' -e 's/chr8/chrVIII/g' -e 's/chr9/chrIX/g' -e 's/chrI0/chrX/g' -e 's/chrI1/chrXI/g' -e 's/chrI2/chrXII/g' -e 's/chrI3/chrXIII/g' -e 's/chrI4/chrXIV/g' -e 's/chrI5/chrXV/g' -e 's/chrI6/chrXVI/g' $outputdir"data"/$name"_filtered".bed
sort -k1,1 -k2,2n $outputdir"data"/$name"_filtered".bed > $outputdir"data"/$name"_filtered_sorted".bed
sed -i '/2-micron/d' $outputdir"data"/$name"_filtered_sorted".bed

# Use bedGraphToBigWig to produce a bigWig file from the .bed file.
echo -e ">>> CREATING bigWig file..."
bedGraphToBigWig $outputdir"data"/$name"_filtered_sorted".bed $annotation"sacCer3.chrom.sizes" $outputdir"data"/$name.bw

# Use deeptools to produce zipped matrix files containing data over your annotation of choice.
# I currently am using the reference-point option and a window of +/- 500bp around the TSS And TES, but this could be modified as needed.
echo -e ">>> COMPUTEMATRIX to produce matrix files..."
computeMatrix reference-point -S $outputdir"data"/$name.bw -R $annotation"$genes" -a 500 -b 500 -o $outputdir"data"/$name"_TSS".gz --missingDataAsZero
computeMatrix reference-point -S $outputdir"data"/$name.bw -R $annotation"$genes" -a 500 -b 500 -o $outputdir"data"/$name"_TES".gz --missingDataAsZero --referencePoint TES
computeMatrix reference-point -S $outputdir"data"/$name.bw -R $annotation"$length" -a 500 -b 500 -o $outputdir"data"/$name"_TSS_len".gz --missingDataAsZero
computeMatrix reference-point -S $outputdir"data"/$name.bw -R $annotation"$length" -a 500 -b 500 -o $outputdir"data"/$name"_TES_len".gz --missingDataAsZero --referencePoint TES
computeMatrix reference-point -S $outputdir"data"/$name.bw -R $annotation"$expression" -a 500 -b 500 -o $outputdir"data"/$name"_TSS_expr".gz --missingDataAsZero
computeMatrix reference-point -S $outputdir"data"/$name.bw -R $annotation"$expression" -a 500 -b 500 -o $outputdir"data"/$name"_TES_expr".gz --missingDataAsZero --referencePoint TES

# Plot the data from matrices produced by deeptools using plotProfile.
echo -e ">>> PLOTTING data as $image files..."
plotProfile -m $outputdir"data"/$name"_TSS".gz -o $outputdir"images"/$name"_TSS".$image --refPointLabel 'Distance from TSS (bp)' --samplesLabel '$name ChIP-exo' --colors black
plotProfile -m $outputdir"data"/$name"_TES".gz -o $outputdir"images"/$name"_TES".$image --refPointLabel 'Distance from TES (bp)' --samplesLabel '$name ChIP-exo' --colors black
plotProfile -m $outputdir"data"/$name"_TSS_len".gz -o $outputdir"images"/$name"_TSS_len".$image --refPointLabel 'Distance from TSS (bp)' --samplesLabel '$name ChIP-exo'
plotProfile -m $outputdir"data"/$name"_TES_len".gz -o $outputdir"images"/$name"_TES_len".$image --refPointLabel 'Distance from TES (bp)' --samplesLabel '$name ChIP-exo'
plotProfile -m $outputdir"data"/$name"_TSS_expr".gz -o $outputdir"images"/$name"_TSS_expr".$image --refPointLabel 'Distance from TSS (bp)' --samplesLabel '$name ChIP-exo'
plotProfile -m $outputdir"data"/$name"_TES_expr".gz -o $outputdir"images"/$name"_TES_expr".$image --refPointLabel 'Distance from TES (bp)' --samplesLabel '$name ChIP-exo'


######## VERSIONS #############
echo -e "\n>>> VERSIONS:"
echo -e "\n>>> SAMTOOLS VERSION:"
samtools --version
echo -e "\n>>> BAMCOVERAGE VERSION:"
bamCoverage --version
echo -e "\n>>> DEEPTOOLS VERSION:"
deeptools --version
echo -e ">>> END: ChIP-exo analysis complete."
