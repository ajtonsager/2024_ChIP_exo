# Workflow for plotting factor occupancy over transcribed regions from published ChIP-exo data
Andrew J. Tonsager
## Introduction
This document contains a bash script and a series of instructions to process BAM sequence alignment files and generate figures over desired genomic regions using published ChIP-exo data published in [Rossi et al., 2021.](https://pubmed.ncbi.nlm.nih.gov/33692541/) The metadata including the BAM files is located at [yeastepigenome.org](http://yeastepigenome.org/), a repository containing data for over 400 yeast factors. These instructions will allow for recaptiulating the figures that I've generated for Chapter 3 of my Dissertation, "RNAPII-associated factors maintain nucleosome features over the yeast genome".
## 1. Setup
### Set up conda environments
I first created a new conda environment `DT-ChIP-exo` with both `deeptools` and `bedGraphtoBigWig` since both are used with this bash script. Conda was installed in Linux using the Anaconda distribution, with instructions that can be found here: [Installing Conda in Linux](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html)
```
conda create -n DT-ChIP-exo
source activate DT-ChIP-exo
conda install bioconda::deeptools
conda install bioconda::ucsc-bedgraphtobigwig
```
### Obtain metadata from ChIP-exo data on yeastepigenome.org
The bash script requires the bam file for a factor of interest to already be downloaded to your working directory. I obtained all data from [yeastepigenome.org](http://yeastepigenome.org/) using the following steps:
1. Search for each factor with the search bar
2. Left click on "META DATA"
3. **Right** click on the direct download arrow, and selecting "Copy link address"
4. In your terminal, use `wget` to download the data to your working directory
An example of the link for a replicate from Spn1 is shown below:
```
wget https://www.datacommons.psu.edu/download/eberly/pughlab/yeast-epigenome-project/12107_YEP.zip
```
5. This data was unzipped using `gunzip *.zip`.
6. I copied the `DT-ChIP-exo.sh` bash script into the unzipped folder.
   
## 2. File organization
As specified in the bash script, this pipeline will produce a new folder titled [factor name] containing folders `/images/` and `/data/`. The factor name is inputted when running the bash script. My folder structure looks like this, but the metadata can be changed as needed if the desired folder structure looks different.
```
── ChIP-exo_data/
  │   └── Rossi-data/                              <- Contains each folder downloaded and unzipped from yeastepigenome.
  |       ├── [SampleID]_YEP/                           <- Folders unzipped are labelled by the SampleID assigned in the Rossi et. al 2021 study.
  |           ├── [SampleID]_AllGenes_ColorBar                <- These folders contain a variety of files not used in the analysis.
  |             ...
  |           ├── [SampleID]_filtered.bam                     <- This folder MUST contain the .bam file which is here by default.
  |             ...
  |           ├── DT-ChIP-exo.sh                              <- Currently the bash script must be run within this folder, so put it here.
  |
  │   └── input_files/                             <- Folder containing annotation .gtf and .txt files.
  │       ├── Length_quartiles.txt                      <- Text file containing TSS annotations ranked into quartiles by gene length.
  |       ├── Pelechano.2013_mRNA.gtf                   <- GTF annotation file containing TSS annotations from Pelechano et al. 2013
  |       ├── Rpb3_quartiles.txt                        <- Text file containing TSS annotations ranked into quartiles by Rpb3 occupancy.
  |
  │   └── [factor name]/                           <- Folder created with the name inputted at the command line
  │       ├── year-month-day_output                     <- Folder created with the date the analysis was run
  │           ├── data/                                      <- Folder created with .bed, .bw and .gz files produced during the analysis
  |           ├── images/                                    <- Folder created with images with file type specified in bash script.
```

## 3. Metadata section to modify
This section should only be modified if the folder structure is different from that listed above and to designate the desired file format for images (I used both .png and .pdf).
```
#Factor name. This pulls the name from the command line
name=$1

#This is where the annotation files live:
annotation="../../input_files/"

#This is the name of the full gene list file:
genes="Pelechano.2013_mRNA.gtf"

#This is the name of your length quartiles file:
length="Length_quartiles.txt"

#This is the name of your expression quartiles file:
expression="Rpb3_quartiles.txt"

#This is the output_directory:
DATE=`date +%Y-%m-%d`
outputdir="../../"$name"/"$DATE"_output/"

#This is your desired file type for images (options: png, pdf, svg, eps, plotly)
image="png"
```

## 4. Usage
Once the bash script is downloaded to the unzipped folder containing the input .bam file for your factor of interest, simply run the following with the corresponding [factor name]:
```
bash DT-ChIP-exo.sh <factor name>
```
Details on how each line of the code works is annotated within the `DT-ChIP-exo.sh` script.
