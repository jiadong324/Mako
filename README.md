

<img src="https://github.com/jiadong324/Mako/blob/master/supports/mako_logo.png" alt="mako_logo" style="zoom:0.01%;" align=center/>

## Mako

Mako is a bottom-up guided model-free CSV detection tool. It first builds a mutational signal graph and utilizes pattern growth to detect maximal subgraphs as CSVs.

## Install and run

Mako requires Java JDK (>=1.8), we provide a prebuilt JAR package **Mako.jar** for directly usage.

**Note:** All results from the paper is made by **Mako.jar**, which can be found under /out/artifacts/Mako_jar/.

### Dependency

- htsjdk (https://github.com/samtools/htsjdk): A Java API for processing high-throughput sequencing (HTS) data.
- Python (V>=3.6): This is required for creating Mako configuration file. 
  - Required package: pysam, pandas, numpy

### Usage

```
git clone https://github.com/jiadong324/Mako.git
```

##### Step 1: create Mako configuration file

```
python process.py config -b sample.bam -n 30000 -w ./working_dir/ -s sampleName
```

**Note:** The BAM file  **must** under your working_dir.  The working_dir is default output directory of Mako results.

The configuration file is *sampleName.mako.cfg*. Meanwhile, a insert size distribution figure *sampleName.config.png* is created.

Example of Mako configuration file for HG00514
```
readlen:126
mean:566
stdev:160
bam:/path/to/sample.bam
workDir:/path/to/working_dir/
name:HG00514
```

##### Step 2: run detection

```
# Get help info of Mako
java -jar /path/to/Mako.jar

============  Mako info  ==================

Version: V1.0
Author: Jiadong Lin (Xi'an JiaoTong University)
Contact: jiadong324@gmail.com
Usage:   java -jar /path/to/Mako.jar fa=/path/to/ref.fa bamCfg=/path/to/bam.cfg

================================================

Required inputs:
fa=    Reference genome
bamCfg=  BAM configuration file

Options:
minAf=    minimum allele frequency of a node (default=0.2)
minWeight=    minimum weight of a node (default=10)
cutStd=  cutoff to determine abnormal insert size read-pairs (default=3)
maxD=    maximum distance to cluster abnormal read-pairs nodes (default=readLen)
minFreq=    minimum frequency to report discovered subgraphs (default=10)
minQ=    minimum mapping quality of read (default=20)
pMax=    maximal range for adjacent edge linked node search (default=4)
chrom=   given a specific region, for whole genome if it is not given. (e.g chr1:1000-2000)
```

Run SV discovery with Mako configuration file.
```
java -jar Mako.jar fa=/path/to/ref.fa bamCfg=/path/to/sampleName.mako.cfg
```

##### Step 3: filter with CXS score (optional)

```
# Usage info
python process.py filter -h 

Usage: process.py [options]

Options:
  -h, --help  show this help message and exit
  -i INPUT    Input Mako callset
  -o OUT      Output of filtered Mako callset by CXS
  -c CXS      CXS threshold
  -f FORMAT   output format (original, bed)

# Filter with cxs>4, and output in BED format 
python process.py filter -i /path/to/sampleName.mako.sites.txt -o /path/to/filtered_output.bed -c 4 -f bed
```

### Run demo data

Note that the demo data is small, large minFreq is not recommended.
```
# Create configuration file
python /path/to/process.py config -b NA19240.30X.chr20.1000K-2000K.bam -N 30000 -w ./working_dir/ -s NA19240

# Run Mako
java -jar /path/to/Mako.jar fa=/path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa bamCfg=/path/to/NA19240.mako.cfg minFreq=1 chrom=chr20
```

### Output file

**sampleName.superitems.txt:** Mako created nodes for mutational signal graph. Each record in the file contains 19 columns of information. Detailed explanation can be found at https://github.com/jiadong324/Mako/tree/master/supports/File_heading.png 

**sampleName.mako.sites.txt:** Mako detect SVs. Additional information of each record can be found in the file heading.

**sampleName.mako.BNDs.txt:** Break-ends are not able to involve in the current mutational signal graph, including inter-chromosomes.

### Parameter settings

**minAf:** calculates the ratio between abnormal reads and all reads, estimating the allele frequency of signal node.

**pMax:** constrains the pattern growth through adjacent edges. The actual distance on genome is calculated as *pMax* * *fragMean*. Typically, SV size is enriched at 300bp and 1000bp. Thus the distance is recommended to be larger than 1000bp. Mako uses *pMax=4* as default value if the estimated average fragment size is 500-600bp.

**minFreq:** for whole genome detection, we would expect one subgraph appears at least on the half number of genomes. If you run only one chromosome, this value is recommended to be 1.

**maxD:** this is used cluster signal nodes produced by abnormal paired-end alignment. And split or clipped is restricted to the same position.

### Known issues

1. Please make sure the reference used for running Mako is identical to the alignment one.
2. ...

## Contact

If you have any questions, please feel free to contact: jiadong324@gmail.com

## License
The manuscript is under review, license details will be specified later.

