## Mako
Mako is a bottom-up guided model-free CSV detection tool. It first builds a mutational signal graph and utilizes pattern growth to detect maximal subgraphs as CSVs.

## Install and run

Mako requires Java JDK (>=1.8), we provide a prebuilt JAR package **Mako.jar** for directly usage.

**Note:** All results from the paper is made by **Mako.jar**, which can be found under /out/artifacts/Mako_jar/.

#### Dependency

- htsjdk (https://github.com/samtools/htsjdk): A Java API for processing high-throughput sequencing (HTS) data.
- Python (V>=3.6): This is required for creating Mako configuration file. 
  - Required package: pysam, pandas, numpy

#### Usage

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
minAf=    given the threshold for nodes to include in the signal graph (default=0.2)
minWeight=    given the threshold for nodes to include in the signal graph (default=10)
cutStd=  given the cutoff to determine abnormal insert size read-pairs (default=3)
maxD=    given the maximum distance to cluster abnormal read-pairs nodes (default=readLen)
minFreq=    given the minimum frequency to report discovered subgraphs (default=10)
minQ=    given the minimum mapping quality of read (default=20)
pMax=    the maximal range for adjacent edge linked node search (default=4)
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

#### Run demo data

```
# Create configuration file
python /path/to/process.py config -b NA19240.30X.chr20.1000K-2000K.bam -N 30000 -w ./working_dir/ -s NA19240

# Run Mako
java -jar /path/to/Mako.jar fa=/path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa bamCfg=/path/to/NA19240.mako.cfg minFreq=1 chrom=chr20
```

#### Output file format

**sampleName.superitems.txt:** Mako created nodes for mutational signal graph. Each record in the file contains 19 columns of informations. See table below for details.

![File_headlings](https://github.com/jiadong324/Mako/blob/master/supports/File_headlings.png)

**sampleName.mako.sites.txt:** Mako detect SVs. Additional information of each record can be found in the file heading.

## Contact
If you have any questions, please feel free to contact: jiadong324@gmail.com

## License
The manuscript is under review, license details will be specified later.
