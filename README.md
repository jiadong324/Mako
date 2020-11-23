

<img src="https://github.com/jiadong324/Mako/blob/master/supports/mako_logo.png" alt="mako_logo" style="zoom:0.01%;" align=center/>

## Mako

Mako is a bottom-up guided model-free CSV detection tool. It first builds a mutational signal graph and utilizes pattern growth to detect maximal subgraphs as CSVs.

Please check the [wiki](https://github.com/jiadong324/Mako/wiki) page for more details.

## Install and run

Mako requires Java JDK (>=1.8), we provide a prebuilt JAR package **Mako.jar** for directly usage. See [Release](https://github.com/jiadong324/Mako/releases).

**Note:** All results from the paper is made by **Mako.jar**, which is a beta version.

### Dependency

- htsjdk (https://github.com/samtools/htsjdk): A Java API for processing high-throughput sequencing (HTS) data.
- Python (V>=3.6): This is required for creating Mako configuration file. 
  - Required package: pysam, pandas, numpy

### Usage

```
git clone https://github.com/jiadong324/Mako.git
```

##### Step 1: create Mako configuration file


**Note:** The BAM file  **must** under your working_dir.  The working_dir is default output directory of Mako results.

```
python process.py config -b sample.bam -n 30000 -w ./working_dir/ -s sampleName -f /path/to/ref.fa.fai
```

The configuration file is *sampleName.mako.cfg*.

##### Step 2: run detection

```
# Get help info of Mako

java -jar Mako.jar

```

Run SV discovery with Mako configuration file.
```
java -jar Mako.jar -R /path/to/ref.fa -F /path/to/sampleName.mako.cfg
```

##### Step 3: Convert Mako output to standard VCF format

```
python ParseMako.py tovcf -m sampleName_mako_calls.txt -o sampleName_mako.vcf -r /path/to/ref.fa -s sampleName
```

### Run demo data

Note that the demo data is small, large minFreq is not recommended.
```
# Create configuration file
python /path/to/process.py config -b NA19240.30X.chr20.1000K-2000K.bam -N 30000 -w ./working_dir/ -s NA19240

# Run Mako
java -jar /path/to/Mako.jar fa=/path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa bamCfg=/path/to/NA19240.mako.cfg minFreq=1 chrom=chr20
```

### Known issues

1. Please make sure the reference used for running Mako is identical to the alignment one.
2. ...

## Contact

If you have any questions, please feel free to contact: jiadong324@gmail.com

## License
The manuscript is under review, license details will be specified later.

