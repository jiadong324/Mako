/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package option;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.List;
import java.util.ArrayList;

/**
 *
 * @author jiadonglin
 */
public class Params {
    
    public int readLen;
    public int fragMean;
    public int fragStd;
    public int cutStd = 3;
    public int minFreq = 10;
    public int clusteringDist;
    public int minMapQ = 20;
    public int minWeight = 10;
    public int patternMaxSpan;
    public String chrom = null;
    public int givenRegionS;
    public int givenRegionE;
    public String workDir;
    public String bamFile;
    public String fastaFile;
    public String fastaIndexFile;
    public String regionMaskFile;
    public String superitemOut;
    public String bndOutPath;
    public String svOut;
    public String indelOut = null;
    public String abnormalSignalOut = null;
    public String mergedPatternOut = null;
    public String frequentPatternOut = null;
    public boolean hasParamInput = true;
    public double minAf = 0.2;
//    public List<String> excludableRegion;

    public Params(){
        
    }
    
    public void loadParam(String[] args) throws IOException{
        
        if (args.length == 0){
            printOptionsHelp();
            hasParamInput = false;
        }else{
            loadTerminalInputParams(args);
        }
    }

    public String getParams() {
        StringBuilder sb = new StringBuilder();
        sb.append("fragMean: ");
        sb.append(fragMean);
        sb.append(" fragStd: ");
        sb.append(fragStd);
        sb.append(" readLen: ");
        sb.append(readLen);
        sb.append(" clustD: ");
        sb.append(clusteringDist);
        sb.append(" patternMaxSpan: ");
        sb.append(patternMaxSpan);
        sb.append(" minAF: ");
        sb.append(minAf);
        sb.append("\n");
        return sb.toString();
    }

    public String createMakoHeader(String libStats) {
        StringBuilder sb = new StringBuilder();
        sb.append("## Software: Mako V1.0\n");
        sb.append("## Running parameters: ");
        sb.append(libStats);
        sb.append("## <LinkType: Arp_span, connection found between sub-graphs (paired edges)>\n");
        sb.append("## <LinkType: Arp_self, connection found within a sub-graph (paired edges)>\n");
        sb.append("## <LinkType: Split, breakpoint estimated with split alignment>\n");
        sb.append("## <LinkType: Cross, breakpoint estimated with read cross alignment>\n");
        sb.append("## <Pattern: the SV graph represented by node type>\n");
        sb.append("## <PatternInfo: attributes of each node in the graph>\n");
        sb.append("## <sa: Number of split reads evidence>\n");
        sb.append("## <rp: Number of abnormal read pair evidence>\n");
        sb.append("## <cr: Number of cross matched reads>\n");
        sb.append("#CHROM\tSTART\tEND\tLinkType\tLinkInfo\tPattern\tPatternInfo\n");

        return sb.toString();
    }

    private void printOptionsHelp(){
        System.out.println("============  Mako info  ==================");
        System.out.println("\nVersion: V1.0\nAuthor: Jiadong Lin (Xi'an JiaoTong University)\nContact: jiadong324@gmail.com");
        System.out.println("Usage:\t java -jar /path/to/Mako.jar fa=/path/to/ref.fa bamCfg=/path/to/bam.cfg\n");
        System.out.println("================================================\n");
        StringBuilder sb = new StringBuilder();


        sb.append("\nRequired inputs:");
        sb.append("\nfa=    Reference genome\n");
        sb.append("bamCfg=  BAM configuration file\n");

        sb.append("\nOptions:");
        sb.append("\nminAf=    given the threshold for nodes to include in the signal graph (default=0.2)\n");
        sb.append("minWeight=    given the threshold for nodes to include in the signal graph (default=10)\n");
        sb.append("cutStd=  given the cutoff to determine abnormal insert size read-pairs (default=3)\n");
        sb.append("maxD=    given the maximum distance to cluster abnormal read-pairs nodes (default=readLen)\n");
        sb.append("minFreq=    given the minimum frequency to report discovered subgraphs (default=10)\n");
        sb.append("minQ=    given the minimum mapping quality of read (default=20)\n");
        sb.append("pMax=    the maximal range for adjacent edge linked node search (default=4)\n");
        sb.append("chrom=   given a specific region, for whole genome if it is not given. (e.g chr1:1000-2000)\n");


//        sb.append("freOut=  output path of frequent raw patterns (optional, default=False)\n");
//        sb.append("sigOut=  output path of abnormal alignments (optional, default=Flase)\n");
//        sb.append("mergeOut=  output path of merged patterns (optional, default=False)\n");
//        sb.append("indelOut=    output path of INDELS (optional, default=False)\n");
//        sb.append("regionMask=   regions to exclude during SV discovery (optional)\n");
        
        System.out.println(sb.toString());
    }

    
    private void loadTerminalInputParams(String[] args) throws IOException{
        for (int i = 0; i < args.length; i++){

            String[] argTokens = args[i].split("=");
            if (argTokens[0].equals("bamCfg")){
                readBamConfigFile(argTokens[1]);                
            }            
            if (argTokens[0].equals("cutStd")){
                cutStd = Integer.parseInt(argTokens[1]);
            }
            if (argTokens[0].equals("maxD")){
                clusteringDist = Integer.parseInt(argTokens[1]);
            }
            if (argTokens[0].equals("minQ")){
                minMapQ = Integer.parseInt(argTokens[1]);
            }
            
            if (argTokens[0].equals("pMax")){
                patternMaxSpan = Integer.parseInt(argTokens[1]) * fragMean;
            }
            
            if (argTokens[0].equals("chrom")){
                String givenRegion = argTokens[1];
                if (!givenRegion.contains(":")){
                    chrom = givenRegion;
                }else{
                    chrom = givenRegion.split(":")[0];
                    givenRegionS = Integer.parseInt(givenRegion.split(":")[1].split("-")[0]);
                    givenRegionE = Integer.parseInt(givenRegion.split(":")[1].split("-")[1]);
                }
            }
            if (argTokens[0].equals("fa")){
                fastaFile = argTokens[1];
                fastaIndexFile = fastaFile + ".fai";
            }
            if (argTokens[0].equals("minAf")){
                minAf = Double.parseDouble(argTokens[1]);
            }
            if (argTokens[0].equals("fm")){
                fragMean = Integer.parseInt(argTokens[1]);
            }
            if (argTokens[0].equals("itemOut")){
                superitemOut = workDir + argTokens[1];
            }
            if (argTokens[0].equals("svOut")){
                svOut = argTokens[1];
            }
            if (argTokens[0].equals("minFreq")){
                minFreq = Integer.parseInt(argTokens[1]);
            }
            if (argTokens[0].equals("minWeight")){
                minWeight = Integer.parseInt(argTokens[1]);
            }
            // Debug usage files
            if (argTokens[0].equals("sigOut")){
                abnormalSignalOut = workDir + "wgs.abnormal.signals.txt";
            }
            if (argTokens[0].equals("mergeOut")){
                mergedPatternOut = workDir + "wgs.merged.patterns.txt";
            }
            if (argTokens[0].equals("freOut")){
                frequentPatternOut = workDir + "wgs.frequent.patterns.txt";
            }
            if (argTokens[0].equals("indelOut")){
                indelOut = workDir + "wgs.indels.txt";
            }
        }
//        loadExcluableRegion();
    }
    
    private void readBamConfigFile(String cfgFile) throws IOException{
        FileInputStream fin = new FileInputStream(new File(cfgFile));
        BufferedReader myInput = new BufferedReader(new InputStreamReader(fin));
        String thisLine;
        String sampleName = "";
//        StringBuilder sb = new StringBuilder();
        while ((thisLine = myInput.readLine()) != null){
            String[] tokens = thisLine.split(":");

            if (tokens[0].equals("bam")){
                bamFile = tokens[1];
            }
            if (tokens[0].equals("mean")){
                fragMean = Integer.parseInt(tokens[1]);
                patternMaxSpan = 4 * fragMean;
            }
            if (tokens[0].equals("stdev")){
                fragStd = Integer.parseInt(tokens[1]);
            }
            if (tokens[0].equals("readlen")){
                readLen = Integer.parseInt(tokens[1]);
                clusteringDist = readLen;
            }
            if (tokens[0].equals("workDir")){
                workDir = tokens[1];

            }
            if (tokens[0].equals("name")){
                sampleName = tokens[1];
            }
        }
        // Default output
        superitemOut = workDir + sampleName + ".superitems.txt";
        svOut = workDir + sampleName + ".mako.sites.txt";
        bndOutPath = workDir + sampleName + ".mako.BNDs.txt";
    }   
    
//    private void loadExcluableRegion() throws IOException{
//        String cenTelRegionFile = "/Users/jiadonglin/SV_data/ref_genome/centromeres.bed";
//        String gapRegionFile = "/Users/jiadonglin/SV_data/ref_genome/hg38.gap.bed";
//        String lowMapQRegionFile = "/Users/jiadonglin/SV_data/ref_genome/hg38.RegionsExcludable.bed";
        
//        String cenTelRegionFile = workDir + "centromeres.bed";
//        String gapRegionFile = workDir + "hg38.gap.bed";
//        String lowMapQRegionFile = workDir + "hg38.RegionsExcludable.bed";
        
//        excludableRegion = new ArrayList<>(3);
//        excludableRegion.add(cenTelRegionFile);
//        excludableRegion.add(gapRegionFile);
//        excludableRegion.add(lowMapQRegionFile);
//    }
}
