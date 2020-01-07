/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mako;
import option.Params;
import structures.SequenceDatabase;
import fspm.CFSPM;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import readers.FileReader;
import readers.SignalReader;


/**
 *
 * @author jiadonglin
 */
public class Mako {

      /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException{
    // TODO code application logic here

        Params paramLoader = new Params();
        paramLoader.loadParam(args);
        int minSup = paramLoader.minFreq;

        if (paramLoader.hasParamInput){
            String makoParams = paramLoader.getParams();
            System.out.println("===================================");
            System.out.println(makoParams);
            System.out.println("Working directory: " + paramLoader.workDir);
            System.out.println("===================================\n");

            SignalReader myReader = new SignalReader(paramLoader.fragMean, paramLoader.fragStd, paramLoader.cutStd,
             paramLoader.readLen, paramLoader.clusteringDist, paramLoader.minMapQ, paramLoader.fastaIndexFile);

            myReader.doWork(paramLoader.bamFile, paramLoader.chrom, paramLoader.givenRegionS, paramLoader.givenRegionE,
             paramLoader.superitemOut, paramLoader.abnormalSignalOut, paramLoader.bndOutPath);


             SequenceDatabase sequenceDatabase = new SequenceDatabase();

             System.out.println("Graph nodes created, start detection!!\nSV output path: " + paramLoader.svOut);

             System.out.println("\nPattern growth parameters " +
             "\nmaxSpan:" + paramLoader.patternMaxSpan + "\nFreq:" + paramLoader.minFreq + "\nminAf:" + paramLoader.minAf + "\nminWeight:" + paramLoader.minWeight);

             BufferedWriter svRegionWriter = new BufferedWriter(new FileWriter(paramLoader.svOut));

             String makoHeader = paramLoader.createMakoHeader(makoParams);
             svRegionWriter.write(makoHeader);

             sequenceDatabase.loadSequencesFromFile(paramLoader.superitemOut, paramLoader.minAf, paramLoader.minWeight);

             CFSPM algoContiguousFSPM = new CFSPM(minSup, paramLoader.patternMaxSpan);

             algoContiguousFSPM.setParams(myReader.chromNameMap, paramLoader.regionMaskFile);
             algoContiguousFSPM.runAlgorithm(sequenceDatabase, paramLoader.frequentPatternOut, paramLoader.mergedPatternOut, svRegionWriter, paramLoader.fastaFile);
             algoContiguousFSPM.printAlgoStatistics();

        }
        /**

        String work_dir = "D:/mako_works/yri_csvs/illumina/mako_debug";
        String svOut = work_dir + "/mako_debug.txt";
        String superitemOut = work_dir + "/NA19240" + ".superitems.txt";
        String mergedPatterns = work_dir + "/merged_patterns.txt";
        String fastaIndexFile = "D:/data/ref_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai";
        String fastaFile = "D:/data/ref_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa";


//        String svOut ="D:/mako_works/skbr3/mako_skbr3.txt";
//        String bndOut = "D:/mako_works/skbr3/mako_skbr3.bnds.txt";
//        String superitemOut = "D:/mako_works/skbr3/skbr3.superitems.txt";
//
//        String bamFile = "D:/data/bams/SKBR3_550bp_pcrFREE_S1_L001_AND_L002_R1_001.101bp.bwamem.ill.mapped.sort.bam";
//        String fastaIndexFile = "D:/data/ref_genome/hs37d5.fa.fai";
//        String fastaFile = "D:/data/ref_genome/hs37d5.fa";

        SignalReader myReader = new SignalReader(569, 156, 3, 126, 126, 20, fastaIndexFile);

//        myReader.doWork(bamFile, null, 0, 0, superitemOut, null, bndOut);

        SequenceDatabase sequenceDatabase = new SequenceDatabase();

        BufferedWriter svRegionWriter = new BufferedWriter(new FileWriter(svOut));

        sequenceDatabase.loadSequencesFromFile(superitemOut, 0.2);

        CFSPM algoContiguousFSPM = new CFSPM(10, 2000);


        algoContiguousFSPM.setParams(myReader.chromNameMap, null);
        algoContiguousFSPM.runAlgorithm(sequenceDatabase, null, null, svRegionWriter, fastaFile);
        algoContiguousFSPM.printAlgoStatistics();
        */

    }
}
