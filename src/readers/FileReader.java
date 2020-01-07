/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package readers;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import htsjdk.samtools.cram.structure.CramCompressionRecord;

import java.io.BufferedInputStream;
import java.io.InputStream;
import java.util.List;


/**
 *
 * @author jiadonglin
 */
public class FileReader {

    
    public FileReader(){
    }
    
    public SamReader openBamFile(String bamFilePath, ValidationStringency stringency, boolean includeFileSource) throws IOException{
        SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().validationStringency(stringency).enable(SamReaderFactory.Option.DONT_MEMORY_MAP_INDEX);
        if(includeFileSource){
            samReaderFactory.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS);            
        }
//        samReaderFactory.samRecordFactory(DefaultSAMRecordFactory.getInstance());
        final SamReader samReader = samReaderFactory.open(new File(bamFilePath));
        return samReader;
    }
    

    
    public ReferenceSequenceFile readFastaFile(String faFile) throws IOException{                        
        File fa = new File(faFile);
        return ReferenceSequenceFileFactory.getReferenceSequenceFile(fa);        
    }


    public String[] getInfoOfChrom(String faIdxFile, String chrom) throws IOException{
        String[] chromNameMap = new String[24];
        FileInputStream fin = new FileInputStream(new File(faIdxFile));
        BufferedReader myInput = new BufferedReader(new InputStreamReader(fin));
        int chrIdx = 0;
        if (chrom != null) {
            String thisLine;
            while ((thisLine = myInput.readLine()) != null){
                String[] tokens = thisLine.split("\t");
                String refSeq = tokens[0];
                // escape "M", "MT", "chrM"
                if (refSeq.equals(chrom)) {
                    chromNameMap[chrIdx] = chrom;
                }
                chrIdx ++;
            }

        }else{

            String thisLine;

            while ((thisLine = myInput.readLine()) != null){
                String[] tokens = thisLine.split("\t");
                String refSeq = tokens[0];
                // escape "M", "MT", "chrM"
                if (refSeq.contains("M") || refSeq.contains("_")){
                    break;
                }

                chromNameMap[chrIdx] = refSeq;
                chrIdx += 1;

            }
        }
        return chromNameMap;
    }
    
}
