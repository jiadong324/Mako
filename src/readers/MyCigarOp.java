/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package readers;

/**
 *
 * @author jiadonglin
 */
public class MyCigarOp {
    String op;
    int opLength;

    public MyCigarOp(String op, int opLength) {
        this.op = op;
        this.opLength = opLength;
    }
    
    public String getOp(){
        return op;
    }
    public int getOpLength(){
        return opLength;
    }
}
