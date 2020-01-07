/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package structures;
import java.util.*;

/**
 *
 * @author jiadonglin
 */
public class Sequence {
    
    /**
     * A sequence is like <(ARP_LARGE_INSERT)(MS)(SM)(ARP_LARGE_INSERT)>
     * There might be several superitems in a set
     */
    private List<List<SuperItem>> itemsets = new ArrayList<List<SuperItem>>();
    /** Each sequence has an identical id */
    private int id;
    
    private List<SuperItem> smallIndelSuperitem = new ArrayList<SuperItem>();
    
    /**
     * Sequence constructor
     * @param id 
     */
    public Sequence(int id){
        this.id = id;
    }
    /**
     * Add a new superitem-set to the sequence.
     * @param itemset 
     */
    public void addItemset(List<SuperItem> itemset){
        itemsets.add(itemset);
    }
    /**
     * Return a string representation of this sequence
     */
    public String sequence2String(){
        StringBuilder r = new StringBuilder("");
        for (List<SuperItem> itemset : itemsets){
            r.append('(');
            for (SuperItem item : itemset){
                r.append(item.type);
                r.append(' ');
            }
            r.append(')');
        }
        return r.append("  ").toString();
    }
    public void printSequence(){
        System.out.println(sequence2String());
    }
    public int getId(){
        return id;
    }
    /** Return all superitem-sets in the sequence */
    public List<List<SuperItem>> getItemsets(){
        return itemsets;
    }
    public Sequence getSubSequence(int itemsetIdx){
        Sequence newSequence = new Sequence(getId());
        List<List<SuperItem>> newItemsets = itemsets.subList(itemsetIdx, itemsets.size());
        newSequence.itemsets = newItemsets;
        return newSequence;
    }
    
    
    /** Get an superitem-sets at a given position in the sequence*/
    public List<SuperItem> get(int idx){
        return itemsets.get(idx);
    }
    /** Get the size of this sequence, number of superitem-sets*/
    public int size(){
        return itemsets.size();
    }
    
    public void setSmallIndelSuperitem(List<SuperItem> alist){
        smallIndelSuperitem = alist;
    }
    public List<SuperItem> getLastSet(){
        return itemsets.get(itemsets.size());
    }
    
    public SuperItem superItemAtPos(int indexItemset, int indexItem){
        List<SuperItem> superitemSet = get(indexItemset);
        return superitemSet.get(indexItem);
    }
    public Sequence removeInFreSuperitemFromSequence(Map<String, Integer> superitemTypeCount, double relativeMinSup){
        Sequence sequence = new Sequence(getId());
        for (List<SuperItem> itemset : itemsets){
            List<SuperItem> newItemset = removeInFreSuperitemFromItemset(itemset, superitemTypeCount, relativeMinSup);
            if(newItemset.size() != 0){
                sequence.addItemset(newItemset);
            }
        }
        return sequence;
    }
    public List<SuperItem> removeInFreSuperitemFromItemset(List<SuperItem> itemset, Map<String, Integer> superitemTypeCount, double relativeMinSup){
        List<SuperItem> newItemset = new ArrayList<SuperItem>();       
        for (SuperItem superitem : itemset){
            String superitemType = superitem.type;
            int support = superitemTypeCount.get(superitemType);
            if (support >= relativeMinSup){
                newItemset.add(superitem);
            }
        }
        return newItemset;
    }
    
}
