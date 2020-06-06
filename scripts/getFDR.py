#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin, Xi'an Jiaotong University, Leiden University

@contact: jiadong324@gmail.com

@time: 2020/5/31
'''

import matplotlib.pylab as plt
import seaborn as sns
import os
import pysam
from Bio import SeqIO, SeqRecord, Seq
from Bio.Alphabet import DNAAlphabet
from intervaltree import IntervalTree

def read_mako(input, output, exclude, bam, seqout):
    exclude_dict = {}

    with open(exclude, 'r') as f:
        for line in f:
            entries = line.strip().split("\t")
            start = int(entries[1])
            end = int(entries[2])
            if entries[0] not in exclude_dict:
                exclude_dict[entries[0]] = IntervalTree()
                exclude_dict[entries[0]][start:end] = (start, end)
            else:
                exclude_dict[entries[0]][start:end] = (start, end)
    counter = 0
    all_events = 0
    writer = open(output, 'w')
    samfile = pysam.AlignmentFile(bam, 'rb')
    with open(input, 'r') as f:
        for line in f:
            if '#' in line:
                continue
            all_events += 1
            entries = line.strip().split("\t")
            info_tokens = entries[4].split(";")
            cxs = 0
            for token in info_tokens:
                if token.split("=")[0] == 'cxs':
                   cxs = int(token.split("=")[1])

            this_chrom_exclude_tree = exclude_dict[entries[0]]
            sv_size = int(entries[2]) - int(entries[1])
            if sv_size >= 50000 or entries[0] == "chrY":
                continue
            if cxs >= 10:
                if this_chrom_exclude_tree.overlap(int(entries[1]), int(entries[2])):
                    continue
                tmpdir = seqout + '{0}_{1}_{2}'.format(entries[0], entries[1], entries[2])

                if not os.path.exists(tmpdir):
                    os.mkdir(tmpdir)
                write_seq_to_fastq(samfile, entries[0], int(entries[1]), int(entries[2]), tmpdir)

                writer.write("{0}\t{1}\t{2}\n".format(entries[0], entries[1], entries[2]))
                counter += 1
            # scores.append(cxs)
    # sns.kdeplot(scores)
    # plt.show()
    # writer.close()

    print("All events: ", all_events)
    print("CXS>=10", counter)

def write_seq_to_fastq(bamfile, chrom, start, end, out_dir):

    for read in bamfile.fetch(chrom, start, end):
        new_read_name = read.query_name.replace("/","_")
        out_fa_name = out_dir + '/' + new_read_name + '.fa'

        seq_record = SeqRecord.SeqRecord(Seq.Seq(read.query_sequence,DNAAlphabet),name=read.query_name,id=read.query_name, description="")
        with open(out_fa_name, 'w') as f:
            SeqIO.write(seq_record, f, 'fasta')


if __name__ == '__main__':
    workdir = '/Users/jiadonglin/SV_Research/onGoingProjects/Mako/'
    hg00514 = '/Users/jiadonglin/Data/HG00514/HG00514.nglmr.srt.bam'
    na19240 = '/Users/jiadonglin/Data/NA19240/NA19240.nglmr.srt.bam'
    # hg00733 = '/Users/jiadonglin/Data/HG00733/HG00733.nglmr.srt.bam'
    # read_mako(workdir + 'NA19240.mako.sites.txt', workdir + 'mako_csv_fdr/NA19240.cxs10.exclude.50K.txt', workdir + 'mako_csv_fdr/grch38.excludeRegions.bed', na19240, workdir + 'mako_csv_fdr/NA19240_plots/')