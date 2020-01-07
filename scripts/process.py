#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin, Xi'an Jiaotong Univeristy, Leiden University

@contact: jiadong324@gmail.com

@time: 2019/11/4
'''


import sys

from optparse import OptionParser
import pysam
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import re

parser = OptionParser()


class Interval:
    def __init__(self, chrom, start, end, pattern, sample, interval_str):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.pattern = pattern
        self.sample = sample
        self.interval = interval_str

    def overlap(self, interval, max_dist, len_prop):
        this_size = self.end - self.start
        inter_size = interval.end - interval.start

        min_len = this_size * len_prop
        max_len = this_size * (2 - len_prop)

        # Two intervals overlap
        if min(self.end, interval.end) >= max(self.start, interval.start) and self.chrom == interval.chrom:
            # Check breakpoint distance
            if abs(self.start - interval.start) <= max_dist and abs(self.end - interval.end) <= max_dist:
                # Check if SV size matches
                if inter_size >= min_len and inter_size <= max_len:
                    return True

        return False

    def toString(self):

        out_str = "{0}\t{1}\t{2}\t".format(self.chrom, self.start, self.end)
        interval_tokens = self.interval.split(";")
        sample_tokens = self.sample.split(";")

        sample_str = ""
        for i in range(len(sample_tokens)):
            sample_str += "{0},{1};".format(sample_tokens[i], interval_tokens[i])

        out_str += sample_str[:-1] + "\t" + self.pattern

        return out_str


def mako_to_vcf(mako, out, ref, sample):
    '''
    Convert Mako raw output to standard VCF format
    :param mako:
    :param out:
    :param ref:
    :param sample:
    :return:
    '''

    print(sample, " convert to VCF ...")

    columns = (0, 1, 2, 3, 4)
    names = ("chr", "start", "end", "filter", "info")

    calls = pd.read_table(mako, header=None, usecols=columns, names=names)
    sorted_calls = calls.sort_values(['chr', 'start'], ascending=[True, True])

    reference = pysam.FastaFile(ref)

    sorted_calls["sample_name"] = sample
    sorted_calls["call_id"] = "."
    sorted_calls["quality"] = "."

    sorted_calls["reference"] = sorted_calls.apply(
        lambda row: reference.fetch(row.chr, row.start, row.start + 1).upper(), axis=1)
    sorted_calls["alt"] = "."

    sorted_calls["info"] = sorted_calls.apply(
        lambda row: "END={0};SAMPLES={1};{2}".format(row.end, row.sample_name, row.info), axis=1)

    simple_calls = sorted_calls[["chr", "start", "call_id", "reference", "alt", "quality", "filter", "info"]].rename(
        {"chr": "#CHROM", "start": "POS", "reference": "REF", "call_id": "ID", "quality": "QUAL", "info": "INFO",
         "alt": "ALT", "filter": "FILTER"}, axis=1)

    with open(out, "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##source=Mako V1.0\n")
        vcf.write(
            '##INFO=<ID=ARP_Span,Number=1,Type=String,Description="SVs are combined from two patterns linked through read-pairs">' + "\n")
        vcf.write('##INFO=<ID=ARP_Self,Number=1,Type=String,Description="SVs supported by read-pairs">' + "\n")
        vcf.write('##INFO=<ID=Split,Number=1,Type=String,Description="SVs supported by split alignment">' + "\n")
        vcf.write(
            '##INFO=<ID=Cross,Number=1,Type=String,Description="SVs supported by local sequence cross match">' + "\n")
        vcf.write(
            '##INFO=<ID=Realign,Number=1,Type=String,Description="SVs discovered from clipped Super-Item realignment">' + "\n")
        vcf.write(
            '##INFO=<ID=SMALL_INSERT,Number=1,Type=String,Description="SVs supported by read pairs of small insert size">' + "\n")
        vcf.write(
            '##INFO=<ID=OEM,Number=1,Type=String,Description="one-end-unmapped reads supported breakpoint">' + "\n")
        vcf.write(
            '##INFO=<ID=Pattern,Number=1,Type=String,Description="Super-Item pattern of this SV">' + "\n")
        vcf.write(
            '##INFO=<ID=Region,Number=1,Type=String,Description="Pattern spanned genome region">' + "\n")
        vcf.write(
            '##INFO=<ID=Weights,Number=1,Type=Integer,Description="Number of reads in each Super-Item">' + "\n")
        vcf.write(
            '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele fraction of each Super-Item (#abnormal reads / #reads)">' + "\n")
        vcf.write(
            '##INFO=<ID=Ori,Number=1,Type=String,Description="Reads orientation of each Super-Item">' + "\n")
        vcf.write('##INFO=<ID=LinkedMapQ,Number=1,Type=Integer,Description="Linked ARP Super-Item quality">' + "\n")
        vcf.write('##INFO=<ID=LinkedItem,Number=1,Type=String,Description="Linked Super-Item type">' + "\n")
        vcf.write(
            '##INFO=<ID=LinkedWeight,Number=1,Type=String,Description="Linked Super-Item allele fraction">' + "\n")
        vcf.write(
            '##INFO=<ID=Cross_read,Number=1,Type=Integer,Description="Number of read support sequence cross match">' + "\n")
        vcf.write(
            '##INFO=<ID=Cross_len,Number=1,Type=Integer,Description="Maximum string match of cross match">' + "\n")
        vcf.write('##INFO=<ID=Split_sup,Number=1,Type=Integer,Description="Number of split aligned reads">' + "\n")
        simple_calls.to_csv(vcf, sep="\t", index=False)


def merge_multiple_makos(sample_files, mako_dir, out_file, max_dist, len_prop):
    '''
    Merge multiple Mako call set
    :param sample_files:
    :param mako_dir:
    :param out_file:
    :param max_dist:
    :param len_prop:
    :return:
    '''

    intervals = []
    for line in open(sample_files, "r"):
        file_name = line.strip()
        mako_file_path = mako_dir + file_name
        sample_name = file_name.split(".")[0]

        cur_sample_sv_num = 0

        for line in open(mako_file_path, "r"):
            tmp = line.strip().split("\t")
            chrom = tmp[0]
            start = int(tmp[1])
            end = int(tmp[2])

            cur_sample_sv_num += 1
            sv_info_tokens = tmp[4].split(';')
            pattern_str = ""
            for token in sv_info_tokens:
                if token.split('=')[0] == 'Pattern':
                    pattern_str = token.split('=')[1]
                    break

            interval_str = "{0},{1},{2}".format(chrom, start, end)

            cur_interval = Interval(chrom, start, end, pattern_str, sample_name, interval_str)

            intervals = add_interval(intervals, cur_interval, max_dist, len_prop)

        # print sample_name + ", " + str(cur_sample_sv_num) + " SVs processed .."
        print("Merge sample: {0} total entries: {1}".format(sample_name, len(intervals)))

    writer = open(out_file, "w")

    for interval in intervals:
        writer.write(interval.toString() + "\n")
    writer.close()


def add_interval(intervals, new_interval, max_dist, len_prop):
    new_intervals = []

    num = len(intervals)

    if num == 0:
        new_intervals.append(new_interval)
        return new_intervals

    if new_interval.end < intervals[0].start or new_interval.start > intervals[num - 1].end:

        if new_interval.end < intervals[0].start:
            new_intervals.append(new_interval)

        new_intervals.extend(intervals)

        if new_interval.start > intervals[num - 1].end:
            new_intervals.append(new_interval)

        return new_intervals

    for i in range(len(intervals)):
        ele = intervals[i]
        overlap = ele.overlap(new_interval, max_dist, len_prop)
        # Overlapped
        if not overlap:
            new_intervals.append(ele)

            # check if given interval lies between two intervals
            if i < num and new_interval.start > intervals[i].end and new_interval.end < intervals[i + 1].start:
                new_intervals.append(new_interval)

            continue

        new_start = min(ele.start, new_interval.start)
        new_pattern = ele.pattern
        new_end = max(ele.end, new_interval.end)
        new_sample = ele.sample
        new_interval_str = ele.interval

        while i < num and overlap:
            new_end = max(intervals[i].end, new_interval.end)
            new_pattern += ";" + new_interval.pattern
            new_sample += ";" + new_interval.sample
            new_interval_str += ";" + new_interval.interval
            if i == num - 1:
                overlap = False
            else:
                overlap = intervals[i + 1].overlap(new_interval, max_dist, len_prop)

            i += 1

        i -= 1

        new_intervals.append(
            Interval(new_interval.chrom, new_start, new_end, new_pattern, new_sample, new_interval_str))

    return new_intervals


def mako_filter(in_file, out_file, cxs, format):
    '''
    Filter mako raw call site with different evidence
    :param in_file:
    :param out_file:
    :return:
    '''
    writer = open(out_file, 'w')
    scores = []
    csvs_num = 0
    all_calls = 0
    for line in open(in_file, 'r'):
        if "#" in line:
            continue
        all_calls += 1
        tmp = line.strip().split("\t")
        sv_info_tokens = tmp[4].split(";")[1:]
        cx_score = int(tmp[4].split(";")[0].split("=")[1])
        # if cx_score not in scores:
        #     scores.append(cx_score)
        if cx_score > cxs:
            csvs_num += 1
            if format == "mako":
                writer.write(line)
            elif format == "bed":
                sv_len = int(tmp[2]) - int(tmp[1])
                out_str = "{0}\t{1}\t{2}\t{3}\n".format(tmp[0], tmp[1], tmp[2], sv_len)
                writer.write(out_str)
    print("Number of calls after filtering: ", csvs_num)
    writer.close()
    # print(np.percentile(scores, 25))

def mako_config(bam, num, mad, out, sample):
    required = 97
    restricted = 3484
    flag_mask = required | restricted

    read_length = 0
    read_counter = 0

    L = []

    bam_file = pysam.AlignmentFile(bam, "r")

    for read in bam_file.fetch():
        if read_counter >= num:
            break

        cigar = read.cigarstring
        if cigar == None:
            continue

        read_length = get_read_length(cigar)
        flag = read.flag
        refname = read.reference_name
        mate_refname = read.next_reference_name
        isize = read.template_length

        valid = mate_refname == refname and flag & flag_mask == required and isize >= 0

        if valid:
            read_counter += 1
            L.append(isize)

    L = np.array(L)
    L.sort()
    med, umad = unscaled_upper_mad(L)
    upper_cutoff = med + mad * umad
    L = L[L < upper_cutoff]

    mean = int(np.mean(L))
    stdev = int(np.std(L))

    print("Estimated mean: {0}, std: {1}, read length: {2}".format(mean, stdev, read_length))

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.hist(L, bins=50)
    ax.set_title("Insert size distribution")

    ax.set_xlabel("insert size (bp)")
    plt.tight_layout()
    # plt.show()

    bam_abs_path = out + bam

    plt.savefig(out + sample + ".config.png")

    outStr = 'mean:' + str(mean) + '\nstdev:' + str(stdev) + '\nreadlen:' + str(
        read_length) + '\nworkDir:' + out + '\nbam:' + bam_abs_path + '\nname:' + sample

    sys.stdout = open(out + sample + '.mako.cfg', 'w')
    print(outStr)


def unscaled_upper_mad(xs):
    """Return a tuple consisting of the median of xs followed by the
    unscaled median absolute deviation of the values in xs that lie
    above the median.
    """
    med = np.median(xs)
    return med, np.median(xs[xs > med] - med)


def get_read_length(cigar):
    cigarPattern = '([0-9]+[MIDNSHP])'
    cigarSearch = re.compile(cigarPattern)
    atomicCigarPattern = '([0-9]+)([MIDNSHP])'
    atomicCigarSearch = re.compile(atomicCigarPattern)

    readLen = 0
    if (cigar != '*'):
        cigarOpStrings = cigarSearch.findall(cigar)

        for opString in cigarOpStrings:
            cigarOpList = atomicCigarSearch.findall(opString)[0]
            readLen += int(cigarOpList[0])

    return readLen

def get_mako_sub(mako, out):

    ind_svs = []

    for line in open(mako, "r"):
        if "#" in line:
            continue

        tmp = line.strip().split("\t")
        chrom = tmp[0]

        sv_info = tmp[4]
        if ";;" in sv_info:
            sv_info_tokens = sv_info.split(";;")
            for info_token in sv_info_tokens:
                tmp_token = info_token.split(",")
                if len(tmp_token) == 2 and "-" in tmp_token[0]:
                    this_start = tmp_token[0].split("-")[0]
                    this_end = tmp_token[0].split("-")[1]
                    this_sv = (chrom, this_start, this_end, tmp_token[1])
                    ind_svs.append(this_sv)
                else:
                    for i in range(2, len(tmp_token)):
                        token = tmp_token[i]
                        if "-" in token:
                            this_start = token.split("-")[0]
                            this_end = token.split("-")[1]
                            this_sv = (chrom, this_start, this_end, tmp_token[0] + "," + tmp_token[1] + "," + tmp_token[i + 1])
                            ind_svs.append(this_sv)
        else:
            sv_info_tokens = sv_info.split(",")
            for i in range(2, len(sv_info_tokens)):
                token = sv_info_tokens[i]
                if "-" in token:
                    this_start = token.split("-")[0]
                    this_end = token.split("-")[1]
                    this_sv = (chrom, this_start, this_end, sv_info_tokens[0] + "," + sv_info_tokens[1] + "," + sv_info_tokens[i + 1])
                    ind_svs.append(this_sv)


    writer = open(out, "w")

    for sv in ind_svs:

        out_str = "{0}\t{1}\t{2}\t{3}\n".format(sv[0], sv[1], sv[2], sv[3])

        writer.write(out_str)

    writer.close()


def print_config_params():
    print("\nCreate the config file of a BAM for Mako input")
    print("\nParameters:")
    print(" -b\tBAM file to process")
    print(" -w\tWorking directory (BAM file under this directory)")
    print(" -n\tNumber of reads used to estimate")
    print(" -s\tsample name")


def print_filter_params():
    print("\nFilter Mako calls support by ARP_Self.\nWe are working on score function for more reliable result filtering")
    print("\nParameters:")
    print(" -i\tMako  callset")
    print(" -o\tfiltered Mako callset")
    print(" -a\tminimum allele fraction (0.2 is recommended)")


def print_merge_params():
    print("\nMerge mutiple Mako calls to a single file")
    print("\nNote: each Mako output has to be sorted")
    print("Calls from different samples are merged if their size are within a range and breakpoint distance is close.")
    print("\nFor example: sv1 = [a, b] and sv2 = [c, d] ")
    print("Similar size range [size(sv1) * len_prop, size(sv1) * (2 - len_prop)]\tbp_dist=abs(c-a) and abs(d-b)")

    print("\nParameters:")
    print(" -s\tList of sample names to merge")
    print(" -c\tThe directory of Mako output")
    print(" -o\tMerged output file")
    print(" -d\tMax distance between breakpoints")
    print(" -l\tLength proportion (recommended 0.5, len_prop)")


def print_tovcf_params():
    print("\nConvert Mako output to a VCF format")
    print("\nParameters:")
    print(" -m\tMako original call set")
    print(" -o\tvcf output")
    print(" -r\treference genome")
    print(" -s\tsample name")


script_name = sys.argv[0]
if len(sys.argv) < 2:
    print('=======================================================')
    print('process.py         Last Update:2019-7-20\n')
    print('This script is used to process Mako raw callset\n')
    print('Usage:')
    print('process.py [options] <parameters>\n')
    print('Options:')
    print('config: Create config file for Mako input ')
    print('filter: filter Mako calls ')
    print('makotovcf: convert Mako calls to standard VCF format')
    # print('merge: merge mutiple Mako calls')

    print("=======================================================")
else:
    option = sys.argv[1]

    if option == "config":
        if len(sys.argv) == 2:
            print_config_params()

        else:
            parser.add_option("-b", dest='bam', help='BAM file to config')
            parser.add_option("-n", type=int, dest="N", help="Number of samples used for estimation")
            parser.add_option("-m", dest="mads", type=int, default=30,
                              help="Outlier cutoff in # of median absolute deviations (unscaled, upper only)")
            parser.add_option("-w", dest='out', help='Working directory')
            parser.add_option("-s", dest="name", help="Name of the sample")

            (options, args) = parser.parse_args()
            mako_config(options.bam, options.N, options.mads, options.out, options.name)

    elif option == "makotovcf":
        if len(sys.argv) == 2:
            print_tovcf_params()
        else:
            parser.add_option("-m", dest="mako")
            parser.add_option("-o", dest="out")
            parser.add_option("-r", dest="ref")
            parser.add_option("-s", dest="sample")
            (options, args) = parser.parse_args()

            if not options.mako:
                parser.error("Mako call not given")

            if not options.out:
                parser.error("VCF output not given")

            if not options.ref:
                parser.error("Please give reference genome")

            if not options.sample:
                parser.error("Please give sample name")

            mako_to_vcf(options.mako, options.out, options.ref, options.sample)

    elif option == "filter":
        if len(sys.argv) == 2:
            print_filter_params()
        else:
            parser.add_option("-i", dest="input", help="Input Mako callset")
            parser.add_option("-o", dest="out", help="Output of filtered Mako callset by CXS")
            parser.add_option("-c", type=int, dest="cxs", help="CXS threshold")
            parser.add_option("-f", dest="format", help="output format (original, BED)")
            (options, args) = parser.parse_args()
            mako_filter(options.input, options.out, options.cxs, options.format)


# if __name__ == '__main__':
    # work_dir = '/mnt/d/mako_works/yri_csvs/illumina/mako_debug/'
    # work_dir = '/mnt/d/mako_works/trios/'
    # work_dir = '/mnt/d/mako_works/skbr3/mako_calls/'

    # work_dir = '/mnt/d/mako_works/visor_sim/known_csv_chr1/version2/p100_c5/'
    # mako_filter(work_dir + "chr1_sim.mako.sites.txt", work_dir + "mako_filtered.txt")
    # mako_filter(work_dir + "HG00514.mako.sites.txt", work_dir + "mako_BP_filtered.txt")
    # mako_filter(work_dir + "mako_skbr3.txt", work_dir + "mako_bp_filtered.cxs4.skbr3.txt")
    # work_dir = '/mnt/d/data/bams/'
    # work_dir = '/mnt/d/mako_works/skbr3/mako_calls/'
    # mako_filter(work_dir + "mako_skbr3.txt", work_dir + "mako_bp_filtered.skbr3.txt")
    # mako_config(work_dir + "SKBR3_550bp_pcrFREE_S1_L001_AND_L002_R1_001.101bp.bwamem.ill.mapped.sort.bam", 30000, 30, work_dir, 'skbr3')





