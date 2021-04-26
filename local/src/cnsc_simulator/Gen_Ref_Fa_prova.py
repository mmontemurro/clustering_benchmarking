# uncompyle6 version 3.7.4
# Python bytecode 2.7 (62211)
# Decompiled from: Python 3.6.10 | packaged by conda-forge | (default, Apr 24 2020, 16:44:11) 
# [GCC 7.3.0]
# Embedded file name: /home/data/work/phylics_sim/cnsc_simulator/Gen_Ref_Fa.py
# Compiled at: 2020-11-30 10:55:47
from CN import CN
import os, numpy as np

def make_fa(ID, tree, ref, chr_name_array, fa_prefix):
    fa_f_prefix = fa_prefix + str(ID) + '_'
    trace = [
     ID]
    visit = ID
    while visit != 0:
        visit = tree[visit].parentID
        trace.append(visit)

    CN = []
    for i in range(len(trace)):
        j = len(trace) - i - 1
        for cn in tree[trace[j]].cn:
            CN.append(cn)

    new_ref = gen_ref(ref, CN)
    write_ref(new_ref, chr_name_array, fa_f_prefix)


def getlen_ref(template):
    chr_name = []
    len_chr = []
    file = open(template, 'r')
    len_chr_ = 0
    line = file.readline().rstrip('\n')
    while line != '':
        if line[0] == '>':
            chr_name.append(line[1:])
            if len_chr_ != 0:
                len_chr.append(len_chr_)
            len_chr_ = 0
        else:
            len_chr_ = len_chr_ + len(line)
        line = file.readline().rstrip('\n')

    len_chr.append(len_chr_)
    file.close()
    return (chr_name, len_chr)


def init_ref_from_npy(prefix):
    len_chr = np.load(prefix + '.leaf_chrlen.npy')
    leaf_index = np.load(prefix + '.leaf_index.npy')
    chr_name = np.load(prefix + '.chr_name_array.npy')
    ref = np.load(prefix + '.ref_array.npy')
    corres = np.load(prefix + '.corres_array.npy')
    print 'Load:'
    print ref
    print chr_name
    print len_chr
    print leaf_index
    for c in corres[0]:
        for ci in c:
            print ci

    for c in corres[1]:
        for ci in c:
            print ci

    return (
     ref, chr_name, len_chr, corres)


def init_ref(template):
    ref = []
    chr_name = []
    len_chr = []
    file = open(template, 'r')
    line = 'tmp'
    str_ = ''
    line = file.readline().rstrip('\n')
    while line != '':
        if line[0] == '>':
            chr_name.append(line[1:])
            if str_ != '':
                ref.append(str_)
                len_chr.append(len(str_))
            str_ = ''
        else:
            str_ = str_ + line
        line = file.readline().rstrip('\n')

    ref.append(str_)
    len_chr.append(len(str_))
    file.close()
    ref_diploid = [ref, ref]
    print 'Load:'
    print ref_diploid
    print chr_name
    print len_chr
    return (
     ref_diploid, chr_name, len_chr)


def gen_ref(ref, CNs):
    ret_ref = [ row[:] for row in ref ]
    for i in range(len(CNs)):
        ale = CNs[i].get_CN_Ale()
        pos1, pos2 = CNs[i].get_CN_position()
        chr_ = CNs[i].get_CN_chromosome()
        if CNs[i].CN_Del == 0:
            amp_num = CNs[i].get_CN_amp_num()
            str_amp = amp_num * ret_ref[ale][chr_][pos1:pos2]
            ret_ref[ale][chr_] = ret_ref[ale][chr_][:pos2] + str_amp + ret_ref[ale][chr_][pos2:]
        else:
            ret_ref[ale][chr_] = ret_ref[ale][chr_][:pos1] + ret_ref[ale][chr_][pos2:]

    return ret_ref


def read_ref(ref):
    ref_a = []
    for i in [1, 2]:
        file_ = ref + str(i) + '.fa'
        if os.path.isfile(file_):
            file = open(file_, 'r')
            ref_ = []
            str_ = ''
            for line in file:
                line = line.strip()
                if line[0] == '>':
                    if str_ != '':
                        ref_.append(str_)
                        str_ = ''
                else:
                    str_ = str_ + line

            if str_ != '':
                ref_.append(str_)
            file.close()
            ref_a.append(ref_)
        else:
            print file_ + ' does not exist and cannot be opened.'

    return ref_a


def write_ref(ref, chr_name, fasta):
    line_len = 60
    for i in [1, 2]:
        file = open(fasta + str(i) + '.fa', 'w')
        for j in range(len(ref[(i - 1)])):
            file.write('>' + chr_name[j] + '\n')
            if len(ref[(i - 1)][j]) < line_len:
                file.write(ref[(i - 1)][j] + '\n')
            else:
                tmp_str = ref[(i - 1)][j]
                a = 0
                while a + line_len < len(tmp_str):
                    file.write(tmp_str[a:a + line_len] + '\n')
                    a = a + line_len

                if len(tmp_str) > 0:
                    file.write(tmp_str[a:] + '\n')

        file.close()


def make_ref(ID, tree, ref):
    trace = [
     ID]
    visit = ID
    while visit != 0:
        visit = tree[visit].parentID
        trace.append(visit)

    CN = []
    for i in range(len(trace)):
        j = len(trace) - i - 1
        for cn in tree[trace[j]].cn:
            CN.append(cn)

    return gen_ref(ref, CN)
# okay decompiling Gen_Ref_Fa.pyc
