import piperine
import numpy as np
import re
import argparse
from pkg_resources import resource_stream
import stickydesign
import os


parser = argparse.ArgumentParser(description='The script creates the inputs for the heyristic measures equations.')
parser.add_argument('-bn','--basename', default=None, type=str, help='Basename - File_Name of .pil, .seqs and .mfe to be evaluated')
parser.add_argument('-sig','--signals', nargs='*', default=None, type=str, help='The signal-complexes of the system, ex. -sig solo-A solo-B solo-C')
parser.add_argument('-g','--dGtarget', default = None, help='Wanted dGTarget')
parser.add_argument('-t','--ted', action='store_true', default='False', help='Additionaly, claculate BN% max, Max Defect Component, BN% avg')

args = parser.parse_args()

ted=args.ted
inter= args.dGtarget
SeqSeriesName = args.basename
signals = args.signals

# ted=False
# inter = None
# SeqSeriesName = 'Solov1_0_1'
# signals = ['solo-IA', 'solo-IB','solo-R','solo-TCC']


pil_file = SeqSeriesName+'.pil'
mfe_file = SeqSeriesName+'.mfe'
seq_file = SeqSeriesName+'.seqs'

TopStranddict = {}      # All toehold nucleotide and some more, postions of a strand to be checked for secondary structure

TopStrandlist = []      # All top strands
BaseStrandlist = []     # All complementary toeholds
BMlist = []             # All branch migration domains 
NotToInteract = {}      # t* : all top strands or top strand parts not containing t* ,
domains_list =[]        # All domains toeholds and their parts
th_strs = []            # All toeholds present in signals

strand_dict = {}        # Strand : its domains
            

complex_names = []      # All complexes names (2 or more strands in the molecule)
seq_dict = {}           # All complexes, strands, domains and their sequences
cmplx_dict = {}         # All complexes (from 1 strand) and their secondary structure
toeholds_arr = []       # All sequences of toehold present in signals
toeholds_ma = [] ## Currently no automation for this

p = open(pil_file,'r')
m = open(mfe_file,'r')
s = open(seq_file,'r')

m_l = m.readlines()
p_l = p.readlines()

for l in p_l:
    l_w = l.split()
    for w in range(len(l_w)):
        if l_w[w] == 'sequence' or l_w[w] == 'sup-sequence':
            if re.search('\-',l_w[w+1]):
                domains_list.append(l_w[w+1])           #creating domain_list from sequence or sup-sequence line
        if l_w[w] == 'strand':
            if not re.search('\*', l):
                TopStrandlist.append(l_w[w+1])          #creating TopStrandlist from strand lines without *
                jj=l_w[w+3:-2]
                strand_dict[l_w[w+1]]=jj                #creating a dict strand_name : its domains
                if len(jj) != 1:
                    for i in range(len(jj)):
                        if not re.search('-t',jj[i]):
                            if jj[i] not in BMlist:
                                BMlist.append(jj[i])    #creating BMlist from not toehold domains of TopStrands
            else:
                for i in range(len(l_w)):
                    if re.search('-t',l_w[i]):
                        if l_w[i] not in BaseStrandlist:
                            BaseStrandlist.append(l_w[i])   #creating BaseStrandlist from -t* in a strand line
        ## signal toeholds

for c_t in BaseStrandlist:
    interactions=[]
    for strand in strand_dict.keys():
        if c_t.replace("*", "") not in strand_dict[strand]:
            if strand not in interactions:
                interactions.append(strand)
        else:
            for dom in strand_dict[strand]:
                if dom not in interactions and dom != c_t.replace("*", ""):
                    interactions.append(dom)
    NotToInteract[c_t]=interactions                     #creating NoToInteract

for l in range(len(m_l)):
    if m_l[l][0].isdigit():
        words=m_l[l].split(':')
        cmplx_name=words[1].rstrip()                    #cmplx_name
        cmplx_scnd_st=m_l[l+2].rstrip()
        cmplx_dict[cmplx_name]=cmplx_scnd_st            # Append cmplx_name in cmplx_dict{}
        if re.search('\(',cmplx_scnd_st):
            complex_names.append(cmplx_name)            # Append cmplx_name with '(' in cmplx_scnd_st in complex_names 

seqs_lines = s.readlines()
for line in seqs_lines:
    seqs_wordList=line.split()
    for i in range(len(seqs_wordList)):
        a = seqs_wordList[i]
        if a == 'sequence': 
            sequence_name = seqs_wordList[i+1]       
            sequence_seq = seqs_wordList[i+3]
            seq_dict[sequence_name] = sequence_seq          
            comp_seq_name = sequence_name+'*'
            comp_seq_nucs=[]
            for pos in range(1,len(sequence_seq)+1):
                n=sequence_seq[-pos]
                if n=='A':
                    n_c = 'T'
                if n=='T':
                    n_c = 'A'
                if n=='G':
                    n_c = 'C'
                if n=='C':
                    n_c = 'G'
                comp_seq_nucs.append(n_c)
                comp_seq_seq=("".join(comp_seq_nucs))
            seq_dict[comp_seq_name] = comp_seq_seq
        if a == 'strand':
            strand_name = seqs_wordList[i+1]
            strand_seq = seqs_wordList[i+3]
            seq_dict[strand_name] = strand_seq
            comp_strand_name = strand_name+'*'
            comp_strand_nucs=[]
            for pos in range(1,len(strand_seq)+1):
                n=strand_seq[-pos]
                if n=='A':
                    n_c = 'T'
                if n=='T':
                    n_c = 'A'
                if n=='G':
                    n_c = 'C'
                if n=='C':
                    n_c = 'G'
                comp_strand_nucs.append(n_c)
                comp_strand_seq=("".join(comp_strand_nucs))
            seq_dict[comp_strand_name] = comp_strand_seq
        if a == 'structure':
            structure_name = seqs_wordList[i+1]       
            structure_seq = seqs_wordList[i+3]
            seq_dict[structure_name] = structure_seq

for sig in signals:
    a=[]
    b=[]
    for d in strand_dict[sig]:
        if re.search('-t', d): 
            a.append(d)
        if a not in th_strs:
            th_strs.append(a)                        
            for toe in a:
                b.append(seq_dict[toe])
            toeholds_ma.append(b)

for t , in th_strs:
    toeholds_arr.append(seq_dict[t])                #creating toeholds_arr

for topstrand in TopStrandlist:
    for dom_po in range(len(strand_dict[topstrand])):
        dom = strand_dict[topstrand][dom_po]
        if re.search('-t',dom):
            po=[]
            s=seq_dict[dom]
            toe_start=seq_dict[topstrand].find(s)+1
            start=toe_start
            toe_end=len(dom) + start -1
            end = toe_end
            if dom_po-1 in range(len(strand_dict[topstrand])):
                start=start-3
            if dom_po+1 in range(len(strand_dict[topstrand])):
                end=end+3
            for i in range(start,end):
                po.append(i)
            TopStranddict[topstrand]=po                                 # TopStranddict Inclues postions of toeholds in a  strand


################## 

def score_toeholds(toeholds,dGTarget):
    tops = lambda s: 4*s[:,:-1]+s[:,1:]
    nt = { 'a': 0, 'c': 1, 'g': 2, 't': 3 }
    taildG = 1.3
    initdG = 0

    dsb = resource_stream('stickydesign', 'params/dnastackingbig.csv')
    nndG_full = -np.loadtxt(dsb ,delimiter=',')
    nndG = nndG_full[np.arange(0,16),15-np.arange(0,16)]
    dsb.close()

    dgl = resource_stream('piperine', 'data/dnadangle.csv')
    dgldG_full = -np.loadtxt(dgl ,delimiter=',')
    dgldG_fixedC = dgldG_full[1, np.arange(4) + 4 * nt['c']]
    dgl.close()

    toeholds = [ th_set[0] for th_set in toeholds]
    toeholds_flanked = [ 'c' + th.lower() + 'c' for th in toeholds]
    ends = stickydesign.endarray(toeholds_flanked, 'TD')

    seqs_len = np.size(ends, 1)
    cols_external = np.arange(seqs_len-1)
    tops_external = tops(ends[:, cols_external])
    nndG_external = np.sum(nndG[tops_external], 1)
    dgldG_external = dgldG_fixedC[ends[:, seqs_len-2]]
    e_vec_ext = nndG_external + dgldG_external - taildG - initdG

    seqs_len = np.size(ends, 1)
    cols_internal = np.concatenate((np.arange(seqs_len-2), [seqs_len-1]))
    tops_internal = tops(ends[:, cols_internal])
    nndG_internal = np.sum(nndG[tops_internal], 1)
    e_vec_int = nndG_internal - taildG - initdG


    e_vec_all = np.concatenate( (e_vec_int, e_vec_ext))
    e_err = np.abs(e_vec_all.mean() - dGTarget)
    e_rng = e_vec_all.max() - e_vec_all.min()
    return (e_err, e_rng)

############  from here evaluation

quiet=True
clean = True
compile_params = ()

cwd =SeqSeriesName[:-2]                                                 
take_sys = '/home/lamphs/graphs/Create_Sequnce/peppercompiler_test/solov/stickydesign-peppercompiler/example/'+cwd

heuristics_inputs = (TopStrandlist, complex_names, BaseStrandlist, TopStranddict, BMlist,NotToInteract)

ssm_scores = piperine.tdm.Spurious_Weighted_Score(take_sys, domains_list, seq_dict, \
    compile_params=compile_params,includes=None, clean=clean)                                               # Seems to run normally.. Resuls ok..
ssm_names = ['WSAS', 'WSIS', 'WSAS-M', 'WSIS-M', 'Verboten', 'Spurious']

css_scores = piperine.tdm.NUPACK_Eval(seq_dict, TopStrandlist, BaseStrandlist, NotToInteract,\
                ComplexSize = 2, T = 37.0, material = 'dna', \
                clean=True, quiet=True)                                                                     # Runs normally.. Acceptable results
css_names = ['TSI avg', 'TSI max', 'TO avg', 'TO max']

if ted==True:
    ted_scores = piperine.tdm.NUPACK_Eval_bad_nucleotide(seq_dict, cmplx_dict, complex_names,\
        prefix='tube_ensemble', clean=False)                                                                  # veeeery slow ... fix it
    ted_names = ['BN% max', 'Max Defect Component','BN% avg']

bm_scores = piperine.tdm.BM_Eval(seq_dict, BMlist, toeholds_arr)
bm_names = ['WS-BM', 'Max-BM']                                                                              # Runs normally.. Acceptable results

ss_scores = piperine.tdm.SS_Eval(seq_dict, TopStranddict, T = 37.0, material = 'dna', clean=clean)          # Runs normally.. Acceptable results
ss_names = ['SSU min', 'SSU avg', 'SSTU min', 'SSTU avg']

if inter!=None:
    dGtarget = inter
else:
    dGtarget = stickydesign.enhist('S',6,energetics=stickydesign.EnergeticsBasic(temperature=37), \
    adjacents=['c','g'],alphabet='h')[2]['emedian'] 

th_scores = score_toeholds(toeholds_ma,dGtarget)                                            # Runs normally.. Acceptable results
th_names = ['dG Error', 'dG Range']

score_list = [ \
    css_scores, \
        bm_scores, \
            ss_scores, \
    #         ted_scores, \
                 ssm_scores,\
                    th_scores]

names_list = [ \
    css_names, \
        bm_names, \
            ss_names, \
                 ssm_names, \
                    th_names]
if ted==True:
    names_list.append(ted_names)
    score_list.append(ted_scores)

scores = [ elem for sub in score_list for elem in sub]
score_names  = [ elem for sub in names_list for elem in sub]

output = (scores, score_names)

scoreslist=[]

scores = [SeqSeriesName] + scores
scoreslist.append(scores)
score_names = ['Set Index'] + score_names
scores = [score_names] + scoreslist

if not os.path.exists(cwd+'_scores.csv'):
    with open(cwd+'_scores.csv', 'w') as f:
        f.write(','.join(score_names))
        f.write('\n')
        f.writelines( [','.join(map(str, l)) + '\n' for l in scoreslist])
else:
    with open(cwd+'_scores.csv', 'a') as f:
        f.write('\n')
        f.writelines( [','.join(map(str, l)) + '\n' for l in scoreslist])

