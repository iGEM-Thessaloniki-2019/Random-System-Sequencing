import re
import stickydesign
from stickydesign.plots import *
import matplotlib.pyplot as plt
import argparse 

parser = argparse.ArgumentParser(description='The script creates a set of toeholds.')
parser.add_argument('-s','--sets', default=1,type=int, help='The number of sets to be created')
parser.add_argument('-n','--number', default=5, help='length of a set')
parser.add_argument('-l','--length', default=6, help='Toehold length')
parser.add_argument('-f','--sys_file', help='Sys file')
parser.add_argument('-v','--verbose', action = 'store_true', help='Print sets and toehold sequence')
parser.add_argument('-g','--dGtarget', default = 0, help='Wanted dGTarget')

args = parser.parse_args()

sets=args.sets
number=args.number
endlength=args.length
sys_file=args.sys_file
interaction=args.dGtarget

for t in range(sets):
    # easyends finds the toehold set
    toeholds=stickydesign.easyends('S',endlength ,number=number, maxspurious=0.4,interaction=interaction, \
 tries=1,energetics=stickydesign.EnergeticsBasic(temperature=37), adjs=['c','g'],alphabet='h')
    # The EnergeticsBasic() model is the most suitable for toeholds, even if it contrasts the README.md Use section.
    #print("nt = {'a': 0, 'c': 1, 'g': 2, 't': 3}")

    ###Convert numbers to DNA as toeholds is an array of numbers insted nucletoids
    toes_list=[] # an array containing toeholds in 'AATTCC' format
    toes_seq=[] ## an array containing toeholds in 'A','A','T','T','C','C' format
    for i in range(len(toeholds)):
        toes=[]
        b=[]
        for j in range(len(toeholds[i])):
            if toeholds[i][j]==0:
                toes.append('A')
            elif toeholds[i][j]==1:
                toes.append('C')
            elif toeholds[i][j]==2:
                toes.append('G')
            elif toeholds[i][j]==3:
                toes.append('T')
        if toes!= []:		
            toes_seq.append(toes)
        toes_list.append("".join(toes))
    if args.verbose:
        print(toes_list)

    file1=sys_file
    fixed_name=file1.split(".")[0]
    fixed_file=fixed_name+'_'+str(t)+".fixed"

    sys = open(file1, 'r')

    sys_lines = sys.readlines()

    for line in sys_lines:
        sys_wordList=line.split()
        for i in range(len(sys_wordList)):
            a = sys_wordList[i]
            if a == 'component':        
                name1 = sys_wordList[i+1]
            if a == 'import':
                comp_file = sys_wordList[i+1] + '.comp'

    comp=open(comp_file,'r')

    comp_lines = comp.readlines()

    toe_names=[]
    for line in comp_lines:
        comp_wordList=line.split()
        for i in range(len(comp_wordList)):
            k=[]
            a = comp_wordList[i]
            for j in a:
                if j=='<':
                    k.append(j)
                    continue
                if j =='t':
                    k.append(j)
                    continue
                if j=='>':
                    k.append(j)
            if len(k) == 3:
                toe_names.append(comp_wordList[i-2])

    f=open(fixed_file,'w')

    for toe in range(len(toe_names)):
        line_template="sequence {}-{} = {}\n".format(name1,toe_names[toe],toes_list[toe])
        f.write(line_template)
        if args.verbose:
            print(line_template)

    f.close()
