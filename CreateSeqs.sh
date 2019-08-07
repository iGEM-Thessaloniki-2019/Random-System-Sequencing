#!/bin/bash

### To specify the number of sequence-sets to be produced, be carefull to change both  -s in line  5 and i<= in line 7
# Carefull, -s starts counting from 1 while i starts from 0. -s=i+1
python write_fixed.py -s 2 -f Solov1_0.sys

for (( i=0; i<=1; i++ ))
do
	pepper-compiler --fixed=Solov1_0_$i.fixed --output Solov1_0_$i.pil --save Solov1_0_$i.save Solov1_0.sys
	pepper-design-spurious Solov1_0_$i.pil
	pepper-finish Solov1_0_$i.mfe

	python EvaluationInputs.py -bn Solov1_0_$i -sig solo-IA solo-IB solo-R solo-TCC
done