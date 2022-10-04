#!/bin/bash
## Date: 4 October, 2022
## Author: Drew Neavin
## Reason: RNA velocity on iPSC village samples for pseudotime following the Allevin tutorial https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/



if [[ $pool == "DRENEA_1" ]] || [[ $pool == "DRENEA_2" ]] || [[ $pool == "DRENEA_3" ]] || [[ $pool == "DRENEA_4" ]] || [[ $pool == "DRENEA_5" ]] || [[ $pool == "DRENEA_6" ]] 
then
	velocyto run10x --samtools-threads $T -vv $TENxDIR/$pool $GTF
else
	velocyto run10x --samtools-threads $T -vv $TENxDIR/$pool/$pool $GTF
fi
