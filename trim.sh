#!/bin/bash

########################### Trim #############################
# cd /home/bif/wash/qual

# Após o uso dos programas FasqQC e Samstat talvez seja necessário
# aplicar algum processo de limpeza, isso irar melhorar o resultados finais e
# diminuir a taxa de erro de analises posteriores.

# 1. Trimar as pontas das sequências por qualidade de base
DynamicTrim.pl -h 20 -i ERR844339.fastq 

# 2. Podemos usar o cutadapt sequencias contaminantes conhecidas.
cutadapt -a TGGAATTCTCGG 10_S5_R1_001.fastq > 10_S5_R1_001.ct.fastq
    
# 3. O programa trim_galore, quando temos um sequenciamento illumina
# e não sabemos os adaptadores usados no processo de sequenciamento
trim_galore 10_S5_R1_001.fastq
