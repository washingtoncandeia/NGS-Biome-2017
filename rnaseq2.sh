#!/bin/bash
# Washington C. Araujo
# IMT - UFRN
# Dia 3 -21/07/2017 - Tutorial RNA-seq 2 miRNA

# Regras para login no servidor:

# Interno à UFRN:
# ssh -p 4422 bif@10.7.5.38

# Externo à UFRN:
# ssh -p 4422 bif@177.20.147.141

########### Tutorial RNA-seq 2 ###############
# /home/bif/wash/rna2

# Criar links simbólicos
# Abaixo, os sequenciamentos das amostras:
ln -s /home/treinamento/NGS/RNAseq/sample_data/SRR326279_R1.fastq .
ln -s /home/treinamento/NGS/RNAseq/sample_data/SRR326280_R1.fastq .
ln -s /home/treinamento/NGS/RNAseq/sample_data/SRR326281_R1.fastq .
ln -s /home/treinamento/NGS/RNAseq/sample_data/SRR326282_R1.fastq .

#################### Arquivos de referência #############################
# Arquivo fasta com o genoma de referência;
# Arquivos de index do genoma de referência;
# Arquivo fasta com os miRNAs referência para a espécie (utilizaremos o miRBase);
# Arquivo fasta com os miRNAs maduros para a espécie (utilizaremos o miRBase);
# Arquivo fasta de predição dos loops das sequências dos miRNAs para a espécie, os hairpins (utilizaremos o miRBase);

# Todos estes arquivos acima, estão aqui:
ln -s /home/treinamento/NGS/RNAseq/small_ref/ .

#################### Preparando o Dado Inicial #############################
# Antes de executar o miRDeep2, os dados devem ser pré-processados para remover adaptadores. 
# Isso pode ser feito usando o cutadapt, para cada amostra .fastq
# Cutadapt - permite encontrar (matches) entre as reads e os adaptadores, retirando-os conforme opções.

cutadapt -b AATCTCGTATGCCGTCTTCTGCTTGC -O 3 -m 17 \ 
         -f fastq SRR326279_R1.fastq > SRR326279_R1.ct.fastq

cutadapt -b AATCTCGTATGCCGTCTTCTGCTTGC -O 3 -m 17 \ 
         -f fastq SRR326280_R1.fastq > SRR326280_R1.ct.fastq 

cutadapt -b AATCTCGTATGCCGTCTTCTGCTTGC -O 3 -m 17 \ 
         -f fastq SRR326281_R1.fastq > SRR326281_R1.ct.fastq 

cutadapt -b AATCTCGTATGCCGTCTTCTGCTTGC -O 3 -m 17 \ 
         -f fastq SRR326282_R1.fastq > SRR326282_R1.ct.fastq 

# Opções:
# -b --> os tipos de adaptadores a serem trimados. Neste caso, indica-se que 5' ou 3' devem ser removidos.
# -O --> --overlap N, usado para reduzir o número de trimagens falsas de bases. 
# -m --> --minimum-length N, joga fora da análise reads menores que este valor N
# -f --> --format, opção de informação sobre que formato de arquivo está sendo usado para input


#################### Mapeamento #############################

# O pré-processamento para uso de miRDeep2 é feito com mapper.pl
# mapper.pl processa as reads e as mapeia contra o genoma.

# Usaremos o script mapper.pl para processar as leituras 
# e mapear-las contra o genoma de referência:

# 1. SRR326279_R1.ct.fastq
mapper.pl SRR326279_R1.ct.fastq \ 
            -e \ 
            -p small_ref/hg19_chr1 \ 
            -s SRR326279.pr.fa \ 
            -t SRR326279.mr.arf \ 
            -h -m -i -j  

# 2. SRR326280_R1.ct.fastq
mapper.pl SRR326280_R1.ct.fastq \ 
            -e \ 
            -p small_ref/hg19_chr1 \ 
            -s SRR326280.pr.fa \ 
            -t SRR326280.mr.arf \ 
            -h -m -i -j

# 3. SRR326281_R1.ct.fastq
mapper.pl SRR326281_R1.ct.fastq \ 
            -e \ 
            -p small_ref/hg19_chr1 \ 
            -s SRR326281.pr.fa \ 
            -t SRR326281.mr.arf \ 
            -h -m -i -j

# 4. SRR326282_R1.ct.fastq
mapper.pl SRR326282_R1.ct.fastq
            -e 
            -p small_ref/hg19_chr1 \ 
            -s SRR326282.pr.fa \ 
            -t SRR326282.mr.arf \ 
            -h -m -i -j

## Opções:
# -e --> input file está no formato .fastq
# -p --> uso: -p genome; mapeia contra o genoma (deve estar indexado por bowtie-build). 
# -s --> imprime as reads processadas para este arquivo de saída (SRR326279.pr.fa, p.ex.)
# -t --> imprime os mapeamentos de reads para este arquivo de saída (SRR326279.mr.arf, p.ex.)
# -h --> parsing para formato .fasta
# -m --> collapse reads
# -i --> converte RNA a alfabeto de DNA (para mapear contra o genoma)
# -j --> remove todas as entradas que contêm uma sequência que contenha letras


#################### miRDeep2 #############################
# miRDeep2.pl é utilizado para correr todos os scripts necessários do pacote miRDeep2.
# Assim, premite miRDeep2 fazer profunda análise do sequenciamento para detecção de microRNA.

##--------------------------------------------- Input: -------------------------------------------------------------------------------------
# 1. Arquivo .fasta com as reads processadas por mapper.pl ( função -s): SRR326279.pr.fa, SRR326280.pr.fa, SRR326281.pr.fa,SRR326282.pr.fa; 
# 2. Arquivo .fasta com o genoma correspondente: small_ref/hg19_chr1.fa
# 3. Arquivo .fasta das reads mapeadas anteriormente por mapper.pl (função -t): SRR326279.mr.arf, SRR326280.mr.arf, SRR326281.mr.arf, SRR326282.mr.arf
# 4. Arquivos (opcional) .fasta de miRNAs conhecidos de espécies analisadas: small_ref/mature.hsa.dna.fa e small_ref/hairpin.hsa.dna.fa, 
# e arquivos (opcional) .fasta de miRNA de espécies relacionadas. 

# Observação: none é utilizado para designar animal sem espécie relacionada na miRBase
##-------------------------------------------------------------------------------------------------------------------------------------------
## Identificação sequenciamento:

# 1. SRR326279
miRDeep2.pl SRR326279.pr.fa small_ref/hg19_chr1.fa \ 
            SRR326279.mr.arf small_ref/mature.hsa.dna.fa \ 
            none             small_ref/hairpin.hsa.dna.fa \ 
            -t Human         2> report.log

# 2. SRR326280
miRDeep2.pl SRR326280.pr.fa small_ref/hg19_chr1.fa \ 
            SRR326280.mr.arf small_ref/mature.hsa.dna.fa \ 
            none             small_ref/hairpin.hsa.dna.fa \ 
            -t Human         2> report.log

# 3. SRR326281
miRDeep2.pl SRR326281.pr.fa small_ref/hg19_chr1.fa \ 
            SRR326281.mr.arf small_ref/mature.hsa.dna.fa \ 
            none             small_ref/hairpin.hsa.dna.fa \ 
            -t Human         2> report.log

# 4. SRR326282
miRDeep2.pl SRR326282.pr.fa small_ref/hg19_chr1.fa \ 
            SRR326282.mr.arf small_ref/mature.hsa.dna.fa \ 
            none             small_ref/hairpin.hsa.dna.fa \ 
            -t Human         2> report.log

## Opções:
# -t --> espécie sendo analisada: Human. Usado para fazer link para browser UCSC apropriado.
# 2> --> Usado para reportar (pipe) todo progresso das saídas para um arquivo report.log

