#!/bin/bash
# Washington C. Araujo
# IMT - UFRN
# Dia 2 -20/07/2017 - Tutorial 03-Aula_pratica_variantes.pdf

# Regras para login no servidor:

# Interno à UFRN:
# ssh -p 4422 bif@10.7.5.38

# Externo à UFRN:
# ssh -p 4422 bif@177.20.147.141

######################### Tutorial Chamada de Variantes ######################
# cd /home/bif/wash/bwa
ln -s /home/treinamento/NGS/NC_012967.1.fa .
ln -s /home/treinamento/NGS/SRR5714077_1_s.1.fastq .
ln -s /home/treinamento/NGS/SRR5714077_2_s.1.fastq .

# Como é o arquivo .fasta
less -S NC_012967.1.fa

# Como é o arquivo .fastq
less -S SRR5714077_1_s.1.fastq


########################### 1. BWA - Mapeamento ###########################
# 1.1. Criar os indexes do genoma de referência
bwa index -a is NC_012967.1.fa

# -a --> algorítmo para construção de index BWT. Possui opções disponíveis:
# is --> algorítmo de tempo linear para construção de uma matriz/array de sufixos. Para genomas menores, como de bacterias.
# bwtsw --> para genomas maiores, como o humano.

# Usar samtools para gerar index de arquivo fasta
samtools faidx NC_012967.1.fa  # Saída: NC_012967.1.fa.fai

# Opções
# faidx --> Sequência index de referência no formato fasta ou extrai subsequência
# da sequência referência indexada. 
# Se nenhuma região for especificada, faidx indexará  o arquivo e crirará <ref.fasta>.fai no disco.

# 1.2. Agora vamos rodar BWA utilizando os arquivos gerados até aqui:
bwa bwasw -t 4 NC_012967.1.fa \ 
    SRR5714077_1_s.1.fastq \ 
    SRR5714077_2_s.1.fastq -f bwa.sam

## Opções:
# bwasw --> invoca algorítmo de alinhamento para sequências mais longas, BWA-SW
# -t --> número de threads usados
# mem --> compartilha características com bwasw, como suporte para reads mais longas e split alignments.
# Porém, mem é recomendado para queries de alta qualidade sendo mais rápida e acurada.

########################### 2. Samtools ###########################
# Converter arquivo .sam em .bam 
# Utilizar Samtools para manipulação e extrair estatísticas básicas

# 2.1. Converter SAM to BAM
samtools view -b -S bwa.sam -o bwa.bam
# Opções:
# -b --> Output in the BAM format. 
# -S --> Indica que o arquivo input é SAM.
# -o --> Especifica o nome do arquivo de saida.

# 2.2. Visualizar
samtools view bwa.bam | less -S

# 2.3. Visualizar apenas sequências não mapeadas
samtools view -f 4 bwa.bam | less -S
# -f --> extrai apenas aquelas reads que batem com uma FLAG sam.
# Neste caso, -f 4 filtra apens aquelas reads com FLAG 4.
# -f 4 --> reads que falharam em mapear no genoma de referência.

# 2.4. Visualizando apenas sequências mapeadas
samtools view -F 4 bwa.bam | less -S
# -F --> extrai apenas aquelas reads que batem com uma FLAG sam.
# Neste caso, -F 4 filtra apens aquelas reads com FLAG 4.
# -F 4 --> extrai reads que mapearam no genoma de referência.

# 2.5a. Quantificando as sequências não mapeadas
samtools view -c -f 4 bwa.bam
# -f --> extrai apenas aquelas reads que batem com uma FLAG sam.
# -c --> A saída é o valor de número de reads que bateram com o critério
# O critério é -f 4 = reads que não mapearam no genoma de referência.

# Não mapeadas: 103.224; 

# 2.5b. Sequências Mapeadas
samtools view -c -F 4 bwa.bam
# -F 4 --> extrai reads que mapearam no genoma de referência.
# -c --> A saída é o valor de número de reads que bateram com o critério
# O critério é -F 4 = reads que mapearam no genoma de referência.

# Mapeadas: 341.432

# 2.6. Quantificando as sequências com qualidade MAPQ superior a 42
samtools view -c -q 42 bwa.bam
# -q --> filtro mínimo de qualidade
# -q 42 = escore de qualidade iniciando a partir de 42
# -c -q 42 = conta o número de reads com escore de qualidade acima de 42

# 325.848

#  Qualidade MAPQ superior a 30
samtools view -c -q 30 bwa.bam
# -q --> filtro mínimo de qualidade
# -q 42 = escore de qualidade iniciando a partir de 30
# -c -q 42 = conta o número de reads com escore de qualidade acima de 30

# 327.294

### Site biostars
# samtools view -F 0x04 -c filename.bam
# == samtools view -c -F 4 bwa.bam
# Count the number of alignments (reads mapping to multiple locations counted multiple times)
# 341.430

#------------- Ordenar BAM ---------------------
# 2.7. Ordenar sequências .bam
samtools sort bwa.bam -o bwa.sort.bam

# 2.8. Remover amplificação de PCR
samtools rmdup bwa.sort.bam bwa.rmd.bam
# Uso de rmdup:
# samtools rmdup [-sS] <input.srt.bam> <out.bam> 
# Remove potenciais duplicatas de PCR

# 2.9. Gerar arquivo mpileup
samtools mpileup -f NC_012967.1.fa bwa.rmd.bam > ecoli.mpileup

# Uso mpileup:
# Gera VCF, BCF ou pileup para um ou múltiplos arquivos BAM.
# Registros de alinhamentos são agrupados por identificadores de amostra (SM), 
# nas linhas de cabeçalho @RG.
# Se os identificadores de amostras estiverem ausentes, cada arquivo de entrada é considerado uma amostra. 

# 2.10. Ver o arquivo mpileup
less -S ecoli.mpileup

########################### 3. Varscan: Fazer Chamada de Variantes ###########################

# 3.1. Chamada de variantes com VarScan
varscan mpileup2snp ecoli.mpileup --output-vcf --strand-filter 0 > ecoli.vcf
# Opções:
# mpileup2snp --> chama SNPs de um arquivo mpileup
# --output-vcf --> opção de mpileup2snp; gera um arquivo de saída .VCF
# --strand-filter --> opção de mpileup2snp; ignora variantes com >90% de suporte sobre uma fita


########################### 4. snpEff ###########################
# O arquivo .vcf contém todas as variantes, porém não estão anotadas todas as informações relevantes
# para extrair significado biológico de cada variante. 
# Para isso, executar o processo de anotação de variantes.
snpEff eff Escherichia_coli_B_REL606_uid58803 ecoli.vcf > ecoli.eff.vcf

## Após isso, teremos o arquivo que contém todas as variantes e as informações
## relevantes para extrair o significado biológico de cada variante











