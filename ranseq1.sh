#!/bin/bash
# Washington C. Araujo
# IMT - UFRN
# Dia 3 -21/07/2017 - Tutorial RNA-seq 1

# Regras para login no servidor:

# Interno à UFRN:
# ssh -p 4422 bif@10.7.5.38

# Externo à UFRN:
# ssh -p 4422 bif@177.20.147.141

########### Tutorial RNA-seq 1 ###############
# cd /home/bif/wash/rna1

ln -s /home/treinamento/NGS/RNAseq/adrenal_1.fastq .
ln -s /home/treinamento/NGS/RNAseq/adrenal_2.fastq .
ln -s /home/treinamento/NGS/RNAseq/brain_1.fastq .
ln -s /home/treinamento/NGS/RNAseq/brain_2.fastq .

# Como é o arquivo
less -S adrenal_1.fastq
less -S adrenal_2.fastq
less -S brain_1.fastq
less -S brain_2.fastq

# Contar sequências totais nos arquivos
wc adrenal_1.fastq
wc adrenal_2.fastq
wc brain_1.fastq
wc brain_2.fastq

# Com grep
grep '@ERR' adrenal_1.fastq | wc 
grep '@ERR' adrenal_2.fastq | wc 
grep '@ERR' brain_1.fastq | wc 
grep '@ERR' brain_2.fastq | wc 

# Filtragem com fasQC
fastqc adrenal_1.fastq -o .
fastqc adrenal_2.fastq -o .
fastqc brain_1.fastq -o .
fastqc brain_2.fastq -o .

# Para visualizar, copiar para web
cp adrenal_1.fastq* /home/bif/public_html/wash/
cp adrenal_2.fastq* /home/bif/public_html/wash/
cp brain_1.fastq* /home/bif/public_html/wash/
cp brain_2.fastq* home/bif/public_html/wash/

# Filtragem com Trimmomatic

# Inicialmente link dos adapters:
ln -s /root/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa .

# Comando Trimmomatic
trimmomatic PE -threads 1 \ 
        adrenal_1.fastq adrenal_2.fastq \ 
        adrenal_1_paired.fastq.gz adrenal_1_unpaired.fastq.gz \ 
        ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \ 
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20 \ 
# Descompactar os arquivos considerados filtrados
# e que mantiveram os pares após a filtragem
gzip -d *_paired.fastq.gz 

##################################### 1. Mapeamento ################################################
# Para criar mapeamento, criar links simbólicos
# do genoma de referência e do arquivo .gtf

# Genoma de referência:
ln -s /home/databases/hg19/ . 
# Arquivo .gtf :
ln -s /home/treinamento/NGS/RNAseq/gene19_annotation.gtf .

################### 1.1. TopHat2 #######################
## Tophat2
# Agora vamos rodar TopHat2 utilizando os arquivos gerados até aqui:

# 1. Amostras de RNA-seq adrenal
tophat2 -p 1 -G gene19_annotation.gtf \ # anotação .gtf
        -o thout_adrenal \              # diretório de saída para thout_adrenal
        hg19/hg19 \                     # genoma de referência
        adrenal_1_paired.fastq \        # mate1.fastq
        adrenal_2_paired.fastq          # mate2.fastq

# 2. Amostras RNA-seq cérebro
tophat2 -p 1 -G gene19_annotation.gtf \ # anotação .gtf
        -o thout_brain \                # diretório de saída para thout_adrenal
        hg19/hg19 \                     # genoma de referência
        brain_1_paired.fastq \          # mate1.fastq
        brain_2_paired.fastq            # mate2.fastq

# O arquivo thout_adrenal é o input para cufflinks

# Opções principais utilizadas:
# -G --> indica o arquivo de anotação .gtf
# -o --> indica o nome do diretório para o qual tophat envia a saída (thout_adrenal)
# hg19 --> genoma index base


################################### 2. Montagem #########################################################

################### 2.1. Montagem com Cufflinks #######################
# Como usar:
# cufflinks [options] <aligned_reads.(sam/bam)>

cufflinks -p 4 -o clout_adrenal thout_adrenal/accepted_hits.bam
# -p ---> threads utilizadas
# -o ---> output. Gera diretório clout_adrenal
# thout_adrenal/accepted_hits.bam --> utiliza o arquivo accepted_hits.bam gerado por Tophat2

cufflinks -p 4 -o clout_brain thout_brain/accepted_hits.bam
# -p ---> threads utilizadas
# -o ---> output. Gera diretório clout_brain
# thout_adrenal/accepted_hits.bam --> utiliza o arquivo accepted_hits.bam gerado por Tophat2

###################### Expressão Diferencial #####################
# Para identificar quais genes estão diferencialmente expressos.

# 1. Criar arquivo de anotação de referência para ser utilizado por cuffmerge:  
# vim assemblies.txt: dentro de vim, digitar:
./clout_adrenal/transcripts.gtf
./clout_brain/transcripts.gtf
# Salvar e sair (usar tecla ESC; digitar ":wq" - não usar aspas - e teclar ENTER).

################### 2.2 Cuffmerge #######################
# Cuffmerge é um script que já vem dentro de Cufflinks
# Utilizado para juntar diversas montagens feitas por Cufflinks

# Uso:
# cuffmerge [options]* <assembly_GTF_list.txt>

cuffmerge -g gene19_annotation.gtf -s hg19/hg19.fa -p 1 assemblies.txt
# gene19_annotation.gtf --> arquivo de entrada para Cuffmerge
# Opções:
# -s --> aponta para o DNA genômico de referência (arquivo de referência)
# -g --> arquivo de anotação .gtf 
# -p --> número de threads utilizadas


################### 2.3 Cuffdiff #######################
# Cuffdiff compara os níveis de expressão de genes e transcritos
# em RNA-seq. 
# Informa não apenas quais genes estão up ou downregulated entre duas condições,
# mas também, que genes sofrem splicing diferencial, 
# ou estão sob outros tipos de regulação de isoformas.

# Uso:
# cuffdiff [options]* <transcripts.gtf> \
# <sample1_replicate1.sam[,…,sample1_replicateM.sam]> \
# <sample2_replicate1.sam[,…,sample2_replicateM.sam]> … \
# [sampleN.sam_replicate1.sam[,…,sample2_replicateM.sam]]

cuffdiff -o diff_out \                          # Define o nome do diretório no qual Cuffdiff irá escrever toda a sua saída.
         -b hg19/hg19.fa \                      
         -p 1 \ 
         -L A,B
         -u merged_asm/merged.gtf \ 
         ./thout_adrenal/accepted_hits.bam \ 
         ./thout_brain/accepted_hits.bam

# Opções:
# -b --> Fornecendo o arquivo multifasta, suas leituras foram mapeadas através desta opção, 
#       instrui-lo a executar nosso algoritmo de detecção e correção de polarização que 
#       pode melhorar significativamente a precisão das estimativas de abundância de transcrição. 

# -o --> Define o nome do diretório no qual Cuffdiff irá escrever toda a sua saída.
# -p --> Uso de threads
# -L --> Especifica a label (prefixo, nome) para cada amostra, que será incluído em vários arquivos de saída produzidos por Cuffdiff.
# -u --> Diz para fazer um procedimento de estimativa inicial para mapear o mapeamento com maior precisão para múltiplos locais no genoma.
# Saída do cuffmerge: 
# ./thout_adrenal/accepted_hits.bam 
# ./thout_brain/accepted_hits.bam


################################## Observar Resultados #############################

# Dentro da pasta diff_out encontramos vários resultados interessantes.
# Com ls podemos checar alguns desses arquivos gerados.
ls diff_out

# Vamos manter o foco no arquivo gene_exp.diff 
more diff_out/gene_exp.diff

# Agora vamos selecionar apenas aqueles genes
# diferencialmente expressos pelos testes estatísticos do cuffdiff
grep 'yes$' diff_out/gene_exp.diff

