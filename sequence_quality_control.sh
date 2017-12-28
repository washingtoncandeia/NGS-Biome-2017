#!/bin/bash
# Washington C. Araujo
# IMT - UFRN
# Dia 2 -20/07/2017 - Tutorial 02-Aula_pratica_SeqQualityControl.pdf

# Regras para login no servidor:

# Interno à UFRN:
# ssh -p 4422 bif@10.7.5.38

# Externo à UFRN:
# ssh -p 4422 bif@177.20.147.141

# cd /home/bif/wash/qual

# Link para amostras:

ln -s /home/treinamento/NGS/ERR844339.fastq .
ln -s /home/treinamento/NGS/10_S5_R1_001.fastq .
ln -s /home/treinamento/NGS/polipo.fastq .


Fastscreen
########################### 1. Fast_screen ###########################


################# 1.1. Procurando contaminações #################

fastq_screen --nohits --subset 0 ERR844339.fastq --outdir .
# Opções:
# --nohits --> Extrair reads não mapeadas em nenhum dos genomas de referência.
# --nohits is equivalent to --tag --filter 0000 (zero for every genome screened).
# --subset 0 > Processar um dataset completo.
# --outdir --> Especifica um diretório no qual salvar os arquivos de saída.

# Podemos escolher os bancos de contaminantes, para tanto
# necessitamos de um arquivo .conf, copie o exemplo para o diretório corrente:

cp /home/treinamento/NGS/fastq_screen.conf .

# Apagar resultado anterior
# Apaga arquivos criados por fastq screen ERR844339_no_hits.fastq e rm ERR844339_screen.png
rm ERR844339_*

################## 1.2. Usando fastq_screen.conf ################## 

fastq_screen --nohits --subset 0 --conf fastq_screen.conf ERR844339.fastq --outdir .
## Opções:
# --nohits --> Escreve para um arquivo as sequências que não mapearam em nenhum dos genomas especificados. 
# --subset 0 > Processar um dataset completo.
# --conf --> Busca o arquivo de configuração fastq_screen.conf
# ERR844339.fastq --> amostra analisada
# --outdir --> Especifica um diretório no qual salvar os arquivos de saída.
# Saída: ERR844339_no_hits.fastq  ERR844339_screen.png  ERR844339_screen.txt

# 1.2.1. Mudar para área de trabalho web e criar o diretório:
mkdir /home/bif/public_html/wash

# 1.2.2. Copiar tudo para lá
cp ERR844339_* /home/bif/public_html/wash

### Entrar lá no servidor web e visualizar -> http://177.20.147.141/~bif/wash/

########################### 2. FastqUtils #############################
# cd /home/bif/wash/qual

# 2.1. Estatísticas básicas utilizando fasqUtils

fastqutils stats ERR844339.fastq | more

fastqutils stats ERR844339_no_hits.fastq | more

# 2.2. Escrever a estatística num arquivo de texto
fasqutils stats ERR844339.fastq > fastqutils.stat.txt
fastqutils stats ERR844339_no_hits.fastq > fastqutils_no_hits.stat.txt


########################### 3. FastQC #############################
# cd /home/bif/wash/qual

# 3.1. Obtenção de estatística
fastqc ERR844339.fastq -o .

# 3.2. Copiar o .html gerado para web
cp ERR844339_fastqc.* /home/bif/public_html/wash/qual

# 3.3. Repitir para os outros arquivos .fastqc
fastqc 10_S5_R1_001.fastq -o .
fastqc polipo.fastq -o .
cp *_fastqc.* /home/bif/public_html/wash/qual

########################### 4. SAMStat #############################
# cd /home/bif/wash/qual

# Utilizado para avaliar sequências já alinhadas. 
# O fastQC é utilizado antes do alinhamento.

# Necessita-se de um arquivo .bam
# 4.1. Criar links
ln -s /home/treinamento/NGS/ERR844339.bam .
ln -s /home/treinamento/NGS/polipo.bam . 

# 4.2. Rodar o SAMstat
samstat ERR844339.bam
samstat polipo.bam

# 4.3. Copiar os resultados para área web
cp *.samstat.* /home/bif/public_html/wash/

########################### 5. TRIM Galore! #############################
# cd /home/bif/wash/qual

# Após o uso dos programas FasqQC e Samstat talvez seja necessário
# aplicar algum processo de limpeza. 
# Isso irar melhorar o resultados finais e diminuir a taxa de erro de analises posteriores.

# 5.1. DynamicTrim - Trimar as pontas das sequências por qualidade de base
 
DynamicTrim.pl -h 20 -i ERR844339.fastq 

# 5.2. Cutadapt - Remover adaptadores ou sequencias contaminantes conhecidas.

cutadapt -a TGGAATTCTCGG 10_S5_R1_001.fastq > 10_S5_R1_001.ct.fastq
    
# 5.3. O programa trim_galore, quando temos um sequenciamento illumina
# e não sabemos os adaptadores usados no processo de sequenciamento

trim_galore 10_S5_R1_001.fastq

# Vale ressaltar que após cada processo de limpeza devemos refazer a
# analise de qualidade novamente, e assim verificar a sua efetiva melhora.













