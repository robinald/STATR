##C. Trim low quality reads and adapter sequences
echo "Starting STATR"

printf "Started read trimming...%s\n" "$(date)"
mkdir -p ./2.Trimmed_reads/

java -jar ./Trimmomatic/trimmomatic-0.39.jar SE -phred33 ./1.Raw_reads/Eco_Ctrl1.fastq.gz ./2.Trimmed_reads/Eco_Ctrl1_trimmed.fastq.gz ILLUMINACLIP:./Trimmomatic/adapters/RiboSeq_adapter_as.fa:2:30:6 SLIDINGWINDOW:4:15 MINLEN:12 -threads 8
java -jar ./Trimmomatic/trimmomatic-0.39.jar SE -phred33 ./1.Raw_reads/Eco_Ctrl2.fastq.gz ./2.Trimmed_reads/Eco_Ctrl2_trimmed.fastq.gz ILLUMINACLIP:./Trimmomatic/adapters/RiboSeq_adapter_as.fa:2:30:6 SLIDINGWINDOW:4:15 MINLEN:12 -threads 8

java -jar ./Trimmomatic/trimmomatic-0.39.jar SE -phred33 ./1.Raw_reads/Eco_Exp1.fastq.gz ./2.Trimmed_reads/Eco_Exp1_trimmed.fastq.gz ILLUMINACLIP:./Trimmomatic/adapters/RiboSeq_adapter_as.fa:2:30:6 SLIDINGWINDOW:4:15 MINLEN:12 -threads 8
java -jar ./Trimmomatic/trimmomatic-0.39.jar SE -phred33 ./1.Raw_reads/Eco_Exp2.fastq.gz ./2.Trimmed_reads/Eco_Exp2_trimmed.fastq.gz ILLUMINACLIP:./Trimmomatic/adapters/RiboSeq_adapter_as.fa:2:30:6 SLIDINGWINDOW:4:15 MINLEN:12 -threads 8



##D. Align reads to the genome
mkdir -p ./3.Aligned/
printf "Started building a reference database...%s\n" "$(date)"
bowtie2-build -f ./ReferenceDB/NC_000913.3.fasta ./ReferenceDB/NC_000913.3

printf "Started read mapping...%s\n" "$(date)"
bowtie2 -q -U ./2.Trimmed_reads/Eco_Ctrl1_trimmed.fastq.gz -x ./ReferenceDB/NC_000913.3 -S ./3.Aligned/Eco_Ctrl1_mapping.sam --local
bowtie2 -q -U ./2.Trimmed_reads/Eco_Ctrl2_trimmed.fastq.gz -x ./ReferenceDB/NC_000913.3 -S ./3.Aligned/Eco_Ctrl2_mapping.sam --local

bowtie2 -q -U ./2.Trimmed_reads/Eco_Exp1_trimmed.fastq.gz -x ./ReferenceDB/NC_000913.3 -S ./3.Aligned/Eco_Exp1_mapping.sam --local
bowtie2 -q -U ./2.Trimmed_reads/Eco_Exp2_trimmed.fastq.gz -x ./ReferenceDB/NC_000913.3 -S ./3.Aligned/Eco_Exp2_mapping.sam --local



##E. Decompile alignment file
mkdir -p ./4.Decompiled/
printf "Started decompiling the mapping data...%s\n" "$(date)"
samtools view -bS ./3.Aligned/Eco_Ctrl1_mapping.sam > ./4.Decompiled/Eco_Ctrl1_mapping.bam
samtools view -bS ./3.Aligned/Eco_Ctrl2_mapping.sam > ./4.Decompiled/Eco_Ctrl2_mapping.bam
samtools view -bS ./3.Aligned/Eco_Exp1_mapping.sam > ./4.Decompiled/Eco_Exp1_mapping.bam
samtools view -bS ./3.Aligned/Eco_Exp2_mapping.sam > ./4.Decompiled/Eco_Exp2_mapping.bam

samtools sort ./4.Decompiled/Eco_Ctrl1_mapping.bam ./4.Decompiled/Eco_Ctrl1_mapping.sorted
samtools sort ./4.Decompiled/Eco_Ctrl2_mapping.bam ./4.Decompiled/Eco_Ctrl2_mapping.sorted
samtools sort ./4.Decompiled/Eco_Exp1_mapping.bam ./4.Decompiled/Eco_Exp1_mapping.sorted
samtools sort ./4.Decompiled/Eco_Exp2_mapping.bam ./4.Decompiled/Eco_Exp2_mapping.sorted

bedtools bamtobed -i ./4.Decompiled/Eco_Ctrl1_mapping.sorted.bam > ./4.Decompiled/Eco_Ctrl1.bed
bedtools bamtobed -i ./4.Decompiled/Eco_Ctrl2_mapping.sorted.bam > ./4.Decompiled/Eco_Ctrl2.bed
bedtools bamtobed -i ./4.Decompiled/Eco_Exp1_mapping.sorted.bam > ./4.Decompiled/Eco_Exp1.bed
bedtools bamtobed -i ./4.Decompiled/Eco_Exp2_mapping.sorted.bam > ./4.Decompiled/Eco_Exp2.bed



##F. Validate experiment and visualize genome-wide ribosome profile
printf "Started meta analysis...%s\n" "$(date)"
mkdir -p ./5.Meta_analysis/
python3 ./Python_scripts/ParseGenomeAnnotation.py -i ./ReferenceDB/NC_000913.3.gff3 -o ./ReferenceDB/NC_000913.3_CDS.gff
python3 ./Python_scripts/CheckPeriodicity.py -i ./4.Decompiled/Eco_Ctrl1.bed -a ./ReferenceDB/NC_000913.3_CDS.gff -o ./5.Meta_analysis/Meta_positional_density_Eco_Ctrl1.txt
python3 ./Python_scripts/CheckPeriodicity.py -i ./4.Decompiled/Eco_Ctrl2.bed -a ./ReferenceDB/NC_000913.3_CDS.gff -o ./5.Meta_analysis/Meta_positional_density_Eco_Ctrl2.txt
python3 ./Python_scripts/CheckPeriodicity.py -i ./4.Decompiled/Eco_Exp1.bed -a ./ReferenceDB/NC_000913.3_CDS.gff -o ./5.Meta_analysis/Meta_positional_density_Eco_Exp1.txt
python3 ./Python_scripts/CheckPeriodicity.py -i ./4.Decompiled/Eco_Exp2.bed -a ./ReferenceDB/NC_000913.3_CDS.gff -o ./5.Meta_analysis/Meta_positional_density_Eco_Exp2.txt