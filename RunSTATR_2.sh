##F. Validate experiment and visualize genome-wide ribosome profile
## The shell script is splitted into two scripts for user to select 5' or 3' end to be used.
python3 ./Python_scripts/GenerateProfile.py -i ./4.Decompiled/Eco_Ctrl1.bed -e 5 -a ./ReferenceDB/NC_000913.3_CDS.gff -n 1000000 -o ./5.Meta_analysis/Eco_Ctrl1_Profile.gff
python3 ./Python_scripts/GenerateProfile.py -i ./4.Decompiled/Eco_Ctrl2.bed -e 5 -a ./ReferenceDB/NC_000913.3_CDS.gff -n 1000000 -o ./5.Meta_analysis/Eco_Ctrl2_Profile.gff
python3 ./Python_scripts/GenerateProfile.py -i ./4.Decompiled/Eco_Exp1.bed -e 5 -a ./ReferenceDB/NC_000913.3_CDS.gff -n 1000000 -o ./5.Meta_analysis/Eco_Exp1_Profile.gff
python3 ./Python_scripts/GenerateProfile.py -i ./4.Decompiled/Eco_Exp2.bed -e 5 -a ./ReferenceDB/NC_000913.3_CDS.gff -n 1000000 -o ./5.Meta_analysis/Eco_Exp2_Profile.gff



##G. Differential expression
mkdir -p ./6.Differential_expression/
bedtools coverage -a ./ReferenceDB/NC_000913.3_CDS.gff -b ./4.Decompiled/Eco_Ctrl1.bed -s > ./6.Differential_expression/Eco_Ctrl1.cov
bedtools coverage -a ./ReferenceDB/NC_000913.3_CDS.gff -b ./4.Decompiled/Eco_Ctrl2.bed -s > ./6.Differential_expression/Eco_Ctrl2.cov
bedtools coverage -a ./ReferenceDB/NC_000913.3_CDS.gff -b ./4.Decompiled/Eco_Exp1.bed -s > ./6.Differential_expression/Eco_Exp1.cov
bedtools coverage -a ./ReferenceDB/NC_000913.3_CDS.gff -b ./4.Decompiled/Eco_Exp2.bed -s > ./6.Differential_expression/Eco_Exp2.cov

python3 ./Python_scripts/FormatDESeqInput.py -i ./6.Differential_expression/Eco_Ctrl1.cov:./6.Differential_expression/Eco_Ctrl2.cov:./6.Differential_expression/Eco_Exp1.cov:./6.Differential_expression/Eco_Exp2.cov -o ./6.Differential_expression/DESeq_Input.txt

sudo Rscript ./R_scripts/RunDESeq.R -i ./6.Differential_expression/DESeq_Input.txt -d ./6.Differential_expression/Design_sheet.txt -o ./6.Differential_expression/
