# vcf_annotation

By Aditi Paranjpe


Script to annotate vcf file


Usage : annotate_vcf.pl -input <input_file_name> -output <output_file_name>

Annotated output file (default file name : annotated.tsv) contains following information:

1.Variant type (e.g. insertion, deletion, etc.). 

2.Variant effect (e.g. missense,synonymous, etc.).

3. Read depth at the site of variation. 

4. Number of reads supporting the variant.

5.Percentage of reads supporting the variant versus those supporting reference reads.

6. Allele frequency of variant
