#!/usr/bin/perl
use strict;
use warnings;
use JSON::PP qw(decode_json);
use Getopt::Long;


########################################################################################
# Parsing input and output file name, checking if the file exists.

my $input;
my $output = "annotated.tsv";
my $option;
GetOptions
("input=s" => \$input,
"output=s"   => \$output, 
)
or die("Error in command line arguments. Usage : annotate_vcf.pl -input <input_file_name> -output <output_file_name>\n");

die("ERROR: Input file path must be specified.\nIf output file name is not specified, output will be stored as annotated.vcf in current folder.") unless defined $input;

if ( ! -e $input) {die "File $input does not exist."};

########################################################################################
# Creating a hash for ranking variant types.  Rankings are according to: https://www.targetvalidation.org/variants

my %variants =(
transcript_ablation => 1,
curator_inference => 2,
trinucleotide_repeat_expansion => 3,
frameshift_variant => 4,
stop_gained => 5,
splice_region_variant => 6,
splice_acceptor_variant => 7,
splice_donor_variant => 8,
coding_sequence_variant => 9,
incomplete_terminal_codon_variant => 10,
stop_lost => 11,
protein_altering_variant => 12,
missense_variant => 13,
initiator_codon_variant => 14,
inframe_deletion => 15,
inframe_insertion => 16,
non_coding_transcript_exon_variant => 17,
NMD_transcript_variant => 18,
intron_variant => 19,
mature_miRNA_variant => 20,
"3_prime_UTR_variant" => 21,
"5_prime_UTR_variant" => 22,
non_coding_transcript_exon_variant => 23,
synonymous_variant => 24,
stop_retained_variant => 25,
regulatory_region_variant => 26,
upstream_gene_variant => 27,
downstream_gene_variant => 28,
TF_binding_site_variant => 29,
transcript_amplification => 30,
regulatory_region_amplification => 31,
TFBS_amplification => 32,
regulatory_region_ablation => 33,
TFBS_ablation => 34,
feature_truncation => 35,
feature_elongation => 36,
Regulatory_nearest_gene_five_prime_end => 37,
Nearest_gene_five_prime_end => 38,
sequence_variant => 39,
conservative_inframe_deletion => 40
);
my %reverse_variants = reverse %variants; #creating reverse reference to fetch variant_type using rank.

########################################################################################

open(FT,$input);
open(FS,">$output");
print FS "CHR\tPOS\tREF\tALT\tVARIANT_TYPE\tVARIANT_EFFECT\tDEPTH\tREADS_VARIANT\tPERCENTAGE_VARIANT\tALLELE_FREQUENCY\n";

while(my $line=<FT>)
{
	if($line !~ /#.*/)
	{
		chomp($line);	
		my @parts=split("\t",$line);
		my $chr=$parts[0];
		my $pos=$parts[1];
		my $ref_allele=$parts[3];
		my $alt_allele=$parts[4];
		my @each_alt = split(/,/,$alt_allele);
		my $all_AF = find_AF($parts[7]);					#subroutine find_AF is called for extracting Allele Frequency
		my @all_AF = @$all_AF;
		my ($read_depth,$var_depth,$var_per) = find_depth($parts[9]); #subroutine find_depth is called for extracting read depth , no.of reads supporting variant and for calculating % for each alternate allele.
		my @var_depth = @$var_depth; my @var_per=@$var_per;				#dereferencing arrays.
		
		for(my $i=0;$i<scalar(@each_alt);$i++)
		{
			my ($variant_type,$variant_effect) = process_variants($chr,$pos,$ref_allele,$each_alt[$i]);		#subroutine process_variants is called for identifying variant type and variant effect
			print FS "$chr\t$pos\t$ref_allele\t$each_alt[$i]\t$variant_type\t$variant_effect\t$read_depth\t$var_depth[$i]\t$var_per[$i]\t$all_AF[$i]\n";
		}
	}
}
close(FT);
close(FS);

########################################################################################

sub find_AF {
	my ($inputLine)= $_[0];
	my @find_AF = split(/;/,$inputLine);
		my @AF_list = split(/=/,$find_AF[3]);
		my @AF_list1 = split(/,/,$AF_list[1]);		#splitting lines to extract AF
		return(\@AF_list1);
}

sub find_depth {
	my ($inputLine)= $_[0];
	my @components=split(":",$inputLine);
	my $read_depth = $components[2];
	my @allele_depth = split(",",$components[3]);
	my $ref_depth = shift @allele_depth;		#removing reference depth from @allele_depth
	my @percent_variant;
	for(my $i=0;$i<scalar(@allele_depth);$i++)
	{
			my $per_variant= sprintf "%.2f" , ($allele_depth[$i]/($ref_depth+$allele_depth[$i])) * 100; #
			push(@percent_variant,$per_variant);
	}
	return($read_depth,\@allele_depth,\@percent_variant);
}

sub process_variants {
	my ($chr,$pos,$ref_allele,$alt_allele) =@_;
	my @chr_all = split("chr",$chr);
	$chr = $chr_all[1];
	my $variant_type; my $variant_effect = "--";
	if(length($ref_allele) eq length($alt_allele)) {$variant_type = "SNP";}		#comparing lengths of strings to identify variant type
	elsif(length($ref_allele) > length($alt_allele)) {$variant_type = "Deletion";}
	elsif(length($ref_allele) < length($alt_allele)) {$variant_type = "Insertion";}
	my $query = "$chr-$pos-$ref_allele-$alt_allele";
	print "Querying ExAC for $query ...\n";
	my $exac_output = `curl -s http://exac.hms.harvard.edu/rest/variant/consequences/$query`;
	if($exac_output ne "null")
	{
		my $json_all = decode_json $exac_output;
		my %json = %$json_all;
		my $rank=40;
		if (scalar keys %json > 1)
		{
			foreach (keys %json)
			{
				if($variants{$_} < $rank) {$rank= $variants{$_};}
					
			}
			$variant_effect = $reverse_variants{$rank};
		}
		else
		{
			foreach (keys %json)
			{
					$variant_effect = $_;
			}
		}	
	}
	else 
	{
		$variant_effect = "--";			#if EXAc does not have information for the variant- print "--" 
	}
	return($variant_type,$variant_effect);
}