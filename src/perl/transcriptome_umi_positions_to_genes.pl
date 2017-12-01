#!/home/linux/perl/bin/perl

use strict;

# ONLY FOR TRANSCRIPTOME ALIGNMENTS.
#  
# Clusterized UMI counts for transcriptomem alignments are like this:
# 
# Pp3c10_10290V3.1:88  4
# Pp3c10_10290V3.1:89  19
# Pp3c10_10290V3.1:90  8
# Pp3c10_10290V3.1:91  1
# Pp3c10_10290V3.1:939  1
# Pp3c10_10290V3.1:969  1
#
# and like this for genome aligmnments:
#
# Chr01:10482857  25
# Chr01:10482858  70
# Chr01:10482859  90
# Chr01:10482860  33
# Chr01:10482861  11
# Chr01:10482862  5
# for genome alignments
#
# For genome alignments you have to look sup a GFF. For transcriptomes, the gene ID
# is in the pointer so easy to just aggregate the UMI counts for each gene.
#


unless ($ARGV[0])
{
   usage();
   exit 1;
}


#
# Specify GFF3 file and a file of positions and UMI counts. The latter is like this:
#
my $umi_file = shift;

my $output_dataref;

open (UMI, "<$umi_file")|| die "Couldn't open the file $!";

print STDERR "Counting UMI lines... \n";
my $linecount = 1;
my $linecount = `wc -l $umi_file`;
($linecount) = $linecount =~ /(\d+)/;

my $i=0;

while (<UMI>)
{  
   $i++;
   
   unless ($i % 10000)
   {
      print STDERR "Processed $i of $linecount UMI records... \n";
   }
   
   chomp;
   
   my ($position,$umi_count) = split /\s+/;
   my ($gene,$remainder) = split /V/,$position;
   
   $output_dataref->{$gene} += $umi_count;
}
 


for my $key (sort {$output_dataref->{$b} <=> $output_dataref->{$a}} (keys %{$output_dataref}))
{
   print "$key\t$output_dataref->{$key}\n";
}



sub usage
{
print<<DONE

Usage: $0 positions_umi_counts_file

Like this:

Pp3c10_10170V3.1:177  1
Pp3c10_10270V3.1:1426  1
Pp3c10_10270V3.1:1428  1
Pp3c10_10290V3.1:1026  1
Pp3c10_10290V3.1:1063  2
Pp3c10_10290V3.1:108  1
Pp3c10_10290V3.1:1087  1
Pp3c10_10290V3.1:1091  1
Pp3c10_10290V3.1:1115  1
Pp3c10_10290V3.1:1116  1

DONE
}
