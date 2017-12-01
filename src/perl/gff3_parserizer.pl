#!/home/linux/perl/bin/perl

use strict;

# DO NOT USE FOR TRANSCRIPTOME ALIGNMENTS.
#  
# Rip a GFF3 file into a structure suitable for fast look ups
# of genes just based on a position. A position is a chromosome
# or transcript and a coordinate:
#
# Chr01:11050
#
# Coordinate is made out of 1st, 4th and 5th columns of the GFF file.
# There is a position made for every coordinate in a gene. This is a hash key.
# The hash value is the last column in the GFF, including the name of the gene. 
# So you can look up any position and find out what gene it is in. Like:
#
# Chr01:63062 => ID=Pp3c1_100;gene_id=Pp3c1_100;ancestorIdentifier=Pp3c1_100.v3.1
#

unless ($ARGV[1])
{
   print STDERR "\nUsage: $0 gff_file positions_umi_counts\n\n";
   exit 1;
}

#
# Specify GFF3 file and a file of positions and UMI counts.
#
my ($gff_file,$umi_file) = @ARGV;

my $lookup_dataref;
my $output_dataref;
my $i = 0;

#
# You can look up different features. The default is gene
#
my $target_feature = "gene";

print STDERR "Counting GFF lines... \n";
my $linecount = 1;
my $linecount = `wc -l $gff_file`;
($linecount) = $linecount =~ /(\d+)/;

open (GFF, "<$gff_file")|| die "Couldn't open the file $!";

<GFF>;<GFF>;<GFF>;

print STDERR "Building the lookup data structure... \n";
while (<GFF>)
{
   $i++;
   
   unless ($i % 100000)
   {
      print STDERR "Processed $i of $linecount GFF records... \n";
   }
   
   chomp;
   
   my ($chr,$method,$feature,$start,$stop,$score,$strand,$phase,$attributes) = split /\t/;
   
   unless($feature eq $target_feature)
   {
      next;
   }
   
   #
   # This builds the big hash for looking up gene belonging to a position, like:
   # Chr01:63062 => ID=Pp3c1_100;gene_id=Pp3c1_100;ancestorIdentifier=Pp3c1_100.v3.1
   #   
   for my $nuc ($start ..$stop )
   {
      my $position = $chr . ":" . $nuc;
      $lookup_dataref->{$position} = $attributes;
   }
}

print STDERR "Finished making the structure. Check memory usage... \n";

# Uses 22 Gbyte for genes only with Ppatens_318_v3.3.gene_exons.gff

# while (1)
# {
#    ;
# }

# Clusterized umi counts
# Chr01:10482857  25
# Chr01:10482858  70
# Chr01:10482859  90
# Chr01:10482860  33
# Chr01:10482861  11
# Chr01:10482862  5
#
# FOR genome alignments you have to do the lookup.
#

open (UMI, "<$umi_file")|| die "Couldn't open the file $!";

$i = 0;

print STDERR "Counting UMI lines... \n";
my $linecount = 1;
my $linecount = `wc -l $umi_file`;
($linecount) = $linecount =~ /(\d+)/;

print STDERR "Doing lookups... \n";

while (<UMI>)
{  
   $i++;
   
   unless ($i % 10000)
   {
      print STDERR "Processed $i of $linecount UMI records... \n";
   }
   
   chomp;
   
   my ($position,$umi_count) = split /\s+/;
   
   if (my $attributes = $lookup_dataref->{$position})
   {
      #
      # This may need to be changed depending on the format of the last field in the GFF file
      #
      my ($gene) = $attributes =~ /gene_id=(.*?);/;
      $output_dataref->{$gene}+=$umi_count;
   }
}
 

# Look at the counts
print STDERR "Printing results... \n";

#for my $key  (keys %{$output_dataref})

for my $key (sort {$output_dataref->{$b} <=> $output_dataref->{$a}} (keys %{$output_dataref}))
{
   print "$key   $output_dataref->{$key}\n";
}


