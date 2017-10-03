#!/home/linux/perl/bin/perl

use strict;

# Rip a GFF3 file into a structure suitable for fast look ups
# of genes just based on a position.



unless ($ARGV[0])
{
   print STDERR "\nUsage: $0 gff_file umi_file\n\n";
   exit 1;
}


my ($gff_file,$umi_file) = @ARGV;

my $lookup_dataref;
my $output_dataref;
my $i = 0;
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

# Open clusterized_umi_counts[2].txt to look up the lookup_dataref

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
      $output_dataref->{$attributes}+=$umi_count;
   }
}
 

# Look at the counts
print STDERR "Printing results... \n";

#for my $key  (keys %{$output_dataref})

for my $key (sort {$output_dataref->{$b} <=> $output_dataref->{$a}} (keys %{$output_dataref}))
{
   print "$key   $output_dataref->{$key}\n";
}


