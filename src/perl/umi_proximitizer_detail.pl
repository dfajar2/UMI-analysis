#!/home/linux/perl/bin/perl

use strict;

# Feed on a list of UMIs and their positions
# Use the output of umi_analyzator.pl
# umi_positions.txt
# TGCAGCGTTG   Chr01:629718   Chr06:10727258 Chr15:12754424 Chr16:7912214


unless ($ARGV[0])
{
   print "\nUsage: $0 umi_positions_file\n\n";
   exit 1;
}

# umi positions file from umi_analyzatior.pl
my ($umifile) = @ARGV;

open (UMI, "<$umifile")|| die "Couldn't open the file $!";

my $spacing_data; 

while (<UMI>)
{
   my ($umi,@positions) = split /\s+/;  
   
   my $umi_data; #hash reference
   
   for my $position (@positions)
   {
      my ($chr,$nuc) = split /:/,$position;
      push @{$umi_data->{$chr}},$nuc;
   }

   CHR: while ( my ($chr,$nuc_positions_arrayref) = each %$umi_data)
   {
      # next chromosome if only one UMI
      next CHR if (scalar @$nuc_positions_arrayref == 1);
            
      for ( my $i = 0 ; $i <= $#$nuc_positions_arrayref-1 ; $i++)
      {
         my $spacing = $nuc_positions_arrayref->[$i+1] - $nuc_positions_arrayref->[$i];
	 $spacing_data->{$umi}->{$spacing}++;
      }
      
   }
   
}


for my $umi (sort keys %$spacing_data)
{
   print "$umi\n";
   my $spacing_hashref = $spacing_data->{$umi};
   for my $spacing (sort {$a<=>$b} keys %$spacing_hashref) 
   {
      print "   $spacing   $spacing_hashref->{$spacing}\n" if($spacing < 5);
   }
}
