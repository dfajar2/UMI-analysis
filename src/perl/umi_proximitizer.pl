#!/home/linux/perl/bin/perl

# Count numbers of UMI with 0,1,2,3,4,5 etc. spacings

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

# umi positions file from umi_analyzer_uthash or its siblings
my ($umifile) = @ARGV;

open (UMI, "<$umifile")|| die "Couldn't open the file $!";

my $spacing_data; 

while (<UMI>)
{
   my ($umi,@positions) = split /\s+/;  
   
   # For this UMI positions, key = chromosome, value = arrayref of coordinates on that chromosome 
   my $umi_data; 
   
   # Fill the $umi_data hashref for the current UMI
   for my $position (@positions)
   {
      my ($chr,$nuc) = split /:/,$position;
      push @{$umi_data->{$chr}},$nuc;
   }
   
   # Go through the position data for this UMI, chromosome by chromosome
   CHR: while ( my ($chr,$nuc_positions_arrayref) = each %$umi_data)
   {
      # If a UMI is found just once, call its separation 0
      if (scalar @$nuc_positions_arrayref == 1)
      {
         $spacing_data->{0}++;
      }
      else
      {
         # Go through the list of coords of this UMI and count the occurrences of each distance 
         for ( my $i = 0 ; $i <= $#$nuc_positions_arrayref-1 ; $i++)
         {
            my $spacing = $nuc_positions_arrayref->[$i+1] - $nuc_positions_arrayref->[$i];
	         $spacing_data->{$spacing}++;
         }
      }      
   }   
}


# fix thus so it prints shifts having no observed frequency

my @shifts = sort {$a<=>$b} keys %$spacing_data;
my $max_shift = $shifts[$#shifts];

for my $shift( 0 .. $max_shift )
{
  #my $value = $spacing_data->{$shift};
  #($value>0) ? $value : 0;
  
  my $value=0;
  if($spacing_data->{$shift})
  {
     $value = $spacing_data->{$shift};
  }

   print "$shift\t$value\n";
}


















