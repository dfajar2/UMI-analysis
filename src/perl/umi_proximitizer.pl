#!/home/linux/perl/bin/perl

# Count numbers of UMI with 0,1,2,3,4,5 etc. spacings

use strict;

# Feed on a list of UMIs and their positions
# Use the output of umi_analyzator.pl
# umi_positions.txt
# TGCAGCGTTG   Chr01:629718   Chr06:10727258 Chr15:12754424 Chr16:7912214

#
# This counts the number of UMI (column 2) having shift size in column 1
#
# There is deliberate redundancy. A UMI can have a shift of 1 and 2. The number
# of UMI having no shift at all is also computed. 
#



unless ($ARGV[0])
{
   print "\nUsage: $0 umi_positions_file\n\n";
   exit 1;
}

# umi positions file from umi_analyzer_uthash or its siblings
my ($umifile) = @ARGV;

open (UMI, "<$umifile")|| die "Couldn't open the file $!";

my $spacing_data; 
my $total_umi_count = 0;
my $shifty_umi_count = 0;
my $no_shifty_umi_count = 0;

while (<UMI>)
{
   my ($umi,@positions) = split /\s+/; 
   
   $total_umi_count++;
   
   my $shifty_umi = 0;
   
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
      # Go through the list of coords of this UMI and count the occurrences of each distance 
      for ( my $i = 0 ; $i <= $#$nuc_positions_arrayref-1 ; $i++)
      {
         my $spacing = $nuc_positions_arrayref->[$i+1] - $nuc_positions_arrayref->[$i];
	      $spacing_data->{$spacing}++;
         
         if($spacing > 0)
         {
            $shifty_umi = 1;  
         }
      }              
   } 
   
   if($shifty_umi)
   {
      $shifty_umi_count++; 
   }  
   else
   {
      $no_shifty_umi_count++;
   }
}


# fix thus so it prints shifts having no observed frequency

my @shifts = sort {$a<=>$b} keys %$spacing_data;
my $max_shift = $shifts[$#shifts];

my $checked_no_shifty_umi_count = $total_umi_count - $shifty_umi_count;

print "----------------------------\n";
print "no_shifty_umi_count: $no_shifty_umi_count\n";
print "checked_no_shifty_umi_count: $checked_no_shifty_umi_count\n";
print "shifty_umi_count: $shifty_umi_count\n";
print "total_umi_count: $total_umi_count\n";
print "----------------------------\n";

for my $shift( 1 .. $max_shift )
{
  
  my $value=0;
  if($spacing_data->{$shift})
  {
     $value = $spacing_data->{$shift};
  }

   print "$shift\t$value\n";
}


















