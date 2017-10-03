#!/home/linux/perl/bin/perl

use strict;

# Feed on a list of UMIs and their positions
# umi_positions.txt
# TGCAGCGTTG   Chr01:629718   Chr06:10727258 Chr15:12754424 Chr16:7912214
# These positions are sorted in ascending order
# Group identical UMIs into clusters if they are within N nucleotides of each other.


unless ($ARGV[1])
{
   print "\nUsage: $0 umi_positions_file max_spacing\n\n";
   exit 1;
}

# umi positions file
my ($umifile,$max_spacing) = @ARGV;

open (UMI, "<$umifile")|| die "Couldn't open the file $!";

my $umi_to_position; #hash reference

# open each line at a time
while (<UMI>)
{
   my ($umi,@positions) = split /\s+/;  
   my $umi_data; #hash reference for each UMI
   
   for my $position (@positions)
   {
      my ($chr,$nuc) = split /:/,$position;   
      push @{$umi_data->{$chr}},$nuc;   # $umi_data->{$chr} is an array reference
   }
   
#    # Data look OK
#    while ( my ($chr,$nuc_positions_arrayref) = each %$umi_data)
#    { 
#       print "$umi   $chr   @$nuc_positions_arrayref\n";
#    }  
#      print "\n";

   
   # Now the data are packed by chromosome
   
   # Go through positions on the current chromosome for the current UMI  
   while ( my ($chr,$nuc_positions_arrayref) = each %$umi_data)
   {
      # If only one position on the chromosome there are no spacings
      # and that is the position of the UMI 
      if (scalar @$nuc_positions_arrayref == 1)
      {
         my $cluster_position = $chr . ":" . $nuc_positions_arrayref->[0];
         push@{$umi_to_position->{$umi}},$cluster_position;
      }
      else # What are the position spacings?
      {      
         #print "$umi   $chr   @$nuc_positions_arrayref\n";
         
         POSITION: for ( my $i = 0 ; $i <= $#$nuc_positions_arrayref-1 ; $i++)
         {
            # Get spacing between positions i and i+1 on the current chr for current UMI
            my $spacing = $nuc_positions_arrayref->[$i+1] - $nuc_positions_arrayref->[$i];
            
            #print "   $i   $nuc_positions_arrayref->[$i]   $nuc_positions_arrayref->[$i+1]   $spacing\n";
            
            # is the spacing <= the max?
            # Is there more to come?
            # Keep going
            if ($spacing <= $max_spacing && $i != $#$nuc_positions_arrayref-1)
            {
               ;
            }
            
            # is the spacing <= the max?
            # Is i+1 the last one?
            # Grab i+1
            if ($spacing <= $max_spacing && $i == $#$nuc_positions_arrayref-1)
            {
                #print "   ->$nuc_positions_arrayref->[$i+1]\n";
                my $cluster_position = $chr . ":" . $nuc_positions_arrayref->[$i+1];
                push@{$umi_to_position->{$umi}},$cluster_position;
            }
            
            # Is the spacing bigger than the max?
            # Is i+1 the last one?
            # Grab i and i+1
            elsif($spacing > $max_spacing && $i == $#$nuc_positions_arrayref-1)
            {
                #print "   ->$nuc_positions_arrayref->[$i]\n";
                my $cluster_position = $chr . ":" . $nuc_positions_arrayref->[$i];
                push@{$umi_to_position->{$umi}},$cluster_position;
                #print "   ->$nuc_positions_arrayref->[$i+1]\n";
                $cluster_position = $chr . ":" . $nuc_positions_arrayref->[$i+1];
                push@{$umi_to_position->{$umi}},$cluster_position;
            }
            
            # Is the spacing bigger than the max?
            # Is there more to come?
            # Grab i
            elsif($spacing > $max_spacing)
            {
                #print "   ->$nuc_positions_arrayref->[$i]\n";
                my $cluster_position = $chr . ":" . $nuc_positions_arrayref->[$i];
                push@{$umi_to_position->{$umi}},$cluster_position;
            }

         }      
      }      
   }   
}


# # Look at what happened
# for my $umi (sort keys %$umi_to_position)
# {
#   my @sorted_positions = sort @{$umi_to_position->{$umi}};  
#   print "$umi @sorted_positions\n";
# }

# AAAAAAACCG Chr04:81039 Chr22:1756765
# AAAAAACAGG Chr12:14544387 Chr22:11243819
# AAAAAACCGG Chr22:1756767
# AAAAAACGAT Chr02:13287863 Chr03:11442465 Chr04:20229245 Chr07:8920599 Chr23:12968463
# AAAAAAGGGG Chr10:8780329 Chr12:3531402 Chr19:2679838 Chr24:3236678
# AAAAAATGCG Chr23:11433528
# AAAAACAACG Chr01:6050233
# AAAAACAGGC Chr20:143854 Chr24:4391667
# AAAAACAGGG Chr01:23369203 Chr02:5395536 Chr12:14544388 Chr13:11617036
# AAAAACCGGG Chr20:15075734 Chr22:1756767
# AAAAACGAGG Chr14:6681933


# Turn the umi_to_position hash into position and UMI counts

my %position_umi_counts;

while ( my($umi,$position_array_ref) = each %$umi_to_position )
{
   for my $position(@$position_array_ref)
   {
      $position_umi_counts{$position}++;
   }  
}


for my $position (sort keys %position_umi_counts)
{
  print "$position  $position_umi_counts{$position}\n";
}
