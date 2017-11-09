#!/home/linux/perl/bin/perl

# Count numbers of UMI with 0,1,2,3,4,5 etc. spacings

use strict;

#
# Args: umi_positions file. From umi_analyzer_uthash or similar
#       samfile             Reads in sam format where the umi data came from
#       max_spacing:        Only analyze umi spacings <= this
#
# umi_positions.txt
# TGCAGCGTTG   Chr01:629718   Chr06:10727258 Chr15:12754424 Chr16:7912214
#

unless ($ARGV[2])
{
   print "\nUsage: $0 umi_positions_file samfile max_spacing\n\n";
   exit 1;
}

my ($umifile,$samfile,$max) = @ARGV;

open (UMI, "<$umifile")|| die "Couldn't open the file $!";
open (SAM, "<$samfile")|| die "Couldn't open the file $!";
my $umi_line_count = `wc -l $umifile`; chomp $umi_line_count;

# Index the sam file
print STDERR "Indexing the SAM file (long)...\n";
my $sam_data = indexify($samfile);

UMI:while (<UMI>)
{
   my ($umi,@positions) = split /\s+/;     
   my $umi_spacings = get_spacings(\@positions,$max);   
   
   my $mean_entropy = get_mean_entropy(\@positions,$sam_data);
   my $umi_entropy = get_entropy($umi);
   
   my $mean_mutual_info = mean_pairwise_mutual_information(\@positions,$sam_data);
   my $umi_mutual_info = pairwise_mutual_information($umi);
      
   my $spacing_count = 0;
   $spacing_count = (scalar @$umi_spacings || 0);
   my $has_spacing = 0;
   $has_spacing = 1 if ($spacing_count > 0);
   
   my @sorted_spacings = sort @$umi_spacings;
   my $max_spacing = 0;
   $max_spacing = ($sorted_spacings[$#sorted_spacings] || 0);
   
   my $has_evens = 0;
   my $has_odds = 0;
   my $nevens = 0;
   my $nodds = 0;
   my $mostly = 0;
   
   for my $spacing (@sorted_spacings)
   {      
      if (!($spacing%2))
      {
         $has_evens = 1; 
         $nevens++;      
      }
      
      if ($spacing == 1 || $spacing%2)
      {
         $has_odds = 1;
         $nodds++;
      }      
   }
   
   if ($nodds > $nevens)
   {
      $mostly = 1;
   }
   elsif($nevens > $nodds)
   {
      $mostly = 2;
   }
   
   print "$umi\t$spacing_count\t$has_spacing\t$max_spacing\t$mean_entropy\t$umi_entropy\t$mean_mutual_info\t$umi_mutual_info\t$has_odds\t$has_evens\t$mostly\n";
   #print "$umi   $has_evens   $mean_mutual_info   $umi_mutual_info\n";
}


#
# Calculate mean pairwise mutual information of a set of strings
#
sub mean_pairwise_mutual_information
{
   my ($positions,$sam_data) = @_;
   my $total = 0;
   my $n = 0; 
      
   for my $position(@{$positions})
   {
      my $read = $sam_data->{$position}->{_seq}; 
      
      #print "$position   $read\n";
          
      my $read_info = pairwise_mutual_information($sam_data->{$position}->{_seq});
      $n++;
      $total+=$read_info;
   }
   
   #return 1;
   
   return $total/$n;      
}


#
# Calculate the pairwise mutual information of a DNA string
#
sub pairwise_mutual_information
{
   my $read = shift;  
   substr($read,0,1,"") if ((length($read))%2);
   
   my $cp_read = $read;
   my $n_chars = length($read);
   my $n_pairs = $n_chars/2;
   my %base_n; my %base_p; my %pair_n; 
   
   my $I;
 
   while (my $base = uc(substr($read,0,1,"")))
   {
      $base_n{$base}++;
   }
   
   while(my ($base,$count) = each %base_n)
   {
      $base_p{$base} = $base_n{$base}/$n_chars;
   }
   
   while (my $pair = uc(substr($cp_read,0,2,"")))
   {
      $pair_n{$pair}++;
   }
   
   while(my ($pair,$count) = each %pair_n)
   {
      my $p_pair = $pair_n{$pair}/$n_pairs;     
      my ($x,$y) = split //,$pair;           
      my $i = (($p_pair) * (log2($p_pair/($base_p{$x}*$base_p{$y}))));
      $I+=$i;
   }
      
   return $I;   
}


#
# Calculate mean Shannon entropy of a set of strings
#
# A position is one of these: Chr06:10727258, i.e. chromosome and coordinate
#
# Get the sequence read corresponding to a position by loking up an indexed SAM file
#
# Calculate the entropy for each sequence read in the set and take the average
#
sub get_mean_entropy
{
   my ($positions,$sam_data) = @_;
   my $total = 0;
   my $n = 0;
   
   for my $position(@{$positions})
   {
      my $read = $sam_data->{$position}->{_seq};
      my $read_entropy = get_entropy($sam_data->{$position}->{_seq});
      $n++;
      $total+=$read_entropy;
   }
   return $total/$n;
}

#
# Calculate the Shannon entropy of a string (sequence read)
#
sub get_entropy
{
   my $read = shift;  
   my $n_chars = length($read);   
   
   my $nA = 0;
   my $nT = 0;
   my $nG = 0;
   my $nC = 0;
   
   # 
   # Count the bases
   #  
   while (my $char = uc(substr($read,0,1,"")))
   {
      if ($char eq "A")
      {
         $nA++;
      }
      elsif($char eq "T")
      {
         $nT++;
      }
      elsif($char eq "G")
      {
         $nG++;
      }
      elsif($char eq "C")
      {
         $nC++;
      }
      else
      {
         ;
      }
   }
   
   #
   # Compute frequencies of bases
   #
   my $fA = $nA/$n_chars;
   my $fT = $nT/$n_chars;
   my $fG = $nG/$n_chars;
   my $fC = $nC/$n_chars;
   
   my $e = 0;
   
   #
   # Calculate shannon entropy
   #
   for ($fA,$fT,$fG,$fC)
   {
      $e += ($_*(log2($_))) if ($_);
   }   
   
   return -$e;
}


#
# Log to base 2 of a number
#
sub log2
{
   my $n = shift;
   return log($n)/log(2);
}


#
# Index in memory of a SAM file. Uses a lot of RAM but look up is fast
#
# The index is a hash reference where the key is a position like contig:coodinate
#
sub indexify
{
   my $samfile = shift;
   
   my $count = `grep -v ^@ $samfile | wc -l`;
   chomp $count;

   open (SAM, "<$samfile")|| die "Couldn't open the file $!";

   my $sam_data;
   my $record_counter = 0;

   while (<SAM>)
   {
      next if(/^@/);
 
      $record_counter++;
      unless ($record_counter%1000000)
      {
          print STDERR "Processed $record_counter of $count SAM records\n";
      }
 
      my (@parts) = split /\t/;
      
      my $position = $parts[2] . ":" . $parts[3];
      
      $sam_data->{$position}->{_read} = $parts[0];
      $sam_data->{$position}->{_flags} = $parts[1];
      $sam_data->{$position}->{_contig} = $parts[2];
      $sam_data->{$position}->{_position} = $parts[3];
      $sam_data->{$position}->{_mapq} = $parts[4];
      $sam_data->{$position}->{_cigar} = $parts[5];
      $sam_data->{$position}->{_seq} = $parts[9];
      $sam_data->{$position}->{_qual} = $parts[10];
      $sam_data->{$position}->{_blob} = $parts[11];
 
   }
   
   return $sam_data;
}


#
# Get spacing data for a UMI. A position is contig:coord, like Chr06:10727258
#
# We want: 
#          Array of adjacent spacing sizes, which contains
#          [Number of spacings]
#          [Biggest_spacing]
#          [How many adjacent spacings]
#
#
sub get_spacings
{
   my ($positions,$max) = @_;
   my $umi_data;
   my @spacings;
   
   # Fill the $umi_data hashref for the current UMI
   # This is a hash reference. Key is a chromosome
   # or contig. Value is an array reference of positions 
   # on that chromosome. We are interested in adjacent
   # positions on the same chromosome.
   for my $position (@{$positions})
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
         my $spacing = ($nuc_positions_arrayref->[$i+1] - $nuc_positions_arrayref->[$i]);
         push @spacings,$spacing if ($spacing <= $max);
      }            
   } 
  
   return \@spacings;
}















