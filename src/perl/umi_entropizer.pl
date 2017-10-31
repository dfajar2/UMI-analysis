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

# umi positions file from umi_analyzer_uthash or its siblings
my ($umifile,$samfile,$max) = @ARGV;

open (UMI, "<$umifile")|| die "Couldn't open the file $!";

# Index the sam file
my $sam_data = indexify($samfile);




my $spacing_data; 

while (<UMI>)
{
   my ($umi,@positions) = split /\s+/;  
   
   my $umi_spacings = get_spacings(\@positions,$max);
   
   print "$umi   @{$umi_spacings}\n";
   
   #my $read_entropies;
   #my $umi_entropy;
   
   
}


sub indexify
{
   my $samfile = shift;
   
   my $count = `grep SRR $samfile | wc -l`;
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
      $sam_data->{$parts[0]}->{_flags} = $parts[1];
      $sam_data->{$parts[0]}->{_contig} = $parts[2];
      $sam_data->{$parts[0]}->{_position} = $parts[3];
      $sam_data->{$parts[0]}->{_mapq} = $parts[4];
      $sam_data->{$parts[0]}->{_cigar} = $parts[5];
      $sam_data->{$parts[0]}->{_seq} = $parts[9];
      $sam_data->{$parts[0]}->{_qual} = $parts[10];
      $sam_data->{$parts[0]}->{_blob} = $parts[11];
 
   }
   
   return $sam_data;
}

sub get_spacings
{
   my ($positions,$max) = @_;
   my $umi_data;
   my $spacing_data;
   
   # Fill the $umi_data hashref for the current UMI
   for my $position (@{$positions})
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
         $spacing_data->{0} = 1;
      }
      else
      {
         # Go through the list of coords of this UMI and count the occurrences of each distance 
         for ( my $i = 0 ; $i <= $#$nuc_positions_arrayref-1 ; $i++)
         {
            my $spacing = $nuc_positions_arrayref->[$i+1] - $nuc_positions_arrayref->[$i];
	         $spacing_data->{$spacing} = 1 if ($spacing <= $max);
         }
      }      
   } 
   
   my @spacings = sort {$a<=>$b} keys %{$spacing_data};  
   
   return \@spacings;
}

















