#!/home/linux/perl/bin/perl

use strict;
use Text::Levenshtein::XS qw/distance/;
use Try::Tiny;

# Feed on a sorted SAM file. 
# Identify UMI duplicates at the same position. 
# Collapse duplicates pairwise into most abundant UMI
# Count the UMIs at each position
# Print to STDOUT


unless ($ARGV[1])
{
   print "\nUsage: $0 sorted_samfile max_edit_distance\n\n";
   exit 1;
}


# sorted sam file and max mismatches to tolerate between two UMIs
my ($samfile,$max_edit_distance) = @ARGV;

if (!$max_edit_distance)
{
   print "Must specify max_edit_distance > 0\n";
   exit 1;
}

#my $linecount = 1;
my $linecount = `wc -l $samfile`;
($linecount) = $linecount =~ /(\d+)/;

open (SAM, "<$samfile")|| die "Couldn't open the file $!";

# place to hold the position of the last record
my $last_position = "start";

#hash reference for current position data
my $position_data; 

# count lines
my $line_counter = 0;

# starting time
my $last_current_time = time;
my $current_time = time;

   
#
# This groups key elements of the sam file at each position
#
while (<SAM>) 
{
   $line_counter++;    
   
   unless ($line_counter%1000000)
   {
       $current_time = time;
       my $elapsed_time = $current_time - $last_current_time;
       print STDERR "Processed $line_counter of $linecount SAM records in $elapsed_time seconds...\n";
   }
      
   chomp;
   
   # skip the headers
   next if (/^\@/);
   
   # split the SAM record
   my @parts = split /\s+/;
      
   # skip the multimaps
   next if ($parts[11] !~ /NH:i:1$/);  #print "$parts[11]\n";
   
   # pull out the UMI
   my ($umi) = $parts[0] =~/#([ATCGN]+)/;
 
   # put chromosome and nucleotide into one scalar
   my $position = $parts[2] . ":" . $parts[3];
      
   # Be clever and group lines from the same position.
   # If this is a new position, send the data off to
   # be collapsed and counted. Without the match ne "start"
   # it would send the first line to the collapser.
   
   if($position ne $last_position && $last_position ne "start")
   {
      # New position, so send the LAST position data off for
      # collapsing and counting, empty the structure, and start
      # it off again with the new position data. 
      
      # take a look at the data pre-collapse
      #print "==================================\n";
      #print "$last_position\n"; 
      #print "Pre-collapso:\n"; 
      #while ( my($k,$v) = each  %$position_data)
      #{
      #  print "   $k->$v\n";
      #}
      
      # send for collapsing
      my $collapsed_position_data;
      try
      {
          $collapsed_position_data = collapse($position_data);
      }
      catch
      {
          print "collapse() failed in main loop\n";
          print "\$_: $_\n";
          while ( my($k,$v) = each  %$position_data)
          { 
            print "   $k->$v\n";
          }
      };
      
      my $num_umis = scalar keys %$collapsed_position_data;
      my @umis = keys %$collapsed_position_data;
 
      # print to STDOUT for the time being
      print "$last_position $num_umis @umis\n";
      
      # empty the structure
      %$position_data = ();
      
      # start adding data for new position
      $position_data->{$umi}++;
   }
   else
   {
      # Current position has not changed.
      # Keep adding data for the current position
      $position_data->{$umi}++;
   }
   
   # set the position for next time around    
   $last_position = $position;
   $last_current_time = $current_time;
}





# Collapse UMIs all within a given edit distance of each other 
# to one UMI. 
sub collapse
{
   my $position_data = shift;
   
   # how many are there to begin with
   my $starting_umi_count = scalar keys %$position_data;
  
   my @sorted_hash_keys = sort keys %$position_data;
   
   # traverse the set pairwise
   for (my $i = 0 ; $i <= $#sorted_hash_keys ; $i++)
   {
      for (my $j = $i+1 ; $j <= $#sorted_hash_keys ; $j++)
      {
          # the traversal looks ok
          #print "$i: $position_data->{$sorted_hash_keys[$i]}    $j: $position_data->{$sorted_hash_keys[$j]}\n";
          
          my $umi_i = $sorted_hash_keys[$i];
          my $umi_j = $sorted_hash_keys[$j];
          
          # get the levenstein distance between the pair of UMIs
          my $d = distance($umi_i,$umi_j);
         
          #print "$sorted_hash_keys[$i]:$position_data->{$sorted_hash_keys[$i]}   $sorted_hash_keys[$j]:$position_data->{$sorted_hash_keys[$j]}   Distance = $d\n";
          
          my $count_i = $position_data->{$umi_i};
          my $count_j = $position_data->{$umi_j};
                   
         # if they are close, merge the UMIs into the most abundant and sum the counts
         if ( $d <= $max_edit_distance )
         {
            my $sum = $count_i + $count_j;
            
            if($count_i >= $count_j)
            {
               $position_data->{$sorted_hash_keys[$i]} = $sum;
               $position_data->{$sorted_hash_keys[$j]} = 0;
            }
            else
            {
               $position_data->{$sorted_hash_keys[$i]} = 0;
               $position_data->{$sorted_hash_keys[$j]} = $sum;              
            }
         }
         #else
         #{
         #   print "No merge...\n";
         #}
      } 
   }

   #print "Collapsed data:\n";
   #while ( my ($k,$v) = each %$position_data)
   #{
   #   print "   $k   $position_data->{$k}\n";
   #}

   # get rid of the zeroes
   while ( my ($k,$v) = each %$position_data)
   {
      delete $position_data->{$k} if($v == 0);
   }
   
   my $ending_umi_count = scalar keys %$position_data;
   
   #print "$starting_umi_count -> $ending_umi_count\n";
   
   # i can't prove that once through the collapser will get to the correct solution
   # recursion may not be necessary but it will do the job. Once nothing changes we 
   # we are at a minimum.   
   if ($starting_umi_count != $ending_umi_count)
   {
      try
      {
         collapse($position_data);
      }
      catch
      { 
          print "collapse() failed in recursion loop\n";
          print "\$_: $_\n";
          while ( my($k,$v) = each  %$position_data)
          { 
            print "   $k->$v\n";
          }
      };
   }
   else
   {
      return($position_data);
   }
  
}

