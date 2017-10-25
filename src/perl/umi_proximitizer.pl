#!/home/linux/perl/bin/perl

use strict;
use Statistics::Descriptive;

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


#my $linecount = 1;
#my $linecount = `wc -l $umifile`;
#($linecount) = $linecount =~ /(\d+)/;

open (UMI, "<$umifile")|| die "Couldn't open the file $!";

# open each line at a time

my $stat = Statistics::Descriptive::Full->new();
my %spacing_data; 

while (<UMI>)
{
   my ($umi,@positions) = split /\s+/;  
   
   #print "$umi\n";
   
   my $umi_data; #hash reference
   
   for my $position (@positions)
   {
      my ($chr,$nuc) = split /:/,$position;
      
      #print "\t$chr\t$nuc\n";
      
      push @{$umi_data->{$chr}},$nuc;
   }

   # unpack the umi data to check
      
   CHR: while ( my ($chr,$nuc_positions_arrayref) = each %$umi_data)
   {
      # next chromosome if only one UMI
      next CHR if (scalar @$nuc_positions_arrayref == 1);
            
      #print "$umi   $chr   @$nuc_positions_arrayref\n";
      
      #print "\t$chr: @$nuc_positions_arrayref\n";
      
      # need some kind of distribution statistics on the 
      # array of positions for the chromosome. What are the 
      # position spacings?
            
      for ( my $i = 0 ; $i <= $#$nuc_positions_arrayref-1 ; $i++)
      {
         my $spacing = $nuc_positions_arrayref->[$i+1] - $nuc_positions_arrayref->[$i];
         $stat->add_data($spacing);
         $spacing_data{$spacing}++;
         
	 if($spacing < 5)
	 {
	    print "$umi   $spacing   $nuc_positions_arrayref->[$i]   $nuc_positions_arrayref->[$i+1]\n";
	 }
      }
      
   }
   
}

print "Total UMIs:\t" . $stat->count()  . "\n";
print "Median spacing:\t" . $stat->median() . "\n";
print "Mean spacing:\t" . $stat->mean() . "\n";
print "Mode spacing:\t" . $stat->mode() . "\n";
print "Min spacing:\t" . $stat->min() . "\n";
print "Max spacing:\t" . $stat->max() . "\n\n";

for my $k (sort {$a <=> $b} keys  %spacing_data)
{
  print "$k   $spacing_data{$k}\n";
}


# Total UMIs:     78207
# Median spacing: 1
# Mean spacing:   544420.354610201
# Mode spacing:   1
# Min spacing:    1
# Max spacing:    29063882


# 1   66354
# 2   1905
# 3   636
# 4   177
# 5   93
# 6   41
# 7   15
# 8   16
# 9   16
# 10   4
# 11   10
# 12   6
# 13   5
# 14   3
# 15   6
# 16   6
# 17   6
# 18   5
# 19   9
# 20   10
# 21   4
# 22   7
# 23   9
# 24   8
# 25   8
# 26   3
# 27   7
# 28   6
# 29   4
# 30   4
# 31   8
# 32   6
# 33   6
# 34   6
# 35   3
# 36   6
# 37   6
# 38   2
# 39   8
# 40   6
# 41   2
# 42   4
# 43   1
# 44   1
# 45   1
# 46   4
# 47   5
# 48   5
# 49   3
# 50   1
# 51   2
# 52   2
# 53   2
# 54   4
# 55   1
# 56   2
# 57   3
# 58   1
# 59   2
# 60   4
# 61   1
# 62   5
# 63   1
# 64   4
# 65   5
# 66   5
# 67   3
# 68   4
# 69   2
# 70   1
# 71   4
# 72   2
# 73   2
# 74   3
# 75   2
# 76   1
# 77   2
# 80   3
# 81   2
# 82   2
# 83   2
# 84   1
# 85   2
# 86   1
# 87   2
# 88   1
# 90   2
# 91   2
# 92   3
# 93   2
# 94   8
# 95   6
# 96   6
# 97   1
# 98   8
# 99   5
# 100   2
# 101   4
# 102   3
# 103   6
# 104   5
# 105   2
# 106   5
# 107   5
# 108   6
# 109   2
# 110   4
# 111   6
# 112   6
# 113   5
# 114   5
# 115   4
# 116   4
# 117   15
# 118   4
# 119   9
# 120   7
# 121   2
# 122   3
# 123   9
# 124   6
# 125   8
# 126   2
# 127   12
# 128   5
# 129   1
# 130   5
# 131   7
# 132   11
# 133   10
# 134   4
# 135   8
# 136   13
# 137   5
# 138   8
# 139   6
# 140   9
# 141   9
# 142   7
# 143   4
# 144   3
# 145   8
# 146   8
# 147   6
# 148   4
# 149   9
# 150   4
# 151   5
# 152   2
# 153   8
# 154   9
# 155   6
# 156   4
# 157   4
# 158   7
# 159   9
# 160   7
# 161   4
# 162   5
# 163   11
# 164   13
# 165   11
# 166   8
# 167   7
# 168   8
# 169   9
# 170   10
# 171   7
# 172   5
# 173   4
# 174   5
# 175   4
# 176   4
# 178   6
# 179   3
# 180   9
# 181   9
# 182   2
# 183   5
# 184   3
# 185   4
# 186   9
# 187   6
# 188   4
# 189   6
# 190   2
# 191   2
# 192   3
# 193   4
# 194   1
# 195   5
# 196   5
# 197   7
# 198   2
# 199   4
# 200   2
# 201   3
# 202   7
# 203   6
# 204   7
# 205   7
# 206   6
# 207   2
# 208   6
# 209   4
# 210   4
# 211   3
# 212   3
# 213   3
# 214   4
# 216   2
# 217   1
# 218   5
# 219   2
# 220   2
# 221   5
# 222   3
# 223   4
# 225   3
# 226   1
# 227   6
# 228   5
# 229   8
# 230   6
# 231   8
# 232   7
# 233   5
# 234   3
# 235   4
# 236   2
# 237   1
# 238   1
# 239   6
# 240   3
# 241   5
# 242   5
# 243   2
# 244   2
# 245   2
# 246   5
# 247   9
# 248   4
# 249   5
# 250   1
# 251   3
# 253   1
# 255   2
# 256   1
# 257   1
# 258   1
# 259   2
# 260   2
# 261   2
# 262   4
# 263   2
# 264   1
# 265   3
# 266   3
# 267   2
# 268   4
# 269   1
# 270   4
# 273   3
# 274   3
# 275   3
# 276   3
# 277   2
# 278   3
# 279   1






