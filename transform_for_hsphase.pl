#!/usr/bin/perl -w
open(I,"$ARGV[0]");
my %c={};
while(<I>)
  {next if(/^##/);
   chomp;
   @t=split;
   
   if(/^#CHROM/)
     {for($i=9;$i<@t;$i++)
        {$s[$i]=$t[$i];
	
	}
     }
   else
     {#next if ($c{$t[0]}==1);
      for($i=9;$i<@t;$i++)
        {$genotype=$1 if($t[$i]=~/(\S+)\:/);
	 if($genotype eq "0/0")
	   {push @{$s[$i]},0;}
	 elsif($genotype eq "0/1")
	   {push @{$s[$i]},1;}
	 elsif($genotype eq "1/1")
	   {push @{$s[$i]},2;}
#	 elsif($genotype eq "./.")
         else
	   {push @{$s[$i]},9;}

	}
     push @markers,$t[0]."_".$t[1];
     #$c{$t[0]}=1;
     }
  }
#open(A,">animal.id");
#open(M,">Marker.id");
print  "ID\t@markers\n";
for($i=9;$i<@s;$i++)
	  {print  "$s[$i]\t";
	   print "@{$s[$i]}\n";}
