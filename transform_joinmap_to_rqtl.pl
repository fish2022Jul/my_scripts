#!/usr/bin/perl -w
open(C,"$ARGV[0]"); #combined file of codom.jm.txt and dom.jm.txt
#less codom_JM|perl -e '{while(<>){@t=split(/\s+/,$_,7);{print "$t[2]_$t[4]\t$t[3]\t$t[6]";}}}' >codom.jm.txt
#less dom_JM|perl -e '{while(<>){@t=split(/\s+/,$_,3);{print "$t[2]";}}}' >dom.jm.txt
open(L,"$ARGV[1]"); #LOG file
while(<L>)
	{chomp;
	 @l=split;	
	}
close L;


while(<C>)
	{chomp;
#	 s/>//;s/hk/H/g;s/hh/A/g;s/kk/B/g;s/--/-/g;s/lm/H/g;s/ll/A/g;s/nn/A/g;s/np/H/g;
         s/>//;s/hk/10/g;s/hh/1/g;s/kk/4/g;s/--/NA/g;s/lm/6/g;s/ll/5/g;s/nn/7/g;s/np/8/g;
	 @t=split(/\s+/,$_);

	 push @all_markers,$t[0];
	 for($i=3;$i<int(@t);$i++)
	   {push @{$l[$i-3]},$t[$i];}
	}
print  "id,",join(",",@all_markers),"\n";
foreach $l(@all_markers)
	{print ",1";}
print "\n";
foreach $l(@l)
	{print "id$l,",join(",",@{$l}),"\n";}

