#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

print "#[1]Chrom\t[2]Start\t[3]End\t[4]State\t[5]MeanQual\n";

my (%beg,%prev,%cur,$qual,$nqual);
while (my $line=<STDIN>)
{
    @cur{qw(chr pos state qual)} = split(/\s+/,$line);
    if ( !%beg )
    {
        %beg   = %cur; 
        %prev  = %cur; 
        $qual += $cur{qual};
        $nqual++;
        next; 
    }
    if ( $beg{chr} eq $cur{chr} && $beg{state} eq $cur{state}) 
    { 
        %prev  = %cur;
        $qual += $cur{qual};
        $nqual++;
        next; 
    }
    $qual /= $nqual;
    printf "$beg{chr}\t$beg{pos}\t$prev{pos}\t$beg{state}\t%.1f\n",$qual;
    %beg   = %cur;
    $qual  = $cur{qual};
    $nqual = 1;
}
if ( %beg && %cur )
{
    $qual /= $nqual;
    printf "$beg{chr}\t$beg{pos}\t$prev{pos}\t$beg{state}\t%.1f\n",$qual;
}