#!/usr/bin/perl
use strict;
use warnings;

open IN, "hapmap.sam";
open OUT, ">results";
open UM, ">unmapped.sam";
open ERROR, ">error.log";
open STAT, ">stat.log";
my $flanking =100;

my $multi = 0;
my $unmapped = 0;
my $mapped =0;

$,="\t";
LOOP:while (<IN>) 
{
	next LOOP if (/^\@/);
	my $line = $_;
	chomp;
	my @data = split "\t";
	my ($id, $flag, $chr, $pos, $mapq, $cigar, @tags) = @data[0..5,  11..$#data];

	my $tag = join "\t", @tags;
	#print ($id, $flag, $chr, $pos, $mapq, $cigar, $x0); 
	if ($flag==4)
	{
		$unmapped ++; 
		print UM $line;
		$data[10] .= "\tunmapped";
			if ($id=~/(.+)_[seht]\d+$/)
			{
				$data[0] = $1;
			}
		print OUT @data, "\n";
		next LOOP;
	}	
	if ($tag =~/X0:i:(\d+)/) 
	{
		my $count = $1;
		
		if ($count!=1) 
		{
			$multi ++;

			$data[10] .= "\tmulti_mapped";
			print OUT @data, "\n";
			next LOOP;
		}
	}
	else
	{
		print "Warning x0 tag: $line";
		exit;
	}
	$cigar =~/^(\d+)/;
	my $cigar5p = $1;
	$cigar =~/(\d+)M$/;
	my $cigar3p = $1;
	my @alldel = $cigar =~/(\d+)D/g;
	my @allins = $cigar=~/(\d+)I/g;
	my $endpos = $pos + ($flanking*2);
	foreach  (@alldel) 
	{
		$endpos += $_;
	}
	foreach  (@allins) 
	{
		$endpos -= $_;
	}
	
	my $strand = "";

	if ($flag == 0) 
	{
		$strand = "+";
		my $oid = $id;
		if ($id=~/(.+)_[seht]\d+$/) 
		{
			$oid = $1;
			#print "ori id", $oid, "\n";
			if ($id=~/.+_([th])(\d+)$/) 
			{
				my $type = $1;
				my $clipsize = $2;
				if (($type eq "h") && (($cigar5p + $clipsize)>=($flanking+1)) )
				{
					$pos = $pos + ($flanking - $clipsize) ;
				}
				elsif (($type eq "t") && (($cigar5p)>=($flanking+1)) )
				{
					$pos = $pos + $flanking;
				}
				else
				{
					$unmapped ++;
					print UM  $line;
					$data[10] .= "\tUnmapped2";
					$data[0] = $oid;

					print OUT @data, "\n";
					next LOOP;
				}			
			}
			elsif (($id=~/.+_([e])(\d+)$/) && ($cigar5p>=($flanking+1)) )
			{
				$pos= $pos + $flanking;
			}
			elsif (($id=~/.+_([s])(\d+)$/) && ($cigar5p>=$2) )
			{
				$pos= $pos + $2;
			}
			else
			{
					$unmapped++;
					print UM  $line;
					$data[10] .= "\tUnmapped2";
					$data[0] = $oid;

					print OUT @data, "\n";
					next LOOP;

			}
		}
		else
		{	
			if ($cigar5p>=($flanking+1)) 
			{
				$pos= $pos + $flanking;
			}
			elsif ($cigar3p>=($flanking+1))
			{
				$pos = $endpos -$flanking;
			}
			else
			{
				$unmapped ++;
				print UM  $line;
				$data[10] .= "\tUnmapped2";
				$data[0] = $oid;

				print OUT @data, "\n";

				next LOOP;
			}
		}
	}
	elsif ($flag ==16)
	{
		$strand = "-";
		my $oid = $id;
		if ($id=~/(.+)_[seht]\d+$/) 
		{
			$oid = $1;
			print $oid, "\n";
			if ($id=~/.+_([th])(\d+)$/) 
			{
				my $type = $1;
				my $clipsize = $2;
				if (($type eq "h") && (($cigar5p)>=($flanking+1)) )
				{
					$pos = $pos + $flanking ;
					#print $id, " 16h\n";
				}
				elsif (($type eq "t") && (($cigar5p + $clipsize)>=($flanking+1)) )
				{
					$pos = $pos + ($flanking - $clipsize) ;
					#print $id, " 16t\n";
				}
				else
				{
					$unmapped ++;
					print UM $line;
					$data[10] .= "\tUnmapped2";
					$data[0] = $oid;
					print OUT @data, "\n";
					next LOOP;
				}			
			}
			else
			{
				$unmapped ++;
				print UM $line;
				$data[10] .= "\tUnmapped2";
				$data[0] = $oid;
				print OUT @data, "\n";
				next LOOP;
			}
		}
		else
		{
			if ($cigar5p>=($flanking+1)) 
			{
				$pos= $pos + $flanking;
			}
			elsif ($cigar3p>=($flanking+1))
			{
				$pos = $endpos -$flanking;
			}
			else
			{
				$unmapped ++;
				print UM $line;
				$data[10] .= "\tUnmapped2";
				$data[0] = $oid;
				print OUT @data, "\n";
				next LOOP;
			}
		}
	}
	else
	{
		print ERROR "unknown flag: $line";
		$data[10] .= "\tunmapped3";
		print OUT @data, "\n";
		next LOOP;
	}
	$mapped ++;
	$data[10] .= "\tR2:$chr:$pos:$strand";
	if ($data[0]=~/(.+)_[seht]\d+$/) 
	{
		$data[0] = $1;
	}
	print OUT @data, "\n";

}

print STAT "multi_mapped\t$multi\n";
print STAT "unmapped\t$unmapped\n";
print STAT "Mapped\t$mapped\n";
