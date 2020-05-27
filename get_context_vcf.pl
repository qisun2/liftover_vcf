#!/usr/bin/perl
use strict;
use warnings;

#####read in two files
#####file 1 is the chromosome file
#####file 2 is the hapmap file, 

my $fastafile = $ARGV[0];
my $vcfdir = $ARGV[1];
my $outdir = ".";

my $flanking = 100;

unless ($vcfdir) 
{
	print "Usage: get_context.pl Fastafile VCF file\n";
	print "You need the first 10 columns of the VCF file\n";
	exit;
}

unless (-e $fastafile) 
{
	print "$fastafile does not exist!\n";
	exit;
}


unless (-d $vcfdir) 
{
	print "$vcfdir does not exist!\n"; 
	exit;
}

use Bio::SeqIO;

my $in  = Bio::SeqIO->new(-file => "$fastafile" ,
      -format => 'Fasta');

my %chr2seq;
my %chr2length;
while ( my $seq = $in->next_seq() ) 
{
	my $id = $seq->display_id();
	if ($id=~/unknown/i) 
	{
		$id=0;
	}
	$chr2seq{uc $id} = $seq->seq();
	$chr2length{uc $id} = length($chr2seq{uc $id}) - 100;
}


opendir(my $dh, $vcfdir) || die "can't opendir $vcfdir: $!";
my @vcffiles = grep { /vcf/ && -f "$vcfdir/$_" } readdir($dh);
closedir $dh;

open OUT, ">$outdir/context.txt";
open OUT2, ">$outdir/hapmap.fastq";

foreach  (@vcffiles) 
{
	print "Start processing file $_ ... \n";
	my $file = "$vcfdir/$_";
	process_vcf($file);
}
close IN;
close OUT;


sub process_vcf
{
	my $vcffile = shift;
	if ($vcffile=~/gz$/) 
	{
		open (IN, '-|', "gunzip -c $vcffile") || die "no vcf file $vcffile";;
	}
	else
	{
		open (IN, "$vcffile") || die "no vcf file $vcffile";
	}
	$,="\t";

	HEADER:while (<IN>) 
	{
		last HEADER if (/^#CHROM/)
	}


	my (@data, $id, $refallele, $altalleles, $alleles, $mychr, $pos, $strand, $seqstr, $s1, $s2, $s3, $lastid);
	$strand = "+";
	$lastid = "";
	my $qstr = "a" x (2*$flanking+1);
	LOOP:while (<IN>) 
	{
		chomp;
		next LOOP unless (/\w/);
		@data = split "\t";
		($mychr, $pos, $refallele, $altalleles) = @data[0, 1, 3, 4 ];
		$alleles = "$refallele,$altalleles";
		$id = "${mychr}_${pos}";
		$mychr = uc $mychr;
		unless (exists $chr2seq{$mychr}) 
		{
			print "$mychr does not exist!\n";
			next LOOP;
		}
		#print $pos; exit;
		if ($id eq $lastid){ next LOOP;}
		$lastid = $id;
		if ($pos<($flanking+1)) 
		{	
			$s1 = substr ($chr2seq{$mychr}, 0, $pos-1);
			$s2 =substr($chr2seq{$mychr}, $pos-1, 1);
			$s3 = substr($chr2seq{$mychr}, $pos, $flanking);
			my $seq = $s1.$s2.$s3;
			my $newqstr = "a" x (length $seq);
			print OUT $id, $alleles, $mychr, $pos, $strand, $s1, $s2, $s3, "\n";
			$id .="_s". ($pos-1);
			print OUT2 "\@${id}\n$seq\n+\n$newqstr\n";
		}
		elsif ($pos > $chr2length{$mychr})
		{
			$s1 = substr ($chr2seq{$mychr}, $pos-($flanking+1), $flanking);
			$s2 =substr($chr2seq{$mychr}, $pos-1, 1);
			$s3 = substr($chr2seq{$mychr}, $pos);
			my $seq = $s1.$s2.$s3;
			my $newqstr = "a" x (length $seq);
			print OUT $id, $alleles, $mychr, $pos, $strand, $s1, $s2, $s3, "\n";
			$id .="_e". (length $s3);
			print OUT2 "\@${id}\n$seq\n+\n$newqstr\n";
		}
		else
		{
			$seqstr = substr ($chr2seq{$mychr}, $pos-($flanking+1), 2*$flanking+1);
			$s1 = substr($seqstr, 0, $flanking);
			$s2 =substr($seqstr, $flanking, 1);
			$s3 = substr($seqstr, $flanking+1, $flanking);
			print OUT $id, $alleles, $mychr, $pos, $strand, $s1, $s2, $s3, "\n";
			my $seq = $s1.$s2.$s3;
			my $newqstr = $qstr; 
			if ($seq=~/^(N*)([^N]+)(N*)$/) 
			{
				my ($header, $middle, $tail) = ($1, $2, $3);

				$seq = $middle;
				$newqstr =  "a" x (length $seq);
				if ($header) 
				{
					$id.="_h".(length $header);
				}
				if ($tail) 
				{
					$id.="_t".(length $tail);
				}
			}
			
			
			print OUT2 "\@${id}\n$seq\n+\n$newqstr\n";
		}

	}
}
