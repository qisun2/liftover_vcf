#!/usr/local/bin/perl
use strict;
use warnings;
use Bio::Seq;

our $vcfdir = $ARGV[0];
our $beddir = $ARGV[1];
my $newgenome = $ARGV[2];



unless (-d $vcfdir) 
{
	print "$vcfdir does not exist!\n";
	exit;
}

unless (-d $beddir) 
{
	print "$beddir does not exist!\n";
	exit;
}

unless ($newgenome) 
{
	print "Usage: assign_newpos2vcf.pl vcfdir genome\n";
	print "You need the genome file\n";
	exit;
}


opendir(my $dh, $vcfdir) || die "can't opendir $vcfdir: $!";
my @vcffiles = grep { /vcf.gz$/ && -f "$vcfdir/$_" } readdir($dh);
closedir $dh;
make_header("$vcfdir/$vcffiles[0]");

my $c =0;

opendir($dh, $beddir) || die "can't opendir $beddir: $!";
my @bedfiles = grep { /bed$/ && -f "$beddir/$_" } readdir($dh);
closedir $dh;


our %id2map;
our %chr2seq;
our %chr2length;
our @cases = ();
our @strands = ();

foreach  (@bedfiles) 
{
	print "load bed file $_\n"; 
	open (IN, "$beddir/$_") || die "file $_ does not exists!\n";
	while (<IN>) 
	{
		chomp;
		my @data = split "\t";
		$id2map{$data[3]} = join ":", @data[0,2,4];
	}
	close IN;
}


print "Import genome fasta\n";
use Bio::SeqIO;
my $in  = Bio::SeqIO->new(-file => "$newgenome" ,
      -format => 'Fasta');
while ( my $seq = $in->next_seq() ) 
{
	my $id = $seq->display_id();
	if ($id=~/unknown/i) 
	{
		$id=0;
	}
	$chr2seq{$id} = $seq->seq();
	$chr2length{$id} = length($chr2seq{$id}) - 100;
}



our %ref =  (
"A" => "T", 
"T" => "A", 
"G" => "C", 
"C" => "G", 
"N" => "N",
"R" => "Y",
"Y" => "R",
"S" => "S",
"W" => "W",
"K" => "M",
"M" => "K",
"B" => "V",
"V" => "B",
"D" => "H",
"H" => "D",
"+" => "+",
"-" => "-",
"0" => "0"
);


print "start processing VCF ...\n";
$,="\t";



open OUT, ">$vcfdir/new.vcf";
open OUT2, ">$vcfdir/site.info";
open ERR, ">$vcfdir/errorlog.txt";

foreach  (@vcffiles) 
{
	print "Start processing file $_ ... \n";
	my $file = "$vcfdir/$_";
	process_vcf($file);
}

close OUT;
close ERR;
open STAT, ">$vcfdir/stat.txt";
print STAT "Ref allele matches\t".$cases[0]."\n";
print STAT "Ref allele matches alt\t".$cases[1]."\n";
print STAT "Ref allele matches one of alt\t".$cases[2]."\n";
print STAT "Ref allele no match\t".$cases[3]."\n";

print STAT "\nstrand +\t".$strands[0]."\n";
print STAT "strand -\t".$strands[1]."\n";

close STAT;

system("gzip $vcfdir/site.info");
sub make_header
{
	my $vcffile = shift;
	open HEAD, ">$vcfdir/tmpheader.txt";
	if ($vcffile=~/gz$/) 
	{
		open (IN, '-|', "gunzip -c $vcffile") || die "no vcf file $vcffile";;
	}
	else
	{
		open (IN, "$vcffile") || die "no vcf file $vcffile";
	}


	HEADER:while (<IN>) 
	{
		print HEAD $_;
		last HEADER if (/^#CHROM/i)
	}
	close IN;
	close HEAD;
}
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


	HEADER:while (<IN>) 
	{
		#print OUT $_;
		last HEADER if (/^#CHROM/i)
	}
	$,="\t";
	LOOP2:while (<IN>) 
	{
		next LOOP2 unless (/\w/);
		chomp;
		my $hmline = $_;
		my @data = split "\t";

		my ($mychr, $mypos, $refallele, $altalleles) = @data[0, 1, 3, 4 ];
		my $ori_refallele =  $refallele;
		my $ori_altalleles = $altalleles;
		my $id = "${mychr}-${mypos}";
		$data[2] = $id;
		unless (exists $id2map{$id})
		{
			print ERR "$id does not exist!\n";
			next LOOP2;
		}
		my $map = $id2map{$id};
		if ($map) 
		{
			$map=~/([^:]+):(\d+):([+\-])/;
			my ($newchr, $newpos, $newstrand) = ($1, $2, $3);
			if ($newstrand=~/\+/)
			{
				$strands[0] ++;
			}
			elsif ($newstrand=~/-/) 
			{
				$strands[1] ++;
				$refallele = $ref{$refallele};
				unless ($refallele=~/\w/) 
				{
					print "Reference allele cannot be revcom: $id, $data[3] \n";
					exit;
				}
				my @alt = split /,/, $altalleles;
				$altalleles = "";
				foreach my $a (@alt) 
				{
					my $b = "";
					if ($a=~/</) 
					{
						$b = $a;
					}
					else
					{
						$b= $ref{$a};
					}
					
					if ($b=~/\w/) 
					{
						$altalleles.=$b.",";
					}
					else
					{
						print "Alt allele cannot be revcom: $id, $data[4] \n";
						exit;
					}
				}
				$altalleles =~s/,$//;
				unless ($altalleles=~/\w/) 
				{
					print "Alt allele cannot be revcom: $id, $data[4] \n";
					exit;
				}
			}

			$data[0] = $newchr;
			$data[1] = $newpos;
			$data[3] = $refallele;
			$data[4] = $altalleles;

			my $c0to1 = 0;
			my $c1to0 = 0;
			my $newrefallele = substr($chr2seq{$newchr}, $newpos-1, 1);
			
			if ($refallele eq $newrefallele) 
			{
				print OUT (@data); print OUT "\n";
				print OUT2 0, $newstrand, $mychr, $mypos, $ori_refallele, $ori_altalleles, $refallele, $altalleles, @data[0, 1, 3,4], $c0to1, $c1to0; print OUT2 "\n";
				$cases[0]++;
			}
			elsif ($altalleles eq $newrefallele) 
			{	
				$cases[1]++;
				$data[3] = $newrefallele;
				$data[4] = $refallele;
				
				for (my $i=9; $i<=$#data; $i++) 
				{
					my @g = split /:/, $data[$i];
					if (@g==3) 
					{
						if ($g[0] eq "0/0") 
						{
							$g[0]="1/1";
							$c0to1 ++;
						}
						elsif ($g[0] eq "1/1")  
						{
							$g[0]="0/0";
							$c1to0 ++;
						}
						my @tmp = split /,/, $g[1]; $g[1]=join ",", (reverse @tmp);
						@tmp = split /,/, $g[2]; $g[2]=join ",", (reverse @tmp);
						$data[$i] = join ":", @g;
					}
				}
				print OUT (@data); print OUT "\n";
				print OUT2 1, $newstrand, $mychr, $mypos, $ori_refallele, $ori_altalleles, $refallele, $altalleles, @data[0, 1, 3,4], $c0to1, $c1to0; print OUT2 "\n";
			}
			elsif ($altalleles =~/$newrefallele/) 
			{
				$cases[2]++;
				my @alt = split /,/, $altalleles;
				my $refalleleIndex;
				III:for (my $i=0; $i<=$#alt; $i++) 
				{
					if ($alt[$i] eq $newrefallele) 
					{
						$refalleleIndex = $i;
						last III;
					}
				}
				$alt[$refalleleIndex] = $refallele;
				$data[3] = $newrefallele;
				$data[4] = join ",", @alt;

				for (my $i=9; $i<=$#data; $i++) 
				{
					my @g = split /:/, $data[$i];
					if ($g[0] eq "0/0") 
					{
						$c0to1 ++;
					}
					if ($g[0] eq "1/1") 
					{
						$c1to0 ++;
					}

					if (@g==3) 
					{
						my @tmp = split /\//, $g[0];
						for  (my $t=0; $t<@tmp; $t++) 
						{
							if ($tmp[$t]==0) 
							{
								$tmp[$t]=$refalleleIndex+1;
							}
							elsif ($tmp[$t]==($refalleleIndex+1)) 
							{
								$tmp[$t] = 0;
							}
						}
						$g[0] = join "/", @tmp;

						@tmp = split /,/, $g[1]; my $tt=$tmp[0]; $tmp[0] = $tmp[$refalleleIndex+1]; $tmp[$refalleleIndex+1] = $tt; $g[1]=join ",", @tmp;
						$data[$i] = join ":", @g;
					}
				}
				print OUT (@data); print OUT "\n";
				print OUT2 2, $newstrand, $mychr, $mypos, $ori_refallele, $ori_altalleles, $refallele, $altalleles, @data[0, 1, 3,4], $c0to1, $c1to0; print OUT2 "\n";
			}
			else
			{
				$cases[3]++;
				$data[3] = $newrefallele;
				$data[4] = $refallele.",".$altalleles;
				for (my $i=9; $i<=$#data; $i++) 
				{
					my @g = split /:/, $data[$i];
					if (@g==3) 
					{
						my @tmp = split /\//, $g[0];
						for  (my $t=0; $t<@tmp; $t++) 
						{
							$tmp[$t] ++;
						}
						$g[0] = join "/", @tmp;
						$g[1] = "0,".$g[1];
						@tmp = split /,/, $g[2]; $g[2]=join ",", (reverse @tmp);
						$data[$i] = join ":", @g;
					}
				}
				print OUT (@data); print OUT "\n";
				print OUT2 3, $newstrand, $mychr, $mypos, $ori_refallele, $ori_altalleles, $refallele, $altalleles, @data[0, 1, 3,4], $c0to1, $c1to0; print OUT2 "\n";
			}
			
		}

	}
	close IN;
}


#print "same $same; diff $diff\n";

