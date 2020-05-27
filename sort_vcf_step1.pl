#!/usr/bin/perl
use strict;
use warnings;

my $dir = $ARGV[0];
$dir=~s/\/$//;

my $genome_version = "agpv4";
my %filehandles;

my $dosplit = 0;

if($dosplit)
{
#make array of 10 file handles
for(my $i=1; $i<=10; $i++)
{
    local *FILE;
    open(FILE, ">$dir/tmpchr$i.txt") || die;
    $filehandles{$i} = *FILE;
}



open IN, "$dir/new.vcf";
while (my $line = <IN>) 
{
		$line=~/^(\S+)/;
		my $chr = $1;
		if (exists $filehandles{$chr}) 
		{
			print {$filehandles{$chr}} $line;
		}
}

for(my $i=1; $i<=10; $i++)
{
    close $filehandles{$i};
}

# For now, just split the files....
exit;
}


for(my $i=1; $i<=10; $i++)
{
	system ("cp $dir/tmpheader.txt $dir/${dir}_${genome_version}_chr$i.vcf");
    system ("LC_ALL=C sort -T ./ --buffer-size=120G -k2,2n $dir/tmpchr$i.txt >> $dir/${dir}_${genome_version}_chr$i.vcf");
	system ("bgzip $dir/${dir}_${genome_version}_chr$i.vcf");
	system ("tabix $dir/${dir}_${genome_version}_chr$i.vcf.gz");
	system ("rm new.vcf");
	system ("rm tmp*");
}
