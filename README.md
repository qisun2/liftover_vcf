#### Liftover of SNPs between two genomes

This tool liftover SNPs between two genomes based on 200bp context sequences.

1. Extract 200bp context sequence for each SNP (100bp up- and down-stream )

  ```
  ./get_context_vcf.pl oldGenome.fasta genotype.vcf
  ```

  

After this, you should see two new files in your current directory: context.txt and hapmap.fastq

2. Run bwa to map the SNP flanking sequence to the new genome. "bwa aln" was used instead of "bwa mem" for higher stringency.

  ```
  bwa index newGenome.fasta
  bwa aln -t 8 newGenome.fasta hapmap.fastq |bwa samse  newGenome.fa - hapmap.fastq > hapmap.sam
  ```

  

3. parse sam file 

  ```
  ./get_pos_from_sam.pl
  ```

  This script will create a file named results

  This script will create a file named results

4. Create new.vcf file with uplifted positions.

  ```
  assign_newpos2vcf.pl outputDir newGenome.fasta
  ```

  This script creates new.vcf file with the new positions.

5. cd vcf

6. Sort the new.vcf files by new positions and create a sorted bcf file

```
sort -S 2G  -k1,1V -k2,2n new.vcf > new2.vcf
cat tmpheader.txt new2.vcf >new3.vcf
bgzip new3.vcf
tabix new3.vcf.gz
bcftools view -o genotype_uplifted.bcf -O b new3.vcf.gz
```

