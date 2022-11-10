#!/bin/bash
#SBATCH --job-name=BBmap_RNA_contigs
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@provider.net
#SBATCH --output=logs/BBmap_RNA2.out
#SBATCH --error=logs/BBmap_RNA2.err
#SBATCH --nice=1000
#SBATCH --partition=himem

#load modules

module load samtools
module load bbmap

DATA_FOLDER=/scratch/dmeier/Avdat_transcriptomes/

#Mapping the transcriptome reads to the contigs with BBmap

cd $TMPDIR/

mkdir reads

scp $DATA_FOLDER/reads/*_no_rRNA_*fastq.gz reads/
scp $DATA_FOLDER/Avdat_FinalBins_plus_Unbinned.fasta ./

mkdir BBmap_transcript_mapping

for SAMPLE in $(ls reads/*read1.fastq.gz | sed 's/reads\///g' | sed 's/_no_rRNA_read1.fastq.gz//g');

do

  bbmap.sh \
    ref=Avdat_FinalBins_plus_Unbinned.fasta \
      in=reads/"$SAMPLE"_no_rRNA_read1.fastq.gz \
      in2=reads/"$SAMPLE"_no_rRNA_read2.fastq.gz \
      outm=BBmap_transcript_mapping/"$SAMPLE"_mapped.sam \
      bamscript=BBmap_transcript_mapping/"$SAMPLE".sh \
      killbadpairs=true \
      pairedonly=true \
      minid=99 \
      idfilter=97 \
      threads=32

  bash BBmap_transcript_mapping/"$SAMPLE".sh
  rm BBmap_transcript_mapping/"$SAMPLE"_mapped.sam
      
done

scp -r BBmap_transcript_mapping $DATA_FOLDER/BBmap_transcript_mapping

module unload bbmap
module load subread

featureCounts \
    -T 32 \
    --tmpDir /tmp/dmeier \
    -a Combined_annotation_file_for_featureCounts.csv \
    -F SAF \
    -o BBmap_transcript_mapping_all/Avdat_featureCounts.tsv \
    BBmap_transcript_mapping/*.bam


