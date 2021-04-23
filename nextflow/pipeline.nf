sample = params.sample
qualityScore = params.qualityScore
target = file(params.outDir)
adapterSeq = params.adapterSeq

target.mkdir()

process downloadRef{
	tag "Ecoli Rel606"
	input:
		path target
	output:
		file "${target}/ref/ecoli_rel606.fasta*" into reference

	"""
	mkdir -p ${target}/ref/
	curl -L -o ${target}/ref/ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
	gunzip ${target}/ref/ecoli_rel606.fasta.gz
        bwa index ${target}/ref/ecoli_rel606.fasta
	"""
}

process download_fastq{
	tag "$sample" 
	input:
		path target	
	
	output:
		file "${target}/raw/*_[1,2].fastq" into readPair

	"""
	mkdir ${target}/raw
	fasterq-dump --split-files -O ${target}/raw/ $sample
	"""
}

process runFastQC{
	cpus { 2 }

	input:
		file '*' from readPair.collect()
		path target		

	output:
		file("${target}/fastqc/*.zip") into fastQC_files

	"""
	mkdir ${target}/fastqc
	fastqc --outdir ${target}/fastqc/ -t1 *.fastq
	"""
}

process trimReads{
	tag "<PHRED${qualityScore}"
	cpus { 4 }

	input:
		file '*' from readPair.collect()
		path target

	output:
		file("${target}/trimmed/*") into trimmedReads		

	"""
	mkdir ${target}/untrimmed
	mkdir ${target}/trimmed
	trimmomatic PE -threads 4 *.fastq ${target}/trimmed/${sample}_1.trimmed.fastq ${target}/untrimmed/${sample}_1.untrimmed.fastq ${target}/trimmed/${sample}_2.trimmed.fastq ${target}/untrimmed/${sample}_2.untrimmed.fastq \
ILLUMINACLIP:${adapterSeq}:2:40:15 MINLEN:25 SLIDINGWINDOW:4:${qualityScore}
	"""
}

process align{
	tag "${sample} to REL606 refseq"
	cpus { 4 }

	input: 
		file '*' from trimmedReads.collect()
		path target		
		file '*' from reference.collect()
	output: 
		file("${target}/sam/*") into alignedSAM

	"""
	mkdir ${target}/sam
	bwa mem -o ${target}/sam/${sample}.aligned.sam *.fasta *.fastq
	"""
}

process samToBAM{
	cpus { 2 }
	input:
		file alignedSAM
		path target	
	output:
		file("${target}/bam/*") into results		

	"""
	mkdir ${target}/bam
	samtools view -S -b *.sam -o ${target}/bam/${sample}.bam
	samtools sort ${target}/bam/${sample}.bam -o ${target}/bam/${sample}.sorted.bam
	rm ${target}/bam/${sample}.bam
	"""
}

