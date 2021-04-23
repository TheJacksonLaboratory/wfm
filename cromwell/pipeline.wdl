task downloadRef{
	String target

	command <<<
	mkdir -p ${target}/ref/
	curl -L -o ${target}/ref/ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
	gunzip ${target}/ref/ecoli_rel606.fasta.gz
        bwa index ${target}/ref/ecoli_rel606.fasta	
	>>>

	output{
		File Ref = "${target}/ref/ecoli_rel606.fasta"
	}

	runtime {
		cpus: "4"
	}				
}

task downloadSRA { 
	String sample
	String target

	command <<<
		mkdir -p ${target}/raw/
		fasterq-dump --split-files -O ${target}/raw/ ${sample}	
	>>>

	output{
		Array[File] rawFastQ = ["${target}/raw/${sample}_1.fastq", "${target}/raw/${sample}_2.fastq"]
	}
}

task runFastQC {
	String sample
	String target
	Array[File] rawFastQ

	command <<<
		mkdir ${target}/fastqc
		fastqc --outdir ${target}/fastqc/ -t1 ${sep=' ' rawFastQ}
	>>>

	runtime {
		cpus: "2"
	}
}

task trimReads{
	String sample
	String target
	Int qualityScore
	Array[File] rawFastQ
	File adapterSeq

	command <<<
		mkdir ${target}/trimmed
		mkdir ${target}/untrimmed
		trimmomatic PE -threads 4 ${sep=' ' rawFastQ} ${target}/trimmed/${sample}_1.trimmed.fastq ${target}/untrimmed/${sample}_1.untrimmed.fastq ${target}/trimmed/${sample}_2.trimmed.fastq ${target}/untrimmed/${sample}_2.untrimmed.fastq \
ILLUMINACLIP:${adapterSeq}:2:40:15 MINLEN:25 SLIDINGWINDOW:4:${qualityScore}
	>>>
	
	runtime { 
		cpus: "4"
	}	

	output{
		Array[File] trimmedFastQ = ["${target}/trimmed/${sample}_1.trimmed.fastq", "${target}/trimmed/${sample}_2.trimmed.fastq"]
		Array[File] untrimmedFastQ = ["${target}/untrimmed/${sample}_1.untrimmed.fastq", "${target}/untrimmed/${sample}_2.untrimmed.fastq"]
	}
}

task align{
	String sample
	String target
	Array[File] trimmedFastQ	
	File Ref

	command <<<
        mkdir ${target}/sam
        bwa mem -o ${target}/sam/${sample}.aligned.sam ${Ref} ${sep=' ' trimmedFastQ}
	>>>

	runtime { 
                cpus: "4"
        }

	output{
		File alignedSAM = "${target}/sam/${sample}.aligned.sam"
	}
} 

task samToBAM {
	String sample
	String target
	File alignedSAM

	command <<<
	mkdir ${target}/bam
	samtools view -S -b ${alignedSAM} -o ${target}/bam/${sample}.bam
	samtools sort ${target}/bam/${sample}.bam -o ${target}/bam/${sample}.sorted.bam
	rm ${target}/bam/${sample}.bam	
	>>>

	runtime { 
                cpus: "2"
        }

	output{
		File results = "${target}/bam/${sample}.sorted.bam" 
	}
} 

workflow pipeline {
	Int qualityScore
	String sample
	String target
	File adapterSeq

	call downloadRef {
		input:
			target=target
		} 

	call downloadSRA {
		input:
			sample=sample,
			target=target
		}

	call runFastQC { 
		input:
			sample=sample,
			target=target,
			rawFastQ=downloadSRA.rawFastQ
		}		

	call trimReads{
		input:
			sample=sample,
			target=target,
			qualityScore=qualityScore,
			adapterSeq=adapterSeq,
			rawFastQ=downloadSRA.rawFastQ
		}

	call align {
		input:
			sample = sample,
			target = target,
			trimmedFastQ = trimReads.trimmedFastQ,
			Ref = downloadRef.Ref
		}

	call samToBAM{
		input:
			sample = sample,
			target = target,
			alignedSAM = align.alignedSAM
		}
}


