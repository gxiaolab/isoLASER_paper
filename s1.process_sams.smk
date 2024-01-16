import pandas as pd

METADATA  = pd.read_table(config["METADATA"])
REF_INDEX = pd.read_table(config["REF"] + ".fai", header = None)

workdir: config['workdir']


rule all:
	input:
		expand("bam/{sample}_labeled.sam", sample = METADATA['SM'] )


rule get_splice_juncs:
	params:
		sge_opts = "-cwd -V -l h_data=16G,h_rt=1:00:00 \
			-N get_splice_juncs \
			-o log/get_splice_juncs.out \
			-e log/get_splice_juncs.err",
		TrCleanDir = config["TranscriptCleanDir"]
	input:
		gtf = config["GTF"],
		ref = config["REF"],
	output:
		jxn = config["annot_version"] + ".SJ.tab"
	shell:
		"""
		mkdir -p log 

		python {params.TrCleanDir}/TranscriptClean-2.0.3/accessory_scripts/get_SJs_from_gtf.py \
			--f {input.gtf} \
			--g {input.ref} \
			--o {output.jxn} 		
		"""


rule sort_index_bam:
	params:
		sge_opts = "-cwd -V -l h_data=8G,h_rt=1:00:00 \
			-N sort_index_bam \
			-o log/sort_index_bam.{sample}.err \
			-e log/sort_index_bam.{sample}.out" 
	input:
		bam = "raw_bam/{sample}.bam",
	output:
		bam = "bam/{sample}.sorted.bam",
		bai = "bam/{sample}.sorted.bam.bai"
	shell:
		"""
		mkdir -p bam 

		samtools sort -o {output.bam} {input.bam}
		
		samtools index {output.bam} 
		"""



rule transcript_clean:
	params:
		sge_opts = "-cwd -V -l h_data=14G,h_rt=7:00:00,highp -pe shared 2 \
			-N transcript_clean \
			-o log/transcript_clean.{sample}.out \
			-e log/transcript_clean.{sample}.err",
		prefix = "bam/{sample}",
		tmpdir = "tmp/{sample}",
		threads = 4,
		TrCleanDir = config["TranscriptCleanDir"]
	input:
		bam = "bam/{sample}.sorted.bam",
		ref = config["REF"],
		fai = config["REF"] + ".fai",
		jxn = config["annot_version"] + ".SJ.tab"
	output:
		sam = "bam/{sample}.sorted.sam",
		tc_cleaned_sam    = 'bam/{sample}_clean.sam',
		tc_cleaned_fa     = 'bam/{sample}_clean.fa',
		tc_cleaned_te_log = 'bam/{sample}_clean.TE.log',
		tc_cleaned_log    = 'bam/{sample}_clean.log'
	shell:
		"""
		samtools view -h {input.bam} $(cut -f1 {input.fai}) | samtools calmd - {input.ref} > {output.sam}

		mkdir -p {params.tmpdir}

		python  {params.TrCleanDir}/TranscriptClean-2.0.4/TranscriptClean.py \
			--correctMismatches False \
			--correctIndels False \
			--threads {params.threads} \
			--sam {output.sam} \
			--genome {input.ref} \
			--spliceJns {input.jxn} \
			--outprefix {params.prefix} \
			--tmpDir {params.tmpdir} 
		"""


rule talon_label_reads:
	params:
		sge_opts = "-cwd -V -l h_data=8G,h_rt=1:00:00 \
			-N talon_label_reads \
			-o log/talon_label_reads.{sample}.out \
			-e log/talon_label_reads.{sample}.err ",
		prefix = "bam/{sample}",
		tmp = "tmp/{sample}"
	input:
		sam = "bam/{sample}_clean.sam",
		ref = config["REF"]
	output:
		sam = 'bam/{sample}_labeled.sam',
		tsv = 'bam/{sample}_read_labels.tsv'
	shell:
		"""
		mkdir -p {params.tmp}

		talon_label_reads \
			--f {input.sam} \
			--g {input.ref} \
			--o {params.prefix} \
			--tmpDir {params.tmp} 
		"""
