import pandas as pd

METADATA  = pd.read_table(config["METADATA"])
REF_INDEX = pd.read_table(config["REF"] + ".fai", header = None)

CHROMS = [chrom for chrom in REF_INDEX[0] if ("_" not in chrom)]

workdir: config['workdir']


rule all:
	input:
		["talon/talon_observedOnly.sorted.gtf.gz"] +
		expand("talon/{chrom}.gtf", chrom = CHROMS) ,

	
rule talon_initiate_db:
	params:
		sge_opts = "-cwd -V -l h_data=16G,h_rt=1:00:00 \
			-N talon_initiate_db \
			-o log/talon_initiate_db.out \
			-e log/talon_initiate_db.err",
		genome_assembly = config['genome_assembly'],
		annot_version   = config['annot_version'],
		dataset_id      = config['dataset_id']
	input:
		gtf = config["GTF"]
	output:
		gtf = expand("talon/{chrm}.gtf", chrm = CHROMS),
		db  = expand("talon/{chrm}.db", chrm = CHROMS),
	shell:
		"""
		mkdir -p talon

		for chrom in {CHROMS}
		do
			less {input.gtf} | grep -w $chrom > talon/$chrom.gtf

			talon_initialize_database \
				--f talon/$chrom.gtf \
				--g {params.genome_assembly} \
				--a {params.annot_version} \
				--idprefix {params.dataset_id} \
				--o talon/$chrom 
		done
		""" 


rule talon:
	params:
		sge_opts = "-cwd -V -l h_data=36G,h_rt=2:00:00,highp -pe shared 1 \
			-N talon.{chrom} \
			-o log/talon.{chrom}.out \
			-e log/talon.{chrom}.err",
		genome_assembly = config['genome_assembly'],
		annot_version   = config['annot_version'],
		dataset_id      = config['dataset_id'],
		prefix = "talon/{chrom}_out",
		tmpdir = "talon/{chrom}_tmp" 
	input:
		db  = "talon/{chrom}.db",
		sam = expand('bam/{sample}_labeled.sam', sample = METADATA['SM']),
	output:
		Config    = "talon/{chrom}.config",
		qc_log    = "talon/{chrom}_out/talon_QC.log",
		read_ann  = "talon/{chrom}_out_talon_read_annot.tsv",
	shell:
		"""
		mkdir -p {params.prefix}
		mkdir -p {params.tmpdir}

		echo -n > {output.Config}
		
		for sam in {input.sam}
		do
			chrom_sam=$sam.{wildcards.chrom}.sam

			# less $sam | grep "^@" > $chrom_sam
			# less $sam | awk '$3=="{wildcards.chrom}"' >> $chrom_sam
			
			sm=${{sam##*/}} 
			echo "$sm,$sm,PacBio,$chrom_sam" >> {output.Config}
		done
		
		talon \
			--f {output.Config} \
			--db {input.db} \
			--build {params.genome_assembly} \
			--threads 24 \
			--tmpDir {params.tmpdir} \
			--o {params.prefix}
		"""



rule talon_filter:
	params:
		sge_opts = "-cwd -V -l h_data=20G,h_rt=1:00:00,highp -pe shared 1 \
			-N talon_filter.{chrom} \
			-o log/talon_filter.{chrom}.out \
			-e log/talon_filter.{chrom}.err",
		genome_assembly = config['genome_assembly'],
		annot_version   = config['annot_version'],
		dataset_id      = config['dataset_id'],
		prefix = "talon/{chrom}_out",
		tmpdir = "talon/{chrom}_tmp" 
	input:
		db        = "talon/{chrom}.db",
		qc_log    = "talon/{chrom}_out/talon_QC.log",
	output:
		whitelist = "talon/{chrom}.whitelist.tsv", 
	shell:
		"""
		talon_filter_transcripts \
			--db {input.db} \
			-a {params.dataset_id} \
			--includeAnnot \
			--minCount 5 \
			--minDatasets 1 \
			--o {output.whitelist}
		"""


rule talon_generate_gtf:
	params:
		sge_opts = "-cwd -V -l h_data=12G,h_rt=3:00:00  \
			-N talon_generate_gtf.{chrom} \
			-o log/talon_generate_gtf.{chrom}.out \
			-e log/talon_generate_gtf.{chrom}.err",
		genome_assembly = config['genome_assembly'],
		annot_version   = config['annot_version'],
		prefix = "talon/{chrom}"
	input:
		db = "talon/{chrom}.db",
		whitelist = "talon/{chrom}.whitelist.tsv"
	output:
		gtf = "talon/{chrom}_talon_observedOnly.gtf",
	shell:
		"""
		talon_create_GTF \
			--db {input.db} \
			-a {params.annot_version} \
			--build {params.genome_assembly} \
			--whitelist {input.whitelist} \
			--observed \
			--o {params.prefix}
		"""


rule BuildExonicParts:
	params:
		sge_opts = "-cwd -V -l h_data=14G,h_rt=1:00:00 \
			-N BuildExonicParts \
			-o log/BuildExonicParts.out \
			-e log/BuildExonicParts.err",
		prefix = "talon/transcript_structures",
		TAlleleDir = config["TAlleleDir"],
	input:
		gtf = expand("talon/{chrm}_talon_observedOnly.gtf", chrm = CHROMS)
	output:
		transcript_str = "talon/transcript_structures/Gene_List.pickle",
		talon_comp_gtf = "talon/talon_observedOnly.sorted.gtf.gz",
		talon_comp_tbi = "talon/talon_observedOnly.sorted.gtf.gz.tbi",
	shell:
		"""
		mkdir -p {params.prefix} 
		
		cat {input.gtf} | \
			awk '$1 ~ /^#/ {{print $0;next}} {{print $0 | "sort -k1,1 -k4,4n -k5,5n"}}' | \
				bgzip -c > {output.talon_comp_gtf} 

		tabix {output.talon_comp_gtf}

		python {params.TAlleleDir}/gqv_process_talondb.py {output.talon_comp_gtf} {params.prefix}
		""" 

