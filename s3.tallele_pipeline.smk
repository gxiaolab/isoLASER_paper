import pandas as pd

METADATA  = pd.read_table(config["METADATA"])
# subset metadata to Platform == S1 
METADATA  = METADATA[METADATA["Platform"] == "S2"]
print(METADATA)

REF_INDEX = pd.read_table(config["REF"] + ".fai", header = None)

CHROMS = [chrom for chrom in REF_INDEX[0] if ("_" not in chrom)]

workdir: config['workdir']

rule all:
	input:
		expand("TAllele_new/{sample}.{suffix}", suffix = ["vcf", "mi_summary.tab"], sample = METADATA['SM'] )


rule MergeSams:
    params:
        sge_opts = "-cwd -V -l h_data=14G,h_rt=1:00:00 -N MergeSams.{sample} \
            -o log/MergeSams.{sample}.out \
            -e log/MergeSams.{sample}.err"
    input:
        sam = expand("bam/{{sample}}_labeled.sam.{chrm}.sam", chrm = CHROMS)
    output:
        tmp = "bam/{sample}.allchr_labeled.bam",
        bam = "bam/{sample}.labeled.bam",
        bai = "bam/{sample}.labeled.bam.bai"
    shell:
        """
        samtools merge -O BAM -f {output.tmp} {input.sam} 

        samtools sort {output.tmp} -o {output.bam} 
        
        samtools index {output.bam} 
        """


rule ReadAnnotator:
    params:
        sge_opts = "-cwd -V -l h_data=16G,h_rt=1:00:00,highp \
            -N ReadAnnotator.{sample} \
            -o log/ReadAnnotator.{sample}.out \
            -e log/ReadAnnotator.{sample}.err",
        TAlleleDir = config["TAlleleDir"],
    input:
        bam = "bam/{sample}.labeled.bam",
        tdb = directory("TALON_db"),
    output:
        tmp = "bam/{sample}.annotated.bam",
        bam = "bam/{sample}.annotated.sorted.bam",
        bai = "bam/{sample}.annotated.sorted.bam.bai"
    shell:
        """
        python {params.TAlleleDir}/filter_transgene_and_talon_annotate.py \
            --input-bam  {input.bam} \
            --output-bam {output.tmp} \
            --talon-db   {input.tdb}  

        samtools sort {output.tmp} -o {output.bam} 

        samtools index {output.bam}
        """

		
rule TAllele:
    params:
        sge_opts = "-cwd -V -l h_data=24G,h_rt=4:00:00 -pe shared 1\
            -N TAllele.{sample} \
            -o log/TAllele.{sample}.out \
            -e log/TAllele.{sample}.err",
        prefix = "TAllele_new/{sample}",
        TAlleleDir = config["TAlleleDir"],
    input:
        bam = "bam/{sample}.annotated.sorted.bam",
        tx  = 'TALON_db/transcript_structures',
        ref = config["REF"]
    output:
        vcf = "TAllele_new/{sample}.vcf",
        mis = "TAllele_new/{sample}.mi_summary.tab"
    shell:
        """
        mkdir -p TAllele_new 

        python {params.TAlleleDir}/gqv.py \
            -b {input.bam} \
            -o {params.prefix} \
            -g {input.ref} \
            -t {input.tx} \
            -n 16

        less {params.prefix}.mi_summary.tab | grep "HAP-intronic_part" > {params.prefix}.mi_summary.hap_only.tab 
        """

                               
