import pandas as pd

METADATA = pd.read_table(config["METADATA"])
METADATA = METADATA[METADATA["Platform"] == "S2"]
METADATA = METADATA[METADATA["Tissue"] == "dorsolateral_prefrontal_cortex"]
METADATA = METADATA.reset_index(drop=True)

SampleGroups = METADATA.groupby("Condition").groups 

print(METADATA, SampleGroups)

workdir: config['workdir']



def getSamplesInGroup(var_group):
	return [METADATA.iloc[i]["SM"] for i in SampleGroups[var_group]]


rule all:
	input:
		expand("TAllele_new/{var_group}.merged.mi_summary.tab", var_group = list(SampleGroups))


rule combine_gvcfs:
	params:
		sge_opts = "-cwd -V -l h_data=16G,h_rt=2:00:00 \
			-N combine_gvcfs.{var_group} \
			-o log/combine_gvcfs.{var_group}.out \
			-e log/combine_gvcfs.{var_group}.err",
		TAlleleDir = config["TAlleleDir"]
	input:
		ref = config["REF"],
		vcf = lambda wildcards: \
			[f"TAllele_new/{sm}.vcf" for sm in getSamplesInGroup(wildcards.var_group)]
	output:
		fofn = "TAllele_new/{var_group}.gvcfs.list",
		vcf  = "TAllele_new/{var_group}.merged.gvcf.gz",
		gvcf = "TAllele_new/{var_group}.merged.genotyped.gvcf.gz"
	shell:
		"""
		echo -n > {output.fofn} 
		
		for vcf in {input.vcf} 
		do
			new_vcf=$(echo $vcf | sed 's/.vcf/.gvcf/g')
			bash {params.TAlleleDir}/convert_vcf.sh $vcf $new_vcf 
			echo $new_vcf >> {output.fofn}
		done

		cat {output.fofn}

		gatk --java-options "-Xmx4g" CombineGVCFs \
			--variant {output.fofn} \
			-R {input.ref} \
			-O {output.vcf} 
		
		gatk --java-options "-Xmx4g" GenotypeGVCFs \
			-R {input.ref} \
			-V {output.vcf} \
			-O {output.gvcf} 

		tabix -p vcf -f {output.gvcf}
		"""


rule tallele_joint:
	params:
		sge_opts = "-cwd -V -l h_data=14G,h_rt=1:00:00 \
			-N tallele_joint.{var_group} \
			-o log/tallele_joint.{var_group}.out \
			-e log/tallele_joint.{var_group}.err",
		sm = lambda wildcards: getSamplesInGroup(wildcards.var_group), 
		TAlleleDir = config["TAlleleDir"]
	input:
		vcf = "TAllele_new/{var_group}.merged.genotyped.gvcf.gz",
	output:
		fofn = "TAllele_new/{var_group}.mi.fofn",
		mi   = "TAllele_new/{var_group}.merged.mi_summary.tab"
	shell:
		"""
		echo -n > {output.fofn}

		for sm in {params.sm}
		do
			echo bam/$sm.annotated.sorted.bam TAllele_new/$sm.mi_summary.tab >> {output.fofn}
		done

		cat {output.fofn}

		python {params.TAlleleDir}/gqv_merge_outputs.py \
			-g {input.vcf} \
			-f {output.fofn} \
			-o {output.mi} 
		"""


