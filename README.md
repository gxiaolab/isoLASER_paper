# **Manuscript Preparation**

[link to publication](https://github.com/gxiaolab/isoLASER/)

## **About**

isoLASER performs gene-level variant calls, phasing, and splicing linkage analysis using third-generation RNA sequencing data.

## **Table of contents**
- [Software](#software)
- [Pipeline](#pipeline)

## **Software**

Documentation for the isoLASER software [here.](https://github.com/gxiaolab/isoLASER)

## **Pipeline**

The *snakemake* files show the workflow used in our manuscript to process and analyze the sequencing data. 

- s1: Preprocessing of bam files
- s2: Isoform prediction using TALON
- s3: Running isoLASER
- s4: Running isoLASER joint
