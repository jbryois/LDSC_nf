#!/usr/bin/env nextflow
// Usage: nextflow run ldsc_pipeline.nf 
 
/* To do
 * 1) Add option to add 500bp window or not
 * 2) Add option for tissue specific analysis (without 500bp window, Rscript to gather results and compute pvalue from coefficient z-scores)
 * 3) Add option for continuous phenotype analysis
 * 4) Test on SLURM cluster
 */
 
// LDSC core files path. This can be changed using the --LDSC_files modifier
params.LDSC_files = "/Users/julienbryios/Documents/Data/Projects/LDSC_nf/"

// LDSC necessary files (European population)
plink = params.LDSC_files + "1000G_EUR_Phase3_plink/1000G.EUR.QC."
frq = params.LDSC_files + "1000G_Phase3_frq/1000G.EUR.QC."
weights = params.LDSC_files +"1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."

// By default, the pipeline will run on all bed files in the folder where the nextflow pipeline is located
// This can be changed using the --bed modifier (e.g. --bed /Users/home/myBedFolder/)
params.bed = "./"
bedfiles = Channel.fromPath(params.bed + "*.bed")

// Output directory. By default LDSC files will be outputed in a folder named "output"
// This can be changed using the --outputDir modifier (e.g. --outputDir LDSC_results)
params.outputDir = "output"

// Phenotype file, tab separated with name of phenotype in first column and path in second column
// This can be changed using the --pheno modifier (e.g. --pheno mysumstats.txt)
// The file should contain an identifier and the location of the corresponding sumstats in a tab delimited files (e.g. scz /path/to/LDSC/sumstats/scz.sumstats.gz)
params.pheno = "pheno.txt"

// By default, the pipeline will use the Finucane et al. 2015 annotation
// The model from Gazal et al. 2018 can also be used.
params.model = "Finucane"
if (params.model == "Finucane"){
	model_path = params.LDSC_files + "1000G_EUR_Phase3_baseline/baseline."
} else if (params.model == "Gazal"){
	model_path = params.LDSC_files + "baselineLD_v1.1/baselineLD."
} else {
	error("Please input a valid baseline model: Finucane or Gazal")
	}

log.info """\
 LDSC - N F   P I P E L I N E
 ===================================
 LDSC path	: ${params.LDSC_files}
 LDSC plink	: ${plink}
 LDSC frq	: ${frq}
 LDSC weights	: ${weights}
 bed path	: ${params.bed}
 outputDir	: ${params.outputDir}
 phenotypes	: ${params.pheno}
 hm_snps	: hm_snp.txt
 """


/* 
 *
 *	Code below should not need to be modified to run the pipeline!
 *
 */

// Creates a channel for the different phenotypes to be tested 
Channel
    .fromPath(params.pheno)
    .splitText()
    .map{ line -> tuple(line.split('\t')[0],line.split('\t')[1]) }
    .set { ch_pheno }


// Splits Input bed files per chromosome (1-22). Outputs the inputname and bed files for each chromosome
process SplitInputBedsPerChr {

	input:
	file bed from bedfiles

	output:
	set inputname, file('chr*') into InputBedsPerChr

	script:
  	inputname = bed.toString() - ~/.bed$/
	"""
	awk '\$1 ~ /^chr(1?[0-9]|2[0-2])/ {print \$0 >> \$1;close(\$1)}' < $bed
	"""
	}
	
// Outputs sample name, bedfiles for each chromosome, annotation files and plink files (adapted from: https://github.com/nextflow-io/nextflow/issues/796)
InputBedsPerChr
    .flatMap { item ->
        inputname = item[0];
        files  = item[1];
        files.collect { onefile -> \
        return [ inputname, \
        onefile, \
        file(model_path + "${onefile.simpleName.replaceFirst(/chr/,'')}" + ".annot.gz"), \
        file(plink + "${onefile.simpleName.replaceFirst(/chr/,'')}" + ".bim"), \
        file(plink + "${onefile.simpleName.replaceFirst(/chr/,'')}" + ".bed"), \
        file(plink + "${onefile.simpleName.replaceFirst(/chr/,'')}" + ".fam"), \
        file('hm_snp.txt')]}
    }
    .set { ch_chr }

/* Gets extended bed annotation for the input bed chromosome file (+-500bp), gets bed file for the annotation file (same chr), gets overlapping SNPs
 * between the input bed chromosome and the annotation, add the overlapping SNPs to the annotation file, compute LD scores for the new annotation
 */
process getLDscores {
    publishDir "${params.outputDir}/$inputname/$params.model/LDscores"
	input:
	set val(inputname), file(input_bed), file(baseline_annot), file(plink_bim), file(plink_bed), file(plink_fam), file('hm_snp.txt') from ch_chr
	
	output:
	set inputname, file('input*') into LDscores
	
	script:
	outname = baseline_annot.toString().replaceAll(/baseline(LD)?/,'input')
	outname = outname.replaceAll(/.annot.gz/,'')
	"""
	awk -v OFS='\t' '{print \$1,\$2-500,\$3+500,\$1":"\$2"-"\$3}' $input_bed > $input_bed".500bp.ext"
    awk -v OFS='\t' '{ if (\$2<0) print \$1,0,\$3,\$4; else print \$0}' $input_bed".500bp.ext" > tmp
	mv tmp $input_bed".500bp.ext"
	zcat < $baseline_annot | awk -v OFS='\t' 'NR>1 {print "chr"\$1,\$2,\$2,\$3}' | gzip > $baseline_annot".bed"
	intersectBed -c -a $baseline_annot".bed" -b $input_bed > $input_bed".1000genomes.intersect"
	intersectBed -c -a $baseline_annot".bed" -b $input_bed".500bp.ext" > $input_bed".1000genomes.500bp.ext.intersect"
	awk '{if(\$5!=0) print \$4}' $input_bed".1000genomes.intersect" > $input_bed".1000genomes.intersect.snp"
	awk '{if(\$5!=0) print \$4}' $input_bed".1000genomes.500bp.ext.intersect" > $input_bed".1000genomes.500bp.ext.intersect.snp"
	fast_match2_minimal.gz.pl $input_bed".1000genomes.intersect.snp" $inputname $baseline_annot > tmp
	fast_match2.gz.pl $input_bed".1000genomes.500bp.ext.intersect.snp" $inputname"_500bp_ext" tmp | gzip > $outname".annot.gz"
	ldsc.py --l2 --bfile $plink_bim.baseName --ld-wind-cm 1 --print-snps hm_snp.txt --annot $outname".annot.gz" --out $outname
	rm *.intersect
	rm *.intersect.snp
	rm $baseline_annot".bed"
	rm tmp
	rm *500bp.ext
	rename s'/.augmented//' *
	"""
}

/* Prepare a channel combining each input bed, all LD scores files and the phenotypes to be tested
 */
LDscores
.groupTuple()
.map { mytuple ->
	 def key = mytuple[0]
     def myfiles = mytuple[1].flatten()
     return tuple(key.toString(), myfiles) 
     }
.combine(ch_pheno)
.set { LDscores_join }

/* Gather all LD files for each inputBed and and phenotypes and run partitioned LD score regression for each input file and each phenotype
 */
process GetPhenotypeEnrichment {
    publishDir "${params.outputDir}/$inputname/$params.model/Results/"
	
	input:
	set inputname, path(inputLDscores),pheno,path(sumstats) from LDscores_join
	
	output:
	set inputname, path('*results') into Results

	"""
	ldsc.py --h2 $sumstats --ref-ld-chr $model_path,input. --w-ld-chr $weights --overlap-annot --frqfile-chr $frq --print-coefficients --out ${pheno}_${inputname}
	"""
}
 