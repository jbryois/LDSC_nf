#!/usr/bin/env nextflow
 
// LDSC necessary files (European population)
plink = params.LDSC_files + "1000G_EUR_Phase3_plink/1000G.EUR.QC."
frq = params.LDSC_files + "1000G_Phase3_frq/1000G.EUR.QC."
weights = params.LDSC_files + "1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC."
hmsnps = params.LDSC_files + "hm_snp.txt"

// By default, the pipeline will run on all bed files in the folder where the nextflow pipeline is located
// This can be changed using the --bed modifier (e.g. --bed /Users/home/myBedFolder/)
params.bed = "./"
bedfiles = Channel.fromPath(params.bed + "*.bed")

// Output directory. By default LDSC files will be outputed in a folder named "output"
// This can be changed using the --out modifier (e.g. --out LDSC_results)
params.out = "output"

// Phenotype file, comma separated with name of phenotype in first column and path in second column, no header!
// (e.g. scz,/path/to/LDSC/sumstats/scz.sumstats.gz)
// Location of sumstats file can be modifed using the --pheno modifier (e.g. --pheno myphenos.csv)
params.pheno = "pheno.csv"
pheno_file = file(params.pheno)

// Indicator for tissue association analysis (false by default)
params.tissue = false

// By default, the pipeline will use the Finucane et al. 2015 annotation
// The model from Gazal et al. 2018 can also be used.
params.model = "Finucane"
if (params.model == "Finucane"){
	model_path = params.LDSC_files + "1000G_EUR_Phase3_baseline/baseline."
	analysis_type = "Finucane"
} else if (params.model == "Gazal"){
	model_path = params.LDSC_files + "baselineLD_v1.1/baselineLD."
	analysis_type = "Gazal"
} else {
	error("Please input a valid baseline model: Finucane or Gazal")
	}

// If the tissue modifier is entered, the name of the model will be updated
if (params.tissue){
	analysis_type = params.model + "_tissue"
}

// Info displayed at runtime
log.info """\
 ===================================

 Usage: 
 nextflow run LDSC_nf
  
 Options:
 model: --model Finucane (default) or Gazal
 tissue: --tissue (to run tissue association analysis) (default off)
 bed: --bed /mybeds/location/ (default, current directory) (bed files should be tab delimited with no header!) 
 output directory: --out output (default, output folder current directory)
 phenotype file: --pheno pheno.csv (default name) (two columns: id in the first column (no underscore), sumstats path in the second) (No header)
 LDSC path: --LDSC_files /path/to/LDSC
 
 ===================================
 L D S C       -   P I P E L I N E
 ===================================
 analysis_type	: ${analysis_type}
 tissue mode	: ${params.tissue}
 bed path	: ${params.bed}
 out		: ${params.out}
 phenotypes	: ${params.pheno}
 LDSC path	: ${params.LDSC_files}
 LDSC plink	: ${plink}
 LDSC frq	: ${frq}
 LDSC weights	: ${weights}
 hm_snps	: ${hmsnps}
 ===================================

 To download the latest version of the pipeline run:
 nextflow pull jbryois/LDSC_nf

 ===================================
 """

/* 
 *
 *	CODE BELOW SHOULD NOT NEED TO BE MODIFIED TO RUN THE PIPELINE!
 *
 */
 
// Creates a channel for the different phenotypes to be tested
// The phenotype name should not contain an underscore
Channel
    .fromPath(pheno_file)
    .splitCsv(header: false)
    .set { ch_pheno }

// Splits Input bed files per chromosome (1-22). Outputs the inputname and bed files for each chromosome
process Split_input_beds_per_chr {

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
process Get_LD_scores {
	
	input:
	set val(inputname), file(input_bed), file(baseline_annot), file(plink_bim), file(plink_bed), file(plink_fam), file('hm_snp.txt') from ch_chr
	
	output:
	set inputname, file('input*') into LDscores
	
	script:
	outname = baseline_annot.toString().replaceAll(/baseline(LD)?/,'input')
	outname = outname.replaceAll(/.annot.gz/,'')
	tissue = params.tissue
	"""
	zcat < $baseline_annot | awk -v OFS='\t' 'NR>1 {print "chr"\$1,\$2,\$2,\$3}' | gzip > $baseline_annot".bed"
	intersectBed -c -a $baseline_annot".bed" -b $input_bed > $input_bed".1000genomes.intersect"
	awk '{if(\$5!=0) print \$4}' $input_bed".1000genomes.intersect" > $input_bed".1000genomes.intersect.snp"

 	if [ $tissue = "false" ]
 	then
		awk -v OFS='\t' '{print \$1,\$2-500,\$3+500,\$1":"\$2"-"\$3}' $input_bed > $input_bed".500bp.ext"
	   	awk -v OFS='\t' '{ if (\$2<0) print \$1,0,\$3,\$4; else print \$0}' $input_bed".500bp.ext" > tmp
		mv tmp $input_bed".500bp.ext"
		intersectBed -c -a $baseline_annot".bed" -b $input_bed".500bp.ext" > $input_bed".1000genomes.500bp.ext.intersect"
		awk '{if(\$5!=0) print \$4}' $input_bed".1000genomes.500bp.ext.intersect" > $input_bed".1000genomes.500bp.ext.intersect.snp"
		fast_match2_minimal.gz.pl $input_bed".1000genomes.intersect.snp" $inputname $baseline_annot > tmp
		fast_match2.gz.pl $input_bed".1000genomes.500bp.ext.intersect.snp" $inputname"_500bp_ext" tmp | gzip > $outname".annot.gz"
		rm tmp
		rm *500bp.ext
	 else
		fast_match2_minimal.gz.pl $input_bed".1000genomes.intersect.snp" $inputname $baseline_annot | gzip > $outname".annot.gz"
	fi
	ldsc.py --l2 --bfile $plink_bim.baseName --ld-wind-cm 1 --print-snps hm_snp.txt --annot $outname".annot.gz" --out $outname
	rm *.intersect
	rm *.intersect.snp
	rm $baseline_annot".bed"
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
 
process Get_pLDSC_Results {
	publishDir "${params.out}/$analysis_type/Results_raw/", mode: 'copy', overwrite: true

	input:
	set inputname, path(inputLDscores),pheno,path(sumstats) from LDscores_join
	
	output:
	file('*results') into (Results_ch_h2,Results_ch_tissue)
	
	"""
	ldsc.py --h2 $sumstats --ref-ld-chr $model_path,input. --w-ld-chr $weights --overlap-annot --frqfile-chr $frq --print-coefficients --out ${pheno}_${inputname}
	"""
}

// For enrichment analysis - Gathers the results files and compute final file

process Clean_results_h2_enrichment {
    publishDir "${params.out}/$analysis_type/", mode: 'copy', overwrite: true
    
    input:
	path(results) from Results_ch_h2.collect()
	
	output:
	file('*.txt') into CleanResultsH2

	when:
	!params.tissue
	
	"""
	CleanH2Results.R
	"""
}


// For tissue specific analysis - Gathers the results files, compute pvalue and produce final results file

process Clean_results_tissue_association {
    publishDir "${params.out}/$analysis_type/", mode: 'copy', overwrite: true
    
	input:
	path(results) from Results_ch_tissue.collect()
	
	output:
	file('*.txt') into CleanResultsZ

	when:
	params.tissue
	
	"""
	CleanTissueResults.R
	"""
}

