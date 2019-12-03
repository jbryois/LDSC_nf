# LDSC_nf

Pipeline to run partitioned LDSC on user specified [bed](https://www.ensembl.org/info/website/upload/bed.html) files with [Nextflow](https://www.nextflow.io).
  
# Requirements

1. Nextflow

[Nextflow](https://www.nextflow.io) requires java 8

Make sure that it's installed using the following command:

```
java -version
```

You can then install nextflow using the following command:

```
curl -s https://get.nextflow.io | bash
```

This creates a file named nextflow in the current directory. It can then be moved to your path.

For example:

```
mv nextflow $HOME
```

2. LDSC files

LDSC files can be installed in a specified directory using the following commands:

```

wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_weights_hm3_no_MHC.tgz 
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_ldscores.tgz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_frq.tgz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baselineLD_v1.1_ldscores.tgz

tar -xvzf 1000G_Phase3_weights_hm3_no_MHC.tgz
tar -xvzf 1000G_Phase3_plinkfiles.tgz
tar -xvzf 1000G_Phase3_baseline_ldscores.tgz
tar -xvzf 1000G_Phase3_frq.tgz
tar -xvzf 1000G_Phase3_baselineLD_v1.1_ldscores.tgz
bunzip2 w_hm3.snplist.bz2

tail -n +2 w_hm3.snplist | cut -f 1 > hm_snp.txt
```

# Usage

First you need to pull the latest version from github

```
nextflow pull jbryois/LDSC_nf
```

You can then use the pipeline with the following options:

```
Usage: 
 nextflow run LDSC_nf --LDSC_files /path/to/LDSC
  
 Options:
 model: --model Finucane (default) or Gazal
 tissue: --tissue (to run tissue association analysis) (default off)
 bed: --bed /mybeds/location/ (default, current directory) (bed files should be tab delimited with no header!) 
 output directory: --out output (default, output folder current directory)
 phenotype file: --pheno pheno.csv (default name) (two columns: id in the first column (no underscore), sumstats path in the second) (No header)
 LDSC path: --LDSC_files /path/to/LDSC
``` 

It might be worth to modify the file nextflow.config to add the path of the LDSC files so that you don't have to specify the path of the LDSC installation files everytime

```
params.LDSC_files = "/path/to/my/LDSC/files/"
```
