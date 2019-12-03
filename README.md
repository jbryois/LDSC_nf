# LDSC_nf

Pipeline to run partitioned LDSC on user specified input files.
 
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

This creates a file named nextflow




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