# PTATO
[![DOI](https://zenodo.org/badge/485723666.svg)](https://zenodo.org/badge/latestdoi/485723666)


[![DOI](https://zenodo.org/badge/485723666.svg)](https://zenodo.org/badge/latestdoi/485723666)

The PTA Analysis TOolbox (PTATO) is a comprehensive pipeline designed to filter somatic single base substitutions (SBS), small insertions and deletions (indels) and structural variants (SVs) from PTA-based single-cell whole genome sequencing (WGS) data. More information about the pipeline can be found in the [manuscript](https://www.biorxiv.org/content/10.1101/2023.02.15.528636v1). *Please cite the manuscript if you use PTATO.*

## Dependencies
PTATO is implemented in Nextflow and required installation of the following dependencies:
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
- [Apptainer](https://apptainer.org/)
- [R/4.1.2](https://cran.r-project.org/index.html)

#### R libraries
- [BSgenome](http://bioconductor.org/packages/release/bioc/html/BSgenome.html)
- [copynumber](https://bioconductor.org/packages/release/bioc/html/copynumber.html)
- [cowplot](https://cran.r-project.org/web/packages/cowplot/index.html)
- [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
- [gtools](https://cran.r-project.org/web/packages/gtools/index.html)
- [MutationalPatterns](https://bioconductor.org/packages/release/bioc/html/MutationalPatterns.html)
- [randomForest](https://cran.r-project.org/web/packages/randomForest/index.html)
- [scales](https://scales.r-lib.org/)
- [StructuralVariantAnnotation](https://www.bioconductor.org/packages/release/bioc/html/StructuralVariantAnnotation.html)
- [VariantAnnotation](https://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html)

## Easy Installation (Linux) (recommended)
Download Singularity image SIF (~4gb) Required
[singularity/Apptainer](https://apptainer.org/docs/admin/main/index.html) v1.1.3
```
# 1. Pull singularity image from Docker bootstrap
singularity pull ptato_1.2.0.sif docker://vanboxtelbioinformatics/ptato:1.2.0

# 2. Clone PTATO repository 
git clone git@github.com:ToolsVanBox/PTATO.git

# 3. Run with singularity exec
singularity exec ptato_1.2.0.sif /ptato/nextflow/nextflow run \
[PTATO_dir]/ptato.nf \
-c [PTATO_dir]configs/run_template.config \
-profile slurm -resume

```

## Installing & Setup

1. [Install Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
2. [Install Singularity](https://sylabs.io/guides/3.5/admin-guide/)
3. Install R and the required R libraries
4. [Pull/Clone PTATO](#pull-or-clone)
5. Download and extract the required resource files


## Pull or Clone 
```
git clone git@github.com:ToolsVanBox/PTATO.git
```

## Demo dataset

A test dataset containing all required input files (e.g. BAM, VCF and config files) to run PTATO is available for [download here](https://surfdrive.surf.nl/files/index.php/s/l3FX6eLnTtuVK1g). This demo dataset contains bulk whole genome sequencing (WGS) data of a clonal cell line and a PTA dataset from a single cell derived from this clone. The bulk WGS data can be used as germline control sample.

## Resources
Most required resource files (for the hg38 reference genome) are already included in the PTATO repository. Only the reference genome and the SHAPEIT resources need to be downloaded seperately. First extract the following resources files:
- `[PTATO_dir]/resources/hg38/gripss/gridss_pon_breakpoint.tar.gz`
- `[PTATO_dir]/resources/hg38/cobalt/COBALT_PTA_Normalized_Full.tar.gz`
- `[PTATO_dir]/resources/hg38/smurf/Mutational_blacklists/Fetal_15x_raw_variants_hg38.tar.gz`
- `[PTATO_dir]/resources/hg38/smurf/Mutational_blacklists/MSC_healthyBM_raw_variants_hg38.tar.gz`

#### Reference genome
Please download the reference genome fasta file. Must have the following files:
- *.dict
- *.fasta
- *.fasta.fai

Recommended files are, otherwise they will be created by the pipeline:
- *.fasta.amb
- *.fasta.ann
- *.fasta.bwt
- *.fasta.dict
- *.fasta.gridsscache
- *.fasta.pac
- *.fasta.sa

And put them in this folder `[PTATO_dir]/resources/hg38/`

#### SHAPEIT resources
Download the following SHAPEIT resource files:
- http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/
- https://github.com/odelaneau/shapeit4/blob/master/maps/genetic_maps.b38.tar.gz

Unzip the genetic_maps.b38.tar.gz by using the following command
```
tar -zxf genetic_maps.b38.tar.gz
```

And put them in this folder, respectively:
- `[PTATO_dir]/resources/hg38/shapeit/Phasing_reference/`
- `[PTATO_dir]/resources/hg38/shapeit/shapeit_maps/`

## Running the workflow
To run the PTATO workflow, the following steps have to be performed:
1. Collect input files 
2. Change the run-template.config
3. Start the pipeline

### 1. Collect input files 
PTATO requires the following input files:
1. BAM files 
2. VCF file (multi-sample) generated by GATK

Optionally, to generate basic quality control plots, PTATO also makes use of the following files:

3. wgs_metrics.txt
4. multiple_metrics.alignment_summary_metrics
- *If you do not have the wgs_metrics and alignment_summary_metrics files, change `QC = true` to `QC = false` in the run.config file.*

#### 1. BAM files
- PTATO requires WGS BAM files of at least one germline control sample and one PTA-based test sample.
- The germline control sample is preferably a bulk WGS sample. In our experience, using PTA-based WGS samples as germline control leads to insufficient filtering of germline variants.
- PTATO will only run for the samples for which a BAM file is included in the BAM directory. Sample names are also derived from the BAM file names.

#### 2. Multisample VCF file generated by GATK
- A pre-generated multi-sample VCF file is used as starting point for filtering by PTATO (raw calling of base substitutions and indels is **not** part of the PTATO pipeline). PTATO has been optimized for multi-sample VCFs that were generated using GATK's [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller).
- Importantly, the input VCF should contain both somatic as well as germline variants. The germline variants are required to calculate allelic imbalance, to perform linked read analysis and to perform B-allele frequency analysis.
- All samples specified in the input directory of the BAM files should be present in the input VCF (with the same sample names). 

#### 3. WGS Metrics files *(Optional)*
- PTATO uses `wgs_metrics.txt` to get an overview of the genome coverage in the input samples. Each sample requires its own seperate txt file.
- WGS metrics files can be generated using [CollectWgsMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037269351-CollectWgsMetrics-Picard-) from Picard.

#### 4. Alignment Summary Metrics *(Optional)*
- PTATO uses `alignment_summary_metrics` files to generate an overview of the number of sequencing reads in the input samples. Each sample requires its own seperate summary file.
- alignment_summary_metrics can be generated using [CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360040507751-CollectAlignmentSummaryMetrics-Picard-) from Picard.

#### Input folder structure
Currently, PTATO requires a strict structure of input directories (eg the bam files should be placed in a subdirectory with the name of the individual/donor/patient). It is possible to use links to the original files (bam/vcf), as long as these links are in the appropriate folder structure. The paths to these directories should be included in the `run.config` file (see below). The input files listed above should be structured as follows:
```
/path/to/vcfs_dir
  ./Donor_1
    ./myfile.vcf(.gz)
/path/to/bams_dir
  ./Donor_1
    ./mycontrol.bam
    ./mysample1.bam
    ./mysample2.bam
    ...
/path/to/wgs_metrics
  ./Donor_1
    ./wgs_metrics1.txt
    ./wgs_metrics2.txt
    ...
/path/to/alignment_summary_metrics
  ./Donor_1
    ./alignment_summary_metrics1
    ./alignment_summary_metrics2
    ...
```


### 2. Change the run-template.config

PTATO uses 4 different config files (templates provided in the `[PTATO_dir]/configs/` directory):
1. `run.config`
2. `process.config`
3. `nextflow.config`
4. `resources.config`

The `run.config` needs to be adjusted for each new PTATO run. The `process.config` may have to be changed if necessary. The `nextflow.config` has to be changed once (and can be reused for later runs), to tailor the settings specific for your compute cluster.

#### 1. run.config
The run.config contains the paths to the input files and therefore needs to be adapted for each PTATO run. 

The first three lines of the config file should contain the paths to the other three config files, as follows:

```
includeConfig "${projectDir}/configs/process.config"
includeConfig "${projectDir}/configs/nextflow.config"
includeConfig "${projectDir}/configs/resources.config"
```
*Change the paths if you use different versions of the config files, stored at different locations (for example if you have a separate process.config file for a specific PTATO run*

All of the parameters in the params section can also be supplied on the commandline or can be pre-filled in the run.config file:
```
includeConfig "${projectDir}/configs/process.config"
includeConfig "${projectDir}/configs/nextflow.config"
includeConfig "${projectDir}/configs/resources.config"

params {

  run {
    snvs = true
    QC = true
    svs = false
    indels = true
    cnvs = false
  }

  // TRAINING
  train {
    version = '2.0.0'
  }
  pta_vcfs_dir = ''
  nopta_vcfs_dir = ''
  // END TRAINING

  // TESTING
  input_vcfs_dir = ''
  bams_dir = ''
  // END TESTING

  out_dir = ''
  bulk_names = [
    ['donor_id', 'sample_id'],
  ]

  snvs {
    rf_rds = "${projectDir}/resources/hg38/snvs/randomforest/randomforest_v1.0.0.rds"
  }

  indels {
    rf_rds = ''
    excludeindellist = "${projectDir}/resources/hg38/indels/excludeindellist/PTA_Indel_ExcludeIndellist_normNoGTrenamed.vcf.gz"
  }
  optional {

    germline_vcfs_dir = ''

    short_variants {
      somatic_vcfs_dir = ''
      walker_vcfs_dir = ''
      phased_vcfs_dir = ''
      ab_tables_dir = ''
      context_beds_dir = ''
      features_beds_dir = ''
    }

    snvs {
      rf_tables_dir = ''
      ptato_vcfs_dir = ''
    }

    indels {
      rf_tables_dir = ''
      ptato_vcfs_dir = ''
    }

    qc {
      wgs_metrics_dir = ''
      alignment_summary_metrics_dir = ''
    }

    svs {
      gridss_driver_vcfs_dir = ''
      gridss_unfiltered_vcfs_dir = ''
      gripss_somatic_filtered_vcfs_dir = ''
      gripss_filtered_files_dir = ''
      integrated_sv_files_dir = ''
    }

    cnvs {
      cobalt_ratio_tsv_dir = ''
      cobalt_filtered_readcounts_dir = ''
      baf_filtered_files_dir = ''
    }
  }

}
```

- Under header `run { }` you can specify which parts of PTATO you would like to run (set to  `= true`). For example, if you don't want to run SV calling and filtering, you can set `svs = false`. Note: the `snvs = true` and `cnvs = true` parts of PTATO are required to run the `svs = true` part.
- Under header `// TESTING` you have to specify the paths to the input directories containing the VCF file (`input_vcfs_dir = '/path/to/vcf/'`) and the BAM files (`bams_dir = '/path/to/bams/'`). Please not that the name of the individual/donor should not be included in the path (eg NOT: `'/path/to/vcfs_dir/donor/'`)
```
  // TESTING
  input_vcfs_dir = '/path/to/vcfs_dir/'
  bams_dir = '/path/to/bams_dir/'
  // END TESTING

```
- Under header "bulk_names" you have to specify the name of the individual/donor/patient and the sample_id of the germline control sample. Mutations in this control sample are used to determine which variants are germline or somatic. Mutations in the control sample are excluded from the somatic variants. Multiple control samples can be specified by adding an additional row:
```
  bulk_names = [
    ['donor_id', 'control_sample1'],
    ['donor_id', 'control_sample2'],
  ]
```
- If you would like to run the sequencing QC, you have to specify the paths to the directories containing the wgs_metrics and alignment_summary_metrics files here:
```
    qc {
      wgs_metrics_dir = '/path/to/wgs_metrics_dir'
      alignment_summary_metrics_dir = '/path/to/alignment_summary_metrics_dir'
    }
```
- All other fields (that are empty in the example run.config file) are optional and can be left empty. These files will be generated by PTATO. If you would like to rerun parts of PTATO later, you can specify the files that were previously generated by PTATO. The old files will then be re-used, which saves time and resources.

#### 2. process.config
The process.config contains the general settings for each type of job. Here you can for example change the time and memory that are reserved for each job. This likely requires some tweaking (and trial and error) for your specific compute cluster setup. The required time and memory also depend on the number of samples you would like to include in your PTATO run. For example, if you would like to run PTATO on 10+ samples, you would likely need to increase the time for the somatic variant filtering (eg change `params.smurf.time = '4h'` to `params.smurf.time = '12h'`).

#### 3. nextflow.config
The nextflow.config has to be changed once to specify the base directory and cache dir for your cluster. Specifically, only this part needs to be changed:
```
singularity {
  enabled = true
  autoMounts = true
  runOptions = '-B /hpc -B $TMPDIR:$TMPDIR'
  cacheDir = '/hpc/local/CentOS7/pmc_vanboxtel/singularity_cache'
}
```
- Change `/hpc` in `runOptions` to the base directory of your cluster
- Change the path in `cacheDir` to a cache directory on your cluster


### 3. Start the pipeline
Once you have collected all the input files and changed the required config files you can start the PTATO pipeline.
To start the pipeline on a Slurm workload manager:

```
/path/to/nextflow run /path/to/ptato.nf -c /path/to/run.config --out_dir /path/to/output_directory -profile slurm -resume
```

## Acknowledgements and References
Also see the references in the [manuscript](https://www.biorxiv.org/content/10.1101/2023.02.15.528636v1). PTATO includes the following external software:
- [GRIDSS2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02423-x) and [GRIPSS](https://github.com/hartwigmedical/hmftools/blob/master/gripss/README.md): Cameron, D.L., Baber, J., Shale, C., Valle-Inclan, J.E., Besselink, N., van Hoeck, A., Janssen, R., Cuppen, E., Priestley, P., and Papenfuss, A.T. (2021). GRIDSS2: comprehensive characterisation of somatic structural variation using single breakend variants and structural variant phasing. Genome Biol 22, 1–25. 10.1186/s13059-021-02423-x
- [SHAPEIT4](https://www.nature.com/articles/s41467-019-13225-y): Delaneau, O., Zagury, J.-F., Robinson, M.R., Marchini, J.L., and Dermitzakis, E.T. (2019). Accurate, scalable and integrative haplotype estimation. Nat Commun 10, 5436. 10.1038/s41467-019-13225-y
- [COBALT](https://github.com/hartwigmedical/hmftools/blob/master/cobalt/README.md): Priestley, P., Baber, J., Lolkema, M.P., Steeghs, N., de Bruijn, E., Shale, C., Duyvesteyn, K., Haidari, S., van Hoeck, A., Onstenk, W., et al. (2019). Pan-cancer whole-genome analyses of metastatic solid tumours. Nature 575, 210–216. 10.1038/s41586-019-1689-y)
- [CIRCOS](https://circos.co/): Krzywinski, M. et al. Circos: an Information Aesthetic for Comparative Genomics. Genome Res (2009) 19:1639-1645


## Known issues and future development
- Currently only tested for the hg38 reference genome. Can in principle be run for other reference genomes as well, as long as the required input files are available (eg ShapeIt maps etc.)
- Currently only tested for slurm
- More documentation and code how to analyze/interpret PTATO output files will be added