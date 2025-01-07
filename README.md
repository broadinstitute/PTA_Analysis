Here's a **README.md** draft for your GitHub repository that describes the analysis pipeline for PTA and WGS data:

---

# PTA-WGS Analysis Pipeline

This repository contains a **workflow description and tools** for analyzing **Primary Template-directed Amplification (PTA)** and **Whole Genome Sequencing (WGS)** data. The pipeline includes all steps necessary for high-quality variant detection, including mapping, duplicate marking, realignment, variant calling, and filtration.

## Overview

This pipeline processes WGS reads mapped against the **GRCh38 reference genome**, focusing on:
- Handling PTA-specific artifacts using advanced machine learning models (PTATO Random Forest Model).
- High-quality variant calling with **GATK**.
- Comprehensive variant filtration for **somatic and germline mutation detection**.

### Key Features
- Support for both **FASTQ** and **CRAM** inputs.
- Modular WDL workflow for running on **Terra** or locally.
- Incorporation of existing GATK tools and custom PTA-specific artifact handling.

---

## Workflow Steps

### 1. **Input Files**
   - **FASTQ or CRAM files:** Raw sequencing data.
   - **Reference Genome:** GRCh38.
   - Additional inputs: BED files, BAM indices, and genome-specific annotations.

### 2. **Pipeline Steps**
#### **Step 1: Mapping Reads** (if FASTQ provided)
   - Tool: `BWA mem`
   - Output: SAM/BAM files.

#### **Step 2: Processing CRAM/BAM Files**
   - Duplicate marking using **Sambamba**.
   - Realignment with **GATK BaseRecalibrator**.

#### **Step 3: Variant Calling**
   - Tool: `GATK HaplotypeCaller`
   - Multi-sample mode with `EMIT_ALL_CONFIDENT_SITES`.

#### **Step 4: Variant Filtration**
   - Tool: `GATK VariantFiltration`
   - Filters:
     - `QD < 2.0`
     - `MQ < 40.0`
     - `FS > 60.0`
     - `HaplotypeScore > 13.0`
     - `SOR > 4.0` and more.

#### **Step 5: PTA Artifact Filtering**
   - **PTATO Random Forest Model** trained with 26 genomic features:
     - Allelic imbalance (most important).
     - DNA replication timing.
     - Distance to the nearest gene.
     - Repeat regions and sequence context.

---

## Installation and Dependencies

### Requirements
- Docker
- Cromwell/WDL runtime (for Terra)
- Tools included:
  - **BWA** (v0.7.17)
  - **Sambamba** (v0.6.8)
  - **GATK** (v4.1.3.0)
  - **Samtools** (v1.10)

### Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/<your-username>/PTA-WGS-Analysis.git
   cd PTA-WGS-Analysis
   ```

2. Install dependencies:
   - Use Docker containers for all tools.

3. Prepare input files (FASTQ or CRAM, Reference Genome).

---

## Running the Pipeline

### **On Terra**
1. Import the WDL workflows into Terra.
2. Set up input JSON files for each step.

### **Locally**
Run the pipeline using Cromwell:
```bash
java -jar cromwell.jar run main.wdl -i inputs.json
```

---

## Repository Contents

- `workflows/`: WDL workflows for each pipeline step.
- `tasks/`: Individual WDL tasks for mapping, duplicate marking, realignment, and variant calling.
- `docker/`: Dockerfiles for each tool used in the pipeline.
- `examples/`: Example input files and configuration JSONs.

---

## Resources

- GATK Documentation: [https://gatk.broadinstitute.org/](https://gatk.broadinstitute.org/)
- PTATO Tool: [https://github.com/ToolsVanBox/PTATO](https://github.com/ToolsVanBox/PTATO)
- Reference Genome (GRCh38): [Ensembl](https://www.ensembl.org/index.html)

---

## Contributing

Contributions are welcome! Please submit a pull request or open an issue to suggest improvements.

---

## License

This repository is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

Let me know if you'd like to add more specific details or adjust the sections!
