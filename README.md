
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]



<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://research.hsr.it/en/divisions/neuroscience/stem-cells-and-neurogenesis.html">
  <h3 align="center"> Sessa Lab </h3>
  </a>
  <p align="center">
    Stem Cells and Neurogenesis DIBIT2 C1, Floor 4, Room 50A
  </p>
</p>

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
  <li>
      <a href="#about-the-lab">About the lab</a>
  </li>
  <li>
      <a href="#people">People</a>
  </li>
    <li>
      <a href="#pipelines">Pipelines</a>
    </li>
    <li>
      <a href="#scripts">Scripts</a>
    </li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>



<!-- ABOUT THE LAB -->
## About The Lab

Our lab has a strong interest in developing novel technologies in stem cells, genetic cell reprogramming and CRISPR/Cas9 gene editing for better modeling and treating neurological disorders. Patient’s derived iPS cells (iPSCs) offer a superior cellular model to recapitulate the key pathophysiological defects underlying the disease. In addition, CRISPR/Cas9 gene-editing provides a fast and efficient system to prove the direct association between a gene mutation and a specific cellular trait. The group has established numerous strategies for direct cell reprogramming to generate induced neuronal and glial cells for accelerating cellular modeling of human disorders. Moreover, we have established iPSCs from patients suffering from various diseases including Dravet, Rett and ASD syndromes, NBIA and Parkinson’s disease. CRISPR/Cas9 gene editing is a crucial tool in the lab to generate isogenic control iPSCs or to introduce targeted gene mutations. Lately, we have conceived and validated new approaches for correcting mutated genes or modulating their expression by CRISPR technology in vitro and in vivo. To vehiculate these tools in the brain and set up strategies of in vivo gene-therapy, this lab is producing new variants of adeno-associated viruses (AAV) that combined high targeting efficiency, tissue spreading and safety.

<!-- people -->
## People
<ul>
  <li> Alessandro Sessa </p>
<p align="center">
<img src="https://lh4.googleusercontent.com/CbldvbQ7_euzbNgWBFJKVrfk_RSMI9stVjr7w6UTs7rjt5tObjcbyeuf2k51n6plsbfVOzqw_4JCgCDzIiarJivWcFc14o4eNiLBFzUdKXvRwjHyW8YWQwINaffJk48W6Q=w1280" width="300" height="300" /> </p> 
</p> AS is within the group since 2005. Under Vania’s supervision, he obtained the PhD in Molecular Medicine in 2008. During the years, he has work on the genetic and epigenetic factors that settle the development of the brain spanning from the importance of specific-cell populations, to the role of molecular features and pathways. More recently, as staff scientist and team leader, he focused his research in (i) modeling rare diseases aiming to unravel their basis, and (ii) developing tools to intercept pathogenic transformations.
</p>
</li>
  <li>Mattia Zaghi <p align="center">
<img src="photo/zaghi.jpg" width="300" height="350" /> </p>  Mattia Zaghi has obtained is BSC and MSC degree in biotechnology at Vita-Salute San Raffaele University. During his master, he has joined Dr. Alessandro Sessa's group in VB Unit focusing his scientific interest in characterizing the epigenetic bases of neuro-developmental disorders such as intellectual disabilities and autism. He then enrolled in the molecular medicine PhD at Vita-Salute San Raffaele university and completed it, with a thesis about the role of SETBP1-SET protein axes in controlling chromatin structure and gene expression during neural development, in the context of rare sindrome Schinzel-Giedion.

As post-doctoral fellow in the lab is continuing his investigation regarding neurodevelopmental disorders and at the same time implementing protocols for several genomics experiments carried out in the group.</p>
</p>
</li>
  <li>Federica Banfi  <p align="center">
<img src="photo/Fede.png" width="300" height="350" /> </p> Federica Banfi is post-doc in Unit of “Molecular NeuroEngineering” at the San Raffaele Scientific Institute in Milan. She received her Bachelor’s degree in Biological Sciences (2013) and Master’s degree in Molecular Biology of the Cell (2015) at the University of Milan. She obtained her Ph.D in Molecular Medicine (Neuroscience curriculum) in 2020 at Vita-Salute San Raffaele University under the supervision of Dr. Alessandro Sessa, focusing on understanding the roles of SETBP1 in brain development and Schinzel-Giedion syndrome. She gained expertise in in vitro neurodevelopmental disease modelling using induced pluripotent stem cells (iPSc) and their differentiation to neuronal lineage and 3D cerebral organoids. </p>
</p> 
  </li>
  <li>Edoardo Bellini </li>
  </ul>



<!-- Pipelines -->
## Pipelines
In pipeline folder are present all our lab piplines in snakemake. Within every folder are present:
* `.sk` file, that contain all rules sequence.
 ```bash 
 head pipeline/Chip_seq/Chipseq.sk

 # lib
from snakemake.io import glob_wildcards, expand
import glob,os
import pathlib
import pandas as pd
#import multiqc

#config
configfile: "config_ATAC_2.yaml"
```

* `.yaml` file, that must be edit according to organism and statistics and user needs.
 ```bash 
 cat pipeline/Chip_seq/config.yaml

# proj parameters
Project: "prova_chip_seq"
skipH: 6
RAWDATA: "/beegfs/scratch/ric.broccoli/ric.broccoli/prova_chip_seq"
RUN_ID: "prova_chip_seq" # "RUN_id"
# genome
genome: "hg38"
ref_genome_fa: "/beegfs/scratch/ric.broccoli/ric.broccoli/Genomes/hg38/fa/hg38.fa"
chrom_sizes: "/beegfs/scratch/ric.broccoli/ric.broccoli/Genomes/hg38/hg38.chrom.sizes"
# trimming
adapters: "/beegfs/scratch/ric.broccoli/ric.broccoli/adapters/NexteraPE-PE.fa"
# blacklist
blacklist_url: "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz"
# peaks
genome_size_bp: 3209286105
peaks_qvalue: 0.01
broad_cut_off: 0.001
```
* `.csv` file conteins all information abount the metadata of the samples such as sample condition, or whatever you need to specify.
 ```bash
cat pipeline/Chip_seq/samplesheet.csv

[Header],,,,,,,,,,,
IEMFileVersion,4,,,,,,,,,,
Date,2019-06-07 08:32:36+00:00,,,,,,,,,,
,,,,,,,,,,,
[Reads],,,,,,,,,,,
,,,,,,,,,,,
,,,,,,,,,,,
[Settings],,,,,,,,,,,
,,,,,,,,,,,
[Data],,,,,,,,,,,
Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description,Container_Label
,CtrlH3K27acChIP-Seq2_S189_R1_001.fastq.gz,CtrlH3K27acChIP-Seq2_S189_R1_001.fastq.gz,,,UDI0073,CAATTAAC,UDI0073,CGAGATAT,prova_chip_seq,,
,CtrlH3K27acChIP-Seq2_S189_R2_001.fastq.gz,CtrlH3K27acChIP-Seq2_S189_R2_001.fastq.gz,,,UDI0073,CAATTAAC,UDI0073,CGAGATAT,prova_chip_seq,,
,InputCtrlChIP-Seqs_S182_R1_001.fastq.gz,InputCtrlChIP-Seqs_S182_R1_001.fastq.gz,,,UDI0073,CAATTAAC,UDI0073,CGAGATAT,prova_chip_seq,,
,InputCtrlChIP-Seqs_S182_R2_001.fastq.gz,InputCtrlChIP-Seqs_S182_R2_001.fastq.gz,,,UDI0073,CAATTAAC,UDI0073,CGAGATAT,prova_chip_seq,,
 ```

<!-- SCRIPTS -->
## Scripts
In scripts folder are present stand alone script used generally to perform downstream analisys like scRNA-seq clustering, Bulk RNA-seq, or custom plot.

<!-- CONTACT -->
## Contacts


Edoardo Bellini - <bellini.edoardo@hsr.it>

Mattia Zaghi - <zaghi.mattia@hsr.it>

Alessandro Sessa - <sessa.alessandro@hsr.it>

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/othneildrew/Best-README-Template.svg?style=for-the-badge
[contributors-url]: https://github.com/othneildrew/Best-README-Template/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/othneildrew/Best-README-Template.svg?style=for-the-badge
[forks-url]: https://github.com/othneildrew/Best-README-Template/network/members
[stars-shield]: https://img.shields.io/github/stars/othneildrew/Best-README-Template.svg?style=for-the-badge
[stars-url]: https://github.com/othneildrew/Best-README-Template/stargazers
[issues-shield]: https://img.shields.io/github/issues/othneildrew/Best-README-Template.svg?style=for-the-badge
[issues-url]: https://github.com/othneildrew/Best-README-Template/issues
[license-shield]: https://img.shields.io/github/license/othneildrew/Best-README-Template.svg?style=for-the-badge
[license-url]: https://github.com/othneildrew/Best-README-Template/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/othneildrew
[product-screenshot]: images/screenshot.png
