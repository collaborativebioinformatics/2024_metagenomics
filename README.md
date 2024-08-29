# MIMIC Metagenome Simulator

The MIMIC Metagenome simulator creates simulated Oxford Nanopore (ONT) reads from any existing metagenomic community
by taking in a fastq file and analyzing the naturally occurring microbial abundances. 

This project was conceptualized and built during the
Baylor College of Medicine Human Genome Sequencing Center (HGSC) 2024 Hackathon. 


![alt text](docs/img/logo_small.png)

## Table of Contents 
1. [Introduction](#introduction) 
2. [Installation](#installation) 
3. [Dependencies](#dependencies)
4. [How to use It](#how-to-use-it) 
6. [MIMIC workflow](#mimic-workflow) 
7. [Contributors](#contributors) 
8. [References](#references) 


## Introduction

Oxford Nanopore (ONT) sequencing is a key technology in metagenomic studies, offering long reads. 
To facilitate the benchmarking and testing of bioinformatics tools intended to ONT data analysis, 
we built The MIMIC Metagenome simulator. 
The MIMIC Metagenome simulator creates simulated Oxford Nanopore (ONT) reads from any existing 
metagenomic community by taking in a fastq file and analyzing the naturally occurring microbial abundances. 
This tool can be beneficial to accurately evaluate and benchmark bioinformatics algorithms, and 
enable to generate artificial fastq files and abundance tables. Cost-effective simulation and controlled 
environment facilitate bioinformatic tool testing, ensuring that tools perform as expected. 



## Installation

## Dependencies

### To generate simulated reads

- Python 3.9
- NCBI Datasets v15.27.1
- Minimap2 v2.24-r1122
- Samtools v1.15.1
- Biopython
- Pandas
- Ete3 v3.1.2
- BWA v0.7.17
- FastANI
- Kraken2
- Bracken
- Nanosim

### To generate truth table

### Other 


## How to use It


### MIMIC Workflow
![alt text](docs/img/flowchart_v3.png)

## Contributors

Hackathon team: Todd	Treangen, Shwetha	Kumar, Ryan	Doughty, Sumaiya	Khan, Iva	Kotásková, Arthur	Shem Kasambula, Mike	Nute

## References 



