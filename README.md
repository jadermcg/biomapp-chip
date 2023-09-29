# BIOMAPP::CHIP: Large-Scale Motif Analysis

## Overview

`Biomapp::chip` is a C++ computational tool designed for the efficient discovery of biological motifs in large datasets, specifically optimized for ChIP-seq data. Utilizing advanced k-mer counting algorithms and data structures, it offers a streamlined, accurate, and fast approach to motif discovery.

## How it Works

The core functionality of `Biomapp::chip` revolves around its innovative k-mer counting method implemented via a specialized suffix tree data structure known as `SMT` (Sparse Motif Tree). The `SMT` ensures both speed and accuracy in the counting process.

### Algorithmic Workflow

1. **Data Input**: The algorithm starts by taking sequence data as input, typically in FASTA or FASTQ format.
  
2. **K-mer Counting**: Employing the `SMT` data structure, k-mer frequencies in the sequence data are accurately and efficiently counted.
   
3.  **Pre-optimized models**: `Biomapp::chip` toolkit includes an advanced feature of pre-optimized models, which are essential for accelerating the motif discovery process. These models are generated after the `SMT` component has successfully counted the kmers.
  
4. **Motif Discovery**: After obtaining k-mer frequencies, the algorithm proceeds to identify statistically significant motifs using a efficient version of `EM` called `FAST-EM` set through predefined metrics and statistical tests.

5. **Output**: The identified motifs, along with their statistical significance and locations within the sequence data, are outputted for further analysis or visualization.

#### How Pre-optimized Models Are Generated
1. **K-mer Counting via smt**: Initially, the `smt` component performs the task of k-mer counting in the given biological sequences.

2. **Sibling K-mer Discovery via kdive**: Following the counting process, the `kdive` algorithm comes into play. It identifies all sibling kmers with up to `d` mutations from the original k-mer.

3. **Model Optimization**: Based on these sibling kmers, pre-optimized models are generated, which serve as heuristic models for guiding the main algorithm in subsequent runs.

#### Benefits of Using Pre-optimized Models

- **Time Efficiency**: Utilizing pre-optimized models significantly speeds up the motif discovery process.
  
- **High Accuracy**: The sibling kmers identified by `kdive` ensures a nuanced and robust model, making the discovery process more accurate.
  
- **Resource Efficiency**: The use of pre-optimized models reduces the computational resources required, making it ideal for lower-end systems as well.

## Why Use biomapp::chip?

- **Efficiency**: Designed to handle large volumes of data, the tool excels in scenarios where traditional motif discovery algorithms falter.
  
- **Accuracy**: Advanced counting algorithms and statistical tests ensure a high rate of true positive motif discoveries.
  
- **Versatility**: Suitable for a wide range of applications in genomics, especially in the study of transcription factors and their binding sites.

## Installation

Instructions to follow...

## Usage

Uso: biomapp -i <fasta> <options>
Options:

-k <size of kmer>
-n <number of models>
-d <number of mutations>
-e <type of EM. Can be oops, zoops or anr>
-r <number of em iterations>
-f <cutoff for convervenge control>
-c <compression: 0 no compression, 1 LF4 compression>

## Dependencies

To run the program, you will need the following libraries installed on your Ubuntu-based system:

- Boost Filesystem
- R
- BLAS (Basic Linear Algebra Subprograms)
- TBB (Threading Building Blocks)
- Standard C++ Library
- GNU C Library
- GNU Multiple Precision Arithmetic Library
- GCC support library
- Other miscellaneous libraries (readline, pcre2, lzma, bz2, z, tirpc, ICU, tinfo, gssapi, krb5)

You can install these dependencies using `apt` with the following commands:

```bash
sudo apt update
sudo apt install libboost-filesystem-dev r-base libblas-dev libtbb-dev libstdc++6 libc6 libgomp1 libgcc1 libreadline8 libpcre2-dev liblzma5 libbz2-1.0 zlib1g libtirpc-dev libicu-dev libtinfo6 libgssapi-krb5-2 libkrb5-3 libk5crypto3 libcom-err2 libkrb5support0 libkeyutils1 libresolv2
```
