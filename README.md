<img src="logo.jpeg" width="350" height="350">

# BIOMAPP::CHIP: Large-Scale Motif Analysis

## Overview

`Biomapp::chip` is a Linux/C++ computational tool designed for the efficient discovery of biological motifs in large datasets, specifically optimized for ChIP-seq data. Utilizing advanced k-mer counting algorithms and data structures, it offers a streamlined, accurate, and fast approach to motif discovery.

## Quick start

### Installation

#### Step 1: Download the Binary Files
First, you need to download all the contents of the `bin` directory from this GitHub repository. You can do this easily with svn. First install subversion: ```sudo apt install subversion``` and then run the following command within the directory in which you want to install Biomapp::chip: ```svn export https://github.com/jadermcg/ biomapp-chip/trunk/bin```.

#### Step 2: Place the Binary Files
After downloading, place the binary files in a directory of your choice on your computer. For example, you could place them in a folder called `biomapp_chip/bin` under your home directory. If you created the biomapp_chip directory and downloaded the files there, run this command line to make files executable ```sudo chmod -Rf u+x bin```.

#### Step 3: Update the PATH Environment Variable
Lastly, you need to update your PATH environment variable to include the directory where you placed the binary files. You can do this using the `export` command in Linux. Open your terminal and run the following command:

```bash
export PATH=$PATH:/path/to/your/biomapp_chip/bin
```
You can place this command inside your local user's .bashrc file so that you don't have to type it every time a command terminal is opened.

```
echo 'export PATH=$PATH:/path/to/your/biomapp_chip/bin' >> ~/.bashrc
source ~/.bashrc
```

#### Dependencies

To run the program, you will need the following libraries installed on your Linux Ubuntu/Debian-based system:
- Armadillo Linear Algebra
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
sudo apt install libarmadillo-dev libboost-filesystem-dev r-base libblas-dev libtbb-dev libstdc++6 libc6 libgomp1 libgcc-s1 libreadline8 libpcre2-dev liblzma5 libbz2-1.0 zlib1g libtirpc-dev libicu-dev libtinfo6 libgssapi-krb5-2 libkrb5-3 libk5crypto3 libcom-err2 libkrb5support0 libkeyutils1 liblz4-dev
```

Finally, you need to install some dependencies in R using these commands:
```
Rscript -e 'if (!require("BiocManager")) install.packages("BiocManager", dependencies=TRUE)'
Rscript -e 'if (!require("seqLogo")) BiocManager::install("seqLogo")'
Rscript -e 'if (!require("evd")) BiocManager::install("evd")'
```

#### Usage
```
Use: biomapp -i <fasta> <options>
Options:

-k <size of kmer>
-n <number of models>
-d <number of mutations>
-e <type of EM. Can be oops, zoops or anr>
-r <number of em iterations>
-f <cutoff for convervenge control>
-c <compression: 0 no compression, 1 LF4 compression>
```
#### Example
To understand how the program works, you can run Biomapp::chip on the example dataset that is provided in the project root.

```
biomapp -i MA0003.4.fasta.masked.dust -k 14 -n 5 -d 2-e zoops -r 1000 -f 0.001
```
This execution of Biomapp::chip will process the dataset MA0003.4.fasta.masked.dust using kmers of length 14. The -n 5 argument specifies that the five most optimal PWM models will be generated. The -d 2 parameter indicates that up to two mutations are allowed within each model. Regarding the Expectation-Maximization (EM) algorithm, the chosen type is "zoops," as denoted by the -e zoops parameter. The stopping criteria for the EM involve two components: the maximum number of iterations, set at 1000 (indicated by -r 1000), and a minimum threshold for improvement between successive solutions, set at 0.001 (indicated by -f 0.001). The convergence process will be reached when either the number of iterations exceeds 1000 or the improvement between solutions falls below 0.001.

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

## Why Use Biomapp::chip?

- **Efficiency**: Designed to handle large volumes of data, the tool excels in scenarios where traditional motif discovery algorithms falter.
  
- **Accuracy**: Advanced counting algorithms and statistical tests ensure a high rate of true positive motif discoveries.
  
- **Versatility**: Suitable for a wide range of applications in genomics, especially in the study of transcription factors and their binding sites.
