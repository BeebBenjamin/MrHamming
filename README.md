MrHamming
==========

<p align="justify">MrHamming is a Python version 3.6 terminal program which can be used to calculate pair-wise Hamming distances (character per character) between genetic or protein sequences contained within a fasta file and is licensed under the MIT license.</p>  
  
<p align="justify">MrHamming, formally known as MrPipeline, was originally written as a part of a more complex pipeline that made use of other bioinformatics software such as RaXML and FastTree, but now only calculates the Hamming distance.  The program as available in this GitHub repository, is a stripped down version of the previously mentioned pipeline that was used in order to generate the data available in our pre-publication article available <a href='https://www.biorxiv.org/content/10.1101/289892v1'>here</a>.  The program uses no model of evolution, so 'divergence' is assessed simply by comparing strings of the same length on a 'character by character' basis e.g. ABC -> ACC would have a reported divergence of 1/3 because one out of three characters is different.  In addition, sometimes, genetic sequences have gaps denoted by a '-' or a 'N', and in order to accurately assess distances, these and the characters that are being compared to should be ignored.  With this in mind MrHamming has two modes, one assuming no gaps which does a straight-forward pair-wise alignment and distance calulation, and another which builds a consensus sequence of all the gaps across sequences in an alignment and uses this to exclude or 'mask' parts of the sequences where gaps are present i.e. it implements complete deletions comparable to other existing software e.g. MEGA X.  For example, A-C -> A--, would have the following mask: A--, and a divergence score of 0.  The example shown here is simplified, but also highlights one of the implications of analysing sequences with too many gaps using this method, i.e. if too much information is missing it is possible to get artificially low divergence values which are therefore not likely to be phylogenetically informative.</p>

Dependencies
============

This project requires pipenv and Python 3.6.
Install and configure a Python virtual environment:

```bash
$> pipenv --python3.6`
$> pipenv install
```
Sample Data
===========
<p align="justify">Two sample fasta files are provided. One containing randomly generated base pair sequences of the same length with no gaps, and one which is the same but with random gaps.</p>

Instructions
============
<p align="justify">The script mrpipeline.py requires an input file and an output file passing as arguments. The input file must be in fasta format, and contain at least two aligned sequences of equal length.  The output file can be called anything.</p>

<p align="justify">Basic use: for calculating the pair-wise distances between sequences (with no gaps), open a terminal in Linux/GNU or Mac OS X, and type the following: </p>
  
 ```bash
 S> sudo python mrpipeline.py input_file output_file"
 ```
<p align="justify">N.B Replace the input file with the file name of the fasta file you wish to generate pair-wise distances for.  Replace output_file with any file name you wish that is allowed by your operating system. Do not use this mode if there are any gaps in the sequences being analysed because this will generate artificial divergence values.</p>

<p align="justify">If the sequences being analysed have gaps in them i.e some of them are partial coverage, MrHamming can generate a consensus mask in order to 'mask out' these regions.  In order to make use of this feature run MrHamming in the terminal with the addition of the 'p' flag as follows:</p> 

```bash
$> sudo python mrpipeline.py input_file output_file p
```  
<p align="justify">N.B Comparing partial sequences can lead to divergence values that are not a true reflection of the evolutionary divergence of two sequences, approach with caution.</p>

<p align="justify">When finished processing the output file from both modes (gaps and no-gaps) will contain a series of lines containing something like the following:</p>

fasta1 fasta2 0.754385964912  
fasta1 fasta3 0.769298245614  
fasta1 fasta4 0.763157894737  
fasta2 fasta3 0.75701754386  
fasta2 fasta4 0.747368421053  
fasta3 fasta4 0.763157894737  

<p align="justify">The above example output shows all the pair-wise comparisons between the sample sequences followed by their divergence value.  These are the raw values and need to be averaged between groups if you wish to compare a group vs another.</p>

To Do
=====
1. Addition of an indicator of % gaps for partial mode will allow a rough indicator of quality of input data and potentially allow this for use in validation experiments i.e. how informative are partial sequences vs complete sequences phylogenetically?

2. Addition of an API interface with GenBank through Biopython to make use of accession numbers in order to retreive species information automatically and group organisms according in intra group comparisons.

3. Add a statistical check to test if distances are normally distributed and aggregate values as appropriate.
