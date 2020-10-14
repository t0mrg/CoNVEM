# CoNVEM

This is a legacy command-line version [CoNVEM - CNV allele frequency estimation by expectation maximisation](http://apps.biocompute.org.uk/convem/). 

## Requirements

The code requires Python 2.7 and will not run under Python 3.x without some edits to the script.

## Usage

The command line version will take interactive input and prints the results directly to screen after completion. To run the script:

```
python2.7 convem_cmdline.py
```

At the prompt enter some CNV counts. The example data (for you to test the script is working) are:

```
0,0,5,17,28,32,25,14,4,2,0
```

The input requirements are: CNV bin counts as comma-delimited list, starting with bin zero, and using a zero to represent any bin for which there are no subjects

## Citation

If you publish results from this tool please cite the original paper:

Gaunt TR, Rodriguez S, Guthrie PA, Day IN. An expectation-maximization program for determining allelic spectrum from CNV data (CoNVEM): insights into population allelic architecture and its mutational history. [Hum Mutat. 2010 31:414-420](https://onlinelibrary.wiley.com/doi/abs/10.1002/humu.21199) (PMID: 20077501)
