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

If you publish results from this tool or produce a derived version please cite the original paper:

Gaunt TR, Rodriguez S, Guthrie PA, Day IN. An expectation-maximization program for determining allelic spectrum from CNV data (CoNVEM): insights into population allelic architecture and its mutational history. [Hum Mutat. 2010 31:414-420](https://onlinelibrary.wiley.com/doi/abs/10.1002/humu.21199) (PMID: 20077501)

## Further notes on re-use/adaptation

This is a basic modification of the web-based ConVEM and is simply provided to enable others to adapt and use the code directly. This script takes interactive input, but could easily be adapted to batch input.

The script uses a multi-threading approach that could be converted to Python multi-processing (v 2.7 as the code is not compatible with Python 3.0). We avoided multi-processing to protect our web-server from excessive load, but on a multi-core system this could give a big speed benefit. Essentially the parallel runs would be handled in the same way.

### Sample session

```
$ python2.7 convem_cmdline.py
Enter CNV bin counts as comma-delimited list, starting with bin zero:
0,0,5,17,28,32,25,14,4,2,0
Possible solutions:


Result 1 :  100 starts (chi^2 = 0.58, mean 125 iterations)


=======================


Result 1
Copy number	0	1	2	3	4	5	6	7	8	9	10
Observed count	0	0	5	17	28	32	25	14	4	2	0
Expected count	0	0	5	16	29	32	25	13	5	1	0
Allele freq	0.0	0.2019	0.317	0.3106	0.1393	0.0278	0.0035	0.0	0.0	0.0	0.0


===========================================
```
