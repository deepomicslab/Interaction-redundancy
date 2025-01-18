# Interaction redundancy
**Interaction redundancy** was proposed to describe the conserved functional roles of microbial interactions within a community. We provided AURORA (quAntify the redUndancy in micRObial inteRactions using pArtial information decomposition) package here to quantify the interaction redundancy within a microbial community from microbial abundance data.


![image](https://github.com/deepomicslab/Interaction-redundancy/blob/main/pipeline.png)

# AURORA package
## Prerequisite
The package was written in Python3. Following Python packages should be installed:
+ numpy 2.1.3
+ pandas 2.2.3

## Install
```shell
git clone https://github.com/deepomicslab/Interaction-redundancy.git
```

## Usage
```shell
python interaction_redundancy_cal.py --ICN_ref ref_ICN.csv --abundance_file example_data/abundance_example.tsv
```
Parameters are shown bellow:

+ --ICN\_ref:	The reference interaction content network file (OTU interaction x KEGG interaction).	
+ --abundance\_file:	The taxon abundance data (OTU x sample).
+ --output\_file:	The output file (default: IR_output.txt).
+ --distance\_measure:	The distance measures used for calculating interaction redundancy, including
                        weighted\_jaccard\_distance, euclidean\_distance, correlation\_distance, and
                        manhattan\_distance (default: weighted\_jaccard\_distance).

For more information, please use python interaction_redundancy_cal.py -h. 

## Output
Each row in the output file represents a sample, with the first column displaying the interaction redundancy value and the second column displaying the interaction diversity value.

## Maintainer
Ruo Han Wang ruohawang2-c@my.cityu.edu.hk
