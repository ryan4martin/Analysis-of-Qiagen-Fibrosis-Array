# Analysis-of-Qiagen-Fibrosis-array

This script can be implemented to analyze data collected using the Qiagen RT2 ProfilerTM PCR array platform. The script will provide the fold change for each gene as well as perform select t-tests between treatments of choice for each gene.

The dataset imported needs to have the first column with the names of each well on the PCR plate and the first row containing 'Gene'. All other columns will contain the Ct values with the respective name in the first row.

