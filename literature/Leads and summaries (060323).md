# Sources and summary for Introduction of report

##1. Microarray technology

###Article: [_Microarray Technology_ ](https://www.genome.gov/genetics-glossary/Microarray-Technology)  

**Summary**:  
Small nucleic acid fragments (DNA or RNA) are bound to a solid surface, this is called a "Chip".
This Chip is then bathed in with fragmented DNA or RNA isolated from a study sample.
This is also fragmented and is able to recognize and bind to the chip-fragments via complementary base pairing, which causes a fluorescent light in that particular cell.
The light is detected and outputted into a visual coloring of the chip.  
The technology allows for genomic analysis without sequencing, this cuts down on costs.
Able to visualize gene expression and gene product (RNA) present in tissues/cells.
Also makes SNPs able to be visualized.

###Article: [_Mammary Gland: Gene Networks Controlling Development and Involution_](https://doi.org/10.1016/B978-0-08-100596-5.00883-0)  

_[Has a small section about Microarray Technology]_  

**Summary**:  
It allows for quantitative and simulataneous monitoring of expression of thousands of genes in tissues or cells.
Important is that these are single-stranded (red labeled) DNA molecules that are able to hybridize to the chips.
The DNA fragments that are bound to the chips can be changed to the DNA of specific genes depending on the gene that is studied.
With this technique, it is possible to measure the abundance of each RNA molecule in a sample.
This is done by reverse transcribing the RNA to cDNA molecules.
During the synthesis and processing of the cDNA, green fluorescent dye is incorporated into the molecules.
The abundance of each RNA molecule can then be measured using the variations in fluorescent light that is detected and put into a numerical value based on the light intensity.
Test samples are always measured against a reference microarray or control sample.
This reference data is then used to normalise the test data. 
The final data for each sample of interest are the average of the signal intensity of the sample divided by the signal intensity of the reference sample across both dye combinations.
This ratio is then used to compare RNA abundance among samples.
After which it goes through statistical analysis to generate a list of differentially expressed genes, meaning genes that are effected by certain test conditions.

###Article: [_Microarray and its applications_](https://journals.lww.com/jpbs/Fulltext/2012/04002/Microarray_and_its_applications.48.aspx)
**Summary**:
The principle of microarrays is that complementary sequences will bind to each other.
The mRNA is converted into a stable cDNA form which is labeld by fluorochrome dyes Cy3 (green) and Cy5 (red).
Unknown DNA molecules are cut into fragments using restriction endonucleases, the fragments are later marked with fluororescent markers.
These are allowed to interact with the probes on the chips. 
The DNA molecules that are able to bind to the chip and can be detected using a laser beam, which activates the fluorochromes

###Article: [_Microarray Data Analysis and Mining Tools_](https://d1wqtxts1xzle7.cloudfront.net/39088129/Paper_1-Saravanakumar_Selvaraj-libre.pdf?1444487431=&response-content-disposition=inline%3B+filename%3DMicroarray_Data_Analysis_and_Mining_Tool.pdf&Expires=1677932600&Signature=XS6toQqDlcZTAyhtsr84na9iOgSDrj6AuYrIjAqpwbxWS~2L7YPI74HJJ-p7RS7F3gmilk266MHh9Qhp0uYPz3ZQKDmYLRPd8PKp38DIuOVoDJB4HurnbuF2LGGYLWfXQfi6WdnWPHrLJ3iYCo-ih5SYvurhfEya6KdFZtffmoCmaGI-VX-Kj3GqacG7~TtrAd0B69q96G0et~OpZ2hV1POppO7QIpykUx42UXCGp1XOcPvFig04hd7hvcxRTr~s8h0oW0asphxqIzahY1MRjuzppzyAMtA9YRWaNVAHIxsysLDW5Su3gs7owNkU-7pk0MWOE~mASxxIfGqSumRZcQ__&Key-Pair-Id=APKAJLOHF5GGSLRBV4ZA)
**Summary**:  
Data mining including classification and clustering used to extract useful knwoledge from microarray data.
This paper highlights such tools and technologies.
Data sets are very large, it's important to reduce the set to those genes that are best distinguished between two cases or classes.
This differential gene expression is the first task of an in depth microarray analysis.
Clustering is an approach to classify data in to groupss of samples that show similar patterns, which are characteristic in that group.
In the past, the fold change in expression was used as an assumption for differentially expressed genes, making them biologically significant when a high fold change was present.
Significance was tested using t-test, f-statistics, models, ANOVA,...  
Clustering is used to find co-regulated and functionally related groups.
This is particularly interesting when we have expression levels in the whole organism.
There are three common types of clustering:
- Hierarchical clustering
- K-means clustering
- self-organising maps

<ins>Hierarchical clustering:</ins>
Iteratively group genes together that highly correlate in their expression. This is then repeated on the groups itself.
It is represented as a dendrogram (a large branching tree), which will later on be clustered visually.  
<ins>K-means clustering:</ins>
Data mining/machine learning algorithm to cluster observations into groups with related observations without prior knowledge of relationships.
A vector space is is partitioned into K parts, and the algorithm calculates the center point in each subspace and adjusts the partition so that each vector is assigned to the cluster the center of which is the closest.  
<ins>Self-organising map:</ins>
Neural network based clustering approach, similar to K-means.

Another part of microarray data analysis is classification of datasets.
A classifier will find a rule to allow for the assignment of different samples to a particular class.
This should be done by allowing an algorithm to be trained on a sufficient amount of sample cases.
Examples of algorithms are k Nearest Neihbors (kNN), Artificial Neural Networks, weighted voting ans dupport vector machines (SVM).


---

##2. K-Mean Clustering
It's a clustering algorithm. 
It will put k points in a vector space and assign all the datapoints of the dataset to the nearest point (centroid) meaning it will be assigned to this particular cluster.
When all the data points are assigned, the average of the vectors of each point in a cluster is calculated and this will be the position of the new centroid.
This process is repeated until convergence.
Maybe good link? Not sure:
https://www.javatpoint.com/k-means-clustering-algorithm-in-machine-learning

###Article: [_Research on k-means Clustering Algorithm_](https://sci-hub.ru/10.1109/IITSI.2010.74)
Clustering is used to classify raw data reasonably and uses data mining to look for hidden patterns in datasets.
It groups data objects into disjointed clusters so that the data in the same cluster are similar, whilst the data in seperate clusters will differ.
K-means clustering is a numerical, unsupervised, non-deterministic iterative method.
It's a very effective way that can produce good clustering results 


---
###3. Paper: Adaptive quality-based clustering of gene expression profiles
**Summary**:  
They analyse drawbacks in the classical algorithms, mainly that they require the predefinition of arbitrary parameters (like amount of clusters) and that they force every gene into a cluster despite a low correlation with other cluster members.
It describes a novel adaptive quality-based clustering algorithm to tackle some of the drawbacks
This is done by finding a sphere in the high-dimensional representation of the data, where the "density" of gene expression profiles reaches the highest local maximum.
Based on a preliminary estimate of the radius of a cluster (quality based approach).
Secondly they derive an optimal radius of the cluster (adaptive approach) so that only significantly co-expressed genes are included in the same cluster.
A model is fitted to estimate this radius using an EM-algorithm or an Expectation-maximization algorithm. 
This two step method was succesfully verified using existing data sets


