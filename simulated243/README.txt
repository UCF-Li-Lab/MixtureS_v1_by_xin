Simulated polymorphic sites for each reference genome
####################################################################################################################
1.File name format:
Example1(group1-group5): mutations_NC_009515.1_0.01_1
(1) NC_009515.1: it is the reference genome
(2) 0.01: it is the mutation rate(0.01%)
(3) 1: It is strain ID. In this study, we at most had 4 strains. So the ID is from 1 to 4. 

Example2(group6_group8): mutations_NC_009515.1_0.01_1_0.1_t1
(1) NC_009515.1: it is the reference genome
(2) 0.01: it is the mutation rate(0.01%)
(3) 1: It is strain ID. In this study, we at most had 4 strains. So the ID is from 1 to 4. 
(4) 0.1: evolution rate which is the percentage of SNPs similar to each other
(5) t1: the evolution type. see our supplementary for detail

2.File content format:
Within each file, each line represents a polymorphic site. It is in the format of 'polymorpihc_site,location'. For example, one line is written as 'C,112419', so its mutation is C at position 112419.
####################################################################################################################
The polymorphic sites corresponding each genome are stored in different group which is same as supplementary Table S1.

Group6, Group7 and Group8:
We defined three types of similarity among strains: (1) Type1: For datasets with three strains, two strains are set to share 10% to 60% polymorphic sites, and the third one is independent to these two strains. The sixth group contained 54 datasets with Type1 evolutionary relationship. (2) Type2: strain1 and strain2 shared 10% to 60% of polymorphic sites, strain3 shared 10% sites with strain1 and strain2, while strain4 did not share polymorphic sites with the other strains. The seventh group contained 54 datasets with Type2 evolutionary relationship (3) Type3:  strain1 and strain2 shared 10% to 60% of polymorphic sites, strain3 and strain4 shared 10% to 60% of polymorphic sites, while strain1 and strain2 shared no polymorphic sites with strain3 and strain4. 

The first, second, third and forth strain above are strain with ID 1, 2, 3 and 4 respectively. 
