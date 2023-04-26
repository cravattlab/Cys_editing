## Abstract
Covalent chemistry represents an attractive strategy for expanding the ligandability of the proteome, and chemical proteomics has revealed numerous electrophile-reactive cysteines on diverse human proteins. Determining which of these covalent binding events impact protein function, however, remains challenging. Here, we describe a base-editing strategy to infer the functionality of cysteines by quantifying the impact of their missense mutation on cancer cell proliferation. The resulting atlas, which covers >13,800 cysteines on >1,750 cancer dependency proteins, confirms the essentiality of cysteines targeted by covalent drugs and, when integrated with chemical proteomic data, identifies essential, ligandable cysteines in >110 cancer dependency proteins. We further show that a stereoselective and site-specific ligand targeting an essential cysteine in TOE1 inhibits the nuclease activity of this protein through an apparent allosteric mechanism. Our findings thus describe a versatile method and valuable resource to prioritize the pursuit of small-molecule probes with high function-perturbing potential.


## Contents
This repository contains code used to 
1) create sgRNA libraries targeting different amino acids including cysteines by base editing.
2) process pool library screen dropout data and make plots

## Base editing sgRNA design
This code is used to design sgRNAs based on the PAM sequences and the base editor type.
It also predicts the associated base editing outcomes.

## Base editing screen data analysis
The input data has been processed by PoolQ (Broad Institute, https://portals.broadinstitute.org/gpp/public/software/poolq) 


## References
1.	Hwang, G.-H. et al. Web-based design and analysis tools for CRISPR base editing. BMC Bioinformatics 19, 542 (2018).
2.	Clement, K. et al. CRISPResso2 provides accurate and rapid genome editing sequence analysis. Nat Biotechnol 37, 224–226 (2019).
3.	Gaudelli, N. M. et al. Programmable base editing of A•T to G•C in genomic DNA without DNA cleavage. Nature 551, 464–471 (2017).
4.	Komor, A. C., Kim, Y. B., Packer, M. S., Zuris, J. A. & Liu, D. R. Programmable editing of a target base in genomic DNA without double-stranded DNA cleavage. Nature 533, 420–424 (2016).
5.	Richter, M. F. et al. Phage-assisted evolution of an adenine base editor with improved Cas domain compatibility and activity. Nat Biotechnol 38, 883–891 (2020).
6.	Thuronyi, B. W. et al. Continuous evolution of base editors with expanded target compatibility and improved activity. Nat Biotechnol 37, 1070–1079 (2019).
