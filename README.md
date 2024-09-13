# R scripts for quality control of NEBNext Immune Sequencing library preparation

## 1. PCR2 quantification
The amount of BCR/TCR template RNA can vary between samples, even when identical input RNA amount is used. For this reason, the optimal PCR2 cyclecount (N° of cycles at which a sufficient amount of product is generated without exhausting the reagents) must be determined for each sample individually. The optimal PCR2 cyclecount is defined as the cycle at which 2/3rds of the maximal fluorescence is reached.

- **Input data**: Gene expression result and Cq result files (.csv format) from CFX96 Real-Time System.
- **Parameters**: Wells to filter from analysis; Samples to highlight; Cycle cutoff values; ...
- **Outputs**: PDF file with relevant metadata, amplification plots, and optimal cyclecount for each sample.

## 2. TapeStation size distribution
Once the libraries have been prepped, their size distribution must be assessed to ensure that they are free from non-specific amplicons. This step also allows to estimate the library concentrations. 

- **Input data**: TapeStation .XML results file, Electropherogram .csv file, and annotations .csv file if non-standard sample loading order.
- **Parameters**: Wells to exclude, dilution factor, Reference sizes of BCR/TCR/IGH peaks.
- **Outputs**: PDF file with relevant metadata, size distribution and concentration estimations.

## 3. qPCR library quantification
The last step of libraryprep is to accurately quantify the libraries. This is acchieved through qPCR with standards.

- **Input data**: Gene expression result and Cq result files (.csv format) from CFX96 Real-Time System.
- **Parameters**: Wells to exclude, dilution factor, Reference sizes of BCR/TCR/IGH peaks.
- **Outputs**: PDF file with relevant metadata, amplification curves, and absolute sample concentrations with SD.
