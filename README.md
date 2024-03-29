# SMDIC
a software package for automated identification of Somatic Mutation-Driven Immune Cells

>A computing tool is developed to automated identify somatic mutation-driven immune cells. The operation modes including: i) inferring the relative abundance matrix of tumor-infiltrating immune cells and integrating it with a particular gene mutation status, ii) detecting differential immune cells with respect to the gene mutation status and converting the abundance matrix of significant differential immune cell into two binary matrices (one for up-regulated and one for down-regulated), iii) identifying somatic mutation-driven immune cells by comparing the gene mutation status with each immune cell in the binary matrices across all samples, and iv) visualization of immune cell abundance of samples in different mutation status.

### how to install

```R
Installation method：

1. library(devtools); 
   install_github("hanjunwei-lab/SMDIC")
2. install.packages("SMDIC")

Use：
library(SMDIC)
```

Please cite the following article when using `SMDIC`:

> [1] Jiang Y , Zheng B , Yang Y ,et al.Identification of Somatic Mutation-Driven Immune Cells by Integrating Genomic and Transcriptome Data[J].Frontiers in Cell and Developmental Biology, 2021, 9:715275.DOI:10.3389/fcell.2021.715275.

