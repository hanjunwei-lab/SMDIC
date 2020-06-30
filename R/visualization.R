#' @title plotwaterfall
#' @description Function‘plotwaterfall‘ plots the waterfall for mutation genes which drive immune cells.
#' @param maffile The name of mutation annotation file (MAF) format data. It must be an absolute path or the name  relatived to the current working directory.
#' @param mutcell.summary The result of `mutcellsummary` function
#' @param cellnumcuoff a threshold value (4 as the default value). The mutation genes which drive at least "cellnumcuoff" cells are retained for drawing an waterfall.
#' @importFrom maftools oncoplot
#' @importFrom maftools read.maf
#' @export
#' @examples
#' file<-"dir" #dir must be an absolute path or the name  relatived to the current working directory.
#' \dontrun{plotwaterfall(maffile = dir,mutcell.summary = summary,cellnumcuoff =0)}
plotwaterfall <- function(maffile,mutcell.summary,cellnumcuoff=4) {
  ## requireNamespace("maftools")|| stop("package maftools is required,please install package maftools")
  gene.top<-mutcell.summary[which(mutcell.summary[,"count"]>=cellnumcuoff),1]
  maf<-read.maf(maf = maffile)
  oncoplot(maf = maf,
           genes = gene.top,
           fill = T,
           fontSize = 0.8 ,
           showTumorSampleBarcodes = F)
}

#' @title plotCoocMutex
#' @description Function `plotCoocMutex` plots the co-occurrence and mutual exclusivity plots for mutation genes which drive immune cells.
#' @param maffile The name of mutation annotation file (MAF) format data. It must be an absolute path or the name  relatived to the current working directory.
#' @param mutcell.summary The result of `mutcellsummary` function
#' @param cellnumcuoff A threshold value (4 as the default value). The mutation genes which drive at least "cellnumcuoff" cells are retained for drawing a co-occurrence and mutual exclusivity plots.
#' @importFrom maftools somaticInteractions
#' @importFrom maftools read.maf
#' @references Gerstung M, Pellagatti A, Malcovati L, et al. Combining gene mutation with gene expression data improves outcome prediction in myelodysplastic syndromes. Nature Communications. 2015;6:5901. doi:10.1038/ncomms6901.
#' @export
#' @examples
#' cellmatrix<-GetExampleData("cellmatrix") # Cell abundance matrix
#' mutmatrix<-GetExampleData("mutmatrix") # A binary mutations matrix
#' mutcell<-GetExampleData("mutcell") # The result of `mutcorcell` funtion
#' summary<-summary<-mutcellsummary(mutcell = mutcell,mutmatrix = mutmatrix,cellmatrix=cellmatrix)
#' file<-"dir" #dir must be an absolute path or the name  relatived to the current working directory.
#' \dontrun{plotCoocMutex(maffile = dir,mutcell.summary = summary,cellnumcuoff =0)}
plotCoocMutex <- function(maffile,mutcell.summary,cellnumcuoff=4) {
  ## requireNamespace("maftools")|| stop("package maftools is required,please install package maftools")
  gene.top<-mutcell.summary[which(mutcell.summary[,"count"]>=cellnumcuoff),1]
  maf<-read.maf(maf = maffile)
  somaticInteractions(maf = maf, genes = gene.top, pvalue = c(0.05, 0.01))
}


#' @title heatmapcell
#' @description A function to draw clustered heatmaps for the cells driven by a somatic mutation.
#' @param gene Somatic mutant gene name
#' @param cellmatrix Cell abundance matrix, cellmatrix is the result of function `exp2cell`.
#' @param mutcell A list, mutcell is the result of function `mutcorcell`.
#' @param mutmatrix A binary mutations matrix, in which 1 represents any mutation occurs in a particular gene in a particular sample, otherwise the element is 0.
#' @importFrom pheatmap pheatmap
#' @export
#' @examples
#' mutcell<-GetExampleData("mutcell") # The result of `mutcorcell` function.
#' cellmatrix<-GetExampleData("cellmatrix") # Cell abundance matrix
#' mutmatrix<-GetExampleData("mutmatrix") # A binary mutations matrix
#' heatmapcell(gene = "TP53",mutcell = mutcell,cellmatrix = cellmatrix,mutmatrix = mutmatrix)
heatmapcell <- function(gene,mutcell,cellmatrix,mutmatrix) {
  mutcell<-mutcell$mut_cell
  intersample<-intersect(colnames(cellmatrix),colnames(mutmatrix))
  mutmatrix<-mutmatrix[,intersample]
  cellmatrix<-cellmatrix[,intersample]
  cell.gene<-names(which(mutcell[gene,]==1))
  order.genemut<-order(as.vector(mutmatrix[gene,]))
  mutwithorder<-mutmatrix[gene,order.genemut]
  cell.gene_ordmut<-cellmatrix[cell.gene,order.genemut]
  cell.gene_ordmut_zscore<-t(apply(cell.gene_ordmut,1,scale))
  colnames(cell.gene_ordmut_zscore)<-colnames(cell.gene_ordmut)
  cell.gene_ordmut_zscore[which(cell.gene_ordmut_zscore>2)]<-2
  cell.gene_ordmut_zscore[which(cell.gene_ordmut_zscore<(-2))]<-(-2)
  ann_colors = list(
    genemut = c("#F9A452","#4B86C5")
  )
  annotation_col = data.frame(
    mutstat = factor(mutwithorder)
  )
  pheatmap(cell.gene_ordmut_zscore,scale ="none",show_rownames = TRUE,show_colnames = FALSE,cluster_cols = FALSE,annotation_col =annotation_col ,annotation_colors = ann_colors)
}









getPrintFile<-function(coxph.fit,gene,GeneCoefResult){
  result<-summary(coxph.fit);#txt<-paste(fileName,"_cox.txt",sep="")  # print(result)
  header=t(c("geneName","coef","exp(coef)[HR]","se(coef)","Z[coef/se]","Pr(>|Z|)[pval]","lower.95","upper.95","Likelihood ratio test","Wald test","Log rank test"))
  count<-rownames(result[[7]])

  if(length(count)==1){
    aa<-result$coefficients
    bb<-result$conf.int
    temp<-t(c(gene,aa,bb[c(3,4)],result[[9]][3],result[[12]][3],result[[10]][3]))
    if(is.null(GeneCoefResult)==TRUE){
      GeneCoefResult<-rbind(GeneCoefResult,temp)
      colnames(GeneCoefResult)<-header
    }else{
      GeneCoefResult<-rbind(GeneCoefResult,temp)
    }

  }else{
    aa<-result$coefficients
    bb<-result$conf.int
    n<-dim(aa)[1]
    GeneCoefResult<-data.frame(gene[1:n],aa,bb[,c(3,4)],rep(result[[9]][3],times=n),rep(result[[12]][3],times=n),rep(result[[10]][3],times=n))
    colnames(GeneCoefResult)<-header
  }
  return(GeneCoefResult)
}
getUniOrMultiCOXAnalysis<-function(subprof,clin,method){

  subprof<-t(subprof);
  samples<-rownames(subprof);
  mm<-match(samples,as.character(clin[,1]))
  newsubprof<-subprof[which(is.na(mm)==FALSE), ,drop=FALSE];
  n1<-dim(newsubprof)[2]
  clinInfor<-clin[mm[which(is.na(mm)==FALSE)], ,drop=FALSE];
  n2<-dim(clinInfor)[2]
  cox.data<-cbind(newsubprof,clinInfor[,-1])[,c((n1+1):(n1+2),(1:n1))]

  numdata<-apply(unname(cox.data),c(1,2),as.numeric)

  gname<-colnames(cox.data)[-c(1,2)]
  if(method=="Univariate"){
    GeneCoefResult<-c()
    for(i in 1:length(gname)){
      x<-numdata[,2+i];
      if(sd(x)==0){next
      }else{
        #attach(cox.data);
        coxph.fit<-coxph(Surv(Survival,Events)~x,data = cox.data)
        GeneCoefResult<-getPrintFile(coxph.fit,gname[i],GeneCoefResult);
        #detach(cox.data)
      }
    }
  }else if(method=="Multivariate"){
    GeneCoefResult<-c();
    data<-numdata;
    #attach(cox.data);
    nn<-dim(numdata)[2];
    coxph.fit<-coxph(Surv(Survival,Events)~data[,3:nn],data =cox.data );
    GeneCoefResult<-getPrintFile(coxph.fit,gname,GeneCoefResult);
    #detach(cox.data)
  }
  #UpGeneCoefResult <- GeneCoefResult [which(GeneCoefResult[,3]>1), ,drop=FALSE]
  #LoGeneCoefResult <- GeneCoefResult [which(GeneCoefResult[,3]<1), , drop=FALSE]
  return(GeneCoefResult)
}
getKMdata<-function(subprof2,clin2,label2){

  subprof2<-t(subprof2);#colnames(subprof2)<-paste("gene.",as.character(subprof2[1,]),sep="");subprof2<-subprof2[-1,]
  colnames(subprof2)<-paste("gene.",as.character(colnames((subprof2))),sep="");#subprof2<-subprof2[-1,]
  samples<-rownames(subprof2);
  inter<-intersect(as.character(clin2[,1]),samples)
  mm<-match(as.character(inter),as.character(clin2[,1]))
  nn<-match(as.character(inter),samples)

  newsubprof<-subprof2[nn[which(is.na(nn)==FALSE)],,drop=FALSE];n1<-dim(newsubprof)[2]
  clinInfor<-clin2[mm[which(is.na(mm)==FALSE)],];n2<-dim(clinInfor)[2]
  KMdata<-data.frame(newsubprof,clinInfor[,-1])[,c((n1+1):(n1+2))]
  interT<-intersect(label2[,1],rownames(KMdata))
  mm2<-match(as.character(interT),rownames(KMdata))
  mm2<-mm2[which(is.na(mm2)==FALSE)]
  nn2<-match(as.character(interT),label2[,1])
  nn2<-nn2[which(is.na(nn2)==FALSE)]
  KMdata<-cbind(label2[nn2,],KMdata[mm2,])[,c(1,3:4,2)]

  rownames(KMdata)<-NULL;colnames(KMdata)<-c("Sample","Survival","Events","group")
  KMdata[,"Events"]<-as.numeric(KMdata[,"Events"])
  KMdata[,"group"]<-as.numeric(KMdata[,"group"])-1
  KMdata[,"Events"]<-as.numeric(KMdata[,"Events"])+1
  KMdata<-as.data.frame(KMdata)
  KMdata<-KMdata[,-1]
  return(KMdata)
}
getCutoffValue<-function(coxRes,subprof){
  coef<-coxRes[,c("geneName","coef"),drop=FALSE]
  inter<-intersect(gsub("gene.","",as.character(coef[,1])),as.character(rownames((subprof))))
  nn<-match(as.character(inter),gsub("gene.","",as.character(coef[,1])))
  nn<-nn[which(is.na(nn)==FALSE)]
  nnn<-match(as.character(inter),as.character(rownames((subprof))))
  nnn<-nnn[which(is.na(nnn)==FALSE)]
  coefExp<-c();numprof<-apply(unname(subprof),c(1,2),as.numeric)
  for(i in 1:dim(numprof)[2]){
    coefExp<-c(coefExp,sum(as.numeric(as.character(coef[nn,2]))*as.numeric(as.character(numprof[nnn,i]))))
  }

  median<-median(coefExp)
  return(median)
}
getCoefExpCluster<-function(coxRes,subprof2,subprof){
  coef<-coxRes[,c("geneName","coef"), drop=FALSE ]
  inter<-intersect(gsub("gene.","",as.character(coef[,1])),as.character(rownames((subprof2))))
  nn<-match(as.character(inter),gsub("gene.","",as.character(coef[,1])))
  nn<-nn[which(is.na(nn)==FALSE)]
  nnn<-match(as.character(inter),as.character(rownames((subprof2))))
  nnn<-nnn[which(is.na(nnn)==FALSE)]
  coefExp2<-c();numprof<-apply(unname(subprof2),c(1,2),as.numeric)
  for(i in 1:dim(numprof)[2]){
    coefExp2<-c(coefExp2,sum(as.numeric(as.character(coef[nn,2]))*as.numeric(as.character(numprof[nnn,i]))))
  }
  cutoff<-getCutoffValue(coxRes,subprof)
  coefExp22<-coefExp2
  coefExp22[which(coefExp2>=cutoff)]<-1
  coefExp22[which(coefExp2<cutoff)]<-0
  return(cbind(colnames(subprof2),coefExp22))
}


#' @title survcell
#' @description Function `survcell` draws Kaplan–Meier curves for survival in the above-median and below-median groups for cell risk score. The cell risk score is calaulated by the weighted mean of cells driven by a gene mutation, where the  weight of cells is estimated by the "Univariate" or "Multivariate" cox.
#' @param gene Somatic mutant gene name
#' @param mutcell The result of `mutcorcell` function
#' @param cellmatrix Cell abundance matrix
#' @param method Method must be one of "Univariate" and "Multivariate". The coefficient of cells for risk score are estimated by "Univariate" or "Multivariate" cox proportional risk regression model on cell abundance matrix and overall survival data..
#' @param surv Surv is the survival data, the first column is the sample name, the second column is the survival time, and the third is the survival event.
#' @return Kaplan–Meier curves
#' @importFrom survival coxph
#' @importFrom survival survfit
#' @importFrom survival survdiff
#' @importFrom survival Surv
#' @importFrom grDevices rgb
#' @importFrom survminer ggsurvplot
#' @export
#' @examples
#' mutcell<-GetExampleData("mutcell") # The result of `mutcorcell` function.
#' cellmatrix<-GetExampleData("cellmatrix") # Cell abundance matrix
#' surv<-GetExampleData("surv") # The survival data
#' survcell(gene ="TP53",mutcell=mutcell,cellmatrix=cellmatrix,surv=surv)
survcell <- function(gene,mutcell,cellmatrix,surv,method="Multivariate") {
  mutcell<-mutcell$mut_cell
  cellgene<-names(which(mutcell[gene,]==1))
  colnames(surv)<-c("Samples","Survival","Events")
  cellmatrix.gene<-cellmatrix[cellgene,]
  if (method=="Multivariate") {
    cox.gene<-getUniOrMultiCOXAnalysis(subprof = cellmatrix.gene,clin = surv,method ="Multivariate" )
  }else if (method=="Univariate") {
    cox.gene<-getUniOrMultiCOXAnalysis(subprof = cellmatrix.gene,clin = surv,method ="Univariate" )
  }else {
    print("method must be one of 'Univariate'and'Multivariate'.")
  }
  label<-getCoefExpCluster(coxRes = cox.gene,subprof2 = cellmatrix,subprof = cellmatrix.gene)
  rownames(label)<-label[,1]
  KMdata<-getKMdata(subprof2 = cellmatrix,clin2 =surv ,label2 = label)
  fit <- survfit(Surv(Survival, Events) ~ group, data = KMdata)
  ggsurvplot(fit,pval = TRUE,data = KMdata)
}