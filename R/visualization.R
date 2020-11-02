#' @title plotwaterfall
#' @description Function ‘plotwaterfall‘ plots the waterfall for mutation genes which drive immune cells.
#' @param maffile The name of mutation annotation file (MAF) format data. It must be an absolute path or the name  relatived to the current working directory.
#' @param mutcell.summary The result of `mutcellsummary` function
#' @param cellnumcuoff a threshold value (4 as the default value). The mutation genes which drive at least "cellnumcuoff" cells are retained for drawing an waterfall.
#' @param fontSize font size for gene names. Default 0.8.
#' @param showTumorSampleBarcodes logical to include sample names.
#' @param showTitle Default TRUE
#' @param colors named vector of colors for each Variant_Classification.
#' @importFrom maftools oncoplot
#' @importFrom maftools read.maf
#' @export
#' @examples
#' # get result of `exp2cell` funtion
#' cellmatrix<-GetExampleData("cellmatrix")
#'
#' #get the binary mutations matrix,
#' mutmatrix<-GetExampleData("mutmatrix")
#'
#' # get the result of `mutcorcell` funtion
#' mutcell<-GetExampleData("mutcell")
#'
#' #perform the function mutcellsummary
#' summary<-mutcellsummary(mutcell = mutcell,mutmatrix = mutmatrix,cellmatrix=cellmatrix)
#'
#' #dir is the name of mutation annotation file (MAF) format data.
#' #It must be an absolute path or the name relatived to the current working directory.
#' # maf<-"dir"
#'
#' # mutcell.summary is the result of function mutcellsummary
#'
#' #plot the waterfall for mutation genes which drive immune cells
#' \donttest{plotwaterfall(maffile = maf,mutcell.summary = summary,cellnumcuoff =0)}
plotwaterfall <- function(maffile,mutcell.summary,cellnumcuoff=3,fontSize = 0.8,showTumorSampleBarcodes=F,showTitle = TRUE,colors = NULL) {
  ## requireNamespace("maftools")|| stop("package maftools is required,please install package maftools")
  gene.top<-mutcell.summary[which(mutcell.summary[,"count"]>=cellnumcuoff),1]
  maf<-read.maf(maf = maffile)
  oncoplot(maf = maf,
           genes = gene.top,
           fontSize = fontSize ,
           showTumorSampleBarcodes = showTumorSampleBarcodes,
           showTitle = showTitle,
           colors=colors)
}

#' @title plotCoocMutex
#' @description Function `plotCoocMutex` plots the co-occurrence and mutual exclusivity plots for mutation genes which drive immune cells.
#' @param maffile The name of mutation annotation file (MAF) format data. It must be an absolute path or the name  relatived to the current working directory.
#' @param mutcell.summary The result of `mutcellsummary` function
#' @param cellnumcuoff A threshold value (4 as the default value). The mutation genes which drive at least "cellnumcuoff" cells are retained for drawing a co-occurrence and mutual exclusivity plots.
#' @param fontSize cex for gene names. Default 0.8
#' @importFrom maftools somaticInteractions
#' @importFrom maftools read.maf
#' @references Gerstung M, Pellagatti A, Malcovati L, et al. Combining gene mutation with gene expression data improves outcome prediction in myelodysplastic syndromes. Nature Communications. 2015;6:5901. doi:10.1038/ncomms6901.
#' @export
#' @examples
#' # get the result of `exp2cell` funtion
#' cellmatrix<-GetExampleData("cellmatrix")
#'
#' #get the binary mutations matrix,
#' mutmatrix<-GetExampleData("mutmatrix")
#'
#' # get the result of `mutcorcell` funtion
#' mutcell<-GetExampleData("mutcell")
#'
#' #perform the function mutcellsummary
#' summary<-mutcellsummary(mutcell = mutcell,mutmatrix = mutmatrix,cellmatrix=cellmatrix)
#'
#' #dir is the name of mutation annotation file (MAF) format data.
#' #It must be an absolute path or the name relatived to the current working directory.
#' maf<-"dir"
#' #plot the co-occurrence and mutual exclusivity plots for mutation genes which drive immune cells.
#' \donttest{plotCoocMutex(maffile = maf,mutcell.summary = summary,cellnumcuoff =0)}
plotCoocMutex <- function(maffile,mutcell.summary,cellnumcuoff=3,fontSize = 0.8) {
  ## requireNamespace("maftools")|| stop("package maftools is required,please install package maftools")
  gene.top<-mutcell.summary[which(mutcell.summary[,"count"]>=cellnumcuoff),1]
  maf<-read.maf(maf = maffile)
  somaticInteractions(maf = maf,
                      genes = gene.top,
                      pvalue = c(0.05, 0.1),
                      fontSize=fontSize
  )
}


#' @title heatmapcell
#' @description A function to draw clustered heatmaps for the cells driven by a somatic mutation.
#' @param gene Somatic mutant gene name
#' @param cellmatrix Cell abundance matrix, cellmatrix is the result of function `exp2cell`.
#' @param mutcell A list, mutcell is the result of function `mutcorcell`.
#' @param mutmatrix A binary mutations matrix, which can not only come from the maf2matrix function, but also any binary mutations matrix, in which 1 represents any mutation occurs in a particular gene in a particular sample, otherwise the element is 0.
#' @param annotation_colors list for specifying annotation_row and annotation_col track colors manually. It is possible to define the colors for only some of the features. Check examples for details.
#' @param annotation_row data frame that specifies the annotations shown on left side of the heatmap. Each row defines the features for a specific row. The rows in the data and in the annotation are matched using corresponding row names. Note that color schemes takes into account if variable is continuous or discrete.
#' @param annotation_col similar to annotation_row, but for columns.
#' @param title The title of the plot
#' @param color vector of colors used in heatmap.
#' @param show_rownames	boolean specifying if column names are be shown.
#' @param show_colnames boolean specifying if column names are be shown.
#' @importFrom pheatmap pheatmap
#' @export
#' @examples
#' #get the result of `mutcorcell` function.
#' mutcell<-GetExampleData("mutcell")
#'
#' #get cell abundance matrix which is the result of exp2cell function
#' cellmatrix<-GetExampleData("cellmatrix")
#'
#' #get the binary mutations matrix
#' mutmatrix<-GetExampleData("mutmatrix")
#'
#' # plot significant up-regulation or down-regulation cells heat map specific for breast cancer
#' heatmapcell(gene = "TP53",mutcell = mutcell,cellmatrix = cellmatrix,mutmatrix = mutmatrix)
heatmapcell <- function(gene,mutcell,cellmatrix,mutmatrix,title = NA,show_rownames=TRUE,show_colnames = FALSE,annotation_colors = NA,annotation_row = NA,annotation_col = NA,color =NA) {
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
  #bk<-c(seq(-2,-0.1,by=0.1),seq(0,2,by=0.1))
  #bkcolor<-c(colorRampPalette(colors = c("#4575B4","white"))(length(bk)/2),colorRampPalette(colors = c("white","#D73027"))(length(bk)/2))
  ann_colors = list(
    genemut = c("#F9A452","#4B86C5")
  )
  annotation_col = data.frame(
    mutstat = factor(mutwithorder)
  )
  pheatmap(cell.gene_ordmut_zscore,scale ="none",show_rownames = show_rownames,show_colnames = show_colnames,cluster_cols = FALSE,annotation_col =annotation_col,annotation_row = annotation_row,annotation_colors = ann_colors,main=title,color=color)
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
        coxph.fit<-survival::coxph(Surv(Survival,Events)~x,data = cox.data)
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
  #nn<-match(gsub("gene.","",as.character(coef[,1])),as.character(subprof[,1]))
  nn<-match(as.character(inter),gsub("gene.","",as.character(coef[,1])))
  nn<-nn[which(is.na(nn)==FALSE)]
  #nnn<-match(as.character(subprof[,1]),gsub("gene.","",as.character(coef[,1])))
  nnn<-match(as.character(inter),as.character(rownames((subprof))))
  nnn<-nnn[which(is.na(nnn)==FALSE)]
  coefExp<-c();numprof<-apply(unname(subprof),c(1,2),as.numeric)
  for(i in 1:dim(numprof)[2]){
    coefExp<-c(coefExp,sum(as.numeric(as.character(coef[nn,2]))*as.numeric(as.character(numprof[nnn,i]))))
  }

  median<-median(coefExp)
  #print(median)
  #median<-quantile(coefExp,  probs = c(0.75))
  return(median)
}

getCoefExpCluster<-function(coxRes,subprof2,subprof){
  coef<-coxRes[,c("geneName","coef"), drop=FALSE ]
  inter<-intersect(gsub("gene.","",as.character(coef[,1])),as.character(rownames((subprof2))))
  #nn<-match(gsub("gene.","",as.character(coef[,1])),as.character(subprof2[,1]))
  nn<-match(as.character(inter),gsub("gene.","",as.character(coef[,1])))
  nn<-nn[which(is.na(nn)==FALSE)]
  #nnn<-match(as.character(subprof2[,1]),gsub("gene.","",as.character(coef[,1])))
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
#' @param palette the color palette to be used. Allowed values include "hue" for the default hue color scale; "grey" for grey color palettes; brewer palettes e.g. "RdBu", "Blues", ...; or custom color palette e.g. c("blue", "red"); and scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty". See details section for more information. Can be also a numeric vector of length(groups); in this case a basic color palette is created using the function palette.
#' @param pval logical value, a numeric or a string. If logical and TRUE, the p-value is added on the plot. If numeric, than the computet p-value is substituted with the one passed with this parameter. If character, then the customized string appears on the plot.
#' @param legend.title legend title.
#' @param legend.labs character vector specifying legend labels. Used to replace the names of the strata from the fit. Should be given in the same order as those strata.
#' @param color color to be used for the survival curves.If the number of strata/group (n.strata) = 1, the expected value is the color name. For example color = "blue".If n.strata > 1, the expected value is the grouping variable name. By default, survival curves are colored by strata using the argument color = "strata", but you can also color survival curves by any other grouping variables used to fit the survival curves. In this case, it's possible to specify a custom color palette by using the argument palette.
#' @param title the title of the survival curve
#' @return Kaplan–Meier curves
#' @importFrom survival coxph
#' @importFrom survival survfit
#' @importFrom survival survdiff
#' @importFrom survival Surv
#' @importFrom grDevices rgb
#' @importFrom survminer ggsurvplot
#' @importFrom stats sd
#' @export
#' @examples
#' # get the result of `mutcorcell` function.
#' mutcell<-GetExampleData("mutcell")
#'
#' # get cell abundance matrix which is the result of exp2cell function
#' cellmatrix<-GetExampleData("cellmatrix")
#'
#' # get survival data
#' surv<-GetExampleData("surv")
#'
#' #draw Kaplan–Meier curves
#' survcell(gene ="TP53",mutcell=mutcell,cellmatrix=cellmatrix,surv=surv)
survcell <-
  function(gene,
           mutcell,
           cellmatrix,
           surv,
           method = "Multivariate",
           legend.title = "Strata",
           legend.labs = c("group=0", "group=1"),
           palette = c("#E7B800", "#2E9FDF"),
           color=NULL,
           pval = TRUE,
           title=NULL) {
    mutcell <- mutcell$mut_cell
    cellgene <- names(which(mutcell[gene, ] == 1))
    colnames(surv) <- c("Samples", "Survival", "Events")
    cellmatrix.gene <- cellmatrix[cellgene, ]
    if (method == "Multivariate") {
      cox.gene <-
        getUniOrMultiCOXAnalysis(subprof = cellmatrix.gene,
                                 clin = surv,
                                 method = "Multivariate")
    } else if (method == "Univariate") {
      cox.gene <-
        getUniOrMultiCOXAnalysis(subprof = cellmatrix.gene,
                                 clin = surv,
                                 method = "Univariate")
    } else {
      print("method must be one of 'Univariate'and'Multivariate'.")
    }
    label <-
      getCoefExpCluster(coxRes = cox.gene,
                        subprof2 = cellmatrix,
                        subprof = cellmatrix.gene)
    rownames(label) <- label[, 1]
    KMdata <-
      getKMdata(subprof2 = cellmatrix,
                clin2 = surv ,
                label2 = label)

    fit <- survfit(Surv(Survival, Events) ~ group, data = KMdata)
    ggsurvplot(fit,
               data = KMdata,
               palette = palette,
               legend.title=legend.title,
               legend.labs=legend.labs,
               color=color,
               pval = pval,
               title=title,ggtheme =theme_light()+ theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(plot.title = element_text(hjust=0.5,vjust=0.5)))
  }

