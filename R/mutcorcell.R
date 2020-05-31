#' @title mutcorcell
#' @description  Function `mutcorcell` identifies somatic mutation-driven immune cells by comparing the cell abundance matrix and binary mutations matrix.
#' @param cellmatrix Cell abundance matrix.
#' @param mutmatrix A binary mutations matrix, in which 1 represents any mutation occurs in a particular gene in a particular sample, otherwise the element is 0.
#' @param samfdr.cutoff False Discovery Rate cutoff for output in significant immune cells
#' @param nperms Number of permutations used by SAM to estimate False Discovery Rates
#' @param fisher.cutoff False Discovery Rate(fisher.adjust=TRUR) or P-Value(fisher.adjust=FALSE) cutoff for Fisher's exact test
#' @param fisher.adjust Logical,tell if corrects p-values
#' @return A list of three matrices: a binary numerical matrix which shows the immune cells driven by somatic mutant gene; two numerical matrix which show the pvalue and fdr of the immune cells driven by somatic mutant gene.
#' @importFrom samr SAM
#' @export
#' @examples
#' cellmatrix<-GetExampleData("cellmatrix") # Cell abundance matrix
#' mutmatrix<-GetExampleData("mutmatrix") # A binary mutations matrix
#' mutcell<-mutcorcell(cellmatrix = cellmatrix,mutmatrix = mutmatrix)
#' # The summary for somatic mutations are produced by function `mutcellsummary`.
#' #summary<-mutcellsummary(mutcell = mutcell,mutmatrix = mutmatrix,cellmatrix=cellmatrix)
#' # The summary of the immune cells driven by a mutation are produced by function `gene2cellsummary`.
#' #genecellsummary<-gene2cellsummary(gene="TP53",mutcell=mutcell)
mutcorcell <- function(cellmatrix=cellmatrix,mutmatrix=mutmatrix,samfdr.cutoff=0.05,nperms=100,fisher.cutoff=0.05,fisher.adjust=T) {
  ## requireNamespace("samr")|| stop("package samr is required,please install package samr")
  intersample<-intersect(colnames(cellmatrix),colnames(mutmatrix))
  mutmatrix<-mutmatrix[,intersample]
  cellmatrix<-cellmatrix[,intersample]
  zscore<-t(apply(cellmatrix, 1, scale))
  mut_cell<-matrix(0,nrow = nrow(mutmatrix),ncol = nrow(cellmatrix))
  rownames(mut_cell)<-rownames(mutmatrix)
  colnames(mut_cell)<-rownames(cellmatrix)
  mut_cell_p<-matrix(0,nrow = nrow(mutmatrix),ncol = nrow(cellmatrix))
  rownames(mut_cell_p)<-rownames(mutmatrix)
  colnames(mut_cell_p)<-rownames(cellmatrix)
  mut_cell_fdr<-matrix(0,nrow = nrow(mutmatrix),ncol = nrow(cellmatrix))
  rownames(mut_cell_fdr)<-rownames(mutmatrix)
  colnames(mut_cell_fdr)<-rownames(cellmatrix)
  for (mut.no in 1:dim(mutmatrix)[1]) {
    print(paste(mut.no,":mut.no is",mut.no))
    label<-mutmatrix[mut.no,]+1
    if (sum(label==2)<=1) {
      next
    }
    samfit<-SAM(as.matrix(cellmatrix),label,resp.type="Two class unpaired",genenames=rownames(cellmatrix),nperms=nperms,fdr.output = samfdr.cutoff)
    if (class(samfit$siggenes.table$genes.up)=="matrix") {
      siggene.up<-cbind(samfit$siggenes.table$genes.up[,1],as.numeric(samfit$siggenes.table$genes.up[,7])/100);
    }else{
      siggene.up<-cbind(samfit$siggenes.table$genes.up[1],as.numeric(samfit$siggenes.table$genes.up[7])/100);
    }
    if (class(samfit$siggenes.table$genes.lo)=="matrix") {
      siggene.down<-cbind(samfit$siggenes.table$genes.lo[,1],as.numeric(samfit$siggenes.table$genes.lo[,7])/100);
    }else{
      siggene.down<-cbind(samfit$siggenes.table$genes.lo[1],as.numeric(samfit$siggenes.table$genes.lo[7])/100);
    }
    sig.up<-siggene.up[siggene.up[,2]<samfdr.cutoff,1];
    sig.down<-siggene.down[siggene.down[,2]<samfdr.cutoff,1];
    sig.all<-c(sig.up,sig.down)
    if (length(sig.all)==0) {
      next
    }
    genesets.sigpath<-cellmatrix[sig.all,]
    zscore.genesets.sigpath<-zscore[sig.all,]
    zscore.sigup<-NULL
    zscore.sigdown<-NULL

    if (length(sig.up)>1) {
      zscore.sigup<-zscore.genesets.sigpath[sig.up,]
      zscore.sigup[zscore.sigup<=2.0]<-0
      zscore.sigup[zscore.sigup>2.0]<-1
    }
    if(length(sig.up)==1){
      if (is.null(nrow(zscore.genesets.sigpath))) {
        zscore.sigup<-zscore.genesets.sigpath
      }else{
        zscore.sigup<-zscore.genesets.sigpath[sig.up,]
      }
      zscore.sigup[zscore.sigup<=2.0]<-0
      zscore.sigup[zscore.sigup>2.0]<-1
    }
    if (length(sig.down)>1) {
      zscore.sigdown<-zscore.genesets.sigpath[sig.down,]
      zscore.sigdown[zscore.sigdown>=(-2.0)]<-0
      zscore.sigdown[zscore.sigdown< (-2.0)]<-1
    }else if(length(sig.down)==1){
      if (is.null(nrow(zscore.genesets.sigpath))) {
        zscore.sigdown<-zscore.genesets.sigpath
      }else{
        zscore.sigdown<-zscore.genesets.sigpath[sig.down,]
      }
      zscore.sigdown[zscore.sigdown>=(-2.0)]<-0
      zscore.sigdown[zscore.sigdown< (-2.0)]<-1
    }
    zscore01<-rbind(zscore.sigup,zscore.sigdown)
    if (!is.null(nrow(zscore01))) {
      rownames(zscore01)<-c(as.character(sig.up),as.character(sig.down))
    }else{
      names(zscore01)<-c(as.character(sig.up),as.character(sig.down))
    }
    mutation<-which(mutmatrix[mut.no,]==1)
    notmutation<-which(mutmatrix[mut.no,]==0)
    p<-c()
    if (nrow(zscore01)!=1) {
      for(i in 1:nrow(zscore01)){
        diffsample<-which(zscore01[i,]==1)
        diffnotsample<-which(zscore01[i,]==0)
        a<-length(intersect(mutation,diffsample))
        b<-length(intersect(mutation,diffnotsample))
        c<-length(intersect(notmutation,diffsample))
        d<-length(intersect(notmutation,diffnotsample))
        Convictions<-matrix(c(a,c,b,d),nrow = 2,dimnames = list(c("mutation", "notmutation"),c("sample diff", "sample not diff")))
        fisher_test<-fisher.test(Convictions)
        p<-c(p,fisher_test$p.value)
      }
    }else{
      diffsample<-which(zscore01==1)
      diffnotsample<-which(zscore01==0)
      a<-length(intersect(mutation,diffsample))
      b<-length(intersect(mutation,diffnotsample))
      c<-length(intersect(notmutation,diffsample))
      d<-length(intersect(notmutation,diffnotsample))
      Convictions<-matrix(c(a,c,b,d),nrow = 2,dimnames = list(c("mutation", "notmutation"),c("sample diff", "sample not diff")))
      fisher_test<-fisher.test(Convictions,alternative = "greater")
      p<-fisher_test$p.value
    }

    if (length(p)==1) {
      if (p<fisher.cutoff) {
        print(as.character(sig.all))
        mut_cell[mut.no,as.character(sig.all)]<-1
        mut_cell_p[mut.no,as.character(sig.all)]<-p
      }else{
        next
      }
    }else{
      if (fisher.adjust==TRUE) {
        fdr<-p.adjust(p,method = "fdr")
        zscore.genesets.cutofffisher<-zscore.genesets.sigpath[which(fdr<fisher.cutoff),]
        fdr<-fdr[which(fdr<fisher.cutoff)]
        print(rownames(zscore.genesets.cutofffisher))
        mut_cell[mut.no,rownames(zscore.genesets.cutofffisher)]<-1
        mut_cell_fdr[mut.no,rownames(zscore.genesets.cutofffisher)]<-fdr
        p<-p[which(fdr<fisher.cutoff)]
        mut_cell_p[mut.no,rownames(zscore.genesets.cutofffisher)]<-p
      }else{
        zscore.genesets.cutofffisher<-zscore.genesets.sigpath[which(p<fisher.cutoff),]
        p<-p[which(p<fisher.cutoff)]
        print(rownames(zscore.genesets.cutofffisher))
        mut_cell[mut.no,rownames(zscore.genesets.cutofffisher)]<-1
        mut_cell_p[mut.no,rownames(zscore.genesets.cutofffisher)]<-p
        fdr<-p.adjust(p,method = "fdr")
        fdr<-fdr[which(p<fisher.cutoff)]
        mut_cell_fdr[mut.no,rownames(zscore.genesets.cutofffisher)]<-fdr
      }

    }
  }
  delrow<-which(rowSums(mut_cell)!=0)
  delcol<-which(colSums(mut_cell)!=0)
  mut_cell<-mut_cell[delrow,]
  mut_cell<-mut_cell[,delcol]
  mut_cell_p<-mut_cell_p[delrow,]
  mut_cell_p<-mut_cell_p[,delcol]
  mut_cell_fdr<-mut_cell_fdr[delrow,]
  mut_cell_fdr<-mut_cell_fdr[,delcol]
  return(list(mut_cell=mut_cell,mut_cell_p=mut_cell_p,mut_cell_fdr=mut_cell_fdr))
}

#' @title mutcellsummary
#' @description Function `mutcellsummary` is a generic function used to produce summaries of the results of `mutcorcell` function.
#' @param mutcell The result of `mutcorcell` funtion.
#' @param mutmatrix A binary mutations,in which 1 represents any mutation occurs in a particular gene in a particular sample, otherwise the element is 0. matrix
#' @param cellmatrix Cell abundance matrix
#' @return The result summaries have four columns. The first column is somatic mutant gene names, the second column is the immune cell names driven by the somatic mutation, the third column is the number of the immune cell, the fourth column is the mutation rate.
#' @export
#' @examples
#' mutcell<-GetExampleData("mutcell") # The result of `mutcorcell` funtion
#' cellmatrix<-GetExampleData("cellmatrix") # Cell abundance matrix
#' mutmatrix<-GetExampleData("mutmatrix") # A binary mutations matrix
#' summary<-mutcellsummary(mutcell = mutcell,mutmatrix = mutmatrix,cellmatrix=cellmatrix)
mutcellsummary <- function(mutcell,mutmatrix,cellmatrix) {
  mutcell<-mutcell$mut_cell
  intersample<-intersect(colnames(cellmatrix),colnames(mutmatrix))
  mutmatrix<-mutmatrix[,intersample]
  mutcell.run<-mutcell
  mutcell.summary<-matrix(0,nrow =nrow(mutcell.run),ncol = 4)
  mutcell.summary<-as.data.frame(mutcell.summary)
  for (i in 1:nrow(mutcell.run)) {
    mutcell.summary[i,1]<-rownames(mutcell.run)[i]
    mutcell.summary[i,2]<-paste0(names(which(mutcell.run[i,]==1)),collapse = ",")
    mutcell.summary[i,3]<-length(names(which(mutcell.run[i,]==1)))
    mutcell.summary[i,4]<-sum(mutmatrix[rownames(mutcell.run)[i],])/dim(mutmatrix)[2]
  }

  colnames(mutcell.summary)<-c("gene","cells","count","mut rate")
  mutcell.summary<-mutcell.summary[order(as.numeric(mutcell.summary[,3]),decreasing = TRUE),]
  return(mutcell.summary)
}

#' @title gene2cellsummary
#' @description Function `gene2cellsummary` is a generic function used to produce result summaries of the immune cells driven by a somatic mutation.
#' @param gene Somatic mutant gene name
#' @param method Method must be one of "xCell","ssGSEA" and "CIBERSORT".
#' @param mutcell The result of `mutcorcell` funtion.
#' @return A matrix shows the short name, full name, pvalue, fdr of the cells driven by a somatic mutation
#' @export
#' @examples
#' mutcell<-GetExampleData("mutcell") # The result of `mutcorcell` funtion.
#' genecellsummary<-gene2cellsummary(gene="TP53",mutcell=mutcell)
gene2cellsummary <- function(gene,method="xCell",mutcell) {
  mutcellp = mutcell$mut_cell_p
  mutcellfdr=mutcell$mut_cell_fdr
  mutcell = mutcell$mut_cell
  cellname<-names(which(mutcell[gene,]==1))
  resmatrix<-matrix(0,length(cellname),5)
  colnames(resmatrix)<-c("gene","cell name","Full name","pvalue","fdr")
  rownames(resmatrix)<-cellname
  resmatrix[,1]<-gene
  resmatrix[,2]<-cellname
  if (method=="xCell") {
    resmatrix[,3]<-cell64[cellname,2]
  }else if(method=="ssGSEA"){
    resmatrix[,3]<-cell24[cellname,2]
  }else {
    resmatrix[,3]<-cellname
  }
  resmatrix[,4]<-mutcellp[gene,cellname]
  resmatrix[,5]<-mutcellfdr[gene,cellname]
  return(resmatrix)
}

