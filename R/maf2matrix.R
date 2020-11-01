#' @title maf2matrix
#' @description Function `maf2matrix` use mutation annotation file (MAF) format data to build a binary mutations matrix.
#' @param maffile The name of mutation annotation file (MAF) format data. It must be an absolute path or the name relatived to the current working directory.
#' @param percent A threshold value(one percent as the default value). The genes with a given mutation frequency equal or greater than the threshold value are retained for the following analysis.
#' @param nonsynonymous Logical, tell if extract the non-silent somatic mutations (nonsense mutation, missense mutation, frame-shif indels, splice site, nonstop mutation, translation start site, inframe indels).
#' @return A binary mutations matrix, in which 1 represents any mutation occurs in a particular gene in a particular sample, otherwise the element is 0.
#' @importFrom utils read.delim
#'
#' @export
#'
#' @examples
#' #get path of the mutation annotation file.
#' maf = system.file('extdata', 'example.maf.gz', package = 'SMDIC')
#'
#' # perform function `maf2matrix`.
#' mutmatrix.example<-maf2matrix(maf)
maf2matrix <- function(maffile,
                       percent = 0.01,
                       nonsynonymous = TRUE) {
  maf <-
    read.delim(
      maffile,
      header = TRUE,
      comment.char = '#',
      stringsAsFactors = FALSE
    )
  mafneed <-
    maf[, c("Hugo_Symbol",
            "Tumor_Sample_Barcode",
            "Variant_Classification")]
  if (nonsynonymous == TRUE) {
    mafmissense <-
      mafneed[which(
        mafneed$Variant_Classification == "Missense_Mutation" |
          mafneed$Variant_Classification == "Frame_Shift_Del" |
          mafneed$Variant_Classification == "Frame_Shift_Ins" |
          mafneed$Variant_Classification == "In_Frame_Del" |
          mafneed$Variant_Classification == "Nonsense_Mutation" |
          mafneed$Variant_Classification == "In_Frame_Ins" |
          mafneed$Variant_Classification == "Splice_Site" |
          mafneed$Variant_Classification == "Nonstop_Mutation" |
          mafneed$Variant_Classification == "Translation_Start_Site"
      ),]
    mafselect <- mafmissense
  } else {
    mafselect <- mafneed
  }
  sample_id <- unique(mafselect$Tumor_Sample_Barcode)
  gene_id <- unique(mafselect$Hugo_Symbol)
  mafselectmatrix <-
    matrix(0, nrow = length(gene_id) , ncol = length(sample_id))
  rownames(mafselectmatrix) <- gene_id
  colnames(mafselectmatrix) <- sample_id
  for (i in 1:nrow(mafselect)) {
    gene_one <- mafselect$Hugo_Symbol[i]
    sample_one <- mafselect$Tumor_Sample_Barcode[i]
    mafselectmatrix[gene_one, sample_one] <- 1
  }
  sample_id <- substr(sample_id, 1, 15)
  sample_id <- gsub("-", ".", sample_id)
  colnames(mafselectmatrix) <- sample_id
  cutmutnum = ceiling(ncol(mafselectmatrix) * percent)
  mafmatrix <-
    mafselectmatrix[rowSums(mafselectmatrix) >= cutmutnum,]
  return(mafmatrix)
}
