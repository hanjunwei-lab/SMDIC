#' @title Get the example data
#' @description Get the example data from SMDIC package.
#' @usage GetExampleData(exampleData)
#' @param exampleData A character, should be one of "exp.example", "cellmatrix", "mutcell", "mutmatrix", "surv".
#' @export
#' @details The function `GetExampleData(ExampleData = "mutmatrix)")` obtains the mutations matrix
#' @references Subramanian, A., Tamayo, P., Mootha, V.K., Mukherjee, S., Ebert, B.L., Gillette, M.A., Paulovich, A., Pomeroy, S.L., Golub, T.R., Lander, E.S. et al. (2005) Gene set enrichment analysis: a knowledgebased approach for interpreting genome-wide expression profiles. Proc Natl Acad Sci U S A, 102, 15545-15550.
GetExampleData<-function(exampleData){

  if(!exists("envData")) {
    envData<-initializeSMDIC()
  }

  if (exampleData=="mutmatrix")
  {
    dataset<- get("mutmatrix",envir=envData)
    return(dataset)
  }

  if (exampleData=="cellmatrix")
  {
    dataset<- get("cellmatrix",envir=envData)
    return(dataset)
  }
  if (exampleData=="exp.example")
  {
    dataset<- get("exp.example",envir=envData)
    return(dataset)
  }

  if (exampleData=="mutcell")
  {
    dataset<- get("mutcell",envir=envData)
    return(dataset)
  }

  if (exampleData=="surv")
  {
    dataset<- get("surv",envir=envData)
    return(dataset)
  }


}
