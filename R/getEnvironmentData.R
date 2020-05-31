initializeSMDIC<-function(){
  utils::data("envData",package="SMDIC")
}

Getenvir<-function(envData){

  if(!exists("envData")) initializeSMDIC()
  return(get(envData,envir=envData))
}
