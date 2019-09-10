#This code is for filtering out non-clonal mutations. 
#We will only keep clonal mutations around which have CCF >= 0.6
VAP.filter.clonal.0.6ccf = function(VAP.ccf){
  
  ccf.samps.bool = grepl(pattern = "ccf$",x = colnames(VAP.ccf)) #which columns end in ccf 
  
  
  #find max CCF across all columns which have CCF at end of name
  max = do.call(pmax, VAP.ccf[ccf.samps.bool])
  
  #Determine if the max value of the row across all CCF columns is >= 0.6 
  max.bool = max >= 0.6
  
  #select only columns which have a clonal mutation (one of the mutations is CCF >= 0.6)
  small.file = VAP.ccf[max.bool,]
  
  return(small.file)
}
