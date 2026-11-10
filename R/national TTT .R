#a secondary version of all TTT functions that allow for the
#age and nativity distribution of each risk population to be
#included when calculating differential screening rates.

#this function returns default values for ttt which would be total population
#'@name def_ttt_nat_ag
#'@export
def_ttt_nat_ag<-function(){
  ttt_list<-vector("list", 9)
  names(ttt_list)<-c("AgeDistUS", "AgeDistNUS", "NRiskGrp", "FrcScrn",
                     "StartYr", "EndYr", "RRprg", "RRmu", "RRPrev")
  ttt_list[[1]]<-rep(1/11,11)
  ttt_list[[2]]<-rep(1/11,11)
  ttt_list[[3]]<-0 #size of population
  ttt_list[[4]]<-0 #fraction screened
  ttt_list[[5]]<-2018
  ttt_list[[6]]<-2050
  ttt_list[[7]]<-1
  ttt_list[[8]]<-1
  ttt_list[[9]]<-1
  return(ttt_list)
}
