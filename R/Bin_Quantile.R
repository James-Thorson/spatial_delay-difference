Bin_Quantile <-
function(Obj, Nregions=4){
  Nbreaks = Nregions - 1
  Breaks = quantile(Obj, prob=seq(0,1,length=Nbreaks+2))
  Region = sapply( Obj, FUN=function(Num){ sum(Num>=Breaks) })
    Region = ifelse(Region==(Nregions+1),Nregions,Region)
    if( is.array(Obj) ) Region = array( Region, dim=dim(Obj))
  Return = list("Nregions"=Nregions, "Region"=Region, "Lwr"=Breaks[-c(Nbreaks+2)], "Upr"=Breaks[-1])
  return(Return)
}
