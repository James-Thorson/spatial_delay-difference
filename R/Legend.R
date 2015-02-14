Legend <-
function( Bin, Col, RowSet=NULL, Pcex=3, Tcex=3, Digits=3 ){
  if( is.null(RowSet) ) RowSet=1:Bin$Nregions
  plot( 1, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", frame.plot=FALSE, xaxt="n", yaxt="n", xaxs="i", yaxs="i", mar=c(0,0,0,0) )
  Y = seq( 0, 1, length=length(RowSet)+2)[-c(length(RowSet)+1:2)] + 1/length(RowSet)/2
  points( x=rep(0.1,length(RowSet)), y=Y, col=Col(Bin$Nregions)[RowSet], pch=20, cex=Pcex)
  text( x=rep(0.2,length(RowSet)), y=Y, labels=paste(formatC(Bin$Lwr[RowSet],format="f",digits=Digits),"to",formatC(Bin$Upr[RowSet],format="f",digits=Digits)), pos=4, cex=Tcex)
}
