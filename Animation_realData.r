
source('Data_real.r')

postate<-read.csv(file="KappaState2019-05-16.txt",sep="",header=FALSE)

end=max(dim(postate) )
burnin = 7000 #30000
mstate<-apply(postate[burnin:end,],2,mean)


#max( ax[-which( is.na(ax) ) ] )
#min( ax[-which( is.na(ax) ) ] )

pl_state = mstate
pl_state[which(pl_state>1.5)] = 17      # BM
pl_state[which(pl_state<1.5)] = 19      # OU
#pl_state[which( is.na(ax) ) ] = NA

pl_state = matrix(c(pl_state),ncol=5)

break_seconds = 1.5

#animation_plot<-function( break_seconds)
#{
     plot(0,0,xlim=c(-1,57),ylim=c(-1,45),xlab="x-dimension",ylab="y-dimension",col='white')

       for(qq in 1:46)
      {
          lines( ax[qq,1],ay[qq,1],type="p",pch= pl_state[qq,1] ,col="red",cex=1.5 )
          lines( ax[qq,2],ay[qq,2],type="p",pch= pl_state[qq,2],col="blue",cex=1.5  )
          lines( ax[qq,3],ay[qq,3],type="p",pch= pl_state[qq,3],col="forestgreen",cex=1.5  )
          lines( ax[qq,4],ay[qq,4],type="p",pch= pl_state[qq,4],col="black" ,cex=1.5 )
          lines( ax[qq,5],ay[qq,5],type="p",pch= pl_state[qq,5],col="purple" ,cex=1.5 )

       Sys.sleep(break_seconds)
          lines( ax[qq,1],ay[qq,1],type="p",pch= pl_state[qq,1],col="white",cex=2  )
          lines( ax[qq,2],ay[qq,2],type="p",pch= pl_state[qq,2],col="white",cex=2  )
          lines( ax[qq,3],ay[qq,3],type="p",pch= pl_state[qq,3],col="white",cex=2  )
          lines( ax[qq,4],ay[qq,4],type="p",pch= pl_state[qq,4],col="white",cex=2  )
          lines( ax[qq,5],ay[qq,5],type="p",pch= pl_state[qq,5],col="white",cex=2  )
       }
#}



