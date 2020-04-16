
source('Setup.r')
source('GenStates.r')
source('lam_lik.r')
set.seed(3809)

iter =  50000

nstate = all_s[,-1] 
nsw_ind = sw_ind
ntime = matrix( all_t, nrow=1)

ostate = nstate
osw_ind = nsw_ind
otime = ntime
osp_state = ostate[which(osw_ind=='SP'),]

## propose a sequence of states which are not the true states
uplen = alen-1 #49
start_time =2
pro_st = GenStates( swlamda,deltat,  uplen,start_time, osp_state,ostate,otime,osw_ind)

ostate = pro_st$pro_state                    ## states at actual switching and sampling points
otime = pro_st$pro_time                      ## time at actual switching and sampling points
osw_ind = pro_st$pro_sw_ind                  ## switch indicators at actual swtiching and sampling points
osp_state = ostate[which(osw_ind=='SP'),]    ## states at sampling points

opon_ind = pro_st$pro_pon_ind                ## switch indicators at all proposed swiching and sampling points 
opon_time = pro_st$pro_pon_time              ## time at all proposed swiching and sampling points 
opon_probAll = pro_st$pro_pon_probAll        ## probability at all proposed swiching and sampling points

 #  sumstate = apply(ostate,1,sum)[-1]
 #  test = opon_ind[which(opon_ind!='pn')]

 #  test[which(test=='sp')] = 0
 #  test[which(test=='OB')] = 1
 #  test[which(test=='BO')] = -1
 #  test = as.numeric(test)

 # if( sum( sumstate-c(cumsum(test)+5) ) !=0 )
 # { break }


ptm <- proc.time()

#### start iteration
for(j in 1:iter)
#for(j in 25001:40000)
{

  for(jj in 1:5){  ## state estimation loop 
 ## propose states list for all animals with switching time
     uplen = 3  #3 is the minimum, in fact we only update the mid point, the start and end states are fixed
     #start_time = sample( seq(2,94,by=2),1)
     start_time = sample( seq(2,(mxsamp - (uplen-1)*deltat),by=deltat ),1)
     end_time = start_time+(uplen-1)*deltat
     
     pro_st = GenStates( swlamda,deltat,  uplen,start_time, osp_state,ostate,otime,osw_ind,opon_ind,opon_time)
      nstate = pro_st$pro_state
      ntime = pro_st$pro_time
      nsw_ind = pro_st$pro_sw_ind
      npon_ind = pro_st$pro_pon_ind
      npon_time = pro_st$pro_pon_time

      iterbreak = pro_st$iterbreak

    if(iterbreak == 1)
    {
      old_par = c(alpha,beta,rho,sigma,theta,Bsigma,swlamda)
      par =  c(alpha,beta,rho,sigma,theta,Bsigma,swlamda)

      resst = Run_KF(par,old_par,osp_state, ostate,otime,osw_ind,opon_ind,opon_time,olik, nstate,ntime,nsw_ind, npon_ind,npon_time ,staccept)
      olik = resst$olik

      osp_state = resst$SPstate
      #swlamda = resst$par[8:9]

      ostate = resst$ostate
      otime = resst$otime
      osw_ind = resst$osw_ind

      opon_ind = resst$opon_ind
      opon_time = resst$opon_time

      staccept = resst$accept
    } 
  }

## propose new switching rates ##########
nlambdaOB <-  rnorm(1,swlamda[1],0.05) 
nlambdaBO <-  rnorm(1,swlamda[2],0.1) 
nswlamda=c(nlambdaOB ,nlambdaBO)

if( 0.02<nlambdaOB & nlambdaOB <nlambdaBO & nai*nlambdaBO<Kappa)
  {
#########################################  switching parameter

      lam_oldlik = lam_lik( swlamda , opon_ind)
      lam_newlik = lam_lik( nswlamda , opon_ind)

      lam_sigHR <- exp(lam_newlik - lam_oldlik ) 
        if(runif(1) < lam_sigHR) 
         {
           swlamda = nswlamda
           lamaccept = lamaccept+1
          } 

#########################################  
   }


#########################################  diffusion parameter

  nalpha<-  rnorm(1,alpha,prop.alpha)   # 0.1
  old_par = c(alpha,beta,rho,sigma,theta,Bsigma)
  par = c(nalpha,beta,rho,sigma,theta,Bsigma)
  resa = Run_KF(par,old_par,osp_state, ostate,otime,osw_ind,opon_ind,opon_time,olik, ostate,otime,osw_ind,npon_ind,npon_time,alaccept)
         
  alpha = resa$par[1]
  olik = resa$olik
  alaccept = resa$accept
  
    nrho<-( rnorm(1,rho,prop.rho) )   
    old_par = c(alpha,beta,rho,sigma,theta,Bsigma)
    par = c(alpha,beta,nrho,sigma,theta,Bsigma)
    resr = Run_KF(par,old_par,osp_state, ostate,otime,osw_ind,opon_ind,opon_time,olik, ostate,otime,osw_ind,npon_ind,npon_time,rhaccept)
          
    rho = resr$par[3]
    olik = resr$olik
    rhaccept = resr$accept

   nsigma<- ( rnorm(1,sigma,prop.sigma) )   
    old_par = c(alpha,beta,rho,sigma,theta,Bsigma)
    par = c(alpha,beta,rho,nsigma,theta,Bsigma)
    ress = Run_KF(par,old_par,osp_state, ostate,otime,osw_ind,opon_ind,opon_time,olik, ostate,otime,osw_ind,npon_ind,npon_time,sigaccept)
    sigma = ress$par[4]
    olik = ress$olik
    sigaccept = ress$accept 

   nBsigma<- ( rnorm(1,Bsigma,prop.Bsigma) )  
   old_par = c(alpha,beta,rho,sigma,theta,Bsigma)
   par = c(alpha,beta,rho,sigma,theta,nBsigma)
   resBs = Run_KF(par,old_par,osp_state, ostate,otime,osw_ind, opon_ind,opon_time,olik, ostate,otime,osw_ind,npon_ind,npon_time,Bsigaccept)
   Bsigma = resBs$par[7]
   olik = resBs$olik
   Bsigaccept = resBs$accept

## print output
     if(j%%2==0)
     {
     cat(file=filetheta4, alpha, "\n", append = TRUE)
     cat(file=filetheta6, sigma, "\n", append = TRUE)
     cat(file=filetheta7, rho, "\n", append = TRUE)
     cat(file=filetheta8, Bsigma, "\n", append = TRUE)
     cat(file=filetheta9, osp_state, "\n", append = TRUE)
     cat(file=filetheta10, swlamda, "\n", append = TRUE)
     }
    
     if(j%%5==0)
     {
     #cat("iteration",j,"alpha",alaccept,"beta",beaccept,"rho",rhaccept,"sigma",sigaccept,"theta",theaccept,"Bsigma",Bsigaccept,"State", staccept,"\n")
     #cat("iteration",j,"State", staccept,"alpha",alaccept,"beta",beaccept,"\n") # "sigma",sigaccept
     cat("iteration",j,"State", staccept,"alpha",alaccept,"rho",rhaccept,"Bsigma",Bsigaccept,"sigma",sigaccept,"lambdas",lamaccept,"iterkappa",pro_st$iterkappa,"\n")
     }
}

proc.time() - ptm


