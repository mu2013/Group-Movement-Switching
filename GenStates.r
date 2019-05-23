
library(mvtnorm)
library(Matrix)

set.seed(3809)

## input
## swlamda  switching rate, ou-bm or bm-ou  predined swlamda=c(0.1,0.5) 
## deltat  time distance between two sampling/observation points
## uplen  the length of updating fragment
## start_time the starting time of updating fragment
## osp_state the old states at all sampling/observation points
## ostate  the old state at all switching and sampling/observation points
## otime   the old time at all switching and sampling/observation points
## osw_ind the old indicator 'SW' or 'SP' for all switching and sampling/observation points
## opon_ind the old ndicator of all potential switching and actual switching 'OB' 'BO'  'pn'
## opon_time the old time at all potential switching and actual switching points.


## output
## pro_state is the proposed states at all switching and sampling/observation points
## pro_time is the proposed time at all switching and sampling/observation points
## pro_sw_ind is the indicator 'SW' or 'SP' for all switching and sampling/observation points
## pro_pon_ind is the indicator of all potential switching and actual switching 'OB' 'BO'  'pn'
## pro_pon_time is the proposed time at all potential switching and actual switching points.

GenStates<- function(swlamda,deltat,  uplen,start_time, osp_state,ostate,otime,osw_ind,  opon_ind,opon_time)
{
    nai = dim(osp_state)[2]

    ## update start from 4 and update at 6 and update end at 8, states for 4 and 8 are fixed
    end_time = start_time+(uplen-1)*deltat
    up_start = which(otime == start_time )
    up_end = which(otime == end_time )
    upNsamp = uplen

    upstate = ostate[ up_start, ]
    endstate = ostate[ up_end, ]

    #Kappa = 3.5

iterkappa = 1  ## count how many iterations to get the right state at the end point
iterbreak = 1  ## if it is iterbreak  we reject the new states and use the old states

repeat
{
    pon_ind = c()  ## 'pn' , 'sw'  , 'BO' , 'OB'
    pon_probAll = c()
    pon_list = c()

    pstate<- upstate  #rep(1,nai)   # state 1 attract to black dot     state 2 unattract to black dot (brownian motion)
    pall_s <- pstate
    psw_ind <- c('SP')
    pall_t <- c(0) ## uptime

    ouindex=which(pstate==1); n_ou = length(ouindex)
    bmindex=which(pstate==2); n_bm = length(bmindex)
          
    psamp = 1 
    pswt_old = 0
    #rest T ,how far from current time to the grid time point
    resT = deltat

    ## propose first switch time
      #Kappa = nai*max(swlamda) #nai*swlamda[1]
      P_ac_sw = ( n_ou*swlamda[1] + n_bm*swlamda[2] )/Kappa
      npsw = 0
      repeat{
        ac_sw = sample( c(0,1), 1, FALSE, prob=c(1-P_ac_sw,P_ac_sw) )
        npsw=npsw+1
        if(ac_sw==1) break
      }
      #pswt = sum( rexp(npsw,Kappa) ) #switch time
      pswt = 0
      po_t = rexp(npsw,Kappa)
      pswt = sum(po_t)

      resT = pswt - deltat
      move2 = pswt
      move1 = deltat
    
      ## another way of propose pontential switching points
      # Kappa = nai*swlamda[2]
      # n_posw = rpois(1, Kappa*deltat)
      # t_po_sw = sort( runif(n_posw,min=0,max=deltat)  )
      # P_ac_sw = ( n_ou*swlamda[1] + n_bm*swlamda[2] )/Kappa
      # ?? = sample(t_po_sw,size=1,prob=c(0.1,0.1,0.9) )

    repeat
    {
    ##case 1, switch animal state

      if(resT<0)     ## switching happen first , before sampling
      {
      
      if(psamp == upNsamp)
        {
          break
        }
      ### move to switch point
        t = move2  
      ### make the move
        ouindex=which(pstate==1)
        bmindex=which(pstate==2)
     
      ## pick the switch animal  according to the switching probability - swlamda[pstate]/de
        de=sum(swlamda[pstate])

        lnai = c(1:nai)
        lppstate = swlamda[pstate]/de
        SwitchAnimal = sample( lnai,1,FALSE, prob = lppstate )

      ## switch animal state 
        pstate[SwitchAnimal] = 3 - pstate[SwitchAnimal]
        

        ######### record all pontential and actual switching time  
        if(pstate[SwitchAnimal] ==1 ){
                 pon_ind = c( pon_ind, rep('pn', npsw-1 ),'BO' )
                 #pon_probAll = c( pon_probAll, rep((1-de/Kappa), npsw-1 ), length(bmindex)*swlamda[2]/Kappa )
                 pon_probAll = c( pon_probAll, rep((1-de/Kappa), npsw-1 ), swlamda[2]/Kappa )
          }else{
                 pon_ind = c( pon_ind, rep('pn', npsw-1 ),'OB' )
                 #pon_probAll = c( pon_probAll, rep((1-de/Kappa), npsw-1 ), length(ouindex)*swlamda[1]/Kappa)
                 pon_probAll = c( pon_probAll, rep((1-de/Kappa), npsw-1 ), swlamda[1]/Kappa)
          }
        pon_list = c( pon_list, po_t )  #pon_swt 
        ######### ######### ######### ######### ######### ######### ######### 


        ## demo state and switch time
        psw_ind= rbind(psw_ind,'SW')
        pall_s = rbind(pall_s,pstate)
        pall_t = cbind(pall_t,pswt+ sum(pswt_old) )
        pswt_old = c(pswt_old,pswt)

        ## update index   
        ouindex=which(pstate==1); n_ou = length(ouindex)
        bmindex=which(pstate==2); n_bm = length(bmindex)


        ## propose next switching time
        #pswt=rexp(1, sum( length(ouindex)*swlamda[1],length(bmindex)*swlamda[2]) ) 
        #Kappa = nai*max(swlamda)#[2]#nai*swlamda[1]
        P_ac_sw = ( n_ou*swlamda[1] + n_bm*swlamda[2] )/Kappa
        npsw = 0
        repeat{
          ac_sw = sample( c(0,1), 1, FALSE, prob=c(1-P_ac_sw,P_ac_sw) )
          npsw=npsw+1
          if(ac_sw==1) break
        }

        po_t = rexp(npsw,Kappa)
        pswt = sum(po_t)

        #pon_swt = cumsum(po_t) 
        #po_resT = pon_swt - abs(resT)

        move1 = abs(resT)
        move2 = pswt
        resT = pswt - abs(resT)

     }

    ###case 2, do not switch animal state
    if(resT > 0)            #sampling happan first  resT>0
     {
      ### move to sampling point
          t = move1  
          ouindex=which(pstate==1)
          bmindex=which(pstate==2)

    ### recording location and state information# 
          if(psamp == upNsamp)
                {
                  break
                } else
                {         
                  # demo state and sampling time
                  pdemot = psamp*2
                  psamp = psamp+1

                  psw_ind= rbind(psw_ind,'SP')
                  pall_s=  rbind(pall_s,pstate)
                  pall_t = cbind(pall_t,pdemot)
                 }

           move2 = abs(resT)
           move1 = deltat
           resT = resT - deltat

         if(psamp == upNsamp)
          {
          ## correct pontential sw time
           pon_listR= cumsum(pon_list)
           last_pon = cumsum(po_t) + pon_listR[ length(pon_listR) ]
           last_pon_ind = which( ( last_pon - (upNsamp-1)*deltat ) < 0 )

           pon_listR = c(pon_listR,last_pon[last_pon_ind] )
           pon_indR = c(pon_ind, rep( 'pn', length(last_pon_ind) ) )
           pon_probAll = c(pon_probAll, rep( (1-sum(swlamda[pstate])/Kappa), length(last_pon_ind) ) )
          }             
      }

    }

iterkappa = iterkappa +1

  if( identical( pall_s[dim(pall_s)[1],] , endstate) ) 
   {
    break
   } 

## if we cannot make the proposed state the same as the end time point in 5000 iterations, 
## we just reject the proposed states and use the old one. 
  if(iterkappa == 5000)   
   {
    iterbreak = 2
    break
   }
##########################################################################################
}
  

if( iterbreak == 1 )
{
  pall_t = pall_t + start_time
  otime_end= dim(otime)[2]


## put proposed fregments of switching time into the old time sequence 
  if(up_end < otime_end)
  {
    pall_time = cbind( t(otime[1:(up_start-1)]),pall_t,t(otime[(up_end+1):dim(otime)[2]] )  )
    pall_state = rbind( ostate[1:(up_start-1),], pall_s ,ostate[(up_end+1):dim(otime)[2],] )
    nsw_ind = c( osw_ind[1:(up_start-1)] , psw_ind  , osw_ind[(up_end+1):dim(otime)[2]] )
   } else if(up_end == otime_end)
   {
   pall_time = cbind( t(otime[1:(up_start-1)]), pall_t )
   pall_state = rbind( ostate[1:(up_start-1),], pall_s )
   nsw_ind = c( osw_ind[1:(up_start-1)] , psw_ind )
   }   
      

   if( missing(opon_ind) ) 
    {  
     ## add sampling/observation time and indicator to the potential switching list 
     ## we need the sampling/obsesrvation indicator to find out which fragment we are updating
     all_pont = c( pon_listR + start_time , seq(start_time,end_time,by=deltat ) )
     ## we need to sort the vector after we combined pontential switching time and sampling/observation time
      ooo = sort(all_pont, index.return=TRUE) 
      pro_pon_time =  ooo$x

     all_pon_ind = c(pon_indR , rep('sp',uplen) )  
     all_pro_pon_probAll = c(pon_probAll, rep('sp',uplen))  
      ## order all indicators as the index we derived from the order index of all potential time point 
      pro_pon_ind = all_pon_ind[ooo$ix]
      pro_pon_probAll = all_pro_pon_probAll[ooo$ix]
    } else{
    ## add sampling/observation time and indicator to the potential switching list 
    ## we need the sampling/obsesrvation indicator to find out which fragment we are updating
        if(up_end < otime_end)
        {
          all_pont = c( pon_listR + start_time , seq(start_time,end_time,by=deltat ) )
          all_pon_ind = c(pon_indR , rep('sp',uplen) ) 
          all_pon_probAll = c(pon_probAll, rep('sp',uplen))
          ooo = sort(all_pont, index.return=TRUE) 

           pon_start = which( opon_time ==start_time); pon_end = which(opon_time==end_time)
           if(pon_start==1){ ## if the updating start from time =2 which are the firt point, we can ignore the old time.
             pro_pon_ind = c( all_pon_ind[ooo$ix], opon_ind[(pon_end+1):length(opon_ind)] )
             pro_pon_time = c(  ooo$x, opon_time[(pon_end+1):length(opon_time)] )
             pro_pon_probAll = c(  all_pon_probAll[ooo$ix], opon_probAll[(pon_end+1):length(opon_probAll)] )
           }else{## if the updating do not start from time =2 which are the firt point, we need to add the old time back to the list.
             pro_pon_ind = c( opon_ind[1:(pon_start-1)], all_pon_ind[ooo$ix], opon_ind[(pon_end+1):length(opon_ind)] )
             pro_pon_time = c( opon_time[1:(pon_start-1)], ooo$x, opon_time[(pon_end+1):length(opon_time)] )
             pro_pon_probAll = c( opon_probAll[1:(pon_start-1)], all_pon_probAll[ooo$ix], opon_probAll[(pon_end+1):length(opon_probAll)] )
           }        
        }else if( up_end == otime_end )
        {
          all_pont = c( pon_listR + start_time , seq(start_time,end_time,by=deltat ) )
          all_pon_ind = c(pon_indR , rep('sp',uplen) ) 
          all_pon_probAll = c(pon_probAll, rep('sp',uplen) )
          ooo = sort(all_pont, index.return=TRUE) 

           pon_start = which( opon_time ==start_time)
           pro_pon_ind = c( opon_ind[1:(pon_start-1)], all_pon_ind[ooo$ix] ) #cbind( opon_ind[1:(pon_start-1)], all_pon_ind[ooo$ix] )
           pro_pon_time = c( opon_time[1:(pon_start-1)], ooo$x ) #cbind( opon_time[1:(pon_start-1)], ooo$x )
           pro_pon_probAll = c( pon_probAll[1:(pon_start-1)], all_pon_probAll[ooo$ix] )
        }
    }
 
 } else {

  pall_state =0 ; pall_time=0; nsw_ind = 0; pro_pon_ind = 0; pro_pon_time =0; pro_pon_probAll=0
 }


 return( list('pro_state'=pall_state,'pro_time'=pall_time,'pro_sw_ind'=nsw_ind,'pro_pon_ind'=pro_pon_ind,'pro_pon_time'=pro_pon_time,'pro_pon_probAll'=pro_pon_probAll,'iterbreak'=iterbreak,'iterkappa'=iterkappa) )
}


  # Kappa = max(swlamda)
  # pro_all_time=c(0)
  # repeat
  # {
  #  # ouindex=which(state==1)
  #  # bmindex=which(state==2)
  #  pon_sw_t = rexp( 1, Kappa*nai )
  #  pro_all_time = cbind( pro_all_time, + pon_sw_t )
  
  #  if(pon_sw_t > mxtime)
  #   {
  #     break
  #   } else {
  #   P_ac_sw = ( n_ou*swlamda[state] + n_bm*swlamda[state] )/Kappa*nai   
  #   de=sum(swlamda[state])
  #   SwitchAnimal =sample( c(1,2,3,4,5,6),1,FALSE, prob=c( swlamda[state[1]]/de,swlamda[state[2]]/de,swlamda[state[3]]/de,swlamda[state[4]]/de,swlamda[state[5]]/de,swlamda[state[6]]/de) )
  #   state[SwitchAnimal] = 3 - state[SwitchAnimal]   
  #   }         
  # }

