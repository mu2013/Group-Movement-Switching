
## input
## lambda_sw switching rate, ou-bm or bm-ou.
## pro_pon_ind  the indicator vector of states, pn - potential switch, OB - ou to bm, BO - bm to ou.

## output
## the loglikelihood given the lambdas.

lam_lik<-function( lambda_sw, pro_pon_ind)
{
  lambda_BO = lambda_sw[2]
  lambda_OB = lambda_sw[1]

  p_BO = lambda_BO /Kappa #
  p_OB = lambda_OB /Kappa # 

  all_pon = pro_pon_ind[ which( pro_pon_ind != 'sp') ]
  pn_inx = which(all_pon=='pn')
  OB_inx = which(all_pon=='OB')
  BO_inx = which(all_pon=='BO')

  num_pon = all_pon
  num_pon[pn_inx] = 0
  num_pon[OB_inx] = -1
  num_pon[BO_inx] = 1
  num_pon = as.numeric(num_pon)

  cum_pon =cumsum(num_pon) + nai  ## if nai  in cum_pon, all animal doing OU at that time.

  num_pon[OB_inx] = p_OB
  num_pon[BO_inx] = p_BO
  num_pon[pn_inx] = 1 - (p_BO*(nai-cum_pon[pn_inx])  + p_OB*(cum_pon[pn_inx]))

  sum(log(num_pon))
}




