clear all 
set more off
use http://fmwww.bc.edu/repec/bocode/t/traindata.dta
set seed 90210

capture ssc install mixlogit 

mixlogit y, rand(price contract local wknown tod seasonal) group(gid) id(pid)

mata: 
real matrix drawb_betaW(beta, W) {
    return( mean(beta) + rnormal(1, cols(beta) ,0 ,1 ) * cholesky(W)' )
}
end

mata:
real matrix drawW_bbeta(beta, b) {
    v  = rnormal( cols(b) + rows(beta), cols(b), 0, 1)
    S1 = variance(beta :- b)
	S=invsym((cols(b)*I(cols(b)) + rows(beta)*S1)/(cols(b) + rows(beta)))
    L  = cholesky(S)
    R = (L*v')*(L*v')' / (cols(b) + rows(beta))
    return(invsym(R))
}
end

mata:
real scalar lncp(real rowvector beta_rn,
                 real rowvector b,
                 real matrix Winv,
                 real matrix ldetW,
                 real matrix y,
                 real matrix Xr,
                 real matrix cid)
{
    
    real scalar i, lnp, lnprior
    real matrix z, Xrp, yp, mus
 
    z   = panelsetup(cid, 1)
    
    lnp = 0
    
    for (i=1; i<=rows(z); i++) {
        Xrp  = panelsubmatrix(Xr, i, z)
         yp  = panelsubmatrix(y,  i, z)
         mus = rowsum(Xrp:*beta_rn)
         max = max(mus)
         sum = max + ln(colsum(exp(mus :- max)))
         lnp = lnp + colsum(yp:*mus) - sum
    }

    lnprior= -1/2*(beta_rn - b)*Winv*(beta_rn - b)' - 1/2*ldetW - cols(b)/2*ln(2*pi())

    return(lnp + lnprior)
}
end

mata:
st_view(y=., .,   "y")
st_view(X=., .,   "price contract local wknown tod seasonal")
st_view(pid=., ., "pid")                                            // Individual identifier
st_view(gid=., ., "gid")                                            // choice occasions
end

mata:
m = panelsetup(pid, 1)

b = J(1, 6, 0)
W = I(6)*6
beta = b :+ sqrt(diagonal(W))':*rnormal(rows(m), cols(b), 0, 1)
end

mata:

its    = 20000
burn   = 10000
nb     = cols(beta)
bvals  = J(0, cols(beta), .)
Wvals  = J(0, cols(rowshape(W, 1)), .)

propVs = J(rows(m), 1, rowshape(W, 1))
propms = J(rows(m), 1, b)
accept = J(rows(m), 1, 0)

atarg  = .25
lam    = J(rows(m), 1, 2.38^2/nb)
damper = 1

for (i=1; i<=its; i++) {
    
    b = drawb_betaW(beta, W/rows(m))
    W = drawW_bbeta(beta, b)
    
    bvals = bvals \ b
    Wvals = Wvals \ rowshape(W, 1)
    
    beta_old = beta
    
    Winv  = invsym(W)
    ldetW = ln(det(W))
    
    for (j=1; j<=rows(m); j++) {
        
        yi   = panelsubmatrix(y, j,m)
        Xi   = panelsubmatrix(X, j, m)
        gidi = panelsubmatrix(gid, j, m)
        
        propV = rowshape(propVs[j, ], nb)
    
        beta_old = beta[j, ]
        beta_hat = beta[j, ] + lam[j]*rnormal(1,nb,0,1)*cholesky(propV)' 
		
        old = lncp(beta_old, b, Winv, ldetW, yi, Xi, gidi)
        pro = lncp(beta_hat, b, Winv, ldetW, yi, Xi, gidi)
        
        if  (pro == . )     alpha = 0
        else if (pro > old) alpha = 1
        else                alpha = exp(pro - old)
        
        if (runiform(1, 1) < alpha) {
		    beta[j, ] = beta_hat
			accept[j] = accept[j] + 1
		}

        lam[j] = lam[j]*exp(1/(i+1)^damper*(alpha - atarg))
        propms[j, ] = propms[j, ] + 1/(i + 1)^damper*(beta[j, ] - propms[j, ])
        propV       = propV + 1/(i + 1)^damper*((beta[j, ] - propms[j, ])'(beta[j, ] - propms[j, ]) - propV)
		_makesymmetric(propV)
        propVs[j, ] = rowshape(propV, 1)    
    }    
}

arates = accept/its       
end

preserve 
clear
getmata (b*) = bvals
sum b*

gen t = _n
tsset t
forvalues i=1/6 {
    quietly tsline b`i', saving(bg`i'.gph, replace)
	local glist `glist' "bg`i'.gph"
}
graph combine `glist'

clear
getmata arates lam

hist arates 
hist lam

     

