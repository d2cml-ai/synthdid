cd "F:\work\synthdid\vignettes\data"

// import delimited "F:\work\synthdid\vignettes\data\quota_example.csv", clear
mata: mata clear
*minimization
mata:
real vector lambda(matrix A, matrix b, matrix x, eta, zeta, maxIter, mindecrease) {    
    row = rows(A)
    col = cols(A)
    vals = J(1, maxIter, .)
    t=0
    dd=1
    
    while (t<maxIter & (t<2 | dd>mindecrease)) {
        t++
        Ax = A * x'	
        hg = (Ax - b)' * A + eta * x
        i = select((1..cols(hg)), colmin(hg :== min(hg)))[1,1]
        dx = -x
        dx[1,i] = 1 - x[1,i]
        v = abs(min(dx))+abs(max(dx))
        if (v==0) {
            x = x
            err = (A, b) * (x' \ -1)
            vals[1,t] = zeta^2 * (x * x') + (err' * err) / row
            if (t>1) {
                dd = vals[1,t-1] - vals[1,t]
            }
        }
        else {
            derr = A[.,i] - Ax
            step = -(hg) * dx' :/ ((derr' * derr) + eta * (dx * dx'))
            conststep = min((1, max((0, step)))) 
            x = x + conststep * dx  
            err = (A, b) * (x' \ -1)
            vals[1,t] = zeta^2 * (x * x') + (err' * err) / row
            if (t>1) {
                dd = vals[1,t-1] - vals[1,t]
            }
        }
    }
    return(x)
}
end

*fw.step 
mata:
real vector fw(matrix A, matrix b, matrix x, eta) {    
    Ax = A * x'	
    hg = (Ax - b)' * A + eta * x
    i = select((1..cols(hg)), colmin(hg :== min(hg)))[1,1]
    dx = -x
    dx[1,i] = 1 - x[1,i]
    v = abs(min(dx))+abs(max(dx))
    if (v==0) {
        x = x
    }
    else {
        derr = A[.,i] - Ax
        step = -(hg) * dx' :/ ((derr' * derr) + eta * (dx * dx'))
        conststep = min((1, max((0, step)))) 
        x = x + conststep * dx  
    }
    return(x)
}
end

*spar function
mata:
    real matrix sspar(matrix V) {
        W = J(1,length(V),.)
        for (i=1; i<=length(V); ++i) {
            W[1,i] = V[1,i] <= max(V)/4 ? 0 : V[1,i]
        }
        W = W :/ sum(W)
        return(W)
    }
end
	
*sum normalize
mata:
    real vector sum_norm(matrix O) {
        sum_o = sum(O)
        if (sum_o!=0) {
            O = O / sum_o
        }
        else {
            O = J(1, cols(O), 1/cols(O)) 
        }
        return(O)
    }
end

*simple merge bt two vectors
mata:
    real matrix smerge(matrix A, matrix B) {
        v = J(rows(A), 1, .)
        A = (A, v)
        for (i=1; i<=rows(A); i++) {
            for (j=1; j<=rows(B); j++) {
                if (A[i,1]==B[j,1]) A[i,2]=B[j,2]
            }
	    }
        A = A[.,2]
        return(A)
}
end	

*projected covariates
mata:
    void projected(Y, Yprojected, Beta) {
        K = cols(Y)
        cdat = Y[selectindex(Y[,6]:==0),(1,2,4,8..K)]
        cdat = select(cdat, rowmissing(cdat):==0)
        X = cdat[.,4::cols(cdat)]
        NX = cols(X)
        yearFEs = uniqrows(cdat[.,3])
        for (fe=1;fe<=rows(yearFEs);fe++) {
            fevar = cdat[.,3]:==yearFEs[fe]
            X = (X,fevar)
        }
        unitFEs = uniqrows(cdat[.,2])
        for (fe=1;fe<=rows(unitFEs);fe++) {
            fevar = cdat[.,2]:==unitFEs[fe]
            X = (X,fevar)
        }            
        y = cdat[.,1]
        XX = quadcross(X,1 , X,1)
        Xy = quadcross(X,1 , y,0)
        beta = invsym(XX)*Xy
        Beta = beta[1::NX]
        X = Y[.,8::K]
        Yprojected = Y[.,1]-X*Beta
}
end

use quota, clear

local 1 = "womparl"
local 2 = "country"
local 3 = "year"
local 4 = "quota"

tempvar touse
//
mark `touse' `if' `in'
**Check if group ID is numeric or string
local clustvar "`2'"

local stringvar=0
cap count if `2'==0
if _rc!=0 {
    local stringvar=1
    local groupvar `2'
    tempvar ID
    egen `ID' = group(`2')
    local varlist `1' `ID' `3' `4'
    tokenize `varlist'
}
else {
    tokenize `varlist'
}

qui xtset `2' `3'

qui count if `1'==.

qui count if `4'==.

qui count if `4'!=0&`4'!=1

qui sum `3'
qui sum `4' if `3'==r(min)

qui sum `4'

tempvar a1 a2
qui gen `a1'=`3' if `4'==1
qui bys `2': egen `a2'=min(`a1')
qui levelsof `a2', local(A2)



tempvar test
qui bys `2' (`3'): gen `test'=`4'-`4'[_n-1] 
qui sum `test'
qui count if `test'!=0&`test'!=1&`test'!=.

tempvar treated ty tyear n
qui gen `ty' = `3' if `4'==1 
qui bys `2': egen `treated' = max(`4') 
qui by  `2': egen `tyear'   = min(`ty') 

sort `3' `treated' `2'

gen `n'=_n
qui sum `3'
local mint=r(min)
qui putmata ori_id=`2' ori_pos=`n' if `3'==`mint' & `tyear'==., replace
if (length("`if'")+length("`in'")>0) {
    restore
}


local control_opt = 0

local m=1

mata: data0 = st_data(.,("`1' `2' `2' `3' `4' `treated' `tyear' `conts'"))
mata: data  = sort(data0, (6,2,4))
mata: units = panelsetup(data,2)
/// units
mata: NT = panelstats(units)[,3]
mata: treat=panelsum(data[.,(2,5)],units)
mata: treat[,1]=treat[,1]/NT
mata: treat[,2]=1*(treat[,2]:>=1) + 0*(treat[,2]:==1)
mata: Nco = sum(data[,7]:==.)/NT
mata: controlID = uniqrows(select(data[.,2],data[,7]:==.))  
mata: jk = 1
mata: inference = 0
mata: controls = 0
mata: mt = 1


mata: uniqID=(uniqrows(select(data[.,2],data[,7]:==.)) \ uniqrows(select(data[.,2],data[,7]:!=.)))
mata: N = panelstats(units)[,1]
mata: Beta = J(1, 1, .)
mata: data_ori=data
mata: trt = select(uniqrows(data[,7]),uniqrows(data[,7]):!=.)
mata: tau    = J(rows(trt),1,.)
mata: tau_wt = J(1,rows(trt),.)

mata: Tcond = J(rows(trt), 1, .)
mata: Ncond = J(rows(trt), 1, .)

mata: Omega = J(Nco, rows(trt),.)
mata: Lambda = J(NT, rows(trt),.)
mata: Beta = J(cols(data)-7, rows(trt),.)
mata: OMEGA  = st_matrix("omega")
mata: LAMBDA = st_matrix("lambda")
mata:
    for(yr=1;yr<=rows(trt);++yr) {
// 	for(yr=1;yr<=1;++yr) {
        cond1 = data[,7]:==trt[yr]
        cond2 = data[,7]:==.
        cond = cond1+cond2
        yNtr = sum(cond1)
        yNco = sum(cond2)
        ydata = select(data,cond)
        yunits = panelsetup(ydata,2)
        yNT = panelstats(yunits)[,3]
        yNG = panelstats(yunits)[,1]
        yN  = panelstats(yunits)[,2]
        yNtr = yNtr/yNT
        yNco = yNco/yNT
        yTpost = max(ydata[,4])-trt[yr]+1
        pretreat = select(uniqrows(data[,4]),uniqrows(data[,4]):<trt[yr])
        Npre  = rows(pretreat)
        Npost = yNT-Npre
        //Calculate Zeta
        ndiff = yNG*(yNT-1)
        ylag = ydata[1..(yN-1),1]
        ylag = (. \ ylag)
        diff = ydata[.,1]-ylag
        first = mod(0::(yN-1),yNT):==0
        postt = ydata[,4]:>=trt[yr]
        dropc = first+postt+ydata[,6]
        prediff = select(diff,dropc:==0)
        sig_t = sqrt(variance(prediff))
        EtaLambda = 1e-6
        
        if (mt!=3) {
            EtaOmega = (yNtr*yTpost)^(1/4)
        }

        
        yZetaOmega  = EtaOmega*sig_t
        yZetaLambda = EtaLambda*sig_t
        //Generate Y matrices
        ids = uniqrows(ydata[.,2])
        ytreat = ydata[,7]:==trt[yr]
        ytreat=panelsum(ytreat,yunits)
        ytreat=ytreat/yNT
        Y = rowshape(ydata[.,1],yNG)
        Y0 = select(Y,ytreat:==0)
        Y1 = select(Y,ytreat:==1)
        //Generate average of outcomes for each unit over time
        promt = mean(Y0[,(Npre+1)::yNT]')'
        //Calculate input matrices (pre-treatment and averages)
        A_l = Y0[,1..Npre]:-mean(Y0[,1..Npre])
        b_l = promt:-mean(promt)
        
        if (mt!=3) {
            A_o = Y0[,1..Npre]':-mean(Y0[,1..Npre]')
            b_o = mean(Y1[.,1..Npre])':-mean(mean(Y1[.,1..Npre])')
        }

        //Calculate Tau for t

        if (mt!=3) {
            lambda_l1 = J(1,cols(A_l),1/cols(A_l))
			lambda_l = lambda_l1
        }
        lambda_o = J(1,cols(A_o),1/cols(A_o))
        mindec = (1e-5*sig_t)^2

        if (controls==0 | controls==1) {
            if (mt!=2) {
                //Find optimal weight matrices
                eta_o = Npre*yZetaOmega^2
                eta_l = yNco*yZetaLambda^2
                if (mt==1) {
                    lambda_l = lambda(A_l,b_l,lambda_l,eta_l,yZetaLambda,100,mindec)
                    lambda_l = sspar(lambda_l)
                    lambda_l = lambda(A_l, b_l, lambda_l,eta_l,yZetaLambda, 10000,mindec)
                }
                lambda_o = lambda(A_o, b_o, lambda_o,eta_o,yZetaOmega, 100,mindec)
                lambda_o = sspar(lambda_o)
                lambda_o = lambda(A_o, b_o, lambda_o,eta_o,yZetaOmega, 10000,mindec)
            }
            if (inference==0) {
                Lambda[.,yr] =  (lambda_l' \ J(Npost,1,.))
                Omega[.,yr] = lambda_o'
            }
            tau[yr] = (-lambda_o, J(1,yNtr,1/yNtr))*Y*(-lambda_l, J(1,Npost,1/Npost))'
            tau_wt[yr] = yNtr*Npost
			Tcond[yr] = yNtr
			Ncond[yr] = Npost
        }
    }
		

        tau_wt = tau_wt/sum(tau_wt')
        ATT = tau_wt*tau

end
  

