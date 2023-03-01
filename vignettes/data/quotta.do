



cd "F:\work\synthdid\vignettes\data"


// import delimited "F:\work\synthdid\vignettes\data\quota_example.csv", clear

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
// mata: data1 = st_data(., ("`4' `treated'" "`1'"))
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
mata:
//         if (jk==1) {
            uniqID=(uniqrows(select(data[.,2],data[,7]:==.)) \ uniqrows(select(data[.,2],data[,7]:!=.)))
            N = panelstats(units)[,1]
//         }
// end
		//		
//         //matrix for beta
//         if (inference==0 & (controls==0)) {
            Beta = J(1, 1, .)
//         }
//		
//         //Adjust for controls in xysnth way
//         //save original data for jackknife
        data_ori=data
//		
		trt = select(uniqrows(data[,7]),uniqrows(data[,7]):!=.)
//		
//         //Iterate over years, calculating each estimate
        tau    = J(rows(trt),1,.)
        tau_wt = J(1,rows(trt),.)
// 		att = tau_wt * tau
//         if (inference==0) {
            Omega = J(Nco, rows(trt),.)
            Lambda = J(NT, rows(trt),.)
//         }
			for(yr=1;yr<=1;++yr) {
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
//             ids = uniqrows(ydata[.,2])
            ytreat = ydata[,7]:==trt[yr]
            ytreat=panelsum(ytreat,yunits)
            ytreat=ytreat/yNT
            Y = rowshape(ydata[.,1],yNG)
            Y0 = select(Y,ytreat:==0)
            Y1 = select(Y,ytreat:==1)
//             //Generate average of outcomes for each unit over time
            promt = mean(Y0[,(Npre+1)::yNT]')'
//             //Calculate input matrices (pre-treatment and averages)
            A_l = Y0[,1..Npre]:-mean(Y0[,1..Npre])
            b_l = promt:-mean(promt)
//			
//             if (mt!=3) {
                A_o = Y0[,1..Npre]':-mean(Y0[,1..Npre]')
                b_o = mean(Y1[.,1..Npre])':-mean(mean(Y1[.,1..Npre])')
//             }
//
//             //Calculate Tau for t
//		
//                if (mt!=3) {
                    lambda_l = J(1,cols(A_l),1/cols(A_l))
//                 }
                lambda_o = J(1,cols(A_o),1/cols(A_o))
//                
            mindec = (1e-5*sig_t)^2
			}
end
//
//
// mata: mata describe
