



cd "F:\work\synthdid\vignettes\data"


import delimited "F:\work\synthdid\vignettes\data\quota_example.csv", clear


local 1 = "womparl"
local 2 = "country"
local 3 = "year"
local 4 = "quota"

tempvar touse
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
// if `"`r(balanced)'"'!="strongly balanced" {
//     dis as error "Panel is unbalanced."
//     exit 451
// }
qui count if `1'==.
if r(N)!=0 {
    dis as error "Missing values found in dependent variable.  A balanced panel without missing observations is required."
    exit 416
}
qui count if `4'==.
if r(N)!=0 {
    dis as error "Missing values found in treatment variable.  A balanced panel without missing observations is required."
    exit 416
}
qui count if `4'!=0&`4'!=1
if r(N)!=0 {
    dis as error "Treatment variable takes values distinct from 0 and 1."
    exit 450
}
qui sum `3'
qui sum `4' if `3'==r(min)
if r(max)!=0 {
    dis as error "Certain units are treated in the first period of the panel."
    dis as error "Any units which are treated the entire panel should not be included in the synthetic DID procedure."
    exit 459
}
qui sum `4'
if (r(min)==0 & r(max)==0)==1 {
    di as error "All units are controls."
    exit 459
}

tempvar a1 a2
qui gen `a1'=`3' if `4'==1
qui bys `2': egen `a2'=min(`a1')
qui levelsof `a2', local(A2)
foreach t of local A2 {
    qui sum `4' if `3'>=`t'
    if r(min)==r(max) {
        di as error "All units are treated in some adoption time."
        exit 459
    }
}


tempvar test
qui bys `2' (`3'): gen `test'=`4'-`4'[_n-1] 
qui sum `test'
qui count if `test'!=0&`test'!=1&`test'!=.
if r(N)!=0 {
    local e1 "to only change from untreated to treated, or remain untreated."
    dis as error "Units are observed to change from treated (earlier) to untreated (later)."
    dis as error "A staggered adoption is assumed in which units are assumed `e1'"
    exit 459
}
drop `test'

local SEN "standard error needs"
if "`vce'"=="jackknife" {
    tempvar t1 t2
    qui gen `t1'=`3' if `4'==1
    qui bys `2': egen `t2'=min(`t1')
    qui levelsof `t2', local(T2)
    foreach t of local T2 {
        qui sum `2' if `4'==1 & `t2'==`t'
        if r(min)==r(max) {
            di as error "Jackknife `SEN' at least two treated units for each treatment period"
            exit 451
        }
    }
}

if (length("`if'")+length("`in'")>0) {
    restore
}

tempvar treated ty tyear n
qui gen `ty' = `3' if `4'==1 
qui bys `2': egen `treated' = max(`4') if `touse'
qui by  `2': egen `tyear'   = min(`ty') if `touse'

if (length("`if'")+length("`in'")>0) {
    preserve
    qui keep if `touse'
}

sort `3' `treated' `2'
if "`mattitles'"!="" {
    if `stringvar'==1 qui levelsof `groupvar' if `treated'==0
    else              qui levelsof `2' if `treated'==0
    local rnames = r(levels)
}

gen `n'=_n
qui sum `3'
local mint=r(min)
qui putmata ori_id=`2' ori_pos=`n' if `3'==`mint' & `tyear'==., replace
if (length("`if'")+length("`in'")>0) {
    restore
}

if length("`graph'")!=0&`stringvar'==1 {
    preserve
    qui sum `3' if `touse'
    **Save original state names for later use with graph
    qui keep if `3' == r(min) & `touse'
    keep `groupvar' `2'
    tempfile stateString
    rename `groupvar' stateName
    rename `2' state
    qui save `stateString'
    restore
}
local control_opt = 0
if "`covariates'"!="" {
    _parse_X `covariates'
    if "`r(ctype)'"=="optimized" local control_opt = 2
    if "`r(ctype)'"=="projected" local control_opt = 1
    local conts = r(controls)
    local contname = r(controls)
	
    foreach var of varlist `contname' {
        qui sum `var'
        if r(sd)==0 {
            dis as error "Covariates were found to be constant in the estimation sample."
            dis as error "Please remove constant covariates from the covariate set."
            exit 416
        }
    }
	
    if `control_opt'==2&length(`"`unstandardized'"')==0 {
        local sconts
        foreach var of varlist `conts' {
            tempvar z`var'
            qui egen `z`var''= std(`var') if `touse'
            local sconts `sconts' `z`var''
        }
        local conts `sconts'
    }
	
    tempvar nmiss
    qui egen `nmiss' = rowmiss(`conts')
    qui sum `nmiss' if `touse'
    if r(mean)!=0 {
        dis as error "Missing values found in covariates."
        dis as error "A balanced panel without missing observations is required."
        exit 416
    }
}
local m=1
display "`m'"

mata: data0 = st_data(.,("`1' `2' `2' `3' `4' `treated' `tyear' `conts'"))
mata: data  = sort(data0, (6,2,4))
mata: units = panelsetup(data,2)
mata: NT = panelstats(units)[,3]
mata: treat=panelsum(data[.,(2,5)],units)
mata: treat[,1]=treat[,1]/NT
mata: treat[,2]=1*(treat[,2]:>=1) + 0*(treat[,2]:==1)
mata: Nco = sum(data[,7]:==.)/NT
mata: controlID = uniqrows(select(data[.,2],data[,7]:==.))  
mata: jk = 1
mata: inference = 0
mata: controls = 0
mata:
        if (jk==1) {
            uniqID=(uniqrows(select(data[.,2],data[,7]:==.)) \ uniqrows(select(data[.,2],data[,7]:!=.)))
            N = panelstats(units)[,1]
        }
		
        //matrix for beta
        if (inference==0 & (controls==0)) {
            Beta = J(1, 1, .)
        }
		
        //Adjust for controls in xysnth way
        //save original data for jackknife
        data_ori=data
		
		trt = select(uniqrows(data[,7]),uniqrows(data[,7]):!=.)
		
        //Iterate over years, calculating each estimate
        tau    = J(rows(trt),1,.)
        tau_wt = J(1,rows(trt),.)
		
        if (inference==0) {
            Omega = J(Nco, rows(trt),.)
            Lambda = J(NT, rows(trt),.)
        }
			for(yr=1;yr<=rows(trt);++yr) {
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
			}
end


mata: mata describe
