capture program drop scpc
program scpc, eclass sortpreserve
	syntax , ///
		[ ///
			cvs	///
			k(real 10) ///
			avc(real -1) ///
			latlong		///
			uncond   ///
		]
		//	syntax [if] [in], [
//		avc(real 1.0)
//		]
	
	matrix scpc_Bread=e(V_modelbased)
	if( scpc_Bread[1,1]==.) {
		disp as text "e(V_modelbased) missing; aborting"
		exit(999)
	}	
	if `avc'==-1 {
		local avc 0.03	// set default
	}
	if `avc'<0.001 | `avc'>0.99 {
		disp as text "average correlation must be between 0.001 and 0.99; aborting"
		exit(999)
	}
	
	// ENSURE COORDINATES ARE UNIQUE WITHIN CLUSTERS
	if( "`e(clustvar)'"!="") {
		unab tmp_svarlist: s_*
		qui foreach tmp_svar of varlist `tmp_svarlist'{
			tempvar unique_coords_tag
			tempvar unique_coords_num
			egen `unique_coords_tag' = tag(`e(clustvar)' `tmp_svar')
			egen `unique_coords_num' = total(`unique_coords_tag'), by(`e(clustvar)')
			sum `unique_coords_num'
			if (`r(min)'!=1 | `r(max)'!=1){
				noi disp as text "coordinates not constant within clusters; aborting"
				exit(999)
			}
		}
	}
	// -----------------------------------------
	
	if( "`e(clustvar)'"!="") {
		tempvar es
		gen `es'=0 if e(sample)
		sort `e(clustvar)' `es'
	}
	tempvar scpc_sel
	scpc_setscores, scpc_sel(`scpc_sel')
	mata setOmsWfin(`avc',"`scpc_sel'","`latlong'")
	
	qui sum `scpc_sel', detail
	tempname neff 
	matrix `neff'=r(sum)
	matrix b=e(b)	
	matrix scpc_Bread0=e(V_modelbased)
	local clist 
	local na : colnames b
	if(e(cmd)=="regress"){
		tokenize `na'
		forval j = 1/`= colsof(b)' {
			if "``j''" != "_cons"  {
			local clist `clist' ``j''
			}
		}
		local clist="qui reg scpc_rx `clist'"
	}	
	else if(e(cmd)=="ivregress"){
		local elist `e(insts)'
		local i=`: word count `e(instd)''+1
		
		tempvar ess
		qui gen `ess'=0
		qui replace `ess'=1 if e(sample)
		_estimates hold oreg, restore
		cap drop scpc_fs_*
		local clist 
		tokenize `na'
		forval j = 1/`= colsof(b)' {
			if(`j'<`i'){			
				qui reg ``j'' `elist' if `ess'
				qui predict scpc_fs_`j', xb
				qui replace scpc_fs_`j'=. if !`ess'
				local clist `clist' scpc_fs_`j'
			}
			else{
				local clist `clist' ``j''
			}
		}
		matrix colnames scpc_Bread0 = `clist'
		_estimates unhold oreg		
		local clist="qui ivregress 2sls scpc_rx `e(exogr)' (`e(instd)' = `e(insts)' )"	
	}
	else if("`uncond'"==""){
		disp as text "conditional critical values only supported for ivregress 2sls and regress, use option uncond for only unconditionally valid critical values; aborting"
		exit(999)
	}
	
	local k=min(colsof(b),`k')
	matrix b=b[1..., 1..`k']
	local na : colnames b
	matrix scpctab =b'*(1,1,1,1,1,1)
	matrix scpc_cvs =b'*(1,1,1,1)

	tokenize `na'
	
	// Loop through covariates
	local i = 1	
	while "``i''" != "" {
		tempvar w
		matrix coef = `neff'*scpc_Bread[`i', 1...]
		qui matrix score `w' = coef if `scpc_sel' == 1		
		qui replace `w' = `w' + b[1,`i'] if `scpc_sel' == 1
		if("`uncond'"==""){	
			capture drop scpc_x 
			capture drop scpc_xs			
			matrix coef = `neff'*scpc_Bread0[`i', 1...]	
			qui matrix score scpc_x = coef if e(sample)						
			if( "`e(clustvar)'"=="") {
				qui gen scpc_xs=sign(scpc_x) if e(sample)				
				mata set_Wx_nocluster("`clist'", "`scpc_sel'")
			}
		else{
			qui by `e(clustvar)': egen scpc_xs=total(scpc_x*scpc_x) if e(sample)
			qui replace scpc_xs=scpc_x/sqrt(scpc_xs) if e(sample) & scpc_xs>0			
			mata set_Wx_cluster("`clist'", "`scpc_sel'")
			}
		}		
		mata set_scpcstats("`w'","`scpc_sel'")
		matrix scpctab[`i', 1] = scpcstats[1, 1..6]
		if("`cvs'"=="cvs"){
			mata set_scpccvs()
			matrix scpc_cvs[`i', 1] = scpc_cvs_c[1, 1..4]
		}
		local ++i
	}
	matrix colnames scpctab = "Coef" "Std_Err" "  t  " "P>|t|" "95% Conf" "Interval"
	local n = rowsof(scpctab)
	local ands = `n'*"&"
	local rs &-`ands'
	local pavc : di %5.3f `avc'
	local tit "SCPC Inference for first `k' coefficients"
	matlist scpctab, border(all) title(`tit') cspec(o2& %12s | %9.0g o2 & %9.0g o2 &o1 %5.2f o1& o2 %6.3f o1 & o2 %9.0g & o1 %9.0g o2&) rspec(`rs')
	if("`cvs'"=="cvs"){
//		mata set_scpccvs()
//		matrix rownames scpc_cvs ="Two-Sided" "One-Sided"
		matrix colnames scpc_cvs ="32%" "10%" "5%" "1%" 
		local n = rowsof(scpc_cvs)
		local ands = `n'*"&"
		local rs &-`ands'
		matlist scpc_cvs, border(all) title("Two-Sided Critical values of SCPC t-tests")  cspec(o2& %12s | %6.3f o2 & %6.3f o2 & %6.3f o2 & %6.3f o2 &) rspec(`rs')
		ereturn matrix scpccvs = scpc_cvs
	}
 	// Return results
	ereturn matrix scpcstats = scpctab
	cap drop scpc_*
end
