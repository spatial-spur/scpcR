void set_Wx_cluster(string scalar Vname, string scalar sel)
{	
	external real scalar condflag
	external real matrix Wfin
	external real matrix Wx
	external real vector permfin
	real vector r
	real scalar i
	
	condflag=1
	r=J(rows(Wfin),1,0)
	Wx=Wfin
	for(i=1;i<=cols(Wfin);i++){
		r[permfin]=Wfin[.,i]	
		stata("cap drop scpc_r")
		stata("qui gen scpc_r=.")
		stata("qui replace scpc_r=0 if e(sample)")
		st_store(.,"scpc_r",sel,r)
		stata("cap drop scpc_rx")
		stata("qui by \`e(clustvar)': egen scpc_rx=total(scpc_r) if e(sample)") // BACKSLASH TO PREVENT EVALUATION WHEN FUNCTION IS FIRST READ
		
		stata("cap drop scpc_r")
		if(i==1){
			stata("qui gen scpc_r=scpc_rx*scpc_xs if e(sample)")
		}
		else{
			stata("qui replace scpc_rx=scpc_rx*scpc_xs")
			stata("_estimates hold oreg, restore")		
			stata(Vname)
			stata("qui predict scpc_r, resid")	
			stata("_estimates unhold oreg")
		}
		stata("cap drop scpc_rx")
		stata("qui by \`e(clustvar)': egen scpc_rx=total(scpc_r*scpc_x) if e(sample)") // BACKSLASH TO PREVENT EVALUATION WHEN FUNCTION IS FIRST READ
		
		Wx[.,i]=st_data(.,"scpc_rx",sel)[permfin]
	}
	set_ccspc()
}
