void setOmsWfin(real scalar avc0, string scalar sel, string scalar latlong)
{
	external real matrix s
	external real scalar c0,cmax
	real scalar cv,n
	external real matrix distmat
	real matrix W
	external real scalar random_t
	external real scalar latlongflag
	struct mats vector Oms
	external real scalar qmax, capM, m, cgridfac, capN, minavc

	external struct mats vector Omsfin
	external real matrix Wfin
	external real vector permfin
	external real scalar cvfin
	external real scalar condflag
	
	// FIX BUG WITH COORDINATE VARIABLE ORDERING --------------------------
	
	stata("unab slist: s_*")  // create local macro slist containing a list of all s_* variables

	slist = st_local("slist") // import to mata
	slist = tokens(slist)' // tokenize list into vector
	slist_nums = strtoreal(substr(slist,3,.)) // extract numbers for sorting
	slist_order = order(slist_nums,1) // extract numerical order (alphabetical order would be incorrect for d>9)
	slist = slist[slist_order] // re-order elements of slist
	slist_len = length(slist) // number of elements

	for (i=1; i<=slist_len; i++) {
		if (substr(slist[i],3,.) != strofreal(i)){ 	// check that i'th element of slist is equal to "s_i", prohibiting e.g. "s_2 s_3" or "s_1 s_3"
			stata("disp as text"+char(34)+"s_* variables not continuously numbered starting from 1 (s_1, s_2, etc.)"+char(34))
			exit(999)
		}
	}
	
	s=st_data(.,slist',sel)		// Import slist variables, now correctly ordered and continuously numbered 
								// expects locations in s_1, s_2, s_3... etc. Dimension d is equal to the number of s_ variables 
	
	// ----------------------------------------------------------------------
	
	n=rows(s)
	if(sum(rowmissing(s) :== 0) < n){
		stata("disp as text"+char(34)+"missing value(s) in s_* variable; aborting"+char(34))
		exit(999)
	}
	if(n<5){
		stata("disp as text"+char(34)+"too few locations found in variables s_1, s_2 etc; aborting"+char(34))
		exit(999)
	}		
	stata("disp as text"+char(34)+"found "+strofreal(rows(s), "%6.0f")+" observations / clusters and "+strofreal(cols(s), "%3.0f")+"-dimensional locations in s_*"+char(34))
	latlongflag=latlong!=""
	if(latlongflag==0) {
		stata("disp as text"+char(34)+"using Euclidan norm to compute distance between locations stored in s_*")
	}
	else
	{
		if(cols(s)!=2) {
			stata("disp as text"+char(34)+"with latlong option, there must only be s_1 and s_2 present")
			exit(999)
		}
		else{
			stata("disp as text"+char(34)+"Computing distances on surface of sphere treating s_1 as latitude and s_2 as longitude")
		}
	}
	condflag=0
	setGQxw()
	minavc=0.00001			// minimal avc value for which size control is checked
	cgridfac=1.2			// factor in c-grid for size control
	
	if(avc0>=0.05){
		qmax=10
	}
	else{
		if(avc0>=0.01){
		qmax=20
		}
		else{
			if(avc0>=0.005){
				qmax=60
			}
			else{
				qmax=120
			}
		}
	}
	while(1){
		qmax=min((n-1,qmax))
		if(n<4500){
	// code for small n
			permfin=(1::n)
			distmat=getdistmat(s)			
			c0=getc0fromavc(lvech(distmat),avc0)		
			cmax=getc0fromavc(lvech(distmat),minavc)	
			W=getW(distmat,c0)
			Oms=getOms(distmat,c0,cmax,W)		
			}
		else{
			capN=20			// number of random sets of locations (of size m) used to approximate eigenvectors
			capM=1000000		// number of location pairs used to approximate Om(c) (in production: capM=1000000 or so)
			m=1000			// size of set of locations for which eigenvector is computed, and then Nystrom extended from m to n (in production: m=1000 or so)
							// note: m cannot be larger than n

	//		stata("display c(current_time)")
			random_t=1
			normalize_s(s)	
			LNsetWc0(s,avc0,W,c0,cmax)
	//		"done with W and c0"
	//		stata("display c(current_time)")

			Oms=LNgetOms(s,c0,cmax,W)	
	//		"done with Oms"
	//		stata("display c(current_time)")
		}
		setfinalW(Oms,W,cv)
		if(cols(W)-1<qmax | qmax==n-1) break
		qmax=round(qmax+qmax/2,1)
	}
	stata("disp as text"+char(34)+"SCPC optimal q = "+strofreal(cols(W)-1, "%3.0f")+" for maximal average pairwise correlation = "+strofreal(avc0, "%5.3f")+`pavc'"")
	
	Wfin=W
	Omsfin=Oms
	cvfin=cv

}
