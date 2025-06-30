clear
import delimited "./toy.csv", clear
rename lat s_1
rename lon s_2
order s_1 s_2
reg y x1 x2, robust       // any supported regression
scpc, latlong avc(0.03) cvs    // match the R call options you’ll use
order s_2 s_1
scpc, latlong avc(0.03) cvs    // match the R call options you'll use


gen touse = 1

mata
s=get_s_matrix("touse","latlong")
latlongflag = 1

D = getdistmat(s)
end
clear
getmata c* = D

matrix list e(scpcstats)            // t-stats Stata stores here
matrix b = e(scpcstats)
clear
svmat double b            // put them in a variable
export delimited using "stata_scpc.csv", replace



clear all
sysuse auto, clear
gen s_1=runiform(-90,90) // "lat"
gen s_2=runiform(-180,180) // "lon"
reg mpg weight length, robust
scpc, latlong
order s_2 s_1 
scpc, latlong



    . sysuse auto, clear
    . gen s_1=rnormal(0,1)
    . gen s_2=rnormal(0,1)
    . regress mpg weight length, robust
    . scpc
    . gen clust=round(rnormal(0,10),1)
    . regress mpg weight length, cluster(clust)
    . scpc
