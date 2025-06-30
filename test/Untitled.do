clear all
sysuse auto, clear
gen s_1=rnormal(0,1)
gen s_2=rnormal(0,1)
regress mpg weight length, robust
di `e(clustvar)'
scpc
scpc2
// scpc ,avc(0.03)      
// scpc ,k(4)
// scpc ,avc(0.01) cvs k(1)
gen clust=round(rnormal(0,10),1)
regress mpg weight length, cluster(clust)

mata: "`e(clustvar)'"
// set sortrngstate 123
// //scpc
// set sortrngstate 123
// //scpc2

bys clust: replace s_1 = s_1[1]
bys clust: replace s_2 = s_2[1]

//scpc
scpc2
