sysname <- Sys.info()["sysname"]
if (sysname == "Windows") {
  options("RStata.StataPath" = "\"C:\\Program Files\\Stata18\\StataSE-64\"")
  options("RStata.StataVersion" = 18)
} else if (sysname == "Darwin") {  # macOS
  stata_path_mp <- "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp"
  stata_path_se <- "/Applications/Stata/StataSE.app/Contents/MacOS/stata-se"
  if (file.exists(stata_path_mp)) {
    options("RStata.StataPath" = stata_path_mp)
  } else if (file.exists(stata_path_se)) {
    options("RStata.StataPath" = stata_path_se)
  } else {
    stop("Stata executable not found in expected macOS locations.")
  }
  options("RStata.StataVersion" = 18)
} else {
  stop("Unknown or unsupported OS for RStata setup")
}


rm(list = ls())
set.seed(123)
n  <- 500
x1 <- rnorm(n)
x2 <- rnorm(n)
y  <- 1 + .8*x1 - .3*x2 + rnorm(n)
coords <- cbind(runif(n,-180,180), runif(n,-90,90))   # fake locations
df <- data.frame(y,x1,x2, lon = coords[,1], lat = coords[,2])
df$cl <- sample(1:100, nrow(df), replace = TRUE)  # add a cluster variable
write.csv(df, "toy.csv", row.names = FALSE)      # exchange format

source("R/scpcR.R")  # the canvas file, or devtools::load_all()

# empty vector
diffs <- c()

testfunc_mata <- function(cmd){
  stata_src <- paste0("
  clear
  set maxvar 30000
  qui do \"./test/scpc.ado\"
  import delimited ./toy.csv, clear
  rename lat s_1
  rename lon s_2
  order s_1 s_2

  ",
  cmd
  )
  out = stata(stata_src, data.out = TRUE)
  return(out)
}

## getdistmat ----------------------------------
out_stata <- testfunc_mata("
              mata
              s = st_data(.,\"s_*\")
              D = getdistmat(s)
              end
              clear
              getmata c* = D
              ")

out_R <- .getdistmat(df[,c("lon","lat")], latlong = TRUE)

diffs <- c(diffs, max(out_R - out_stata))

## lvech ----------------------------------
out_stata <- testfunc_mata("
              mata
              s = st_data(.,\"s_*\")
              D = getdistmat(s)
              lvech = lvech(D)
              end
              clear
              getmata c* = lvech
              ")

out_R <- .lvech(.getdistmat(df[,c("lon","lat")], latlong = TRUE))

diffs <- c(diffs, max(out_R - out_stata))

## getavc ----------------------------------
out_stata <- testfunc_mata("
              mata
              s = st_data(.,\"s_*\")
              D = getdistmat(s)
              lvech = lvech(D)
              avc = getavc(0.03,lvech)
              end
              clear
              getmata c* = avc
              ")

out_R <- .getavc(0.03, .lvech(.getdistmat(df[,c("lon","lat")], latlong = TRUE)))

diffs <- c(diffs, max(abs(out_R - out_stata)))

## demeanmat ----------------------------------
out_stata <- testfunc_mata("
              mata
              s = st_data(.,\"s_*\")
              D = getdistmat(s)
              D2 = demeanmat(D)
              end
              clear
              getmata c* = D2
              ")

out_R <- .demeanmat(.getdistmat(df[,c("lon","lat")], latlong = TRUE))

diffs <- c(diffs, max(abs(out_R - out_stata)))

## getc0fromavc ----------------------------------
out_stata <- testfunc_mata("
              mata
              s = st_data(.,\"s_*\")
              D = getdistmat(s)
              lvech = lvech(D)
              c0 = getc0fromavc(lvech,0.03)
              end
              clear
              getmata c* = c0
              ")

out_R <- .getc0fromavc(.lvech(.getdistmat(df[,c("lon","lat")], latlong = TRUE)), 0.03)

diffs <- c(diffs, max(abs(out_R - out_stata)))

## getW ----------------------------------
out_stata <- testfunc_mata("
              mata
              s = st_data(.,\"s_*\")
              D = getdistmat(s)
              qmax = 13
              W = getW(D,0.03)
              end
              clear
              getmata c* = W
              ")

out_R <- .getW(.getdistmat(df[,c("lon","lat")], latlong = TRUE), 0.03, 13)

diffs <- c(diffs, max(abs(abs(out_R) - abs(out_stata))))

## getnc ----------------------------------
out_stata <- testfunc_mata("
              mata
              cgridfac = 1.2
              nc = getnc(1.3, 3.5)
              end
              clear
              getmata c* = nc
              ")

out_R <- .getnc(1.3, 3.5, 1.2)

diffs <- c(diffs, max(abs(out_R - out_stata)))

diffs

## getOms
out_stata <- testfunc_mata("
              mata
              s = st_data(.,\"s_*\")
              D = getdistmat(s)
              qmax = 13
              W = getW(D,0.03)
              cgridfac = 1.2
              Oms = getOms(D, 0.01, 5, W)
              Oms
              Om = Oms[2]
              Omm = Om.mat
              end
              clear
              getmata c* = Omm
              ")

out_R <- .getOms(.getdistmat(df[,c("lon","lat")], latlong = TRUE), 
                 0.01, 
                 5, 
                 .getW(.getdistmat(df[,c("lon","lat")], latlong = TRUE), 0.03, 13), 
                 1.2)[[2]]

diffs <- c(diffs, max(abs(abs(out_R) - abs(out_stata))))
diffs

## gettau
out_stata <- testfunc_mata("
              mata
              s = st_data(.,\"s_*\")
              D = getdistmat(s)
              qmax = 13
              W = getW(D,0.03)
              W
              y = 1::rows(W)
              tau = gettau(y,W)
              end
              clear
              getmata c* = tau
              ")

out_R <- .gettau(1:nrow(df), .getW(.getdistmat(df[,c("lon","lat")], latlong = TRUE), 0.03, 13))

diffs <- c(diffs, max(abs(out_R - out_stata)))

diffs

## getrp ----------------------------------
out_stata <- testfunc_mata("
              mata
              s = st_data(.,\"s_*\")
              D = getdistmat(s)
              qmax = 13
              W = getW(D,0.03)
              cgridfac = 1.2
              Oms = getOms(D, 0.01, 5, W)
              Om = Oms[2]
              Omm = Om.mat
              setGQxw()
              rp = getrp(Omm, 3)
              end
              clear
              getmata c* = rp
              ")

out_R <- .getrp(.getOms(.getdistmat(df[,c("lon","lat")], latlong = TRUE), 
                        0.01, 
                        5, 
                        .getW(.getdistmat(df[,c("lon","lat")], latlong = TRUE), 0.03, 13), 
                        1.2)[[2]], 
                 3)

diffs <- c(diffs, max(abs(out_R - out_stata)))

diffs

## maxrp ----------------------------------
out_stata <- testfunc_mata("
              mata
              s = st_data(.,\"s_*\")
              D = getdistmat(s)
              qmax = 13
              W = getW(D,0.03)
              cgridfac = 1.2
              Oms = getOms(D, 0.01, 5, W)
              setGQxw()
              maxrp = maxrp(Oms, 13, 3)
              end
              clear
              getmata c* = maxrp
              ")

out_R <- .maxrp(.getOms(.getdistmat(df[,c("lon","lat")], latlong = TRUE), 
                        0.01, 
                        5, 
                        .getW(.getdistmat(df[,c("lon","lat")], latlong = TRUE), 0.03, 13), 
                        1.2), 
                 13, 
                 3)$max

diffs <- c(diffs, max(abs(out_R - out_stata)))
diffs

## getcv ----------------------------------

out_stata <- testfunc_mata("
              mata
              s = st_data(.,\"s_*\")
              D = getdistmat(s)
              qmax = 13
              W = getW(D,0.03)
              cgridfac = 1.2
              Oms = getOms(D, 0.01, 5, W)
              setGQxw()
              cv = getcv(Oms, 11, 0.05)
              end
              clear
              getmata c* = cv
              ")

out_R <- .getcv(.getOms(.getdistmat(df[,c("lon","lat")], latlong = TRUE), 
                        0.01, 
                        5, 
                        .getW(.getdistmat(df[,c("lon","lat")], latlong = TRUE), 0.03, 13), 
                        1.2), 
                 11, 
                 0.05)

diffs <- c(diffs, max(abs(out_R - out_stata)))
diffs

## setfinalW ----------------------------------

out_stata <- testfunc_mata("
              mata
              s = st_data(.,\"s_*\")
              D = getdistmat(s)
              qmax = 13
              W = getW(D,0.03)
              cgridfac = 1.2
              Oms = getOms(D, 0.01, 5, W)
              setGQxw()
              cv = .
              setfinalW(Oms, W, cv)
              end
              clear
              getmata c* = W
              ")

out_R <- .setfinalW(.getOms(.getdistmat(df[,c("lon","lat")], latlong = TRUE), 
                            0.01, 
                            5, 
                            .getW(.getdistmat(df[,c("lon","lat")], latlong = TRUE), 0.03, 13), 
                            1.2), 
                     .getW(.getdistmat(df[,c("lon","lat")], latlong = TRUE), 0.03, 13), 
                     13)$W

diffs <- c(diffs, max(abs(abs(out_R) - abs(out_stata))))
diffs

## setOmsWfin ----------------------------------
out_stata <- testfunc_mata("
              mata
              setOmsWfin(0.03, \"\", \"latlong\")
              end
              clear
              getmata c* = Wfin
              ")

out_R <- .setOmsWfin(.getdistmat(df[,c("lon","lat")], latlong = TRUE), 
                     0.03)$Wfin

diffs <- c(diffs, max(abs(abs(out_R) - abs(out_stata))))
diffs




m   <- feols(y ~ x1 + x2, data = df)           # same formula

start.time <- Sys.time()

Rout <- scpc(m,
             df,
             c("lat","lon"),  # names of the coordinates
             cluster = NULL,
             avc     = 0.03,
             latlong = TRUE,    # must match the Stata run
             cvs     = TRUE,
             uncond = FALSE)

end.time <- Sys.time()
R_time <- end.time - start.time
R_time


bla <- testfunc_mata(
  "
  reg y x1 x2, vce(cluster cl)
  timer clear
  timer on 1
  scpc, latlong avc(0.03)
  timer off 1
  timer list 1
    reg y x1 x2, vce(robust)
  timer clear
  timer on 1
  scpc, latlong avc(0.03)
  timer off 1
  timer list 1
  "
)

Rout$scpcstats[c(2,3,1),]
