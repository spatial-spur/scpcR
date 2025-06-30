## R – create tiny toy data (or use any real one)
rm(list = ls())
set.seed(123)
n  <- 500
x1 <- rnorm(n)
x2 <- rnorm(n)
y  <- 1 + .8*x1 - .3*x2 + rnorm(n)
coords <- cbind(runif(n,-180,180), runif(n,-90,90))   # fake locations
df <- data.frame(y,x1,x2, lon = coords[,1], lat = coords[,2])
write.csv(df, "toy.csv", row.names = FALSE)      # exchange format

library(fixest)   # or lm()
source("R/scpcR.R")  # the canvas file, or devtools::load_all()

df  <- read.csv("toy.csv")

df[1,1] <- NA
m   <- feols(y ~ x1 + x2, data = df, data.save = TRUE)           # same formula
m2   <- lm(y ~ x1 + x2, data = df)           # same formula

Rout <- scpc(m,
             df,
             c("lat","lon"),  # names of the coordinates
             avc     = 0.03,
             latlong = TRUE,    # must match the Stata run
             cvs     = TRUE,
             uncond = FALSE,
             k = NULL)

print(Rout)

Rout$scpcstats    # vector of R t-stats
Rout$scpccvs
