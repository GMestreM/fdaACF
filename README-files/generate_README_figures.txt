#### Script for generating the figures shown in README.md ####

rm(list = ls())
library(fdaACF)

path_img <- "README-files/figure-html/"

# autocov_surfc
name_fig <- "autocov_surfc.png"


# b_bridge
name_fig <- "b_bridge.png"

N <- 400
v <- seq(from = 0, to = 1, length.out = 24)
sig <- 1
set.seed(10) # For replication
Y <- simulate_iid_brownian_bridge(N, v, sig)

png(filename = file.path(getwd(),path_img,name_fig),
    width = 540,
    height = 540,
    units = "px")
matplot(t(Y), type = "l",lty = 1,xlab = "v", ylab = "", main = "Functional Brownian Bridge")
dev.off()

# elec_prices
name_fig <- "elec_prices.png"

png(filename = file.path(getwd(),path_img,name_fig),
    width = 540,
    height = 540,
    units = "px")
data(elec_prices)
matplot(t(elec_prices), type = "l",lty = 1,xlab = "Hours", ylab = "Price (€/MWh)")
dev.off()


# FACF_FPACF_Prices
name_fig <- "FACF_FPACF_Prices.png"

v <- seq(from = 1, to = 24)
nlags <- 30
ci <- 0.95


# Autocorrelation function for functional data
FACF <- obtain_FACF(Y = as.matrix(elec_prices), 
                    v = v,
                    nlags = nlags,
                    ci = ci,
                    figure = TRUE)

# Partial autocorrelation function for functional data
FPACF <- obtain_FPACF(Y = as.matrix(elec_prices), 
                      v = v,
                      nlags = nlags,
                      n_harm = 5, 
                      ci = ci,
                      figure = FALSE)

png(filename = file.path(getwd(),path_img,name_fig),
    width = 789,
    height = 405,
    units = "px")
par(mfrow = c(1,2))
plot_FACF(rho = FACF$rho,
          Blueline = FACF$Blueline,
          ci = ci,
          main = "FACF")
plot_FACF(rho = FPACF$rho,
          Blueline = FPACF$Blueline,
          ci = ci,
          main = "FPACF")
par(mfrow = c(1,1))
dev.off()


# FACF_iid
name_fig <- "FACF_iid.png"


N <- 400
v <- seq(from = 0, to = 1, length.out = 24)
sig <- 1
set.seed(10) # For replication
Y <- simulate_iid_brownian_bridge(N, v, sig)

png(filename = file.path(getwd(),path_img,name_fig),
    width = 540,
    height = 540,
    units = "px")
FACF_iid <- obtain_FACF(Y = Y, 
                        v = v,
                        nlags = 30,
                        ci = ci,
                        figure = TRUE)    
dev.off()

# FACF_prices
name_fig <- "FACF_prices.png"