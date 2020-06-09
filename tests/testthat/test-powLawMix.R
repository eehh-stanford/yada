# Test the functions in powLawMix

# From random.org between 1 and 1,000,000
set.seed(937974)

# Simulate data with two ordinal and two continuous variables with no
# correlations

modSpec_homo <- list(meanSpec='powLaw')
modSpec_homo$J <- 2
modSpec_homo$K <- 2
modSpec_homo$M <- c(2,2)
modSpec_homo$hetSpec = 'none'
modSpec_homo$cdepSpec = 'indep'

modSpec_hetero <- modSpec_homo
modSpec_hetero$hetSpec = 'sd_x'
modSpec_hetero$hetGroups <- rep(1,4)

th_y_sim_homo <- list(modSpec=modSpec_homo)
th_y_sim_homo$rho <- c(.75,.5)
th_y_sim_homo$tau <- list()
th_y_sim_homo$tau[[1]] <- c(2,5)
th_y_sim_homo$tau[[2]] <- c(1,3)
th_y_sim_homo$a <- c(2,3)
th_y_sim_homo$r <- c(.45,.10)
th_y_sim_homo$b <- c(1.2,-.5)
th_y_sim_homo$s <- c(.01,.02,.05,.04)

th_y_sim_hetero <- th_y_sim_homo
th_y_sim_hetero$modSpec <- modSpec_hetero
th_y_sim_hetero$kappa <- rep(0.02,4)

xmin <- 0
xmax <- 20
th_x_sim <- list(fitType='uniform',xmin=xmin,xmax=xmax)

#sim_homo <- simMixedCumProbit(th_y_sim_homo,100,th_x_sim,hp_homo)
sim_homo <- simPowLawMix(th_y_sim_homo,th_x_sim,100,modSpec_homo)
sim_hetero <- simPowLawMix(th_y_sim_hetero,th_x_sim,100,modSpec_hetero)


th_y_sim_vect_homo <- theta_y_list2vect(th_y_sim_homo)
th_y_sim_vect_hetero <- theta_y_list2vect(th_y_sim_hetero)

th_y_sim_vect_homo_bar <- theta_y_constr2unconstr(th_y_sim_vect_homo,modSpec_homo)
th_y_sim_vect_hetero_bar <- theta_y_constr2unconstr(th_y_sim_vect_hetero,modSpec_hetero)
