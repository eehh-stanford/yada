# Test functions in theta_y_transforms.R

# Make an ordinal model with three variables for which all the parameters are
# different values
rho   <- 0.5
tau   <- c(-2.0,2.5)
s     <- 0.9
kap   <- 0.02

hp_homo <- list(paramModel='powLawOrdHomo')
hp_homo$J <- 1
hp_homo$M <- 2

th_y_vect_homo <- c(rho,tau,s)

hp_hetero <- list(paramModel='powLawOrdHetero')
hp_hetero$J <- 1
hp_hetero$M <- 2

th_y_vect_hetero <- c(rho,tau,s,kap)

hp_bad <- list(paramModel='not_a_model')
hp_bad$J <- 1
hp_bad$M <- 2

# Test check_model
expect_error(
  check_model(hp_homo$paramModel),
  NA
)

expect_error(
  check_model(hp_hetero$paramModel),
  NA
)

expect_error(
  check_model(hp_bad$paramModel),
  'Unrecognized model not_a_model'
)

# Test is_hetero
expect_equal(
  is_hetero(hp_homo$paramModel),
  F
)

expect_equal(
  is_hetero(hp_hetero$paramModel),
  T
)

# Test is_corr
expect_equal(
  is_corr('powLawHomo'),
  F
)

expect_equal(
  is_corr('powLawOrdHomo'),
  F
)

expect_equal(
  is_corr('powLawOrdHetero'),
  F
)

expect_equal(
  is_corr('powLawMixUncorrHomo'),
  F
)

expect_equal(
  is_corr('powLawMixUncorrHetero'),
  F
)

expect_equal(
  is_corr('powLawMixCorrHomo'),
  T
)

expect_equal(
  is_corr('powLawMixCorrHetero'),
  T
)

expect_error(
  is_corr('not_a_model'),
  'Unrecognized model not_a_model'
)

# Test theta_y_constr2unconstr
expect_error(
  th_y_vect_bar_homo <- theta_y_constr2unconstr(th_y_vect_homo,hp_homo),
  NA
)

expect_error(
  th_y_vect_bar_hetero <- theta_y_constr2unconstr(th_y_vect_hetero,hp_hetero),
  NA
)

# Do the computation directly and make sure the results are the same
th_y_vect_bar_homo2 <- c(log(rho),tau[1],log(tau[2]-tau[1]),log(s))

expect_equal(
  th_y_vect_bar_homo,
  th_y_vect_bar_homo2
)

th_y_vect_bar_hetero2 <- c(log(rho),tau[1],log(tau[2]-tau[1]),log(s),log(kap))

expect_equal(
  th_y_vect_bar_hetero,
  th_y_vect_bar_hetero2
)

# Test theta_y_2unconstr2
expect_error(
  th_y_vect_homo2 <- theta_y_unconstr2constr(th_y_vect_bar_homo,hp_homo),
  NA
)

expect_error(
  th_y_vect_hetero2 <- theta_y_unconstr2constr(th_y_vect_bar_hetero,hp_hetero),
  NA
)

expect_equal(
  th_y_vect_homo,
  th_y_vect_homo2
)

expect_equal(
  th_y_vect_hetero,
  th_y_vect_hetero2
)

# Test theta_y_vect2list
expect_error(
  th_y_list_homo <- theta_y_vect2list(th_y_vect_homo,hp_homo),
  NA
)

expect_error(
  th_y_list_hetero <- theta_y_vect2list(th_y_vect_hetero,hp_hetero),
  NA
)

# Test theta_y_list2vect
expect_error(
  th_y_vect_homo3 <- theta_y_list2vect(th_y_list_homo),
  NA
)

expect_error(
  th_y_vect_hetero3 <- theta_y_list2vect(th_y_list_hetero),
  NA
)

expect_equal(
  th_y_vect_homo,
  th_y_vect_homo3
)

expect_equal(
  th_y_vect_hetero,
  th_y_vect_hetero3
)

# Test get_var_index
hp_uo <- list(paramModel='powLawMixUncorrHomo')
hp_uo$J <- 1
hp_uo$K <- 1
hp_uo$M <- 2

hp_ue <- list(paramModel='powLawMixUncorrHetero')
hp_ue$J <- 1
hp_ue$K <- 1
hp_ue$M <- 2

hp_co <- list(paramModel='powLawMixCorrHomo')
hp_co$J <- 1
hp_co$K <- 1
hp_co$M <- 2

hp_ce <- list(paramModel='powLawMixCorrHetero')
hp_ce$J <- 1
hp_ce$K <- 1
hp_ce$M <- 2

hpList <- list()
hpList[[1]] <- hp_uo
hpList[[2]] <- hp_ue
hpList[[3]] <- hp_co
hpList[[4]] <- hp_ce

# rho
for(hp in hpList) {
  expect_equal(
    get_var_index('rho',hp,j=1),
    1
  )
}

for(hp in hpList) {
  expect_error(
    get_var_index('rho',hp),
    'Neither j nor k is specified'
  )
}

for(hp in hpList) {
  expect_error(
    get_var_index('rho',hp,k=1),
    'Unsupported variable for k being specified'
  )
}

for(hp in hpList) {
  expect_error(
    get_var_index('rho',hp,j=2),
    'j is not between 1 and J'
  )
}

# tau
for(hp in hpList) {
  expect_equal(
    get_var_index('tau',hp,j=1),
    c(2,3)
  )
}

for(hp in hpList) {
  expect_error(
    get_var_index('tau',hp),
    'Neither j nor k is specified'
  )
}

for(hp in hpList) {
  expect_error(
    get_var_index('tau',hp,k=1),
    'Unsupported variable for k being specified'
  )
}

for(hp in hpList) {
  expect_error(
    get_var_index('tau',hp,j=2),
    'j is not between 1 and J'
  )
}

# a, r, b
contVar <- c('a','r','b')
for(vv in 1:length(contVar)) {
  for(hp in hpList) {
    expect_equal(
      get_var_index(contVar[vv],hp,k=1),
      4 + vv-1
    )
  }

  for(hp in hpList) {
    expect_error(
      get_var_index(contVar[vv],hp),
      'Neither j nor k is specified'
    )
  }

  for(hp in hpList) {
    expect_error(
      get_var_index(contVar[vv],hp,j=1),
      'Unsupported variable for j being specified'
    )
  }

  for(hp in hpList) {
    expect_error(
      get_var_index(contVar[vv],hp,k=2),
      'k is not between 1 and K'
    )
  }
}

# test s
for(hp in hpList) {
  expect_equal(
    get_var_index('s',hp,j=1),
    7
  )
}

for(hp in hpList) {
  expect_equal(
    get_var_index('s',hp,k=1),
    8
  )
}

for(hp in hpList) {
  expect_error(
    get_var_index('s',hp),
    'Neither j nor k is specified'
  )
}

for(hp in hpList) {
  expect_error(
    get_var_index('s',hp,j=2),
    'j is not between 1 and J'
  )
}

for(hp in hpList) {
  expect_error(
    get_var_index('s',hp,k=2),
    'k is not between 1 and K'
  )
}

# test kappa
expect_error(
  get_var_index('kappa',hp_uo),
  'kappa requested but model is not heteroskedastic'
)

expect_equal(
  get_var_index('kappa',hp_ue),
  9
)

expect_error(
  get_var_index('kappa',hp_co),
  'kappa requested but model is not heteroskedastic'
)

expect_equal(
  get_var_index('kappa',hp_ce),
  10
)

# test z
expect_error(
  get_var_index('z',hp_uo),
  'z requested but model is not correlated'
)

expect_error(
  get_var_index('z',hp_ue),
  'z requested but model is not correlated'
)

expect_equal(
  get_var_index('z',hp_co),
  c(7,8,9)
)

expect_equal(
  get_var_index('z',hp_ce),
  c(7,8,9)
)
