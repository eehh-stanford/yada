

# Test functions in theta_y_transforms.R

# Test functions that depend only on modSpec and not on the parameter vector:
# check_model
# get_J
# get_K
# get_Gkappa
# get_Gz
# is_hetero
# is_cdep
# get_z_length
# get_num_var

# A one variable, ordinal model that is homoskedastic
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 1
modSpec$M <- 2
modSpec$hetSpec <- 'none' # homoskedastic

expect_error(
  check_model(modSpec),
  NA
)

#if(F) {
expect_equal(
  get_J(modSpec),
  1
)

expect_equal(
  get_K(modSpec),
  0
)

expect_equal(
  get_Gkappa(modSpec),
  0
)

expect_equal(
  get_Gz(modSpec),
  0
)

expect_equal(
  is_hetero(modSpec),
  F
)

expect_equal(
  is_cdep(modSpec),
  F
)

expect_equal(
  get_z_length(modSpec),
  0
)

#expect_error(
#  get_kappa_full(modSpec),
#  'get_kappa_full called, but model is not heteroskedastic'
#)

expect_equal(
  get_num_var('rho',modSpec),
  1
)

expect_equal(
  get_num_var('tau',modSpec),
  2
)

expect_equal(
  get_num_var('a',modSpec),
  0
)

expect_equal(
  get_num_var('r',modSpec),
  0
)

expect_equal(
  get_num_var('b',modSpec),
  0
)

expect_equal(
  get_num_var('s',modSpec),
  1
)

expect_equal(
  get_num_var('z',modSpec),
  0
)

expect_equal(
  get_num_var('kappa',modSpec),
  0
)

expect_equal(
  get_num_var('rho',modSpec,preceding=T),
  0
)

expect_equal(
  get_num_var('tau',modSpec,preceding=T),
  1
)

expect_equal(
  get_num_var('a',modSpec,preceding=T),
  3
)

expect_equal(
  get_num_var('r',modSpec,preceding=T),
  3
)

expect_equal(
  get_num_var('b',modSpec,preceding=T),
  3
)

expect_equal(
  get_num_var('s',modSpec,preceding=T),
  3
)

expect_equal(
  get_num_var('z',modSpec,preceding=T),
  4
)

expect_equal(
  get_num_var('kappa',modSpec,preceding=T),
  4
)

# A one variable, ordinal model that is heteroskedastic
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 1
modSpec$M <- 2
modSpec$hetSpec <- 'sd_x' # heteroskedastic
modSpec$hetGroups <- 1

expect_error(
  check_model(modSpec),
  NA
)

expect_equal(
  get_J(modSpec),
  1
)

expect_equal(
  get_K(modSpec),
  0
)

expect_equal(
  get_Gkappa(modSpec),
  1
)

expect_equal(
  get_Gz(modSpec),
  0
)

expect_equal(
  is_hetero(modSpec),
  T
)

expect_equal(
  is_cdep(modSpec),
  F
)

expect_equal(
  get_z_length(modSpec),
  0
)

expect_equal(
  get_num_var('rho',modSpec),
  1
)

expect_equal(
  get_num_var('tau',modSpec),
  2
)

expect_equal(
  get_num_var('a',modSpec),
  0
)

expect_equal(
  get_num_var('r',modSpec),
  0
)

expect_equal(
  get_num_var('b',modSpec),
  0
)

expect_equal(
  get_num_var('s',modSpec),
  1
)

expect_equal(
  get_num_var('z',modSpec),
  0
)

expect_equal(
  get_num_var('kappa',modSpec),
  1
)

expect_equal(
  get_num_var('rho',modSpec,preceding=T),
  0
)

expect_equal(
  get_num_var('tau',modSpec,preceding=T),
  1
)

expect_equal(
  get_num_var('a',modSpec,preceding=T),
  3
)

expect_equal(
  get_num_var('r',modSpec,preceding=T),
  3
)

expect_equal(
  get_num_var('b',modSpec,preceding=T),
  3
)

expect_equal(
  get_num_var('s',modSpec,preceding=T),
  3
)

expect_equal(
  get_num_var('z',modSpec,preceding=T),
  4
)

expect_equal(
  get_num_var('kappa',modSpec,preceding=T),
  4
)

# A one variable, continuous model that is homoskedastic
modSpec <- list(meanSpec='powLaw')
modSpec$K <- 1
modSpec$hetSpec <- 'none' # homoskedastic

expect_error(
  check_model(modSpec),
  NA
)

expect_equal(
  get_J(modSpec),
  0
)

expect_equal(
  get_K(modSpec),
  1
)

expect_equal(
  get_Gkappa(modSpec),
  0
)

expect_equal(
  get_Gz(modSpec),
  0
)

expect_equal(
  is_hetero(modSpec),
  F
)

expect_equal(
  is_cdep(modSpec),
  F
)

expect_equal(
  get_z_length(modSpec),
  0
)

expect_equal(
  get_num_var('rho',modSpec),
  0
)

expect_equal(
  get_num_var('tau',modSpec),
  0
)

expect_equal(
  get_num_var('a',modSpec),
  1
)

expect_equal(
  get_num_var('r',modSpec),
  1
)

expect_equal(
  get_num_var('b',modSpec),
  1
)

expect_equal(
  get_num_var('s',modSpec),
  1
)

expect_equal(
  get_num_var('z',modSpec),
  0
)

expect_equal(
  get_num_var('kappa',modSpec),
  0
)

expect_equal(
  get_num_var('rho',modSpec,preceding=T),
  0
)

expect_equal(
  get_num_var('tau',modSpec,preceding=T),
  0
)

expect_equal(
  get_num_var('a',modSpec,preceding=T),
  0
)

expect_equal(
  get_num_var('r',modSpec,preceding=T),
  1
)

expect_equal(
  get_num_var('b',modSpec,preceding=T),
  2
)

expect_equal(
  get_num_var('s',modSpec,preceding=T),
  3
)

expect_equal(
  get_num_var('z',modSpec,preceding=T),
  4
)

expect_equal(
  get_num_var('kappa',modSpec,preceding=T),
  4
)

# A one variable, continuous model that is heteroskedastic
modSpec <- list(meanSpec='powLaw')
modSpec$K <- 1
modSpec$hetSpec <- 'sd_x' # heteroskedastic
modSpec$hetGroups <- 1

expect_error(
  check_model(modSpec),
  NA
)

expect_equal(
  get_J(modSpec),
  0
)

expect_equal(
  get_K(modSpec),
  1
)

expect_equal(
  get_Gkappa(modSpec),
  1
)

expect_equal(
  get_Gz(modSpec),
  0
)

expect_equal(
  is_hetero(modSpec),
  T
)

expect_equal(
  is_cdep(modSpec),
  F
)

expect_equal(
  get_z_length(modSpec),
  0
)

expect_equal(
  get_num_var('rho',modSpec),
  0
)

expect_equal(
  get_num_var('tau',modSpec),
  0
)

expect_equal(
  get_num_var('a',modSpec),
  1
)

expect_equal(
  get_num_var('r',modSpec),
  1
)

expect_equal(
  get_num_var('b',modSpec),
  1
)

expect_equal(
  get_num_var('s',modSpec),
  1
)

expect_equal(
  get_num_var('z',modSpec),
  0
)

expect_equal(
  get_num_var('kappa',modSpec),
  1
)

expect_equal(
  get_num_var('rho',modSpec,preceding=T),
  0
)

expect_equal(
  get_num_var('tau',modSpec,preceding=T),
  0
)

expect_equal(
  get_num_var('a',modSpec,preceding=T),
  0
)

expect_equal(
  get_num_var('r',modSpec,preceding=T),
  1
)

expect_equal(
  get_num_var('b',modSpec,preceding=T),
  2
)

expect_equal(
  get_num_var('s',modSpec,preceding=T),
  3
)

expect_equal(
  get_num_var('z',modSpec,preceding=T),
  4
)

expect_equal(
  get_num_var('kappa',modSpec,preceding=T),
  4
)

# A two variable, mixed model that is homoskedastic and conditionally
# independent
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 1
modSpec$K <- 1
modSpec$M <- 2
modSpec$hetSpec <- 'none' # homoskedastic
modSpec$cdepSpec <- 'indep' # conditionally independent

expect_error(
  check_model(modSpec),
  NA
)

expect_equal(
  get_J(modSpec),
  1
)

expect_equal(
  get_K(modSpec),
  1
)

expect_equal(
  get_Gkappa(modSpec),
  0
)

expect_equal(
  get_Gz(modSpec),
  0
)

expect_equal(
  is_hetero(modSpec),
  F
)

expect_equal(
  is_cdep(modSpec),
  F
)

expect_equal(
  get_z_length(modSpec),
  0
)

expect_equal(
  get_non_singleton_groups(modSpec$cdepGroups),
  numeric()
)

expect_equal(
  get_num_var('rho',modSpec),
  1
)

expect_equal(
  get_num_var('tau',modSpec),
  2
)

expect_equal(
  get_num_var('a',modSpec),
  1
)

expect_equal(
  get_num_var('r',modSpec),
  1
)

expect_equal(
  get_num_var('b',modSpec),
  1
)

expect_equal(
  get_num_var('s',modSpec),
  2
)

expect_equal(
  get_num_var('z',modSpec),
  0
)

expect_equal(
  get_num_var('kappa',modSpec),
  0
)

expect_equal(
  get_num_var('rho',modSpec,preceding=T),
  0
)

expect_equal(
  get_num_var('tau',modSpec,preceding=T),
  1
)

expect_equal(
  get_num_var('a',modSpec,preceding=T),
  3
)

expect_equal(
  get_num_var('r',modSpec,preceding=T),
  4
)

expect_equal(
  get_num_var('b',modSpec,preceding=T),
  5
)

expect_equal(
  get_num_var('s',modSpec,preceding=T),
  6
)

expect_equal(
  get_num_var('z',modSpec,preceding=T),
  8
)

expect_equal(
  get_num_var('kappa',modSpec,preceding=T),
  8
)

# A two variable, mixed model that is heteroskedastic and conditionally
# independent
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 1
modSpec$K <- 1
modSpec$M <- 2
modSpec$hetSpec <- 'sd_x' # heteroskedastic
modSpec$hetGroups <- c(1,2)
modSpec$cdepSpec <- 'indep' # conditionally independent

expect_error(
  check_model(modSpec),
  NA
)

expect_equal(
  get_J(modSpec),
  1
)

expect_equal(
  get_K(modSpec),
  1
)

expect_equal(
  get_Gkappa(modSpec),
  2
)

expect_equal(
  get_Gz(modSpec),
  0
)

expect_equal(
  is_hetero(modSpec),
  T
)

expect_equal(
  is_cdep(modSpec),
  F
)

expect_equal(
  get_z_length(modSpec),
  0
)

expect_equal(
  get_num_var('rho',modSpec),
  1
)

expect_equal(
  get_num_var('tau',modSpec),
  2
)

expect_equal(
  get_num_var('a',modSpec),
  1
)

expect_equal(
  get_num_var('r',modSpec),
  1
)

expect_equal(
  get_num_var('b',modSpec),
  1
)

expect_equal(
  get_num_var('s',modSpec),
  2
)

expect_equal(
  get_num_var('z',modSpec),
  0
)

expect_equal(
  get_num_var('kappa',modSpec),
  2
)

expect_equal(
  get_num_var('rho',modSpec,preceding=T),
  0
)

expect_equal(
  get_num_var('tau',modSpec,preceding=T),
  1
)

expect_equal(
  get_num_var('a',modSpec,preceding=T),
  3
)

expect_equal(
  get_num_var('r',modSpec,preceding=T),
  4
)

expect_equal(
  get_num_var('b',modSpec,preceding=T),
  5
)

expect_equal(
  get_num_var('s',modSpec,preceding=T),
  6
)

expect_equal(
  get_num_var('z',modSpec,preceding=T),
  8
)

expect_equal(
  get_num_var('kappa',modSpec,preceding=T),
  8
)

# A two variable, mixed model that is homoskedastic and conditionally dependent
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 1
modSpec$K <- 1
modSpec$M <- 2
modSpec$hetSpec <- 'none' # homoskedastic
modSpec$cdepSpec <- 'dep' # conditionally dependent
modSpec$cdepGroups <- c(1,2)

expect_error(
  check_model(modSpec),
  NA
)

expect_equal(
  get_J(modSpec),
  1
)

expect_equal(
  get_K(modSpec),
  1
)

expect_equal(
  get_Gkappa(modSpec),
  0
)

expect_equal(
  get_Gz(modSpec),
  2
)

expect_equal(
  is_hetero(modSpec),
  F
)

expect_equal(
  is_cdep(modSpec),
  T
)

expect_equal(
  get_z_length(modSpec),
  1
)

expect_equal(
  get_non_singleton_groups(modSpec$cdepGroups),
  numeric()
)

expect_equal(
  get_num_var('rho',modSpec),
  1
)

expect_equal(
  get_num_var('tau',modSpec),
  2
)

expect_equal(
  get_num_var('a',modSpec),
  1
)

expect_equal(
  get_num_var('r',modSpec),
  1
)

expect_equal(
  get_num_var('b',modSpec),
  1
)

expect_equal(
  get_num_var('s',modSpec),
  2
)

expect_equal(
  get_num_var('z',modSpec),
  1
)

expect_equal(
  get_num_var('kappa',modSpec),
  0
)

expect_equal(
  get_num_var('rho',modSpec,preceding=T),
  0
)

expect_equal(
  get_num_var('tau',modSpec,preceding=T),
  1
)

expect_equal(
  get_num_var('a',modSpec,preceding=T),
  3
)

expect_equal(
  get_num_var('r',modSpec,preceding=T),
  4
)

expect_equal(
  get_num_var('b',modSpec,preceding=T),
  5
)

expect_equal(
  get_num_var('s',modSpec,preceding=T),
  6
)

expect_equal(
  get_num_var('z',modSpec,preceding=T),
  8
)

expect_equal(
  get_num_var('kappa',modSpec,preceding=T),
  9
)

# A two variable, mixed model that is heteroskedastic and conditionally
# dependent
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 1
modSpec$K <- 1
modSpec$M <- 2
modSpec$hetSpec <- 'sd_x' # heteroskedastic
modSpec$hetGroups <- c(1,2)
modSpec$cdepSpec <- 'dep' # conditionally dependent
modSpec$cdepGroups <- c(1,2)

expect_error(
  check_model(modSpec),
  NA
)

expect_equal(
  get_J(modSpec),
  1
)

expect_equal(
  get_K(modSpec),
  1
)

expect_equal(
  get_Gkappa(modSpec),
  2
)

expect_equal(
  get_Gz(modSpec),
  2
)

expect_equal(
  is_hetero(modSpec),
  T
)

expect_equal(
  is_cdep(modSpec),
  T
)

expect_equal(
  get_z_length(modSpec),
  1
)

expect_equal(
  get_non_singleton_groups(modSpec$cdepGroups),
  numeric()
)

expect_equal(
  get_num_var('rho',modSpec),
  1
)

expect_equal(
  get_num_var('tau',modSpec),
  2
)

expect_equal(
  get_num_var('a',modSpec),
  1
)

expect_equal(
  get_num_var('r',modSpec),
  1
)

expect_equal(
  get_num_var('b',modSpec),
  1
)

expect_equal(
  get_num_var('s',modSpec),
  2
)

expect_equal(
  get_num_var('z',modSpec),
  1
)

expect_equal(
  get_num_var('kappa',modSpec),
  2
)

expect_equal(
  get_num_var('rho',modSpec,preceding=T),
  0
)

expect_equal(
  get_num_var('tau',modSpec,preceding=T),
  1
)

expect_equal(
  get_num_var('a',modSpec,preceding=T),
  3
)

expect_equal(
  get_num_var('r',modSpec,preceding=T),
  4
)

expect_equal(
  get_num_var('b',modSpec,preceding=T),
  5
)

expect_equal(
  get_num_var('s',modSpec,preceding=T),
  6
)

expect_equal(
  get_num_var('z',modSpec,preceding=T),
  8
)

expect_equal(
  get_num_var('kappa',modSpec,preceding=T),
  9
)

# A four variable, mixed model that is heteroskedastic and conditionally
# dependent, with one variable being NA for each.
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 2
modSpec$K <- 2
modSpec$M <- c(2,3)
modSpec$hetSpec <- 'sd_x' # heteroskedastic
modSpec$hetGroups <- c(1,2,NA,1)
modSpec$cdepSpec <- 'dep' # conditionally dependent
modSpec$cdepGroups <- c(1,NA,2,2)

expect_error(
  check_model(modSpec),
  NA
)

expect_equal(
  get_J(modSpec),
  2
)

expect_equal(
  get_K(modSpec),
  2
)

expect_equal(
  get_Gkappa(modSpec),
  2
)

expect_equal(
  get_Gz(modSpec),
  2
)

expect_equal(
  is_hetero(modSpec),
  T
)

expect_equal(
  is_cdep(modSpec),
  T
)

expect_equal(
  get_z_length(modSpec),
  2
)

expect_equal(
  get_non_singleton_groups(modSpec$cdepGroups),
  2
)

expect_equal(
  get_num_var('rho',modSpec),
  2
)

expect_equal(
  get_num_var('tau',modSpec),
  5
)

expect_equal(
  get_num_var('a',modSpec),
  2
)

expect_equal(
  get_num_var('r',modSpec),
  2
)

expect_equal(
  get_num_var('b',modSpec),
  2
)

expect_equal(
  get_num_var('s',modSpec),
  4
)

expect_equal(
  get_num_var('z',modSpec),
  2
)

expect_equal(
  get_num_var('kappa',modSpec),
  2
)

expect_equal(
  get_num_var('rho',modSpec,preceding=T),
  0
)

expect_equal(
  get_num_var('tau',modSpec,preceding=T),
  2
)

expect_equal(
  get_num_var('a',modSpec,preceding=T),
  7
)

expect_equal(
  get_num_var('r',modSpec,preceding=T),
  9
)

expect_equal(
  get_num_var('b',modSpec,preceding=T),
  11
)

expect_equal(
  get_num_var('s',modSpec,preceding=T),
  13
)

expect_equal(
  get_num_var('z',modSpec,preceding=T),
  17
)

expect_equal(
  get_num_var('kappa',modSpec,preceding=T),
  19
)

# A six variable, mixed model
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 3
modSpec$K <- 3
modSpec$M <- c(2,3,2)
modSpec$hetSpec <- 'sd_x' # heteroskedastic
modSpec$hetGroups <- c(NA,1,2,2,NA,3)
modSpec$cdepSpec <- 'dep' # conditionally dependent
modSpec$cdepGroups <- c(1,2,1,3,NA,2)

expect_error(
  check_model(modSpec),
  NA
)

expect_equal(
  get_J(modSpec),
  3
)

expect_equal(
  get_K(modSpec),
  3
)

expect_equal(
  get_Gkappa(modSpec),
  3
)

expect_equal(
  get_Gz(modSpec),
  3
)

expect_equal(
  is_hetero(modSpec),
  T
)

expect_equal(
  is_cdep(modSpec),
  T
)

expect_equal(
  get_z_length(modSpec),
  5
)

expect_equal(
  get_non_singleton_groups(modSpec$cdepGroups),
  c(1,2)
)

expect_equal(
  get_num_var('rho',modSpec),
  3
)

expect_equal(
  get_num_var('tau',modSpec),
  7
)

expect_equal(
  get_num_var('a',modSpec),
  3
)

expect_equal(
  get_num_var('r',modSpec),
  3
)

expect_equal(
  get_num_var('b',modSpec),
  3
)

expect_equal(
  get_num_var('s',modSpec),
  6
)

expect_equal(
  get_num_var('z',modSpec),
  5
)

expect_equal(
  get_num_var('kappa',modSpec),
  3
)

expect_equal(
  get_num_var('rho',modSpec,preceding=T),
  0
)

expect_equal(
  get_num_var('tau',modSpec,preceding=T),
  3
)

expect_equal(
  get_num_var('a',modSpec,preceding=T),
  10
)

expect_equal(
  get_num_var('r',modSpec,preceding=T),
  13
)

expect_equal(
  get_num_var('b',modSpec,preceding=T),
  16
)

expect_equal(
  get_num_var('s',modSpec,preceding=T),
  19
)

expect_equal(
  get_num_var('z',modSpec,preceding=T),
  25
)

expect_equal(
  get_num_var('kappa',modSpec,preceding=T),
  30
)

# Test get_var_index
# Conditionally independent / homoskedastic
modSpec_io <- list(meanSpec='powLaw')
modSpec_io$J <- 1
modSpec_io$K <- 1
modSpec_io$M <- 2
modSpec_io$cdepSpec <- 'indep'
modSpec_io$hetSpec  <- 'none'

# Conditionally independent / heteroskedastic
modSpec_ie <- list(meanSpec='powLaw')
modSpec_ie$J <- 1
modSpec_ie$K <- 1
modSpec_ie$M <- 2
modSpec_ie$cdepSpec <- 'indep'
modSpec_ie$hetSpec  <- 'sd_x'
modSpec_ie$hetGroups  <- c(1,2)

# Conditionally dependent / homoskedastic
modSpec_do <- list(meanSpec='powLaw')
modSpec_do$J <- 1
modSpec_do$K <- 1
modSpec_do$M <- 2
modSpec_do$cdepSpec <- 'dep'
modSpec_do$cdepGroups <- c(1,2)
modSpec_do$hetSpec  <- 'none'

# Conditionally dependent / heteroskedastic
modSpec_de <- list(meanSpec='powLaw')
modSpec_de$J <- 1
modSpec_de$K <- 1
modSpec_de$M <- 2
modSpec_de$cdepSpec <- 'dep'
modSpec_de$cdepGroups <- c(1,2)
modSpec_de$hetSpec  <- 'sd_x'
modSpec_de$hetGroups  <- c(1,2)

modSpecList <- list()
modSpecList[[1]] <- modSpec_io
modSpecList[[2]] <- modSpec_ie
modSpecList[[3]] <- modSpec_do
modSpecList[[4]] <- modSpec_de

# rho
for(modSpec in modSpecList) {
  expect_equal(
    get_var_index('rho',modSpec,j=1),
    1
  )
}

for(modSpec in modSpecList) {
  expect_equal(
    get_var_index('rho',modSpec),
    1
  )
}

for(modSpec in modSpecList) {
  expect_error(
    get_var_index('rho',modSpec,k=1),
    'Unsupported variable for k being specified'
  )
}

for(modSpec in modSpecList) {
  expect_error(
    get_var_index('rho',modSpec,j=2),
    'j is not between 1 and J'
  )
}

# tau
for(modSpec in modSpecList) {
  expect_equal(
    get_var_index('tau',modSpec,j=1),
    c(2,3)
  )
}

for(modSpec in modSpecList) {
  expect_equal(
    get_var_index('tau',modSpec),
    c(2,3)
  )
}

for(modSpec in modSpecList) {
  expect_error(
    get_var_index('tau',modSpec,k=1),
    'Unsupported variable for k being specified'
  )
}

for(modSpec in modSpecList) {
  expect_error(
    get_var_index('tau',modSpec,j=2),
    'j is not between 1 and J'
  )
}

# a, r, b
contVar <- c('a','r','b')
for(vv in 1:length(contVar)) {
  for(modSpec in modSpecList) {
    expect_equal(
      get_var_index(contVar[vv],modSpec,k=1),
      4 + vv-1
    )
  }

  for(modSpec in modSpecList) {
    expect_equal(
      get_var_index(contVar[vv],modSpec),
      4 + vv-1
    )
  }

  for(modSpec in modSpecList) {
    expect_error(
      get_var_index(contVar[vv],modSpec,j=1),
      'Unsupported variable for j being specified'
    )
  }

  for(modSpec in modSpecList) {
    expect_error(
      get_var_index(contVar[vv],modSpec,k=2),
      'k is not between 1 and K'
    )
  }
}

# test s
for(modSpec in modSpecList) {
  expect_equal(
    get_var_index('s',modSpec,j=1),
    7
  )
}

for(modSpec in modSpecList) {
  expect_equal(
    get_var_index('s',modSpec,k=1),
    8
  )
}

for(modSpec in modSpecList) {
  expect_equal(
    get_var_index('s',modSpec),
    c(7,8)
  )
}

for(modSpec in modSpecList) {
  expect_error(
    get_var_index('s',modSpec,j=2),
    'j is not between 1 and J'
  )
}

for(modSpec in modSpecList) {
  expect_error(
    get_var_index('s',modSpec,k=2),
    'k is not between 1 and K'
  )
}

# test kappa
expect_error(
  get_var_index('kappa',modSpec_io,j=1),
  'kappa requested but model is not heteroskedastic'
)

expect_error(
  get_var_index('kappa',modSpec_io,k=1),
  'kappa requested but model is not heteroskedastic'
)

expect_equal(
  get_var_index('kappa',modSpec_ie,j=1),
  9
)

expect_equal(
  get_var_index('kappa',modSpec_ie,k=1),
  10
)

expect_error(
  get_var_index('kappa',modSpec_do,j=1),
  'kappa requested but model is not heteroskedastic'
)

expect_error(
  get_var_index('kappa',modSpec_do,k=1),
  'kappa requested but model is not heteroskedastic'
)

expect_equal(
  get_var_index('kappa',modSpec_de,j=1),
  10
)

expect_equal(
  get_var_index('kappa',modSpec_de,k=1),
  11
)

# test z
expect_error(
  get_var_index('z',modSpec_io,i1=1,i2=2),
  'z requested but model is not conditionally dependent'
)

expect_error(
  get_var_index('z',modSpec_ie,i1=1,i2=2),
  'z requested but model is not conditionally dependent'
)

expect_equal(
  get_var_index('z',modSpec_do,i1=1,i2=2),
  9
)

expect_equal(
  get_var_index('z',modSpec_de,i1=1,i2=2),
  9
)

expect_equal(
  get_var_index('z',modSpec_do),
  9
)

expect_equal(
  get_var_index('z',modSpec_de),
  9
)

# Test the functions that depend on both modSpec and the parameter vector. These
# functions are:
#
# theta_y_list2vect
# theta_y_vect2list
# theta_y_vect2list
# theta_y_constr2unconstr
# theta_y_unconstr2constr
# get_kappa_full
# get_Sigma0

# Create a model and parameter vector with two ordinal variables and two
# continuous variables, for which there are two groups of conditionally
# dependent variables and all variables have separate heteroskedastic
# parameters.
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 2
modSpec$M <- c(2,2)
modSpec$K <- 2
modSpec$cdepSpec <- 'dep' # conditionally dependent
modSpec$cdepGroups <- c(1,2,1,2)
modSpec$hetSpec <- 'sd_x' # homoskedastic
modSpec$hetGroups <- c(1,2,3,4)

rho      <- c(0.5,0.75)
tau      <- list()
tau[[1]] <- c(-2.0,2.5)
tau[[2]] <- c(-1.0,3.5)
a        <- c(0.1,1.2)
r        <- c(0.3,0.8)
b        <- c(-.25,.25)
s        <- c(.2,.4,.3,.5)
zns      <- c(-.1,.2)
zcr      <- c(.05)
z        <- c(zns,zcr)
kappa    <- c(.02,.01,.03,.01)

th_y_vect <- c(rho,unlist(tau),a,r,b,s,z,kappa)
th_y_list <- theta_y_vect2list(th_y_vect,modSpec)

expect_equal(
  th_y_list$rho,
  rho
)

expect_equal(
  th_y_list$tau,
  tau
)

expect_equal(
  th_y_list$a,
  a
)

expect_equal(
  th_y_list$r,
  r
)

expect_equal(
  th_y_list$b,
  b
)


expect_equal(
  th_y_list$z,
  z
)

expect_equal(
  th_y_list$kappa,
  kappa
)

# Test that theta_y_list2vect yields the original vector
th_y_vect2 <- theta_y_list2vect(th_y_list)
expect_equal(
  th_y_vect,
  th_y_vect2
)

# Test get_kappa_full by building it "from scratch". For this example,
# kappa_full equals kappa
kappa_full <- kappa
kappa_full_matrix <- matrix(kappa_full) %*% t(matrix(kappa_full))

expect_equal(
  get_kappa_full(th_y_list),
  kappa_full
)

expect_equal(
  get_kappa_full(th_y_list,asMatrix=T),
  kappa_full_matrix
)

expect_equal(
  get_kappa_full(th_y_vect,modSpec,asMatrix=T),
  kappa_full_matrix
)

# Test get_Sigma0 by building Sigma0 "from scratch"
Sigma0 <- diag(s^2)
Sigma0[1,2] <- s[1]*s[2]*z[3] # cross group (index 3)
Sigma0[1,3] <- s[1]*s[3]*z[1] # within group 1 (index 1)
Sigma0[1,4] <- s[1]*s[4]*z[3] # cross group (index 3)
Sigma0[2,3] <- s[2]*s[3]*z[3] # cross group (index 3)
Sigma0[2,4] <- s[2]*s[4]*z[2] # within group 2 (index 2)
Sigma0[3,4] <- s[3]*s[4]*z[3] # cross group (index 3)
Sigma0[2,1] <- Sigma0[1,2]
Sigma0[3,1] <- Sigma0[1,3]
Sigma0[4,1] <- Sigma0[1,4]
Sigma0[3,2] <- Sigma0[2,3]
Sigma0[4,2] <- Sigma0[2,4]
Sigma0[4,3] <- Sigma0[3,4]

expect_equal(
  get_Sigma0(th_y_list),
  Sigma0
)
expect_equal(
  get_Sigma0(th_y_vect,modSpec),
  Sigma0
)

x <- c(1.5,2.5)
heteroTerm1 <- matrix(1 + kappa_full*x[1]) # a column vector
heteroTerm2 <- matrix(1 + kappa_full*x[2]) # a column vector
Sigma1 <- Sigma0 * (heteroTerm1 %*% t(heteroTerm1))
Sigma2 <- Sigma0 * (heteroTerm2 %*% t(heteroTerm2))
Sigma <- array(NA,c(2,nrow(Sigma0),ncol(Sigma0)))
Sigma[1,,] <- Sigma1
Sigma[2,,] <- Sigma2

expect_equal(
  get_Sigma(th_y_vect,x[1],modSpec),
  Sigma1
)

expect_equal(
  get_Sigma(th_y_vect,x[2],modSpec),
  Sigma2
)

expect_equal(
  get_Sigma(th_y_vect,x,modSpec),
  Sigma
)

# Test theta_y_constr2unconstr
# Do the calculation directly, then check the function
th_y_vect_unconstr_direct <- c(log(rho),tau[[1]][1],log(diff(tau[[1]])),tau[[2]][1],log(diff(tau[[2]])),log(a),log(r),b,log(s),gtools::logit((1+z)/2),log(kappa))

th_y_vect_unconstr <- theta_y_constr2unconstr(th_y_vect,modSpec)

expect_equal(
  th_y_vect_unconstr_direct,
  th_y_vect_unconstr
)

# Test theta_y_constr2unconstr by ensuring that the reverse transform yields the
# original vector
th_y_vect2 <- theta_y_unconstr2constr(th_y_vect_unconstr,modSpec)

expect_equal(
  th_y_vect2,
  th_y_vect
)

# Create a model and parameter vector with two ordinal variables and two
# continuous variables, for which there are two groups of conditionally
# dependent variables and only one variable is heteroskedastic
# parameters.
modSpec <- list(meanSpec='powLaw')
modSpec$J <- 2
modSpec$M <- c(2,2)
modSpec$K <- 2
modSpec$cdepSpec <- 'dep' # conditionally dependent
modSpec$cdepGroups <- c(1,2,1,2)
modSpec$hetSpec <- 'sd_x' # homoskedastic
modSpec$hetGroups <- c(NA,1,NA,NA)

rho      <- c(0.5,0.75)
tau      <- list()
tau[[1]] <- c(-2.0,2.5)
tau[[2]] <- c(-1.0,3.5)
a        <- c(0.1,1.2)
r        <- c(0.3,0.8)
b        <- c(-.25,.25)
s        <- c(.2,.4,.3,.5)
zns      <- c(-.1,.2)
zcr      <- c(.05)
z        <- c(zns,zcr)
kappa    <- c(.02)

th_y_vect <- c(rho,unlist(tau),a,r,b,s,z,kappa)
th_y_list <- theta_y_vect2list(th_y_vect,modSpec)

expect_equal(
  th_y_list$rho,
  rho
)

expect_equal(
  th_y_list$tau,
  tau
)

expect_equal(
  th_y_list$a,
  a
)

expect_equal(
  th_y_list$r,
  r
)

expect_equal(
  th_y_list$b,
  b
)


expect_equal(
  th_y_list$z,
  z
)

expect_equal(
  th_y_list$kappa,
  kappa
)

# Test that theta_y_list2vect yields the original vector
th_y_vect2 <- theta_y_list2vect(th_y_list)
expect_equal(
  th_y_vect,
  th_y_vect2
)

# Test get_kappa_full by building it "from scratch".
kappa_full <- c(0,kappa,0,0)
kappa_full_matrix <- matrix(kappa_full) %*% t(matrix(kappa_full))

expect_equal(
  get_kappa_full(th_y_list),
  kappa_full
)

expect_equal(
  get_kappa_full(th_y_list,asMatrix=T),
  kappa_full_matrix
)

expect_equal(
  get_kappa_full(th_y_vect,modSpec,asMatrix=T),
  kappa_full_matrix
)

# Test get_Sigma0 by building Sigma0 "from scratch"
Sigma0 <- diag(s^2)
Sigma0[1,2] <- s[1]*s[2]*z[3] # cross group (index 3)
Sigma0[1,3] <- s[1]*s[3]*z[1] # within group 1 (index 1)
Sigma0[1,4] <- s[1]*s[4]*z[3] # cross group (index 3)
Sigma0[2,3] <- s[2]*s[3]*z[3] # cross group (index 3)
Sigma0[2,4] <- s[2]*s[4]*z[2] # within group 2 (index 2)
Sigma0[3,4] <- s[3]*s[4]*z[3] # cross group (index 3)
Sigma0[2,1] <- Sigma0[1,2]
Sigma0[3,1] <- Sigma0[1,3]
Sigma0[4,1] <- Sigma0[1,4]
Sigma0[3,2] <- Sigma0[2,3]
Sigma0[4,2] <- Sigma0[2,4]
Sigma0[4,3] <- Sigma0[3,4]

expect_equal(
  get_Sigma0(th_y_list),
  Sigma0
)

expect_equal(
  get_Sigma0(th_y_vect,modSpec),
  Sigma0
)

x <- c(1.5,2.5)
heteroTerm1 <- matrix(1 + kappa_full*x[1]) # a column vector
heteroTerm2 <- matrix(1 + kappa_full*x[2]) # a column vector
Sigma1 <- Sigma0 * (heteroTerm1 %*% t(heteroTerm1))
Sigma2 <- Sigma0 * (heteroTerm2 %*% t(heteroTerm2))
Sigma <- array(NA,c(2,nrow(Sigma0),ncol(Sigma0)))
Sigma[1,,] <- Sigma1
Sigma[2,,] <- Sigma2

expect_equal(
  get_Sigma(th_y_vect,x[1],modSpec),
  Sigma1
)

expect_equal(
  get_Sigma(th_y_vect,x[2],modSpec),
  Sigma2
)

expect_equal(
  get_Sigma(th_y_vect,x,modSpec),
  Sigma
)

# Test theta_y_constr2unconstr
# Do the calculation directly, then check the function
th_y_vect_unconstr_direct <- c(log(rho),tau[[1]][1],log(diff(tau[[1]])),tau[[2]][1],log(diff(tau[[2]])),log(a),log(r),b,log(s),gtools::logit((1+z)/2),log(kappa))

th_y_vect_unconstr <- theta_y_constr2unconstr(th_y_vect,modSpec)

expect_equal(
  th_y_vect_unconstr_direct,
  th_y_vect_unconstr
)

# Test theta_y_constr2unconstr by ensuring that the reverse transform yields the
# original vector
th_y_vect2 <- theta_y_unconstr2constr(th_y_vect_unconstr,modSpec)

expect_equal(
  th_y_vect2,
  th_y_vect
)
