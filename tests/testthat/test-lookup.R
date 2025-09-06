
test_that("CE lookup returns reasonable functions", {
  set.seed(2)
  lookup <- build_ce_lookup(
    n1=80, n2=80,
    muX_C=3, sigmaX=0.2, gamma0=2, gamma1=5.5,
    mu0=22, sigmaY=1, theta=1, eta=0, piZ=0.5,
    pi_fixed=0.5, error_type="normal", rho_XY=0.3,
    alpha_one_sided=0.05, B_ref=20000, batch_size=2000,
    z1_grid=seq(-3,3,by=0.1), min_in_bin=100
  )
  expect_true(is.numeric(lookup$b_ref) && is.finite(lookup$b_ref))
  expect_true(is.function(lookup$e_fun))
  expect_true(is.function(lookup$c_fun))
  # e(t1) should be small near very negative t1, increasing towards b_ref
  e_low  <- lookup$e_fun(min(lookup$z1_grid))
  e_high <- lookup$e_fun(max(lookup$z1_grid[lookup$z1_grid < lookup$b_ref]))
  expect_true(e_low < e_high)
})
