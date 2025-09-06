
test_that("Reference boundary attains approx alpha under H0", {
  skip_on_cran()
  
  set.seed(1)
  n1 <- 100; n2 <- 100
  pars <- list(muX_C=3, sigmaX=0.15, gamma0=2.0, gamma1=5.5,
               mu0=22, sigmaY=1.0, theta=1.0, eta=0.0,
               piZ=0.5, pi_fixed=0.5, rho_XY=0.3)
  
  lookup <- build_ce_lookup(n1, n2,
                            muX_C=pars$muX_C, sigmaX=pars$sigmaX, gamma0=pars$gamma0, gamma1=pars$gamma1,
                            mu0=pars$mu0, sigmaY=pars$sigmaY, theta=pars$theta, eta=pars$eta, piZ=pars$piZ,
                            pi_fixed=pars$pi_fixed, error_type="normal", rho_XY=pars$rho_XY,
                            alpha_one_sided=0.05, B_ref=30000, batch_size=3000,
                            z1_grid=seq(-3, 3, by=0.1), min_in_bin=200, h0=0.05
  )
  
  # Monte Carlo estimate of overall rejection under H0
  R <- 2000
  rej <- logical(R)
  for (r in seq_len(R)) {
    # Under H0: delta=0 when eta=0
    muX_T <- pars$muX_C + 0 / pars$gamma1
    s1 <- gen_stage1_and_T1(n1, pars$muX_C, muX_T, pars$sigmaX, pars$gamma0, pars$gamma1)
    if (s1$T1 >= lookup$b_ref) {
      rej[r] <- TRUE
    } else {
      df1 <- data.frame(stage="S1",
                        T=c(rep(0,n1),rep(1,n1)),
                        Z=0,
                        Y=c(pars$mu0 + rnorm(n1,0,pars$sigmaY),
                            pars$mu0 + 0 + rnorm(n1,0,pars$sigmaY)))
      # S2 null
      Zc2 <- rbinom(n2,1,pars$piZ); Zt2 <- rbinom(n2,1,pars$piZ)
      Yc2 <- pars$mu0 + pars$theta*Zc2 + rnorm(n2,0,pars$sigmaY)
      Yt2 <- pars$mu0 + 0 + pars$theta*Zt2 + 0*Zt2 + rnorm(n2,0,pars$sigmaY)
      df2 <- data.frame(stage="S2",
                        T=c(rep(0,n2),rep(1,n2)),
                        Z=c(Zc2,Zt2),
                        Y=c(Yc2,Yt2))
      df <- rbind(df1,df2)
      fa <- final_analysis(df, pi_fixed=pars$pi_fixed)
      cz1 <- lookup$c_fun(s1$T1)
      rej[r] <- isTRUE(fa$Zf >= cz1)
    }
  }
  alpha_hat <- mean(rej)
  expect_true(abs(alpha_hat - 0.05) < 0.015)
})
