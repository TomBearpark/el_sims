## A little bit of timing
pacman::p_load('microbenchmark')



est <- gmm(g = moments, x = data, t0 = init, wmatrix = "ident")
est2 <- gmm(g = moments_mom_alt, x = data, t0 = init, wmatrix = "ident")