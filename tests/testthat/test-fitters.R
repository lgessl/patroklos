test_that("zeroSum_with_pheno() works", {

  n_samples <- 10
  n_genes <- 2
  n_pheno_features <- 2
  lambda <- 1

  x <- matrix(rnorm(n_samples * (n_genes+n_pheno_features)), nrow = n_samples)
  colnames(x) <- c(paste0("gene_", 1:n_genes), paste0("pheno_", 1:n_pheno_features, "++"))
  y_binary <- sample(c(0, 1), n_samples, replace = TRUE)
  y_survival <- cbind(rnorm(n_samples), y_binary)
  zs_weights <- rep(c(1, 0), c(n_genes, n_pheno_features))

  fit_obj <- zeroSumEI(
    x = x, y = y_binary, family = "binomial", nFold = 1, lambda = lambda
  )
  expect_equal(fit_obj$zeroSum.weights, zs_weights)

  fit_obj <- zeroSumEI(
    x = x, y = y_survival, family = "cox", nFold = 1, lambda = lambda
  )
  expect_equal(fit_obj$zeroSum.weights, zs_weights)

  colnames(x) <- as.character(seq(n_genes + n_pheno_features))
  zs_weights <- rep(1, ncol(x))
  
  fit_obj <- zeroSumEI(
    x = x, y = y_binary, family = "binomial", nFold = 1, lambda = lambda
  )
  expect_equal(fit_obj$zeroSum.weights, zs_weights)

})
