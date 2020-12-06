test_that("a small example for the table parser", {
    lines_models=c("1 2 3 4 5","6 7 8  9 10","11 12 13\t14 15")
    restab=tabparse(lines_models,c(1,3),c(2,4))
    exptab=as.data.frame(matrix(c(2,12,4,14),nrow=2))
    expect_equivalent(restab,exptab)
})
