test_that("rownames of data frame are the name of list", {
    inputlist=list("1.34"=c(1,2,3),"8.5"=c(3,2,3))
    outdf=as_data_frame_rowname(inputlist)
    expect_equal(names(inputlist),rownames(outdf))
})
