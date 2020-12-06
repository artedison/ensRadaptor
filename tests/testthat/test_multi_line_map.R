test_that("rownames of data frame are the name of list", {
    string=c("aabb","cc","eeab","ce")
    pattern="b\\nc"
    expect_equal(c(2,4),multi_line_map(string=string,pattern=pattern))
    pattern="eea"
    expect_equal(c(3),multi_line_map(string=string,pattern=pattern))
    pattern="eea\\n"
    expect_null(multi_line_map(string=string,pattern=pattern))
})
