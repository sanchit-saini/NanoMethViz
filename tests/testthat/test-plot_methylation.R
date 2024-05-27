test_that("Plotting gene works", {
    # setup
    nmr <- load_example_nanomethresult()

    # test
    p_gene <- expect_no_warning(plot_gene(nmr, "Peg3"))
    p_gene2 <- expect_no_warning(plot_gene(nmr, "Peg3", spaghetti = TRUE))
    expect_s3_class(p_gene, "patchwork")
    expect_s3_class(p_gene, "ggplot")

    expect_s3_class(p_gene2, "patchwork")
    expect_s3_class(p_gene2, "ggplot")
})
