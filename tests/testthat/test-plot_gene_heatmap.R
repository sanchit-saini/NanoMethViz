test_that("Plotting gene methylation heatmap works", {
    # setup
    nmr <- load_example_nanomethresult()

    # test
    p <- expect_no_warning(plot_gene_heatmap(nmr, "Peg3"))
    expect_s3_class(p, "ggplot")

    # test for a bug whereby samples not present in data cause function to hang
    nmr_extra_sample <- load_example_nanomethresult()
    samples(nmr_extra_sample) <- bind_rows(
        samples(nmr_extra_sample),
        c(sample = "foo", group = "bar")
    )
    p <- expect_no_warning(plot_gene_heatmap(nmr_extra_sample, "Peg3"))
    expect_s3_class(p, "ggplot")
})
