test_that("Plotting gene aggregates works", {
    nmr <- load_example_nanomethresult()
    expect_no_error(plot_agg_genes(nmr))
    expect_no_error(plot_agg_tss(nmr))
    expect_no_error(plot_agg_tes(nmr))

    # empty exons should throw error
    exons(nmr) <- tibble::tibble(
        gene_id = character(),
        chr = character(),
        strand = character(),
        start = integer(),
        end = integer(),
        transcript_id = character(),
        symbol = character()
    )
    expect_error(plot_agg_genes(nmr))
    expect_error(plot_agg_tss(nmr))
    expect_error(plot_agg_tes(nmr))
})
