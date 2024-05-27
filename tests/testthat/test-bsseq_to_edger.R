test_that("Conversion functions work", {
    # setup
    nmr <- load_example_nanomethresult()

    # tests
    bsseq <- expect_no_warning(methy_to_bsseq(nmr))
    expect_no_warning(methy_to_bsseq(methy(nmr)))
    regions <- expect_no_warning(exons_to_genes(nmr))
    expect_no_warning(exons_to_genes(NanoMethViz::exons(nmr)))

    # convert via bsseq
    edger_site <- expect_no_warning(bsseq_to_edger(bsseq))
    expect_ncol(edger_site, 2*ncol(bsseq))
    expect_nrow(edger_site, nrow(bsseq))

    # convert directly from nmr
    edger_direct_site <- expect_no_warning(methy_to_edger(nmr))
    expect_ncol(edger_direct_site, 2*ncol(bsseq))
    expect_nrow(edger_direct_site, nrow(bsseq))

    # should be identical
    expect_identical(edger_site, edger_direct_site)

    # convert to edgeR
    edger_region <- expect_no_warning(bsseq_to_edger(bsseq, regions))
    expect_ncol(edger_region, 2*ncol(bsseq))
    expect_nrow(edger_region, nrow(regions))

    # convert to log-methylation-ratio
    lmr_regions <- expect_no_warning(bsseq_to_log_methy_ratio(bsseq, regions))
    expect_ncol(lmr_regions, ncol(bsseq))
    expect_nrow(lmr_regions, nrow(regions))

    expect_warning(
        bsseq_to_log_methy_ratio(bsseq, regions, prior_count = 0.5),
        "prior_count of 1 or higher is recommended"
   )

    expect_error(bsseq_to_log_methy_ratio(bsseq, regions, prior_count = -1))
})
