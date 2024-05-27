test_that("methy_to_bsseq works", {
    # setup
    nmr <- load_example_nanomethresult()
    bss <- methy_to_bsseq(nmr)

    # test
    expect_no_error(methy_to_bsseq(methy(nmr)))
    expect_s4_class(methy_to_bsseq(nmr), "BSseq")
    expect_equal(ncol(bss), 6)
    expect_equal(
        nrow(SummarizedExperiment::colData(bss)),
        nrow(NanoMethViz::samples(nmr))
    )
    expect_equal(
        ncol(SummarizedExperiment::colData(bss)),
        ncol(NanoMethViz::samples(nmr))
    )
    expect_equal(
        colnames(SummarizedExperiment::colData(bss)),
        colnames(NanoMethViz::samples(nmr))
    )
})

test_that("bsseq_to_* works", {
    nmr <- load_example_nanomethresult()
    bss <- methy_to_bsseq(NanoMethViz::methy(nmr))

    # test
    edger_counts <- expect_no_warning(bsseq_to_edger(bss))
    expect_ncol(edger_counts, 12)
    expect_nrow(edger_counts, 4778)
    lmr <- expect_no_warning(bsseq_to_log_methy_ratio(bss))
    expect_ncol(lmr, 6)
    expect_nrow(lmr, 4778)

    expect_equal(ncol(edger_counts), 2 * ncol(bss))
    expect_equal(ncol(lmr), ncol(bss))
})
