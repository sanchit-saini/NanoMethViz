test_that("methy_to_bsseq works", {
    # setup
    nmr <- load_example_nanomethresult()
    bss <- methy_to_bsseq(nmr)

    # test
    expect_no_error(methy_to_bsseq(methy(nmr)))
    expect_true(is(methy_to_bsseq(nmr), "BSseq"))
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
    expect_true(is(edger_counts <- bsseq_to_edger(bss), "matrix"))
    expect_true(is(lmr <- bsseq_to_log_methy_ratio(bss), "matrix"))

    expect_equal(ncol(edger_counts), 2 * ncol(bss))
    expect_equal(ncol(lmr), ncol(bss))
})
