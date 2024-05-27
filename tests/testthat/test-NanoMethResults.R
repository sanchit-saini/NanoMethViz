test_that("NanoMethResults getters and setters work", {
    # setup
    nmr <- load_example_nanomethresult()

    # test
    expect_true(fs::file_exists(NanoMethViz::methy(nmr)))
    expect_type(NanoMethViz::methy(nmr), "character")
    expect_s3_class(NanoMethViz::exons(nmr), "data.frame")
    expect_s3_class(NanoMethViz::samples(nmr), "data.frame")

    expect_no_warning(
        NanoMethResult(
            NanoMethViz::methy(nmr),
            NanoMethViz::samples(nmr)
        )
    )

    expect_no_warning(methy(nmr) <- methy(nmr))
    expect_no_warning(samples(nmr) <- samples(nmr))
    expect_no_warning(exons(nmr) <- exons(nmr))

    expect_error(methy(nmr) <- "invalid_path")
    expect_error(exons(nmr) <- exons(nmr)[, -"strand"])
})
