test_that("Plotting feature clustering", {
    nmr <- load_example_nanomethresult()
    gene_anno <- exons_to_genes(NanoMethViz::exons(nmr))

    expect_no_warning(cluster_regions(nmr, gene_anno, centers = 2, grid_method = "uniform"))
    expect_no_warning(cluster_regions(nmr, gene_anno, centers = 3, grid_method = "uniform"))
    expect_no_warning(cluster_regions(nmr, gene_anno, centers = 2, grid_method = "density"))
})
