# bin	611	smallint(6)	range	Indexing field to speed chromosome range queries.
# chrom	chr1	varchar(255)	values	Reference sequence chromosome or scaffold
# chromStart	3531624	int(10) unsigned	range	Start position in chromosome
# chromEnd	3531843	int(10) unsigned	range	End position in chromosome
# name	CpG: 27	varchar(255)	values	CpG Island
# length	219	int(10) unsigned	range	Island Length
# cpgNum	27	int(10) unsigned	range	Number of CpGs in island
# gcNum	167	int(10) unsigned	range	Number of C and G in island
# perCpg	24.7	float	range	Percentage of island that is CpG
# perGc	76.3	float	range	Percentage of island that is C or G
# obsExp	0.86	float	range	Ratio of observed(cpgNum) to expected(numC*numG/length) CpG in island

# exons
# the data.frame of exon information containing at least columns gene_id, chr, strand, start, end, transcript_id and symbol.

col_names <- c(
    "bin",
    "chrom",
    "chromStart",
    "chromEnd",
    "name",
    "length",
    "cpgNum",
    "gcNum",
    "perCpg",
    "perGc",
    "obsExp"
)

read_cgi_anno <- function(x) {
    x %>%
        read_tsv(col_names = col_names) %>%
        dplyr::rename(
            gene_id = name,
            chr = chrom,
            start = chromStart,
            end = chromEnd
        ) %>%
        mutate(
            transcript_id = gene_id,
            strand = "*",
            symbol = gene_id
        )
}

download_parse_and_save <- function(genome_name, url) {
    temp_path <- tempfile()
    download.file(url, temp_path)

    anno_name <- paste0("inst/cgi_", genome_name, ".rds")
    saveRDS(read_cgi_anno(temp_path), anno_name, compress = "xz")

    fs::file_delete(temp_path)
}

# mm10 ----
download_parse_and_save(
    "mm10",
    "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/cpgIslandExt.txt.gz"
)

# GRCm39 ----
download_parse_and_save(
    "GRCm39",
    "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/cpgIslandExt.txt.gz"
)

# hg19 ----
download_parse_and_save(
    "hg19",
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz"
)

# hg38 ----
download_parse_and_save(
    "hg38",
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz"
)
