#!/usr/bin/env Rscript
#!/usr/bin/env Rscript

if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    install.packages(
        "SeuratObject",
        repos = "https://cloud.r-project.org",
        Ncpus = 42
    )
}

#!/usr/bin/env Rscript

if (!requireNamespace("optparse", quietly = TRUE)) {
    install.packages(
        "optparse",
        repos = "https://cloud.r-project.org",
        Ncpus = 42
    )
}

library(optparse)
library(SeuratObject)

###function for reading cellranger output -> implemented in Seurat
ReadMtx <- function (mtx,
                     cells,
                     features,
                     cell.column = 1,
                     feature.column = 2,
                     cell.sep = "\t",
                     feature.sep = "\t",
                     skip.cell = 0,
                     skip.feature = 0,
                     mtx.transpose = FALSE,
                     unique.features = TRUE,
                     strip.suffix = FALSE){
    
    all.files <- list(
        `expression matrix` = mtx,
        `barcode list` = cells,
        `feature list` = features
    )
    for (i in seq_along(along.with = all.files)) {
        uri <- tryCatch(
            expr = {
                con <- url(description = all.files[[i]])
                close(con = con)
                all.files[[i]]
            },
            error = function(...) {
                return(normalizePath(path = all.files[[i]], winslash = "/"))
            }
        )
        err <- paste("Cannot find", names(x = all.files)[i], "at", uri)
        uri <- httr::build_url(url = httr::parse_url(url = uri))
        if (grepl(pattern = "^[A-Z]?:///", x = uri)) {
            uri <- gsub(pattern = "^://",
                        replacement = "",
                        x = uri)
            if (!file.exists(uri)) {
                stop(err, call. = FALSE)
            }
        }
        else {
            if (!Online(url = uri, seconds = 2L)) {
                stop(err, call. = FALSE)
            }
            if (file_ext(uri) == "gz") {
                con <- url(description = uri)
                uri <- gzcon(con = con, text = TRUE)
            }
        }
        all.files[[i]] <- uri
    }
    cell.barcodes <- read.table(
        file = all.files[["barcode list"]],
        header = FALSE,
        sep = cell.sep,
        row.names = NULL,
        skip = skip.cell
    )
    feature.names <- read.table(
        file = all.files[["feature list"]],
        header = FALSE,
        sep = feature.sep,
        row.names = NULL,
        skip = skip.feature
    )
    bcols <- ncol(x = cell.barcodes)
    if (bcols < cell.column) {
        stop(
            "cell.column was set to ",
            cell.column,
            " but ",
            cells,
            " only has ",
            bcols,
            " columns.",
            " Try setting the cell.column argument to a value <= to ",
            bcols,
            "."
        )
    }
    cell.names <- cell.barcodes[, cell.column]
    if (all(grepl(pattern = "\\-1$", x = cell.names)) &
        strip.suffix) {
        cell.names <- as.vector(x = as.character(
            x = sapply(
                X = cell.names,
                FUN = ExtractField,
                field = 1,
                delim = "-"
            )
        ))
    }
    fcols <- ncol(x = feature.names)
    if (fcols < feature.column) {
        stop(
            "feature.column was set to ",
            feature.column,
            " but ",
            features,
            " only has ",
            fcols,
            " column(s).",
            " Try setting the feature.column argument to a value <= to ",
            fcols,
            "."
        )
    }
    if (any(is.na(x = feature.names[, feature.column]))) {
        na.features <- which(x = is.na(x = feature.names[, feature.column]))
        replacement.column <- ifelse(test = feature.column ==
                                         2,
                                     yes = 1,
                                     no = 2)
        if (replacement.column > fcols) {
            stop(
                "Some features names are NA in column ",
                feature.column,
                ". Try specifiying a different column.",
                call. = FALSE
            )
        }
        else {
            warning(
                "Some features names are NA in column ",
                feature.column,
                ". Replacing NA names with ID from column ",
                replacement.column,
                ".",
                call. = FALSE
            )
        }
        feature.names[na.features, feature.column] <- feature.names[na.features, replacement.column]
    }
    feature.names <- feature.names[, feature.column]
    if (unique.features) {
        feature.names <- make.unique(names = feature.names)
    }
    data <- Matrix::readMM(file = all.files[["expression matrix"]])
    if (mtx.transpose) {
        data <- t(x = data)
    }
    if (length(x = cell.names) != ncol(x = data)) {
        stop(
            "Matrix has ",
            ncol(data),
            " columns but found ",
            length(cell.names),
            " barcodes. ",
            ifelse(
                test = length(x = cell.names) >
                    ncol(x = data),
                yes = "Try increasing `skip.cell`. ",
                no = ""
            ),
            call. = FALSE
        )
    }
    if (length(x = feature.names) != nrow(x = data)) {
        stop(
            "Matrix has ",
            nrow(data),
            " rows but found ",
            length(feature.names),
            " features. ",
            ifelse(
                test = length(x = feature.names) >
                    nrow(x = data),
                yes = "Try increasing `skip.feature`. ",
                no = ""
            ),
            call. = FALSE
        )
    }
    colnames(x = data) <- cell.names
    rownames(x = data) <- feature.names
    data <- as.sparse(x = data)
    return(data)
}

option_list <- list(
    make_option("--solo", type = "character"),
    make_option("--metadata", type = "character"),
    make_option("--output", type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

filtered_dir <- file.path(opt$solo, "Gene", "filtered")

if (!dir.exists(filtered_dir)) {
    stop("STARsolo filtered directory does not exist: ", filtered_dir)
}

barcodes_file <- file.path(filtered_dir, "barcodes.tsv")
features_file <- file.path(filtered_dir, "features.tsv")
matrix_file <- file.path(filtered_dir, "matrix.mtx")

if (!file.exists(barcodes_file)) {
    stop("Missing barcodes.tsv: ", barcodes_file)
}

if (!file.exists(features_file)) {
    stop("Missing features.tsv: ", features_file)
}

if (!file.exists(matrix_file)) {
    stop("Missing matrix.mtx: ", matrix_file)
}

message("Reading STARsolo filtered matrix...")

counts <- ReadMtx(
    mtx = matrix_file,
    features = features_file,
    cells = barcodes_file,
    feature.column = 2
)

message("Creating Seurat object...")

obj <- CreateSeuratObject(
    counts = counts
)

message("Reading Xescape cell-level metadata...")

metadata <- read.delim(
    opt$metadata,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

if (!"Barcode" %in% colnames(metadata)) {
    stop("Metadata file does not contain a 'Barcode' column")
}

if (anyDuplicated(metadata$Barcode)) {
    stop("Duplicate barcodes detected in Xescape metadata")
}

rownames(metadata) <- metadata$Barcode
metadata$Barcode <- NULL

cell_barcodes <- colnames(obj)

missing_metadata <- setdiff(cell_barcodes, rownames(metadata))

if (length(missing_metadata) > 0) {
    warning(
        length(missing_metadata),
        " STARsolo cells have no corresponding Xescape metadata"
    )
}

extra_metadata <- setdiff(rownames(metadata), cell_barcodes)

if (length(extra_metadata) > 0) {
    message(
        length(extra_metadata),
        " metadata barcodes are not present in the STARsolo filtered matrix"
    )
}

metadata <- metadata[
    intersect(cell_barcodes, rownames(metadata)),
    ,
    drop = FALSE
]

metadata <- metadata[
    match(intersect(cell_barcodes, rownames(metadata)), rownames(metadata)),
    ,
    drop = FALSE
]

metadata <- metadata[match(
    cell_barcodes,
    rownames(metadata)
), , drop = FALSE]

rownames(metadata) <- cell_barcodes

message("Adding metadata...")

obj <- AddMetaData(
    object = obj,
    metadata = metadata
)

#Set the correct sample name from given output name
obj$orig.ident <- sub("_Xes.*$", "", opt$output)

message("Cells in Seurat object: ", ncol(obj))
message("Metadata columns: ", ncol(obj[[]]))

message("Saving RDS: ", opt$output)

saveRDS(
    obj,
    file = opt$output
)

message("Done.")
