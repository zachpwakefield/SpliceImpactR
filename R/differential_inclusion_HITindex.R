#' Extract differentially included AFE, ALE from HIT Index data
#'
#' @param test_names character vector of the path of test_names (eg:
#' '/path/to/test1', with the sample located '/path/to/test1.AFEPSI')
#' @param control_names character vector of the path of control_names (eg:
#' '/path/to/control1', with the sample located '/path/to/control1.AFEPSI')
#' @param et string of the event type (eg: 'AFE')
#' @param outlier_threshold the thresholding of the cooks distance,
#' no outlier removal is "Inf"
#' @param minReads threshold to count an exon as present
#' @param min_prop_samples minimum proportion of samples in
#' either phenotype an exon has to be present in to be kept
#' @param chosen_method stats model abbreviation
#' @param geneIndependent bool to (TRUE) remove non expressed genes
#' @return a dataframe with differential inclusion information
#' @importFrom data.table := data.table fread fifelse rbindlist uniqueN
#' @importFrom dplyr arrange filter group_by mutate ungroup
#' @importFrom stats cooks.distance lm logLik pchisq setNames start
#' @importFrom utils data read.csv read.delim read.table write.csv
#' @examples
#'
#' pdir <- system.file("extdata", package="SpliceImpactR")
#' dataDirectory <- paste0(pdir, "/rawData/")
#' test_group <- paste0(dataDirectory, c("test1","test2", "test3"))
#' control_group <- paste0(dataDirectory, c("control1", "control2", "control3"))
#'
#' result <- differential_inclusion_HITindex(test_names = test_group,
#'                                           control_names = control_group,
#'                                           et = "AFE",
#'                                           outlier_threshold = "Inf",
#'                                           minReads = 10,
#'                                           min_prop_samples = 0,
#'                                           chosen_method = "qbGLM"
#'                                           )
#'
#' @export
differential_inclusion_HITindex <- function(test_names,
                                            control_names,
                                            et,
                                            outlier_threshold = c("4/n", "1",
                                                                  "Inf")[3],
                                            minReads = 10,
                                            min_prop_samples = .5,
                                            chosen_method = 'qbGLM',
                                            geneIndependent = TRUE) {

    # Create sample type vector efficiently
    sample_types <- data.table::data.table(
        sample_name = c(test_names, control_names),
        type = rep(c("test", "control"),
                   c(length(test_names), length(control_names))))

    # Load PSI values for each sample and splicing event type
    psi_data_list <- lapply(paste0(sample_types$sample_name, ".", et, "PSI"),
                            data.table::fread)
    names(psi_data_list) <- sample_types$sample_name

    # Bind all data tables together
    psi_data <- data.table::rbindlist(psi_data_list, idcol = "sample_name")
    psi_data <- merge(psi_data, sample_types, by = "sample_name", all.x = TRUE)
    colnames(psi_data)[grep("PSI", colnames(psi_data))] <- 'psi'

    ## remove less than valid reads
    psi_data[, valid_reads := (nUP + nDOWN >= minReads)]
    psi_data <- psi_data[psi_data$valid_reads,]
    ## remove non terminal exons
    if (et == 'ALE' | et == 'HLE') {
        filtered_data <- psi_data %>%
            dplyr::filter(nLE != 0) %>%
            dplyr::group_by(gene, sample_name) %>%
            dplyr::mutate(psi = psi / sum(psi)) %>%
            dplyr::ungroup()
        psi_data <- data.table(filtered_data)
  } else {
        filtered_data <- psi_data %>%
            dplyr::filter(nFE != 0) %>%
            dplyr::group_by(gene, sample_name) %>%
            dplyr::mutate(psi = psi / sum(psi)) %>%
            dplyr::ungroup()
        psi_data <- data.table::data.table(filtered_data)
  }

    psi_data$id <- paste0(psi_data$gene, ";", psi_data$exon)

    keep_ids <- union(
        names(table(psi_data$id[psi_data$type == 'control'])[
            table(psi_data$id[psi_data$type == 'control']) >=
                min_prop_samples*sum(sample_types$type == 'control')]),
        names(table(psi_data$id[psi_data$type == 'test'])[
            table(psi_data$id[psi_data$type == 'test']) >=
                min_prop_samples*sum(sample_types$type == 'test')]))

    psi_data <- psi_data[psi_data$id %in% keep_ids,]


    ## add back in 0s
    if (et %in% c('HFE', 'AFE')) {
        all_gene_exons <- unique(psi_data[, .(gene, exon, strand, nFE)])
        # Ensure each sample has all gene/exon combinations that appear in  sample
        expanded_data <- all_gene_exons[, .(sample_name = sample_types$sample_name),
                                        by = .(gene, exon, strand, nFE)]
        expanded_data <- merge(expanded_data, psi_data,
                               by = c("sample_name", "gene", "exon",
                                      "strand", 'nFE'), all.x = TRUE)
        expanded_data <- merge(expanded_data[,-c('type')], sample_types,
                           by = "sample_name", all.x = TRUE)
  } else {
        all_gene_exons <- unique(psi_data[, .(gene, exon, strand, nLE)])
        # Ensure each sample has all gene/exon combinations that appear in  sample
        expanded_data <- all_gene_exons[, .(sample_name = sample_types$sample_name),
                                        by = .(gene, exon, strand, nLE)]
        expanded_data <- merge(expanded_data, psi_data,
                               by = c("sample_name", "gene", "exon",
                                      "strand", 'nLE'), all.x = TRUE)
        expanded_data <- merge(expanded_data[,-c('type')], sample_types,
                               by = "sample_name", all.x = TRUE)
  }
    expanded_data[is.na(psi), `:=` (
        valid_reads = TRUE,
        psi = 0,
        `sumL-R` = 0,
        nDOWN = 0,
        nUP = 0)]
    expanded_data[, nDiff := abs(nDOWN - nUP)]
    psi_data <- expanded_data

    ## option to add remove 0s for genes not expressed in samples
    if (geneIndependent) {
        psi_data[, has_nonzero_psi := any(psi > 0), by = .(sample_name, gene)]
        psi_data <- psi_data[psi_data$has_nonzero_psi]
        if (nrow(psi_data) == 0) return(NULL)
    }


    # modeling and Cook's distance for outlier detection if enabled
    psi_data[, valid_group := .N > 1 && data.table::uniqueN(type) > 1,
             by = .(gene, exon)]

    threshold <- switch(outlier_threshold,
                        "4/n" = 4 / nrow(sample_types),
                        "1" = 1,
                        as.numeric(outlier_threshold))

    psi_data <- psi_data[, total := sum(nDiff), by = .(sample_name, gene)]
    psi_data[, inclusion := nDiff]
    psi_data[, exclusion := total-inclusion]


    psi_data[, valid_group := .N > 1 && uniqueN(type) > 1, by = .(gene, exon)]
    psi_data <- psi_data[psi_data$valid_group]
    if (nrow(psi_data) == 0) return(NULL)

    psi_data <- getSignificance(psi_data, chosen_method)


    psi_data[, `:=` (
        delta.psi = mean(psi[type == "test" & cooks_d <= threshold]) -
            mean(psi[type == "control" & cooks_d <= threshold]),
        test_average_psi = mean(psi[type == "test" & cooks_d <= threshold]),
        control_average_psi = mean(psi[type == "control" & cooks_d <= threshold])
    ), by = .(gene, exon)]

    # Add zero_test, count_test, and count_control columns and finalize cols
    psi_data[, zero_count := sum(psi == 0), by = .(gene, exon)]
    final_data <- psi_data[, .(
      gene, exon, p.val, delta.psi, test_average_psi, control_average_psi,
      count_test = sum(type == "test"), count_control = sum(type == "control"),
      zero_count
    ), by = .(gene, exon)]

    final_data <- final_data[!duplicated(final_data),-c(1, 2)] %>%
        dplyr::arrange(.data$gene)

    # Adjust p-values for multiple testing
    final_data$p.adj <- p.adjust(final_data$p.val, method = "fdr")
    final_data$p.adj[is.na(final_data$p.adj)] <- 1
    final_data$p.val[is.na(final_data$p.val)] <- 1

    ## Annotate
    final_data$add_inf <- "none"
    final_data$type <- et

    return(data.frame(final_data))

}
