#' gets the changes in tti interactions
#'
#' @param paired_foreground from getPaired
#' @param background proBed from getBackground
#' @param steps the number of steps for viz
#' @param max_vertices_for_viz max number of vertices to plot, saves time and space
#' @param fdr fdr to cut off for geneset enrichment
#' @param ppidm_class threshold of ppidm sets to use
#' @param output_location location to make output
#' @param tti_location if init_tti already performed, location of output or ""
#' @param tgp tgp_biomart from setup_gtf output
#' @param init_edgelist the edgelist out of initTTI if not using saved location
#' @param write_igraphs bool whether to write the graphs out or not (if large runs, takes up a lot of memory)
#' @param enrichHypeR bool set default to FALSE for enrichment of interactions 
#' with changing PPI -- need to install hypeR for this and current version broken
#' @return differences between each tti pair and the overall results
#' @importFrom igraph graph_from_edgelist V make_ego_graph write_graph simplify E layout.fruchterman.reingold
#' @importFrom tidyr crossing
#' @importFrom readr read_csv
#' @importFrom parallel mclapply
#' @importFrom dplyr select relocate
#'
#' @examples
#'
#' pdir <- system.file("extdata", package="SpliceImpactR")
#' dataDirectory <- paste0(pdir, "/")
#' test_group <- paste0(dataDirectory, "rawData/", c("test1","test2", "test3"))
#' control_group <- paste0(dataDirectory, "rawData/", c("control1", "control2", "control3"))
#' data_df <- data.frame(
#'     sample_names = c(control_group, test_group),
#'     phenotype_names = c(
#'       rep("control", length(control_group)),
#'       rep("test", length(test_group))
#'      ),
#'    stringsAsFactors = FALSE
#'   )
#' data_df$utc <- "control"
#' data_df$utc[data_df$phenotype_names == unique(data_df$phenotype_names)[2]] <- "test"
#'
#' transDF <- readr::read_csv(paste0(dataDirectory, "transcripts_limited_transDF.csv"))
#' c_trans <- readr::read_lines(paste0(dataDirectory, "transcripts_limited_c_trans.csv"))
#'
#' transcripts_sample <- list(transDF = transDF,
#'                            c_trans = c_trans)
#'
#' gtf_sample <- list(gtf = readr::read_csv(paste0(dataDirectory, "gtf_limited.csv")),
#'             transcript_gtf = readr::read_csv(paste0(dataDirectory, "transcript_gtf_limited.csv")),
#'             tgp_biomart = readr::read_csv(paste0(dataDirectory, "tgp_biomart_limited"))
#'             )
#' translations_sample <- readr::read_lines(paste0(dataDirectory, "translations_limited.csv"))
#'
#' ip <- readr::read_csv(paste0(dataDirectory, "biomart_ip.csv"))
#' code_regions <- readr::read_csv(paste0(dataDirectory, "biomart_code_regions.csv"))
#' pfam_exon_level <- readr::read_csv(paste0(dataDirectory, "biomart_pfam_exon_level.csv"))
#' fsd_exon_data <- readr::read_csv(paste0(dataDirectory, "biomart_data_sample.csv"))
#' pfam_data = readr::read_csv(paste0(dataDirectory, "biomart_pfam_exon.csv"))
#' biomart_data_sample <- list(ip = ip,
#'                      code_regions = code_regions,
#'                      fsd_exon_data = fsd_exon_data,
#'                      pfam_exon_level = pfam_exon_level,
#'                      pfam_data = pfam_data)
#'
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
#' fg <- getForeground(input = result,
#'                             test_names = test_group,
#'                             control_names = control_group,
#'                             thresh = .1,
#'                             fdr = .05,
#'                             mOverlap = .1,
#'                             exon_type = "AFE",
#'                             output_location = NULL,
#'                             cores = 1,
#'                             gtf = gtf_sample,
#'                             max_zero_prop = 1,
#'                             min_prop_samples = 0,
#'                             translations = translations_sample)
#'
#' bg <- getBackground(input=c(test_group, control_group),
#'                     mOverlap = 0.1,
#'                     cores = 1,
#'                     exon_type = "AFE",
#'                     output_location = NULL, gtf_sample, translations_sample)
#' library(msa)
#' pfg <- getPaired(foreground = fg$proBed,
#'           et = "AFE",
#'           nucleotides = transcripts_sample,
#'           newGTF = gtf_sample,
#'           cores = 1,
#'           output_location = NULL,
#'           saveAlignments = FALSE,
#'           exon_data = biomart_data_sample$fsd_exon_data)
#'
#' pfamData <- getPfam(background = bg,
#'                     foreground = fg,
#'                     pdir,
#'                     output_location = NULL,
#'                     cores = 1,
#'                     biomart_data_sample)
#'
#' domain_data <- getDomainData(fg,
#'                              bg,
#'                              pfg,
#'                              pfamData,
#'                              cores = 1,
#'                              output_location = NULL,
#'                              fdr_use = .25,
#'                              min_sample_success = 1,
#'                              engine = "Pfam",
#'                              repeatingDomains = FALSE,
#'                              topViz = 15)
#'
#' initDDI <- init_ddi(pdir = dataDirectory,
#'                     output_location = NULL,
#'                     ppidm_class = c("Gold_Standard", "Gold", "Silver", "Bronze")[1],
#'                     removeDups = TRUE,
#'                     cores = 1,
#'                     pfam_data = biomart_data_sample$pfam_data)
#'
#' tti <- getTTI(paired_foreground = pfg$paired_proBed,
#'               background = bg$proBed,
#'               steps = 1,
#'               max_vertices_for_viz = 300,
#'               fdr = .05,
#'               ppidm_class = c("Gold", "Silver", "Bronze")[1],
#'               write_igraphs = FALSE,
#'               output_location = NULL,
#'               tti_location = NULL,
#'               tgp = gtf_sample$tgp_biomart,
#'               init_edgelist = initDDI$edgelist,
#'               enrichHypeR = FALSE)
#' @export
getTTI <- function(paired_foreground,
                   background,
                   steps = 1,
                   max_vertices_for_viz = 5000,
                   fdr = .05,
                   ppidm_class = c("Gold", "Silver", "Bronze")[1],
                   write_igraphs = FALSE,
                   output_location = NULL,
                   tti_location = NULL,
                   tgp,
                   init_edgelist,
                   enrichHypeR = FALSE) {
  # Create a directory for storing plots if plot_bool is TRUE
  if (!is.null(output_location)) {
    system(paste0("mkdir ", output_location, "tti"))
    if (write_igraphs) {
      system(paste0("mkdir ", output_location, "tti/transcript_igraph_edgelists"))
    }
  }

  if (!is.null(tti_location)) {
    # Read edgelist output from initTTI -- using the ppidm class used previously
    edgeList <- read.table(paste0(tti_location,  "tti_igraph_edgelist_", paste(ppidm_class, collapse = ""), "_removeDups"),
                           sep = " ", row.names = NULL)
    # Convert the edgelist to a matrix format
    edgeList_Matrix <- matrix(c(edgeList$V1, edgeList$V2), ncol = 2)
  } else {
    edgeList_Matrix <- init_edgelist
  }

  # Filter for transcripts in the sample
  genes_in_sample <- unique(background$gene)
  possible_transcripts <- tgp$transcript_id[tgp$transcript_id %in% background$transcript]


  # Reduce the edgelist to include only edges between possible transcripts
  el_reduce <- edgeList_Matrix[(edgeList_Matrix[,1] %in% possible_transcripts & edgeList_Matrix[,2] %in% possible_transcripts),]

  # Create an undirected graph from the reduced edgelist
  g <- igraph::graph_from_edgelist(el_reduce, directed = FALSE)


  # Analyze differences between paired transcripts in the foreground group
  differences <- lapply(seq(1, (length(paired_foreground$transcript)-1), by = 2), function(tr) {
    # Check if both transcripts in the pair are in the graph
    if (sum(paired_foreground$transcript[c(tr, tr+1)] %in% igraph::V(g)$name) == 2) {
      # Generate ego graphs for the transcripts
      eg <- igraph::make_ego_graph(g, order = steps, nodes = paired_foreground$transcript[c(tr, tr+1)],
                                   mode = c("all", "out", "in")[1], mindist = 0)
      # Optionally write the ego graphs to files
      if (write_igraphs & !is.null(output_location)) {

        igraph::write_graph(eg[[1]], paste0(output_location, "tti/transcript_igraph_edgelists/",
                                            paired_foreground$gene[tr], "_",
                                            paired_foreground$transcript[tr], "_igraph"),
                            format = "ncol")
        igraph::write_graph(eg[[2]], paste0(output_location, "tti/transcript_igraph_edgelists/",
                                            paired_foreground$gene[tr+1], "_",
                                            paired_foreground$transcript[tr+1], "_igraph"),
                            format = "ncol")
      }

      # Identify vertices unique to each ego graph
      graph1_setDiff <- list(paired_foreground$transcript[tr], setdiff(igraph::V(eg[[1]])$name, igraph::V(eg[[2]])$name))
      graph2_setDiff <- list(paired_foreground$transcript[tr+1], setdiff(igraph::V(eg[[2]])$name, igraph::V(eg[[1]])$name))

      # Check whether difference between graphs is just the original seed/transcript -- if not, get TTI plots for output
      g1_check <- !(graph1_setDiff[[1]] %in% graph1_setDiff[[2]] & length(graph1_setDiff[[2]]) == 1)
      g2_check <- !(graph2_setDiff[[1]] %in% graph2_setDiff[[2]] & length(graph2_setDiff[[2]]) == 1)

      # Perform further analysis if there are unique vertices
      if ((length(graph1_setDiff[[2]]) > 0 | length(graph2_setDiff[[2]]) > 0) & (g1_check | g2_check)) {

        # Get iGraph plots for the transcripts
        # if (g1_check | g2_check) {
        tti_igraph <- getTTIiGraphPlot(paired_foreground$transcript[c(tr, tr+1)],
                                       gene = paired_foreground$gene[tr],
                                       full_graph = g,
                                       steps = steps,
                                       max_vertices_for_viz = max_vertices_for_viz,
                                       output_location = output_location)
        # }

        # Perform enrichment analysis for unique vertices
        internal_loop <- lapply(list(graph1_setDiff, graph2_setDiff), function(x) {
          if (sum(paired_foreground$transcript[c(tr, tr+1)] %in% x[[2]]) > 0 & length(x[[2]]) == 1) {
            list(0, NA, NA)
          } else {
            if (length(x[[2]]) > 0) {
              if (enrichHypeR) {
                if (!requireNamespace("hypeR", quietly = TRUE)) {
                  stop("Please install hypeR to use nbGLM.")
                }
                list(x[[2]], getEnrichmentTTI(current_transcript = x[[1]], t_impacts = x[[2]], fdr = fdr, transGeneProt = tgp,
                                              backgroundGenes = genes_in_sample, steps = steps,
                                              output_location = output_location))
              } else {list(x[[2]], NA)}
              
            } else {
              list(x[[2]], NA)
            }
          }

        })

        # Store the iGraph plots in the result
        internal_loop[[1]][[3]] <- tti_igraph[[1]]
        internal_loop[[2]][[3]] <- tti_igraph[[2]]
        internal_loop

      } else {
        # Return NA if there are no unique vertices
        list(list(0, NA, NA), list(0, NA, NA))
      }
      # Return -1 if one or both transcripts are not in the graph
    } else {list(list(-1, NA, NA), list(-1, NA, NA))}
  })

  # Name the list elements with the paired transcript and gene information
  names(differences) <- paste(paired_foreground$transcript[seq(1, length(paired_foreground$transcript), by=2)], ';',
                              paired_foreground$transcript[seq(2, length(paired_foreground$transcript), by=2)], ';',
                              paired_foreground$gene[seq(1, length(paired_foreground$transcript), by=2)], sep = "")


  results <- data.frame(transcript1 = paired_foreground$transcript[seq(1, length(paired_foreground$transcript), by=2)],
                        transcript2 = paired_foreground$transcript[seq(2, length(paired_foreground$transcript), by=2)],
                        gene = paired_foreground$gene[seq(1, length(paired_foreground$transcript), by=2)],
                        transcript1_setdiff = as.numeric(lapply(lapply(differences, function(yy) {
                          unlist(lapply(yy, function(yy2) {length(yy2[[1]])}))}), "[[", 1)),
                        transcript2_setdiff = as.numeric(lapply(lapply(differences, function(yy) {
                          unlist(lapply(yy, function(yy2) {length(yy2[[1]])}))}), "[[", 2))
                        )
  results <- results[results$transcript1_setdiff > 1 | results$transcript2_setdiff > 1,]
  if (!is.null(output_location)) {
    write_csv(results, paste0(output_location, "tti/tti_change_results.csv"))
  }
  # Return the list of differences
  return(list(differences = differences,
              results = results))
}

#' graph helper function
#' @return tti igraphs
#' @keywords internal
getTTIiGraphPlot <- function(paired_transcript, gene, steps, full_graph, max_vertices_for_viz, output_location) {
  # Create ego graphs for each of the paired transcripts
  eg <- igraph::make_ego_graph(
    full_graph,  # The full graph from which ego graphs are derived
    order = steps,  # The number of steps to consider in the ego graph
    nodes = paired_transcript,  # Nodes for which ego graphs are to be created
    mode = "all",  # Consider all directions of edges
    mindist = 0  # Include the node itself in the ego graph
  )


  # Extract ego graphs for each transcript
  e1g <- eg[[1]]  # Ego graph for the first transcript
  e2g <- eg[[2]]  # Ego graph for the second transcript

  # Simplify the ego graphs by removing loops and keeping multiple edges
  e1g <- igraph::simplify(e1g, remove.multiple = FALSE, remove.loops = TRUE)
  e2g <- igraph::simplify(e2g, remove.multiple = FALSE, remove.loops = TRUE)

  # Set edge curvature to FALSE for both graphs
  igraph::E(e1g)$curved <- FALSE
  igraph::E(e2g)$curved <- FALSE

  # Set default vertex and edge colors to 'azure4' for both graphs
  igraph::V(e1g)$color <- rep("azure4", length(igraph::V(e1g)))
  igraph::E(e1g)$color <- rep("azure4", length(igraph::E(e1g)))

  igraph::V(e2g)$color <- rep("azure4", length(igraph::V(e2g)))
  igraph::E(e2g)$color <- rep("azure4", length(igraph::E(e2g)))

  # Highlight unique vertices and edges in e1g not present in e2g
  if (length(setdiff(unique(igraph::V(e1g)), unique(igraph::V(e2g)))) > 0){

    igraph::V(e1g)$color[setdiff(unique(igraph::V(e1g)),
                                 unique(igraph::V(e2g)))] <- rep("red", length(setdiff(unique(igraph::V(e1g)),
                                                                                       unique(igraph::V(e2g)))))
    igraph::E(e1g)$color[setdiff(unique(igraph::E(e1g)),
                                 unique(igraph::E(e2g)))] <- rep("brown", length(setdiff(unique(igraph::E(e1g)),
                                                                                         unique(igraph::E(e2g)))))

    igraph::E(e1g)$color[igraph::E(e1g) %in% igraph::E(e1g)[from(igraph::V(e1g)[setdiff(unique(igraph::V(e1g)),
                                                                                        unique(igraph::V(e2g)))])]] <- rep("brown",
                                                                                                                           length(igraph::E(e1g)[from(igraph::V(e1g)[setdiff(unique(igraph::V(e1g)),
                                                                                                                                                                             unique(igraph::V(e2g)))])]))
  }

  # Highlight unique vertices and edges in e2g not present in e1g
  if (length(setdiff(unique(igraph::V(e2g)), unique(igraph::V(e1g)))) > 0){

    igraph::V(e2g)$color[setdiff(unique(igraph::V(e2g)),
                                 unique(igraph::V(e1g)))] <- rep("chartreuse4", length(setdiff(unique(igraph::V(e2g)),
                                                                                               unique(igraph::V(e1g)))))
    igraph::E(e2g)$color[setdiff(unique(igraph::E(e2g)),
                                 unique(igraph::E(e1g)))] <- rep("blue", length(setdiff(unique(igraph::E(e2g)),
                                                                                        unique(igraph::E(e1g)))))
    igraph::E(e2g)$color[igraph::E(e2g) %in% igraph::E(e2g)[from(igraph::V(e2g)[setdiff(unique(igraph::V(e2g)),
                                                                                        unique(igraph::V(e1g)))])]] <- rep("blue",
                                                                                                                           length(igraph::E(e2g)[from(igraph::V(e2g)[setdiff(unique(igraph::V(e2g)),
                                                                                                                                                                             unique(igraph::V(e1g)))])]))
  }

  # Highlight the nodes representing the paired transcripts in gold
  igraph::V(e1g)$color[V(e1g)$name == paired_transcript[1]] <- "gold"
  igraph::V(e2g)$color[V(e2g)$name == paired_transcript[2]] <- "gold"

  # Plot the ego graphs if plot_bool is TRUE and the number of vertices is within the specified limit
  if (!is.null(output_location)) {
    if (length(igraph::V(e1g)) <= max_vertices_for_viz) {
      pdf(paste0(output_location, "tti/", gene, "_", paired_transcript[1], "_", steps, 'steps_tti_graph.pdf'))
      print(igraph::plot.igraph(e1g,vertex.size=3,vertex.label=NA,main=paired_transcript[1],
                                layout=igraph::layout.fruchterman.reingold(e1g, niter=10000)))
      dev.off()
    } else {message("number of vertices exceeds max vertices for plotting set (max_vertices_for_viz)")}
    if (length(igraph::V(e2g)) <= max_vertices_for_viz) {
      pdf(paste0(output_location, "tti/", gene, "_", paired_transcript[2], "_", steps, 'steps_tti_graph.pdf'))
      print(igraph::plot.igraph(e2g,vertex.size=3,vertex.label=NA,main=paired_transcript[2],
                                layout=igraph::layout.fruchterman.reingold(e2g, niter=10000)))
      dev.off()
    } else {message("number of vertices exceeds max vertices for plotting set (max_vertices_for_viz)")}
  }

  # Return the list of ego graphs for the paired transcripts
  tti_graphs <- list(e1g, e2g)
  names(tti_graphs) <- paired_transcript
  return(tti_graphs)
}
#' enrichment helper function
#' @return geneset enrichment from hypeR
#' @keywords internal
getEnrichmentTTI <- function(current_transcript, t_impacts, fdr, transGeneProt,
                             backgroundGenes, steps, output_location) {
  if (!requireNamespace("hypeR", quietly = TRUE)) {
    stop("Please install hypeR to use nbGLM.")
  }

  # Extract gene names associated with the input transcript impacts
  enrichment_list <- transGeneProt$gene_name[transGeneProt$transcript_id %in% t_impacts]

  # Retrieve gene sets for different GO categories, KEGG, Biocarta, and Hallmark pathways using the hypeR package
  GO.cc <- hypeR::msigdb_gsets("Homo sapiens", "C5", "GO:CC", clean=TRUE)
  GO.mf <- hypeR::msigdb_gsets("Homo sapiens", "C5", "GO:MF", clean=TRUE)
  GO.bp <- hypeR::msigdb_gsets("Homo sapiens", "C5", "GO:BP", clean=TRUE)
  genesetsC2 <- hypeR::msigdb_gsets("Homo sapiens", "C2", "CP:KEGG", clean=TRUE)
  genesetsH <- hypeR::msigdb_gsets("Homo sapiens", "H", clean=TRUE)


  # Perform enrichment analysis using hypergeometric test for each gene set category
  # And generate dot plots for visualizing the enrichment results
  cc_table <- hypeR::hypeR(enrichment_list, GO.cc, background = length(backgroundGenes), test="hypergeometric")
  cc_dots <- hypeR::hyp_dots(cc_table, fdr = fdr, title = "GO Cell Comp Enrichment", merge = TRUE,
                             top = 10, abrv = 150)

  mf_table <- hypeR::hypeR(enrichment_list, GO.mf, background = length(backgroundGenes), test="hypergeometric")
  mf_dots <- hypeR::hyp_dots(mf_table, fdr = fdr, title = "GO Mol FXN Enrichment", merge = TRUE,
                             top = 10, abrv = 150)

  bp_table <- hypeR::hypeR(enrichment_list, GO.bp, background = length(backgroundGenes), test="hypergeometric")
  bp_dots <- hypeR::hyp_dots(bp_table, fdr = fdr, title = "GO Biol Proc Enrichment", merge = TRUE,
                             top = 10, abrv = 150)

  kg_table <- hypeR::hypeR(enrichment_list, genesetsC2, background = length(backgroundGenes), test="hypergeometric")
  kg_dots <- hypeR::hyp_dots(kg_table, fdr = fdr, title = "KEGG Enrichment", merge = TRUE,
                             top = 10, abrv = 150)


  hm_table <- hypeR::hypeR(enrichment_list, genesetsH, background = length(backgroundGenes), test="hypergeometric")
  hm_dots <- hypeR::hyp_dots(hm_table, fdr = fdr, title = "Hallmark Enrichment", merge = TRUE,
                             top = 10, abrv = 150)

  # If plot_bool is TRUE, save the plots to a PDF file
  if (!is.null(output_location)) {
    pdf(paste0(output_location, "tti/", current_transcript, '_unique_tti_', steps, '_steps_enrichment.pdf'))
    print(cc_dots)
    print(mf_dots)
    print(bp_dots)
    print(kg_dots)
    print(hm_dots)
    dev.off()
  }

  # Organize the enrichment tables and plots into a list for each category
  enrichment_tables_plots <- list("cellularComponent" = list(table = cc_table, dots = cc_dots),
                                  "molecularFunction" = list(table = mf_table, dots = mf_dots),
                                  "biologicalProcess" = list(table = bp_table, dots = bp_dots),
                                  "kegg" = list(table = kg_table, dots = kg_dots),
                                  "hallmark" = list(table = hm_table, dots = hm_dots))

  # Return the list containing the enrichment tables and plots
  return(list("cellularComponent" = list(table = cc_table, dots = cc_dots),
              "molecularFunction" = list(table = mf_table, dots = mf_dots),
              "biologicalProcess" = list(table = bp_table, dots = bp_dots),
              "kegg" = list(table = kg_table, dots = kg_dots),
              "hallmark" = list(table = hm_table, dots = hm_dots)))
}


#' initiates tti network for entire genome, can be very time consumng
#'
#' @param pdir the directory of the package
#' @param ppidm_class threshold of ppidm sets to use
#' @param output_location location to make output
#' @param cores number of requested cores
#' @param removeDups takes a longer time to init, but faster down the line and saves a lot of memory to remove duplcate vertices
#' @param pfam_data from biomart_data, pfam_data
#' @return overall tti network in both graph and edgelist form
#' @importFrom igraph graph_from_edgelist V make_ego_graph write_graph simplify E layout.fruchterman.reingold as_edgelist
#' @importFrom tidyr crossing
#' @importFrom readr read_csv
#' @importFrom parallel mclapply
#' @importFrom dplyr select relocate
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom data.table fread
#'
#' @examples
#'
#' pdir <- system.file("extdata", package="SpliceImpactR")
#' dataDirectory <- paste0(pdir, "/")
#'
#' ip <- readr::read_csv(paste0(dataDirectory, "biomart_ip.csv"))
#' code_regions <- readr::read_csv(paste0(dataDirectory, "biomart_code_regions.csv"))
#' pfam_exon_level <- readr::read_csv(paste0(dataDirectory, "biomart_pfam_exon_level.csv"))
#' fsd_exon_data <- readr::read_csv(paste0(dataDirectory, "biomart_data_sample.csv"))
#' pfam_data = readr::read_csv(paste0(dataDirectory, "biomart_pfam_exon.csv"))
#' biomart_data_sample <- list(ip = ip,
#'                      code_regions = code_regions,
#'                      fsd_exon_data = fsd_exon_data,
#'                      pfam_exon_level = pfam_exon_level,
#'                      pfam_data = pfam_data)
#'
#' initDDI <- init_ddi(pdir = dataDirectory,
#'                     output_location = NULL,
#'                     ppidm_class = c("Gold_Standard", "Gold", "Silver", "Bronze")[1],
#'                     removeDups = TRUE,
#'                     cores = 1,
#'                     pfam_data = biomart_data_sample$pfam_data)
#'
#' @export
#'

init_ddi <- function(pdir,
                     output_location = NULL,
                     ppidm_class = c("Gold_Standard", "Gold", "Silver", "Bronze")[1],
                     removeDups = TRUE,
                     cores = 1,
                     pfam_data) {
  pfam_data <- pfam_data[pfam_data$pfam != "",]
  pfam_in <- data.frame(transcriptID = pfam_data$ensembl_transcript_id,
                        geneID = pfam_data$ensembl_gene_id,
                        pfamID = pfam_data$pfam)

  pfam_in <- pfam_in[!is.na(pfam_in$transcriptID),]

  # Select and deduplicate transcriptID and geneID pairs
  gt_df <- pfam_in %>% dplyr::select(transcriptID, geneID)
  gt_df <- gt_df[!duplicated(gt_df), ]

  # Create a list of transcripts associated with each unique domain (V5)
  d_transcripts <- lapply(unique(pfam_in$pfam), function(x) {
    unique(pfam_in$transcriptID[pfam_in$pfam == x])
  })

  # Name the list elements with the unique domain names
  names(d_transcripts) <- unique(pfam_in$pfam)

  # Remove duplicate list elements
  d_transcripts <- d_transcripts[!duplicated(d_transcripts)]

  pdm1 <- data.frame(data.table::fread(paste0(pdir, "PPIDM_GoldDDIs.csv.gz")))

  # Filter interactions based on the specified PPIDM class (e.g., Gold, Silver, Bronze)
  if ("Gold_Standard" == ppidm_class) {
    d6 <- pdm1[pdm1$IN_GOLDSTANDARD == "yes",]
  } else if ("Gold_Standard" %in% ppidm_class) {
    d6 <- pdm1[pdm1$IN_GOLDSTANDARD == "yes",]
    d6 <- pdm1[pdm1$CLASS %in% ppidm_class, c(1, 2)]
  } else {
    d6 <- pdm1[pdm1$CLASS %in% ppidm_class, c(1, 2)]
  }

  # Rename columns for clarity
  colnames(d6) <- c("n1", "n2")

  # Generate a dataframe of transcript pairs from the filtered PPIDM interactions
  tti <- do.call(rbind, parallel::mclapply(seq_along(d6$n1), function(x) {
    # Check if both interacting proteins are in the domain-transcript list
    if (sum(names(d_transcripts) == d6$n1[x]) > 0 & sum(names(d_transcripts) ==
                                                        d6$n2[x]) > 0) {
      a <- d_transcripts[names(d_transcripts) == d6$n1[x]][[1]]
      b <- d_transcripts[names(d_transcripts) == d6$n2[x]][[1]]
      data.frame(tidyr::crossing(a,b)) # Create all combinations of transcripts for the interacting proteins
    }
  }, mc.cores = cores))

  if (removeDups) {
    # Sort and deduplicate the transcript pairs [This can be time consuming]
    list_noDups <- parallel::mclapply(seq_len(nrow(tti)), mc.cores = cores, function(x) {
      sort(as.character(tti[x,]))
    })
    list_noDuplicates <- list_noDups[!duplicated(list_noDups)]
    tti <- data.frame(a = unlist(lapply(list_noDuplicates, "[[", 1)),
                             b = unlist(lapply(list_noDuplicates, "[[", 2)))
  }
  # Create an undirected graph from the deduplicated transcript pairs
  g <- igraph::graph_from_edgelist(as.matrix(tti), directed = FALSE)
  edgelist <- igraph::as_edgelist(g)

  if (!is.null(output_location)) {
    write_graph(
      g,
      paste0(output_location, 'tti_igraph_edgelist_', paste(ppidm_class, collapse = ""), '_removeDups'),
      format = "ncol"
    )
  }

  return(list(graph = g,
              edgelist = edgelist))
}
