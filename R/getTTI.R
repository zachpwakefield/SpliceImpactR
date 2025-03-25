#' gets the tti interactions with various helper functions to do so
#'
#' @param paired_foreground from getPaired
#' @param background proBed from getBackground
#' @param pdir the directory of the package
#' @param steps the number of steps for viz
#' @param max_vertices_for_viz max number of vertices to plot, saves time and space
#' @param fdr fdr to cut off for geneset enrichment
#' @param plot_bool bool to make output plots of not
#' @param ppidm_class threshold of ppidm sets to use
#' @param write_graphs save full graph interactions for each with signfiicant changes
#' @param output_location location to make output
#' @param tti_location if init_tti already performed, location of output or ""
#' @param minOverlap minimum overlap to classify as matched to annotation
#' @param tgp tgp_biomart from setup_gtf output
#' @param cores number of requested cores
#' @return differences between each tti pair and the overall results
#' @importFrom igraph graph_from_edgelist V make_ego_graph write_graph simplify E layout.fruchterman.reingold
#' @importFrom tidyr crossing
#' @importFrom readr read_csv
#' @importFrom parallel mclapply
#' @importFrom dplyr select relocate
#' @importFrom hypeR msigdb_gsets hypeR hyp_dots
#' @export
getTTI <- function(paired_foreground, background, pdir, steps = 1, max_vertices_for_viz = 5000,
                   fdr = .05, plot_bool = TRUE, ppidm_class = c("Gold", "Silver", "Bronze")[1],
                   write_igraphs = TRUE,
                   output_location, tti_location, tgp) {
  # Create a directory for storing plots if plot_bool is TRUE
  if (plot_bool) {
    system(paste0("mkdir ", output_location, "tti"))
    if (write_igraphs) {
      system(paste0("mkdir ", output_location, "tti/transcript_igraph_edgelists"))
    }
  }

  # Read the transcript-gene-protein mapping data for current genome
  # tgp <- read.csv(paste0(pdir, '/gencodev42_transcriptGeneProtein.csv'))

  # Read edgelist output from initTTI -- using the ppidm class used previously
  edgeList <- read.table(paste0(tti_location,  "tti_igraph_edgelist_", paste(ppidm_class, collapse = ""), "_removeDups"),
                         sep = " ", row.names = NULL)

  # Convert the edgelist to a matrix format
  edgeList_Matrix <- matrix(c(edgeList$V1, edgeList$V2), ncol = 2)

  # # Extract genes present in the control and test groups
  # genes_in_sample <- unlist(lapply(strsplit(unique(unlist(lapply(c(control_group, test_group), function(f) {
  #   exon_df <- read.delim(paste0(f, '.exon'))
  #   unique(exon_df$gene)
  # }))), split = "[.]"), "[[", 1))

  # # Filter for transcripts in the sample
  # possible_transcripts <- tgp$transcript_id[tgp$gene_id %in% genes_in_sample]

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
      if (write_igraphs) {

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
                                       plot_bool = TRUE,
                                       output_location = output_location)
        # }

        # Perform enrichment analysis for unique vertices
        internal_loop <- lapply(list(graph1_setDiff, graph2_setDiff), function(x) {
          if (sum(paired_foreground$transcript[c(tr, tr+1)] %in% x[[2]]) > 0 & length(x[[2]]) == 1) {
            list(0, NA, NA)
          } else {
            if (length(x[[2]]) > 0) {
              list(x[[2]], getEnrichmentTTI(current_transcript = x[[1]], t_impacts = x[[2]], fdr = fdr, transGeneProt = tgp,
                                            backgroundGenes = genes_in_sample, steps = steps, plot_bool = plot_bool,
                                            output_location = output_location))
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
  write_csv(results, paste0(output_location, "tti/tti_change_results.csv"))
  # Return the list of differences
  return(list(differences = differences,
              results = results))
}

#' graph helper function
#' @return tti igraphs
#' @keywords internal
getTTIiGraphPlot <- function(paired_transcript, gene, steps, full_graph, max_vertices_for_viz, plot_bool, output_location) {
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
  if (plot_bool) {
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
                             backgroundGenes, steps, plot_bool, output_location) {

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
  if (plot_bool) {
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
#' @return overall tti network
#' @importFrom igraph graph_from_edgelist V make_ego_graph write_graph simplify E layout.fruchterman.reingold
#' @importFrom tidyr crossing
#' @importFrom readr read_csv
#' @importFrom parallel mclapply
#' @importFrom dplyr select relocate
#' @importFrom hypeR msigdb_gsets hypeR hyp_dots
#' @importFrom R.utils gunzip
#' @importFrom biomaRt useEnsembl getBM
#' @export
init_ddi <- function(pdir, output_location, ppidm_class = c("Gold_Standard", "Gold", "Silver", "Bronze")[1], removeDups = TRUE, cores = 1, pfam_data) {
  pfam_data <- pfam_data[pfam_data$pfam != "",]
  pfam_in <- data.frame(transcriptID = pfam_data$ensembl_transcript_id,
                        geneID = pfam_data$ensembl_gene_id,
                        pfamID = pfam_data$pfam)

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


  # Read in the protein-protein interaction domain mapping (PPIDM) dataset
  if (!(file.exists(paste0(pdir, "/PPIDM_GoldDDIs.csv")))) {
    R.utils::gunzip(paste0(pdir, "/PPIDM_GoldDDIs.csv.gz"), remove=FALSE, overwrite=TRUE)
  }

  pdm1 <- readr::read_delim(paste(pdir, "/PPIDM_GoldDDIs.csv",
                                sep = ""), delim = ";")

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
  tti <- do.call(rbind, parallel::mclapply(seq_alone(d6$n1), function(x) {
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

  # Write the graph edgelist to a file
  write_graph(
    g,
    paste0(output_location, 'tti_igraph_edgelist_', paste(ppidm_class, collapse = ""), '_removeDups'),
    format = "ncol"
  )

  return(g)
}
