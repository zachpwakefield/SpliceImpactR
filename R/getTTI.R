getTTI <- function(paired_foreground, pdir = pdir, steps = 1, max_vertices_for_viz = 5000,
                   fdr = .05, plot_bool = T, ppidm_class = c("Gold", "Silver", "Bronze")[1],
                   write_igraphs = T, ddi = "Gold", ddi_type = "pdm", output_location = output_location) {
  # Create a directory for storing plots if plot_bool is TRUE
  if (plot_bool) {
    system(paste0("mkdir ", output_location, "tti"))
  }

  # Read the transcript-gene-protein mapping data for current genome
  tgp <- read.csv(paste0(pdir, '/gencodev42_transcriptGeneProtein.csv'))

  # Read edgelist output from initTTI -- using the ppidm class used previously
  edgeList <- read.table(paste0(output_location,  "tti_igraph_edgelist_", paste(ppidm_class, collapse = ""), "_removeDups"), sep = " ", row.names = NULL)

  # Convert the edgelist to a matrix format
  edgeList_Matrix <- matrix(c(edgeList$V1, edgeList$V2), ncol = 2)

  # Extract genes present in the control and test groups
  genes_in_sample <- unlist(lapply(strsplit(unique(unlist(lapply(c(control_group, test_group), function(f) {
    exon_df <- read.delim(paste0(f, '.exon'))
    unique(exon_df$gene)
  }))), split = "[.]"), "[[", 1))

  # Filter for transcripts in the sample
  possible_transcripts <- tgp$transcript_id[tgp$gene_id %in% genes_in_sample]

  # Reduce the edgelist to include only edges between possible transcripts
  el_reduce <- edgeList_Matrix[(edgeList_Matrix[,1] %in% possible_transcripts & edgeList_Matrix[,2] %in% possible_transcripts),]

  # Create an undirected graph from the reduced edgelist
  g <- igraph::graph_from_edgelist(el_reduce, directed = F)


  # Analyze differences between paired transcripts in the foreground group
  differences <- lapply(seq(1, (length(paired_foreground$transcript)-1), by = 2), function(tr) {
    # Check if both transcripts in the pair are in the graph
    if (sum(paired_foreground$transcript[c(tr, tr+1)] %in% igraph::V(g)$name) == 2) {
      # Generate ego graphs for the transcripts
      eg <- igraph::make_ego_graph(g, order = 1, nodes = paired_foreground$transcript[c(tr, tr+1)],
                                   mode = c("all", "out", "in")[1], mindist = 0)
      # Optionally write the ego graphs to files
      if (write_igraphs) {
        system(paste0("mkdir ", output_location, "tti/transcript_igraph_edgelists"))

        igraph::write_graph(eg[[1]], paste0(output_location, "tti/transcript_igraph_edgelists/", paired_foreground$gene[tr], paired_foreground$transcript[tr], "_igraph"),
                            format = "ncol")
        igraph::write_graph(eg[[2]], paste0(output_location, "tti/transcript_igraph_edgelists/", paired_foreground$gene[tr+1], paired_foreground$transcript[tr+1], "_igraph"),
                            format = "ncol")
      }

      # Identify vertices unique to each ego graph
      graph1_setDiff <- list(paired_foreground$transcript[tr], setdiff(igraph::V(eg[[1]])$name, igraph::V(eg[[2]])$name))
      graph2_setDiff <- list(paired_foreground$transcript[tr+1], setdiff(igraph::V(eg[[2]])$name, igraph::V(eg[[1]])$name))

      # Perform further analysis if there are unique vertices
      if (length(graph1_setDiff[[2]]) > 0 | length(graph2_setDiff[[2]]) > 0) {


        # Check whether difference between graphs is just the original seed/transcript -- if not, get TTI plots for output
        g1_check <- !(graph1_setDiff[[1]] %in% graph1_setDiff[[2]] & length(graph1_setDiff[[2]]) == 1)
        g2_check <- !(graph2_setDiff[[1]] %in% graph2_setDiff[[2]] & length(graph2_setDiff[[2]]) == 1)

        # Get iGraph plots for the transcripts
        if (g1_check | g1_check) {
          tti_igraph <- getTTIiGraphPlot(paired_foreground$transcript[c(tr, tr+1)], full_graph = g, steps = 1, max_vertices_for_viz = 5000, plot_bool = T)
        }

        # Perform enrichment analysis for unique vertices
        internal_loop <- lapply(list(graph1_setDiff, graph2_setDiff), function(x) {
          if (sum(paired_foreground$transcript[c(tr, tr+1)] %in% x[[2]]) > 0 & length(x[[2]]) == 1) {
            list(0, NA, NA)
          } else {
            if (length(x[[2]]) > 0) {
              list(x[[2]], getEnrichmentTTI(current_transcript = x[[1]], t_impacts = x[[2]], fdr = fdr, transGeneProt = tgp,
                                            backgroundGenes = genes_in_sample, steps = steps, plot_bool = plot_bool))
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
                        transcript1_setdiff = unlist(lapply(differences, function(yy) {length(yy[[1]])})),
                        transcript2_setdiff = unlist(lapply(differences, function(yy) {length(yy[[1]])}))
                        )
  # Return the list of differences
  return(list(differences = differences,
              results = results))
}



getTTIiGraphPlot <- function(paired_transcript, steps, full_graph, max_vertices_for_viz, plot_bool) {
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
  e1g <- igraph::simplify(e1g, remove.multiple = F, remove.loops = T)
  e2g <- igraph::simplify(e2g, remove.multiple = F, remove.loops = T)

  # Set edge curvature to FALSE for both graphs
  igraph::E(e1g)$curved <- F
  igraph::E(e2g)$curved <- F

  # Set default vertex and edge colors to 'azure4' for both graphs
  igraph::V(e1g)$color <- rep("azure4", length(igraph::V(e1g)))
  igraph::E(e1g)$color <- rep("azure4", length(igraph::E(e1g)))

  igraph::V(e2g)$color <- rep("azure4", length(igraph::V(e2g)))
  igraph::E(e2g)$color <- rep("azure4", length(igraph::E(e2g)))

  # Highlight unique vertices and edges in e1g not present in e2g
  if (length(setdiff(unique(igraph::V(e1g)), unique(igraph::V(e2g)))) > 0){

    igraph::V(e1g)$color[setdiff(unique(igraph::V(e1g)), unique(igraph::V(e2g)))] <- rep("red", length(setdiff(unique(igraph::V(e1g)), unique(igraph::V(e2g)))))
    igraph::E(e1g)$color[setdiff(unique(igraph::E(e1g)), unique(igraph::E(e2g)))] <- rep("brown", length(setdiff(unique(igraph::E(e1g)), unique(igraph::E(e2g)))))
    igraph::E(e1g)$color[igraph::E(e1g) %in% igraph::E(e1g)[from(igraph::V(e1g)[setdiff(unique(igraph::V(e1g)), unique(igraph::V(e2g)))])]] <- rep("brown", length(igraph::E(e1g)[from(igraph::V(e1g)[setdiff(unique(igraph::V(e1g)), unique(igraph::V(e2g)))])]))
  }

  # Highlight unique vertices and edges in e2g not present in e1g
  if (length(setdiff(unique(igraph::V(e2g)), unique(igraph::V(e1g)))) > 0){

    igraph::V(e2g)$color[setdiff(unique(igraph::V(e2g)), unique(igraph::V(e1g)))] <- rep("chartreuse4", length(setdiff(unique(igraph::V(e2g)), unique(igraph::V(e1g)))))
    igraph::E(e2g)$color[setdiff(unique(igraph::E(e2g)), unique(igraph::E(e1g)))] <- rep("blue", length(setdiff(unique(igraph::E(e2g)), unique(igraph::E(e1g)))))
    igraph::E(e2g)$color[igraph::E(e2g) %in% igraph::E(e2g)[from(igraph::V(e2g)[setdiff(unique(igraph::V(e2g)), unique(igraph::V(e1g)))])]] <- rep("blue", length(igraph::E(e2g)[from(igraph::V(e2g)[setdiff(unique(igraph::V(e2g)), unique(igraph::V(e1g)))])]))
  }

  # Highlight the nodes representing the paired transcripts in gold
  igraph::V(e1g)$color[V(e1g)$name == paired_transcript[1]] <- "gold"
  igraph::V(e2g)$color[V(e2g)$name == paired_transcript[2]] <- "gold"

  # Plot the ego graphs if plot_bool is TRUE and the number of vertices is within the specified limit
  if (plot_bool) {
    if (length(igraph::V(e1g)) <= max_vertices_for_viz) {
      pdf(paste0(output_location, "tti/", paired_transcript[1], "_", steps, 'steps_tti_graph.pdf'))
      print(igraph::plot.igraph(e1g,vertex.size=3,vertex.label=NA,main=paired_transcript[1],
                                layout=igraph::layout.fruchterman.reingold(e1g, niter=10000)))
      dev.off()
    } else {print("number of vertices exceeds max vertices for plotting set (max_vertices_for_viz)")}
    if (length(igraph::V(e2g)) <= max_vertices_for_viz) {
      pdf(paste0(output_location, "tti/", paired_transcript[2], "_", steps, 'steps_tti_graph.pdf'))
      print(igraph::plot.igraph(e2g,vertex.size=3,vertex.label=NA,main=paired_transcript[2],
                                layout=igraph::layout.fruchterman.reingold(e2g, niter=10000)))
      dev.off()
    } else {print("number of vertices exceeds max vertices for plotting set (max_vertices_for_viz)")}
  }

  # Return the list of ego graphs for the paired transcripts
  tti_graphs <- list(e1g, e2g)
  names(tti_graphs) <- paired_transcript
  return(tti_graphs)
}

getEnrichmentTTI <- function(current_transcript, t_impacts, fdr, transGeneProt,
                             backgroundGenes, steps, plot_bool) {

  # Extract gene names associated with the input transcript impacts
  enrichment_list <- transGeneProt$gene_name[transGeneProt$transcript_id %in% t_impacts]

  # Retrieve gene sets for different GO categories, KEGG, Biocarta, and Hallmark pathways using the hypeR package
  GO.cc <- hypeR::msigdb_gsets("Homo sapiens", "C5", "GO:CC", clean=TRUE)
  GO.mf <- hypeR::msigdb_gsets("Homo sapiens", "C5", "GO:MF", clean=TRUE)
  GO.bp <- hypeR::msigdb_gsets("Homo sapiens", "C5", "GO:BP", clean=TRUE)
  genesetsC2 <- hypeR::msigdb_gsets("Homo sapiens", "C2", "CP:KEGG", clean=TRUE)
  geneset_BIOCARTA <- hypeR::msigdb_gsets("Homo sapiens", "C2", "CP:BIOCARTA", clean=TRUE)
  genesetsH <- hypeR::msigdb_gsets("Homo sapiens", "H", clean=TRUE)


  # Perform enrichment analysis using hypergeometric test for each gene set category
  # And generate dot plots for visualizing the enrichment results
  cc_table <- hypeR::hypeR(enrichment_list, GO.cc, background = length(backgroundGenes), test="hypergeometric")
  cc_dots <- hypeR::hyp_dots(cc_table, fdr = fdr, title = "GO Cell Comp Enrichment", merge = TRUE)

  mf_table <- hypeR::hypeR(enrichment_list, GO.mf, background = length(backgroundGenes), test="hypergeometric")
  mf_dots <- hypeR::hyp_dots(mf_table, fdr = fdr, title = "GO Mol FXN Enrichment", merge = TRUE)

  bp_table <- hypeR::hypeR(enrichment_list, GO.bp, background = length(backgroundGenes), test="hypergeometric")
  bp_dots <- hypeR::hyp_dots(bp_table, fdr = fdr, title = "GO Biol Proc Enrichment", merge = TRUE)

  kg_table <- hypeR::hypeR(enrichment_list, genesetsC2, background = length(backgroundGenes), test="hypergeometric")
  kg_dots <- hypeR::hyp_dots(kg_table, fdr = fdr, title = "KEGG Enrichment", merge = TRUE)

  bc_table <- hypeR::hypeR(enrichment_list, geneset_BIOCARTA, background = length(backgroundGenes), test="hypergeometric")
  bc_dots <- hypeR::hyp_dots(bc_table, fdr = fdr, title = "Biocarta Enrichment", merge = TRUE)

  hm_table <- hypeR::hypeR(enrichment_list, genesetsH, background = length(backgroundGenes), test="hypergeometric")
  hm_dots <- hypeR::hyp_dots(hm_table, fdr = fdr, title = "Hallmark Enrichment", merge = TRUE)

  # If plot_bool is TRUE, save the plots to a PDF file
  if (plot_bool) {
    pdf(paste0(output_location, "tti/", current_transcript, '_unique_tti_', steps, '_steps_enrichment.pdf'))
    print(cc_dots)
    print(mf_dots)
    print(bp_dots)
    print(kg_dots)
    print(bc_dots)
    print(hm_dots)
    dev.off()
  }

  # Organize the enrichment tables and plots into a list for each category
  enrichment_tables_plots <- list("cellularComponent" = list(table = cc_table, dots = cc_dots),
                                  "molecularFunction" = list(table = mf_table, dots = mf_dots),
                                  "biologicalProcess" = list(table = bp_table, dots = bp_dots),
                                  "kegg" = list(table = kg_table, dots = kg_dots),
                                  "biocarta" = list(table = bc_table, dots = bc_dots),
                                  "hallmark" = list(table = hm_table, dots = hm_dots))

  # Return the list containing the enrichment tables and plots
  return(list("cellularComponent" = list(table = cc_table, dots = cc_dots),
              "molecularFunction" = list(table = mf_table, dots = mf_dots),
              "biologicalProcess" = list(table = bp_table, dots = bp_dots),
              "kegg" = list(table = kg_table, dots = kg_dots),
              "biocarta" = list(table = bc_table, dots = bc_dots),
              "hallmark" = list(table = hm_table, dots = hm_dots)))
}


# Initiate the TTI network using this function. This can be very time consuming
init_ddi <- function(pdir, output_location, ppidm_class = c("Gold", "Silver", "Bronze")[1], removeDups = T) {
  # Read in the protein coding data from the package directory
  pfam_in <- read.delim(paste0(pdir, '/protein_code_from_gencodev43_headerFix.txt.tsv'),
                        header = F)

  # Extract the transcript ID from the first column
  pfam_in$transcriptID <- unlist(lapply(strsplit(pfam_in$V1, split = "#"),
                                        "[[", 1))
  # Relocate the transcriptID column to the first position
  pfam_in <- pfam_in %>% dplyr::relocate(transcriptID)

  # Extract the gene ID from the first column
  pfam_in$geneID <- unlist(lapply(strsplit(unlist(lapply(strsplit(pfam_in$V1,
                                                                  split = "#"), "[[", 2)), split = ";"), "[[", 1))

  # Make a copy of pfam_in before filtering for Pfam domains only
  pfam_in_i <- pfam_in
  pfam_in <- pfam_in %>% dplyr::filter(V4 %in% c("Pfam"))

  # Select and deduplicate transcriptID and geneID pairs
  gt_df <- pfam_in_i %>% dplyr::select(transcriptID, geneID)
  gt_df <- gt_df[!duplicated(gt_df), ]

  # Create a list of transcripts associated with each unique domain (V5)
  d_transcripts <- lapply(unique(pfam_in$V5), function(x) {
    unique(pfam_in$transcriptID[pfam_in$V5 == x])
  })

  # Name the list elements with the unique domain names
  names(d_transcripts) <- unique(pfam_in$V5)

  # Remove duplicate list elements
  d_transcripts <- d_transcripts[!duplicated(d_transcripts)]


  # Read in the protein-protein interaction domain mapping (PPIDM) dataset
  pdm1 <- readr::read_csv(paste(pdir, "/PPIDM_FullSortedDataset_84K_GSB.csv",
                                sep = ""))

  # Filter interactions based on the specified PPIDM class (e.g., Gold, Silver, Bronze)
  d6 <- pdm1[pdm1$CLASS %in% ppidm_class, c(1, 2)]

  # Rename columns for clarity
  colnames(d6) <- c("n1", "n2")

  # Generate a dataframe of transcript pairs from the filtered PPIDM interactions
  tti <- do.call(rbind, parallel::mclapply(1:length(d6$n1), function(x) {
    # Check if both interacting proteins are in the domain-transcript list
    if (sum(names(d_transcripts) == d6$n1[x]) > 0 & sum(names(d_transcripts) ==
                                                        d6$n2[x]) > 0) {
      a <- d_transcripts[names(d_transcripts) == d6$n1[x]][[1]]
      b <- d_transcripts[names(d_transcripts) == d6$n2[x]][[1]]
      data.frame(tidyr::crossing(a,b)) # Create all combinations of transcripts for the interacting proteins
    }
  }, mc.cores = 20))

  # If removeDups == T (F would save time), sort out duplicated rows through sorting rows
  if (removeDups) {
    # Sort and deduplicate the transcript pairs [This can be time consuming]
    list_noDups <- mclapply(1:nrow(tti), mc.cores = 20, function(x) {
      sort(as.character(tti[x,]))
    })
    list_noDuplicates <- list_noDups[!duplicated(list_noDups)]
    tti <- data.frame(a = unlist(lapply(list_noDuplicates, "[[", 1)),
                             b = unlist(lapply(list_noDuplicates, "[[", 2)))
  }
  # Create an undirected graph from the deduplicated transcript pairs
  g <- igraph::graph_from_edgelist(as.matrix(tti), directed = F)

  # Write the graph edgelist to a file
  write_graph(
    g,
    paste0(output_location, 'tti_igraph_edgelist_', paste(ppidm_class, collapse = ""), '_removeDups'),
    format = "ncol"
  )

  return(g)
}
