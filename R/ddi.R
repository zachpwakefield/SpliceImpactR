getTTI <- function(paired_foreground, pdir = pdir, steps = 1, max_vertices_for_viz = 5000, fdr = .05, plot_bool = T, write_igraphs = T, ddi = "Gold", ddi_type = "pdm", output_location = output_location) {
  if (plot_bool) {
    system(paste0("mkdir ", output_location, "tti"))
  }
  tgp <- read.csv(paste0(pdir, '/gencodev42_transcriptGeneProtein.csv'))

  edgeList <- do.call(rbind, lapply(list.files('/projectnb2/evolution/zwakefield/proteinImpacts/graph_goldStandard_brokendown'), function(x) {
    read.table(paste0('/projectnb2/evolution/zwakefield/proteinImpacts/graph_goldStandard_brokendown/', x), sep = " ", row.names = NULL)
  }))
  edgeList_Matrix <- matrix(c(edgeList$V1, edgeList$V2), ncol = 2)


  genes_in_sample <- unlist(lapply(strsplit(unique(unlist(lapply(c(control_group, test_group), function(f) {
    exon_df <- read.delim(paste0(f, '.exon'))
    unique(exon_df$gene)
  }))), split = "[.]"), "[[", 1))
  possible_transcripts <- tgp$transcript_id[tgp$gene_id %in% genes_in_sample]

  el_reduce <- edgeList_Matrix[(edgeList_Matrix[,1] %in% possible_transcripts & edgeList_Matrix[,2] %in% possible_transcripts),]
  g <- igraph::graph_from_edgelist(el_reduce, directed = F)

  differences <- lapply(seq(1, (length(paired_foreground$transcript)-1), by = 2), function(tr) {
    if (sum(paired_foreground$transcript[c(tr, tr+1)] %in% igraph::V(g)$name) == 2) {
      eg <- igraph::make_ego_graph(g, order = 1, nodes = paired_foreground$transcript[c(tr, tr+1)],
                                   mode = c("all", "out", "in")[1], mindist = 0)
      if (write_graphs) {
        system(paste0("mkdir ", output_location, "tti/transcript_igraph_edgelists"))

        igraph::write_graph(eg[[1]], paste0(output_location, "tti/transcript_igraph_edgelists/", paired_foreground$transcript[tr], "_igraph"),
                            format = "ncol")
        igraph::write_graph(eg[[2]], paste0(output_location, "tti/transcript_igraph_edgelists/", paired_foreground$transcript[tr+1], "_igraph"),
                            format = "ncol")
      }
      graph1_setDiff <- list(paired_foreground$transcript[tr], setdiff(igraph::V(eg[[1]])$name, igraph::V(eg[[2]])$name))
      graph2_setDiff <- list(paired_foreground$transcript[tr+1], setdiff(igraph::V(eg[[2]])$name, igraph::V(eg[[1]])$name))

      if (length(graph1_setDiff) > 0 | length(graph2_setDiff) > 0) {

        tti_igraph <- getTTIiGraphPlot(ttt[c(tr, tr+1)], full_graph = g, steps = 1, max_vertices_for_viz = 5000, plot_bool = T)

        internal_loop <- lapply(list(graph1_setDiff, graph2_setDiff), function(x) {
          if (sum(ttt[c(tr, tr+1)] %in% x[[2]]) > 0 & length(x[[2]]) == 1) {
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

        internal_loop[[1]][[3]] <- tti_igraph[[1]]
        internal_loop[[2]][[3]] <- tti_igraph[[2]]
        internal_loop

      } else {
        list(list(0, NA, NA), list(0, NA, NA))
      }
    } else {list(list(-1, NA, NA), list(-1, NA, NA))}
  })


  names(differences) <- paste(paired_foreground$transcript[seq(1, length(paired_foreground$transcript), by=2)], ';',
                              paired_foreground$transcript[seq(2, length(paired_foreground$transcript), by=2)], ';',
                              paired_foreground$gene[seq(1, length(paired_foreground$transcript), by=2)], sep = "")
  return(differences)
}



getTTIiGraphPlot <- function(paired_transcript, steps, full_graph, max_vertices_for_viz, plot_bool) {
  eg <- igraph::make_ego_graph(
    full_graph,
    order = steps,
    nodes = paired_transcript,
    mode = c("all", "out", "in")[1],
    mindist = 0
  )

  e1g <- eg[[1]]
  e2g <- eg[[2]]

  e1g <- igraph::simplify(e1g, remove.multiple = F, remove.loops = T)
  e2g <- igraph::simplify(e2g, remove.multiple = F, remove.loops = T)

  igraph::E(e1g)$curved <- F
  igraph::E(e2g)$curved <- F

  igraph::V(e1g)$color <- rep("azure4", length(igraph::V(e1g)))
  igraph::E(e1g)$color <- rep("azure4", length(igraph::E(e1g)))


  igraph::V(e2g)$color <- rep("azure4", length(igraph::V(e2g)))
  igraph::E(e2g)$color <- rep("azure4", length(igraph::E(e2g)))


  if (length(setdiff(unique(igraph::V(e1g)), unique(igraph::V(e2g)))) > 0){

    igraph::V(e1g)$color[setdiff(unique(igraph::V(e1g)), unique(igraph::V(e2g)))] <- rep("red", length(setdiff(unique(igraph::V(e1g)), unique(igraph::V(e2g)))))
    igraph::E(e1g)$color[setdiff(unique(igraph::E(e1g)), unique(igraph::E(e2g)))] <- rep("brown", length(setdiff(unique(igraph::E(e1g)), unique(igraph::E(e2g)))))
    igraph::E(e1g)$color[igraph::E(e1g) %in% igraph::E(e1g)[from(igraph::V(e1g)[setdiff(unique(igraph::V(e1g)), unique(igraph::V(e2g)))])]] <- rep("brown", length(igraph::E(e1g)[from(igraph::V(e1g)[setdiff(unique(igraph::V(e1g)), unique(igraph::V(e2g)))])]))
  }
  if (length(setdiff(unique(igraph::V(e2g)), unique(igraph::V(e1g)))) > 0){

    igraph::V(e2g)$color[setdiff(unique(igraph::V(e2g)), unique(igraph::V(e1g)))] <- rep("chartreuse4", length(setdiff(unique(igraph::V(e2g)), unique(igraph::V(e1g)))))
    igraph::E(e2g)$color[setdiff(unique(igraph::E(e2g)), unique(igraph::E(e1g)))] <- rep("blue", length(setdiff(unique(igraph::E(e2g)), unique(igraph::E(e1g)))))
    igraph::E(e2g)$color[igraph::E(e2g) %in% igraph::E(e2g)[from(igraph::V(e2g)[setdiff(unique(igraph::V(e2g)), unique(igraph::V(e1g)))])]] <- rep("blue", length(igraph::E(e2g)[from(igraph::V(e2g)[setdiff(unique(igraph::V(e2g)), unique(igraph::V(e1g)))])]))
  }

  igraph::V(e1g)$color[V(e1g)$name == paired_transcript[1]] <- "gold"
    igraph::V(e2g)$color[V(e2g)$name == paired_transcript[2]] <- "gold"

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
    tti_graphs <- list(e1g, e2g)
    names(tti_graphs) <- paired_transcript
    return(tti_graphs)
}

getEnrichmentTTI <- function(current_transcript, t_impacts, fdr, transGeneProt,
                             backgroundGenes, steps, plot_bool) {

  enrichment_list <- transGeneProt$gene_name[transGeneProt$transcript_id %in% t_impacts]


  GO.cc <- hypeR::msigdb_gsets("Homo sapiens", "C5", "GO:CC", clean=TRUE)
  GO.mf <- hypeR::msigdb_gsets("Homo sapiens", "C5", "GO:MF", clean=TRUE)
  GO.bp <- hypeR::msigdb_gsets("Homo sapiens", "C5", "GO:BP", clean=TRUE)
  genesetsC2 <- hypeR::msigdb_gsets("Homo sapiens", "C2", "CP:KEGG", clean=TRUE)
  geneset_BIOCARTA <- hypeR::msigdb_gsets("Homo sapiens", "C2", "CP:BIOCARTA", clean=TRUE)
  genesetsH <- hypeR::msigdb_gsets("Homo sapiens", "H", clean=TRUE)

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

  if (plot_bool == T) {
    pdf(paste0(output_location, "tti/", current_transcript, '_unique_tti_', steps, '_steps_enrichment.pdf'))
    print(cc_dots)
    print(mf_dots)
    print(bp_dots)
    print(kg_dots)
    print(bc_dots)
    print(hm_dots)
    dev.off()
  }
  enrichment_tables_plots <- list("cellularComponent" = list(table = cc_table, dots = cc_dots),
                                  "molecularFunction" = list(table = mf_table, dots = mf_dots),
                                  "biologicalProcess" = list(table = bp_table, dots = bp_dots),
                                  "kegg" = list(table = kg_table, dots = kg_dots),
                                  "biocarta" = list(table = bc_table, dots = bc_dots),
                                  "hallmark" = list(table = hm_table, dots = hm_dots))
  return(list("cellularComponent" = list(table = cc_table, dots = cc_dots),
              "molecularFunction" = list(table = mf_table, dots = mf_dots),
              "biologicalProcess" = list(table = bp_table, dots = bp_dots),
              "kegg" = list(table = kg_table, dots = kg_dots),
              "biocarta" = list(table = bc_table, dots = bc_dots),
              "hallmark" = list(table = hm_table, dots = hm_dots)))
}
