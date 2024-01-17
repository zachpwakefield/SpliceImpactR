tti <- function(location, steps = 2,
                max_vertices_for_viz = 1000,
                fdr = .05, plot_bool = T,
                ddi = c("Gold", "Silver", "Bronze")[1],
                ddi_type = c("pdm", "3did")[1],
                output_location) {
  system(paste0("mkdir ", output_location, "tti_output"))
  gtt <- readRDS(paste0(location,'/hg38_geneRef_conv.RDS'))
  colnames(gtt) <- c("geneID", "geneName")
  gttf <- gtt[!duplicated(gtt),]

  GO.cc <- hypeR::msigdb_gsets("Homo sapiens", "C5", "GO:CC", clean=TRUE)
  GO.mf <- hypeR::msigdb_gsets("Homo sapiens", "C5", "GO:MF", clean=TRUE)
  GO.bp <- hypeR::msigdb_gsets("Homo sapiens", "C5", "GO:BP", clean=TRUE)
  genesetsC2 <- hypeR::msigdb_gsets("Homo sapiens", "C2", "CP:KEGG", clean=TRUE)
  geneset_BIOCARTA <- hypeR::msigdb_gsets("Homo sapiens", "C2", "CP:BIOCARTA", clean=TRUE)
  genesetsH <- hypeR::msigdb_gsets("Homo sapiens", "H", clean=TRUE)

  i_fgtsv <- read_csv(paste0(output_location, 'paired_fgoutBed.csv'))
  i_fgtsv <- i_fgtsv[i_fgtsv$prot != "none",]
  t_fgtsv <- i_fgtsv
  t_fgtsv$rn <- 1:length(t_fgtsv$transcript)
  keep <- list()
  for (j in 1:(length(i_fgtsv$transcript)-1)) {
    if (i_fgtsv$gene[j] == i_fgtsv$gene[j+1]) {
      keep[[length(keep)+1]] <- c(j, j+1)
      t_fgtsv <- t_fgtsv[!(t_fgtsv$rn %in% c(j, j+1)),]
    } else {
      t_fgtsv <- t_fgtsv[!(t_fgtsv$rn %in% c(j)),]
    }

  }
  fgtsv <- i_fgtsv[unlist(keep),]


  bgtsv <- read.delim(paste0(output_location, 'bgoutFast.fa.tsv'), header  = F)
  bgtsv$transcriptID <- unlist(lapply(strsplit(bgtsv$V1, split = "#"), "[[", 1))
  bgtsv <- bgtsv %>% dplyr::relocate(transcriptID)
  bgtsv$geneID <- unlist(lapply(strsplit(unlist(lapply(strsplit(bgtsv$V1, split = "#"), "[[", 2)), split = ";"), "[[", 1))
  bgtsv_i <- bgtsv
  bgtsv <- bgtsv %>% dplyr::filter(V4 %in% c("Pfam"))

  gt_df <- bgtsv_i %>% dplyr::select(transcriptID, geneID)
  gt_df <- gt_df[!duplicated(gt_df),]


  t_domains <- lapply(unique(bgtsv$transcriptID), function(x) {
    unique(bgtsv$V5[bgtsv$transcriptID == x])
  })

  names(t_domains) <- unique(bgtsv$transcriptID)

  d_transcripts <- lapply(unique(bgtsv$V5), function(x) {
    unique(bgtsv$transcriptID[bgtsv$V5 == x])
  })
  names(d_transcripts) <- unique(bgtsv$V5)
  d_transcripts <- d_transcripts[!duplicated(d_transcripts)]

  if (ddi_type == "3did") {


    d1 <- read_lines('/projectnb2/evolution/zwakefield/proteinImpacts/3did_flat')
    d2 <- gsub("\t", " ", d1[grepl("#=ID", d1) | grepl("#=3D", d1)])

    iter <- c(grep("#=ID", d2), length(d2))
    strong <- unlist(lapply(1:(length(iter)-1), function(x) {
      i2 <- d2[(iter[x]+1):(iter[x+1]-1)]
      i3 <- as.numeric(unlist(lapply(strsplit(i2, split = " "), "[[", 5)))
      if (sum(i3 >= 2.3) >= 1) {
        iter[x]
      }
    }))

    d2 <- d2[strong]
    d3 <- lapply(d2, function(x) gsub("@Pfam", "", gsub("[)]", "", gsub("[(]", "", strsplit(x, split = " ")[[1]][c(2, 3, 5, 6)]))))
    d4 <- lapply(d3, function(x) data.frame(t(data.frame(x))))
    d5 <- do.call(rbind, d4)
    rownames(d5) <- NULL
    d6 <- data.frame(n1 = unlist(lapply(strsplit(d5$X3, split = "[.]"), "[[", 1)),
                     n2 = unlist(lapply(strsplit(d5$X4, split = "[.]"), "[[", 1)))
  } else {
    pdm1 <- read_csv(paste(location,"/PPIDM_FullSortedDataset_84K_GSB.csv",sep=""))
    d6 <- pdm1[pdm1$CLASS == ddi,c(1, 2)]
    d6 <- pdm1[pdm1$IN_GOLDSTANDARD == "yes",c(1, 2)]
    colnames(d6) <-c("n1", "n2")
  }

  t_d_df <- do.call(rbind, lapply(1:length(d6$n1), function(x) {
    if (sum(names(d_transcripts) == d6$n1[x]) > 0 &
        sum(names(d_transcripts) == d6$n2[x]) > 0) {
      a <- d_transcripts[names(d_transcripts) == d6$n1[x]][[1]]
      b <- d_transcripts[names(d_transcripts) == d6$n2[x]][[1]]

      comb <- unlist(lapply(unique(a), function(x) {
        lapply(unique(b), function(y) {
          sort(c(x, y))
        })
      }), recursive = F)
      comb <- comb[!duplicated(comb)]

      df_trans <- data.frame(n1 = unlist(lapply(comb, "[[", 1)),
                             n2 = unlist(lapply(comb, "[[", 2)))
      df_trans$d1 <- d6$n1[x]
      df_trans$d2 <- d6$n2[x]
      df_trans[!duplicated(df_trans),]
    }

  }))
  if (file.exists(paste0(output_location, "tti_output/adjacency_stats.txt"))) {
    system(paste0("rm -f ", output_location, "tti_output/adjacency_stats.txt"))
  }
  domainDiff <- lapply(seq(1, length(fgtsv$transcript)-1, by = 2), function(x) {
    t1 <- fgtsv$transcript[(x)]
    t2 <- fgtsv$transcript[(x+1)]
    t1_con <- unique(c(t_d_df$n2[t_d_df$n1 == t1], t_d_df$n1[t_d_df$n2 == t1]))
    t2_con <- unique(c(t_d_df$n2[t_d_df$n1 == t2], t_d_df$n1[t_d_df$n2 == t2]))
    if (
      (length(setdiff(t1_con,
                      t2_con)) != 0 | length(setdiff(t1_con,
                                                     t2_con)) != 0) & (length(t1_con) != 0 & length(t2_con) != 0)

    )
    {
      readr::write_lines(paste(fgtsv$transcript[(x)], " : ", fgtsv$transcript[(x+1)], "[", unique(gt_df$geneID[gt_df$transcriptID %in% c(fgtsv$transcript[(x)], fgtsv$transcript[(x+1)])]), "]", '\n',
                               " iso1 trans #: ", length(t1_con),'\t', " iso1 genes #: ", length(unique(gt_df$geneID[gt_df$transcriptID %in% t1_con])), '\n',
                               " iso2 trans #: ", length(t2_con),'\t', " iso2 genes #: ", length(unique(gt_df$geneID[gt_df$transcriptID %in% t2_con])), '\n',
                               " intersect trans #: ", length(intersect(t1_con,
                                                                        t2_con)),'\t', " intersect genes #: ", length(unique(intersect(gt_df$geneID[gt_df$transcriptID %in% t1_con],
                                                                                                                                       gt_df$geneID[gt_df$transcriptID %in% t2_con]))), '\n',
                               " iso1 trans sd #: ", length(setdiff(t1_con,
                                                                    t2_con)),'\t', " iso1 genes unique #: ", length(unique(setdiff(gt_df$geneID[gt_df$transcriptID %in% t1_con],
                                                                                                                                   gt_df$geneID[gt_df$transcriptID %in% t2_con]))), '\n',
                               " iso2 trans sd #: ", length(setdiff(t2_con,
                                                                    t1_con)), '\t', " iso2 genes unique #: ", length(unique(setdiff(gt_df$geneID[gt_df$transcriptID %in% t2_con],
                                                                                                                                    gt_df$geneID[gt_df$transcriptID %in% t1_con]))),
                               '\n','\n','\n'
                               , sep = ""), path = paste0(output_location, "tti_output/adjacency_stats.txt"),
                         append = T)

    }


    c(length(t1_con),
      length(t2_con),
      length(intersect(t1_con,
                       t2_con)),
      length(setdiff(t1_con,
                     t2_con)),
      length(setdiff(t2_con,
                     t1_con)),

      length(unique(gt_df$geneID[gt_df$transcriptID %in% t1_con])),
      length(unique(gt_df$geneID[gt_df$transcriptID %in% t2_con])),
      length(unique(intersect(gt_df$geneID[gt_df$transcriptID %in% t1_con],
                              gt_df$geneID[gt_df$transcriptID %in% t2_con]))),
      length(unique(setdiff(gt_df$geneID[gt_df$transcriptID %in% t1_con],
                            gt_df$geneID[gt_df$transcriptID %in% t2_con]))),
      length(unique(setdiff(gt_df$geneID[gt_df$transcriptID %in% t2_con],
                            gt_df$geneID[gt_df$transcriptID %in% t1_con])))
    )
  })


  dtr <- fgtsv$transcript[sort(c(
    (which((unlist(lapply(domainDiff, "[[", 4)) != 0 | unlist(lapply(domainDiff, "[[", 5)) != 0) & (unlist(lapply(domainDiff, "[[", 1)) != 0 & unlist(lapply(domainDiff, "[[", 2)) != 0))*2)-1,
    which((unlist(lapply(domainDiff, "[[", 4)) != 0 | unlist(lapply(domainDiff, "[[", 5)) != 0) & (unlist(lapply(domainDiff, "[[", 1)) != 0 & unlist(lapply(domainDiff, "[[", 2)) != 0))*2))]

  ge <- igraph::graph_from_edgelist(as.matrix(t_d_df[!duplicated(t_d_df[,c(1,2)]),c(1,2)]), directed = F)
  d.ge <- igraph::degree(
    ge,
    v = igraph::V(ge),
    mode = c("all", "out", "in", "total")[1],
    loops = TRUE,
    normalized = FALSE
  )
  print(paste0("degree stats... Min: ", round(as.numeric(summary(d.ge))[1], 3),
               ", Median: ", round(as.numeric(summary(d.ge))[3], 3),
               ", Mean: ", round(as.numeric(summary(d.ge))[4], 3),
               ", Max: ", round(as.numeric(summary(d.ge))[6], 3)))


  ## vars: order, viz max,
  eg <- igraph::make_ego_graph(
    ge,
    order = steps,
    nodes = dtr,
    mode = c("all", "out", "in")[1],
    mindist = 0
  )


  graph_output <- list()
  for (i in seq(1, length(dtr)-1, by = 2)) {
    if (i == 1) {
      graph_output <- list()
    }
    e1g <- eg[[i]]
    e2g <- eg[[i+1]]

    e1g <- igraph::simplify(e1g, remove.multiple = F, remove.loops = T)
    #
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

    igraph::V(e1g)$color[V(e1g)$name == dtr[i]] <- "gold"
    igraph::V(e2g)$color[V(e2g)$name == dtr[i+1]] <- "gold"

    if (plot_bool)
    {
      if (sum(c(length(igraph::V(e1g)),
                length(igraph::V(e2g))) <= max_vertices_for_viz) == 2) {

        pdf(paste0(output_location, "tti_output/", dtr[i], "_", steps, 'steps_tti_graph.pdf'))
        print(igraph::plot.igraph(e1g,vertex.size=3,vertex.label=NA,main=dtr[i],
                                  layout=igraph::layout.fruchterman.reingold(e1g, niter=10000)))
        dev.off()

        pdf(paste0(output_location, "tti_output/", dtr[i+1], "_", steps, 'steps_tti_graph.pdf'))
        print(igraph::plot.igraph(e2g,vertex.size=3,vertex.label=NA,main=dtr[i+1],
                                  layout=igraph::layout.fruchterman.reingold(e2g, niter=10000)))
        dev.off()
      }
    }

    ##get nodes unique to one or both graphs

    ext1_edge <- igraph::V(e1g)$name[igraph::V(e1g)$color == "red"]
    ext2_edge <- igraph::V(e2g)$name[igraph::V(e2g)$color == "chartreuse4"]
    g1_unique <- gt_df$geneID[gt_df$transcriptID %in% unique(ext1_edge)]
    g2_unique <- gt_df$geneID[gt_df$transcriptID %in% unique(ext2_edge)]

    backgroundGenes <- unique(gttf$geneName[gttf$geneID %in% bgtsv$geneID])

    geneInput1 <- gttf$geneName[gttf$geneID %in% g1_unique]
    enrList1 <- unique(geneInput1[geneInput1 != gt_df$geneID[gt_df$transcriptID %in% dtr[i]]])

    geneInput2 <- gttf$geneName[gttf$geneID %in% g2_unique]
    enrList2 <- unique(geneInput2[geneInput2 != gt_df$geneID[gt_df$transcriptID %in% dtr[i+1]]])

    fdrUse <- fdr

    if (length(enrList1) > 0) {
      cc_1 <- hypeR::hyp_dots(hypeR::hypeR(enrList1, GO.cc, background = length(backgroundGenes), test="hypergeometric"),
                              fdr = fdrUse, title = "GO Cell Comp Enrichment", merge = TRUE)
      mf_1 <- hypeR::hyp_dots(hypeR::hypeR(enrList1, GO.mf, background = length(backgroundGenes), test="hypergeometric"),
                              fdr = fdrUse, title = "GO Mol FXN Enrichment", merge = TRUE)
      bp_1 <- hypeR::hyp_dots(hypeR::hypeR(enrList1, GO.bp, background = length(backgroundGenes), test="hypergeometric"),
                              fdr = fdrUse, title = "GO Biol Proc Enrichment", merge = TRUE)
      kg_1 <- hypeR::hyp_dots(hypeR::hypeR(enrList1, genesetsC2, background = length(backgroundGenes), test="hypergeometric"),
                              fdr = fdrUse, title = "KEGG Enrichment", merge = TRUE)
      bc_1 <- hypeR::hyp_dots(hypeR::hypeR(enrList1, geneset_BIOCARTA, background = length(backgroundGenes), test="hypergeometric"),
                              fdr = fdrUse, title = "Biocarta Enrichment", merge = TRUE)
      hm_1 <- hypeR::hyp_dots(hypeR::hypeR(enrList1, genesetsH, background = length(backgroundGenes), test="hypergeometric"),
                              fdr = fdrUse, title = "Hallmark Enrichment", merge = TRUE)

      pdf(paste0(output_location, "tti_output/", dtr[i], '_unique_tti_', steps, '_steps_enrichment.pdf'))
      print(cc_1)
      print(mf_1)
      print(bp_1)
      print(kg_1)
      print(bc_1)
      print(hm_1)
      dev.off()
    }
    if (length(enrList2) > 0) {
      cc_2 <- hypeR::hyp_dots(hypeR::hypeR(enrList2, GO.cc, background = length(backgroundGenes), test="hypergeometric"),
                              fdr = fdrUse, title = "GO Cell Comp Enrichment", merge = TRUE)
      mf_2 <- hypeR::hyp_dots(hypeR::hypeR(enrList2, GO.mf, background = length(backgroundGenes), test="hypergeometric"),
                              fdr = fdrUse, title = "GO Mol FXN Enrichment", merge = TRUE)
      bp_2 <- hypeR::hyp_dots(hypeR::hypeR(enrList2, GO.bp, background = length(backgroundGenes), test="hypergeometric"),
                              fdr = fdrUse, title = "GO Biol Proc Enrichment", merge = TRUE)
      kg_2 <- hypeR::hyp_dots(hypeR::hypeR(enrList2, genesetsC2, background = length(backgroundGenes), test="hypergeometric"),
                              fdr = fdrUse, title = "KEGG Enrichment", merge = TRUE)
      bc_2 <- hypeR::hyp_dots(hypeR::hypeR(enrList2, geneset_BIOCARTA, background = length(backgroundGenes), test="hypergeometric"),
                              fdr = fdrUse, title = "Biocarta Enrichment", merge = TRUE)
      hm_2 <- hypeR::hyp_dots(hypeR::hypeR(enrList2, genesetsH, background = length(backgroundGenes), test="hypergeometric"),
                              fdr = fdrUse, title = "Hallmark Enrichment", merge = TRUE)

      pdf(paste0(output_location, "tti_output/", dtr[i+1], '_unique_tti_', steps, '_steps_enrichment.pdf'))
      print(cc_2)
      print(mf_2)
      print(bp_2)
      print(kg_2)
      print(bc_2)
      print(hm_2)
      dev.off()
    }
    graph_output[[i]] <- e1g
    graph_output[[i+1]] <- e2g
  }
  names(graph_output) <- dtr
  saveRDS(graph_output, paste0(output_location, 'tti_output/igraph_objects.RDS'))
  return(list(dtr = dtr,
              graph_output = graph_output))
}
