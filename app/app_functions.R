
#Query folder
directorio_qr <- "Query"
directorio_db<- "/data"
path_to_db_p <- paste(directorio_db,"AA/rep_proteins.faa", sep="/")
path_to_db <- paste(directorio_db,"DNA/rep_genes.fna", sep="/")

#Internal Blast parameters
Max_num_query <- 100 #Default
C_pus <- 8

if (!exists("big_tbl")) big_tbl <- readRDS(paste(directorio_db,"big_tbl.rds", sep="/"))
if (!exists("ref_coords")) ref_coords<-readRDS(paste(directorio_db, "ref_coords.rds", sep="/"))
if (!exists("cluster_df") ) cluster_df<-readRDS(paste(directorio_db, "cluster_df.rds", sep="/"))

#Functions

cluster <- new_cluster(C_pus)

dash_split <- function(x){
  unlist(strsplit(x,"::"))[1]
}

run_blast<-function(bl, query,cpus, hits, evalue, minIden,minalg, type) {
  lista=query@ranges@NAMES
  if (type == "blastn") {blast_arg=paste("-num_threads", cpus,"-max_target_seqs", hits, "-evalue", evalue, "-soft_masking false","-dust no") }
  if (type == "blastp") {blast_arg=paste("-num_threads", cpus,"-max_target_seqs", hits, "-evalue", evalue, "-soft_masking false") }
  blast_table <- vector("list", length =length(query))
  names(blast_table) <- lista
  for (i in 1:length(query)) {
    nam=lista[i]
    blast_table[[nam]] <- predict(bl, query[i,], BLAST_args = blast_arg)
    LEN=query[i,]@ranges@width
    blast_table[[nam]] <- blast_table[[nam]] %>% mutate(Perc.Query.Coverage = round((Alignment.Length-Gap.Openings)*100/LEN,2)) %>% filter(Perc.Ident > minIden & Perc.Query.Coverage > minalg)  
  }
  return(blast_table)
}

Hit_blast <- function(B_type, query, Cpus, N_hits, Evalue, minPer_Iden, minPer_alg) {
  
  if (B_type == "blastp") {
    blp <- blast(db=path_to_db_p, type = "blastp")
    RESULTS=run_blast(blp, query, Cpus, N_hits, Evalue, minPer_Iden, minPer_alg, B_type)
  }
  
  
  if (B_type == "blastn") {
    bln <- blast(db=path_to_db, type = "blastn")
    RESULTS=run_blast(bln, query, Cpus, N_hits, Evalue, minPer_Iden, minPer_alg, B_type)
  }
  return(RESULTS)
}

get_annotation <- function(selected_hit) {
  Annotation_and_taxonomy_result <- big_tbl %>% filter(GeneID %in% selected_hit) %>% select(GeneID,
                                                                                            Preferred_name,
                                                                                            dbCAN_family, RFAM_description,
                                                                                            RFAM_accession,
                                                                                            PFAM_family,
                                                                                            PFAM_accession,
                                                                                            EC,
                                                                                            COG,
                                                                                            COG_cat,
                                                                                            KEGG,
                                                                                            Description,
                                                                                            best_tax_level,
                                                                                            CAT_assigned_taxonomy,
                                                                                            MMseq2_assigned_taxonomy
  ) %>% rename("best_tax_level"="Eggnog_best_tax_level")
  
}

deconv_clusters <- function(x, co_patern) {

  ind_lista=unlist(strsplit(x,","))
  co_hits=ind_lista[grep(co_patern,ind_lista)]
  ind_lista=ind_lista[!ind_lista %in% co_hits] 
  return(ind_lista)}

get_localisation <- function(selected_hit) {

  map<-NULL
  co_patern="CO::"
  co_hits=selected_hit[grep(co_patern,selected_hit)]
  eq_ind_genes=c()
  if (length(co_hits) > 0) {

    for (coa in co_hits) {
      ind_genes=cluster_df %>% filter(Rep == coa)
      eq_ind_genes=c(eq_ind_genes,
                     unique(
                       sapply(ind_genes$genes_in_cluster, function(x) deconv_clusters(x, co_patern))
                     )
      )
    }
  }
    ind_hits=selected_hit[!selected_hit %in% co_hits]

    ia_genes=c()
    if (length(ind_hits) > 0) {
      
      for (ia in ind_hits) {
        ind_genes=cluster_df %>% filter(Rep == ia)
        ia_genes=c(ia_genes,
                   unique(
                     sapply(ind_genes$genes_in_cluster, function(x) deconv_clusters(x, co_patern))
                   )
        )
      }
    }  

      Sel_hits=unique(c(ia_genes,eq_ind_genes))

  selectes=unique(sapply(Sel_hits, dash_split))
  if (!is.null(selectes[[1]])) {

    hit_seqs <- ref_coords %>% filter(ref_id %in% selectes) %>% partition(cluster) %>% collect()
    if (nrow(hit_seqs) >0) {
      coordinates(hit_seqs) <- ~Lon + Lat
      
      map<-leaflet(hit_seqs) %>% addTiles() %>% addCircles() %>%
        addMarkers(data=hit_seqs,
                   clusterOptions = markerClusterOptions(zoomToBoundsOnClick = F),
                   popup = ~paste(
                     paste('<b>', 'Ref.:', '</b>', hit_seqs$ref_id),
                     paste('<b>',  'Sal.:', '</b>', hit_seqs$Sal),
                     paste('<b>', 'Temp.:','</b>', hit_seqs$Temp),
                     paste('<b>', 'Depth:', '</b>',hit_seqs$Depth),
                     sep = '<br/>'),popupOptions = popupOptions(closeButton = F)
        )
    }
  }
  return(map)
}

expand_kegg <- function(kegg_idin){
  
  df=list("BRITE"=NULL, "PATHWAY" = NULL)
  if (length(grep(",", kegg_idin)) > 0){
    my_lista <- str_split(kegg_idin, ",")%>% unlist()} else {my_lista = list(kegg_idin)}
  for (kegg_id in my_lista) {
    kegg_info <- keggGet(kegg_id)
    try(kegg_BRITE <- kegg_info[[1]]$BRITE %>% as_tibble(), silent=T)
    try(kegg_PATHWAY <-kegg_info[[1]]$PATHWAY %>% as_tibble(), silent = T)
    
    if (exists("kegg_BRITE")) { 
      
      df[["BRITE"]] = list(df[["BRITE"]], kegg_BRITE)  %>% bind_rows(.id = kegg_id) }
    if (exists("kegg_PATHWAY")) { 
      
      df[["PATHWAY"]] = list(df[["PATHWAY"]], kegg_PATHWAY)  %>% bind_rows(.id = kegg_id) } }
  
  return(df)
}


