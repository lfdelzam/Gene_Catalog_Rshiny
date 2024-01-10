rm(list = ls())
library("arrow")
C_pus <- 8


directorio_db= "data"

path_to_db_p <- paste(directorio_db,"AA/rep_proteins.faa", sep="/")
path_to_db <- paste(directorio_db,"DNA/rep_genes.fna", sep="/")
path_to_taxonomy_db <- paste(directorio_db,"Mmseqs2_rep_genes_taxonomy.tsv", sep="/")
path_to_annot_db <- paste(directorio_db,"rep_annotation.tsv", sep="/")
path_to_metadata <- paste(directorio_db,"BS_MAGv2_sample_metadata.tsv", sep="/")
path_to_clusters <- paste(directorio_db,"Rep_clusters_all.tsv", sep="/")
path_to_taxonomy_db_cat <- paste(directorio_db,"CAT_rep_genes_taxonomy.tsv", sep="/")

KEEG_list<-function(x ,sep){
  my_list <- str_split(x, sep)%>% unlist()
  element<-my_list[grep("[K]\\d",  my_list)]
  text <- gsub("ko:","",element) 
  if (length(text) > 0) {
    result <- paste(as.vector(unique(text)), collapse = ',')
  }
  else {
    result <- NA
  }
  
  return(result)
}

COG_list<-function(x ,sep){
  my_list <- str_split(x, sep)%>% unlist()
  element<-my_list[grep("COG\\d+.@",  my_list, fixed = F)]
  text<-gsub("@.*","",element)
  if (length(text) > 0) {
    result <- paste(as.vector(unique(text)), collapse = ',')
  }
  else {
    result <- NA
  }
  return(result)
}

cluster <- new_cluster(C_pus)

if ( !file.exists(paste(directorio_db, "big_tbl.rds", sep="/") )) {
  
  if ( !file.exists(paste(directorio_db, "tax.rds", sep="/") )) {
    tax <- read_tsv(path_to_taxonomy_db, show_col_types = FALSE, col_names=F)
    tax <- tax[-1,]
    colnames(tax) <- c("GeneID", "MMseq2_assigned_taxonomy")
    if (file.exists(path_to_taxonomy_db_cat)) {
      tax_cat <- read_tsv(path_to_taxonomy_db_cat, show_col_types = FALSE, col_names=F)
      tax_cat <- tax_cat[-1,]  
      colnames(tax_cat) <- c("GeneID", "CAT_assigned_taxonomy")
      tax_b <- merge(tax, tax_cat, by= "GeneID")
      tax_b$MMseq2_assigned_taxonomy=gsub("-_cellular organisms;", "",tax_b$MMseq2_assigned_taxonomy)
      saveRDS(tax_b, file=paste(directorio_db, "tax.rds", sep="/"))
    } else { 
      saveRDS(tax, file=paste(directorio_db, "tax.rds", sep="/")) }
    
    
  } else { tax <- readRDS(paste(directorio_db,"tax.rds", sep="/")) }
  
  if ( !file.exists(paste(directorio_db, "ann.rds", sep="/") )) {
    ann <- read_tsv(path_to_annot_db, show_col_types = FALSE, col_names=T)
    colnames(ann)[1] <- "GeneID"
    
    ann <- ann %>% rowwise() %>% #so each command is performed separatly on each row
      mutate(KEGG = KEEG_list(KEGG_ko, ","), .before =6) %>%
      mutate(COG = COG_list(eggNOG_OGs, ","), .before =6) %>% partition(cluster) %>% collect()
    
    saveRDS(ann, file=paste(directorio_db, "ann.rds", sep="/"))
  } else { ann <- readRDS(paste(directorio_db,"ann.rds", sep="/")) }
  
  if (file.exists(path_to_taxonomy_db_cat)) { big_tbl <- tax_b %>%  full_join(ann, by = "GeneID") } else {
    big_tbl <- tax %>%  full_join(ann, by = "GeneID") }
  
  
  saveRDS(big_tbl, file=paste(directorio_db,"big_tbl.rds", sep="/"))
  write_parquet(big_tbl, paste(directorio_db,"big_tbl.parquet", sep="/"))
  
}

if ( !file.exists(paste(path_to_db, ".ndb", sep=""))) {
  makeblastdb(path_to_db, dbtype = "nucl")
}

if ( !file.exists(paste(path_to_db_p, ".pto", sep=""))) {
  makeblastdb(path_to_db_p, dbtype = "prot")
}


#get the coordinates
if (!file.exists(paste(directorio_db, "ref_coords.rds", sep="/"))) {
  ref_coords<-as.data.frame(t(read.table(path_to_metadata)))%>%
    rownames_to_column(.,var="ref_id")%>%
    select(ref_id, Lon, Lat,Sal,Depth,Temp ) %>%
    mutate(Lon=as.numeric(Lon), Lat=as.numeric(Lat)) %>% partition(cluster) %>% collect()
  
  saveRDS(ref_coords, file=paste(directorio_db, "ref_coords.rds", sep="/"))
}


if (file.exists(path_to_clusters) & !exists("cluster_df") ) {
  cluster_df<-read_tsv(path_to_clusters, col_names = c("Rep", "genes_in_cluster"), num_threads = C_pus, show_col_types = FALSE)
  
  saveRDS(cluster_df, file=paste(directorio_db, "cluster_df.rds", sep="/"))
  write_parquet(cluster_df, paste(directorio_db, "cluster_df.parquet", sep="/"))
  
}
