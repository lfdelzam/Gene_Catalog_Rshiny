rm(list = ls())
list_of_packages <- c("tidyverse","multidplyr","arrow")

for (s in list_of_packages) { suppressPackageStartupMessages(library(s, character.only = TRUE)) }

C_pus <- 8


directorio_db= "/data"

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

if ( !file.exists(paste(directorio_db, "big_tbl_red.parquet", sep="/") )) {
  
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
    
    ann <- ann %>% rowwise() %>% #so each command is performed separately on each row
      mutate(KEGG = KEEG_list(KEGG_ko, ","), .before =6) %>%
      mutate(COG = COG_list(eggNOG_OGs, ","), .before =6) %>% partition(cluster) %>% collect()
    
    saveRDS(ann, file=paste(directorio_db, "ann.rds", sep="/"))
  } else { ann <- readRDS(paste(directorio_db,"ann.rds", sep="/")) }
  
  if (file.exists(path_to_taxonomy_db_cat)) { big_tbl <- tax_b %>%  full_join(ann, by = "GeneID") } else {
    big_tbl <- tax %>%  full_join(ann, by = "GeneID") }
  
  big_tbl <- big_tbl %>% select(GeneID,Preferred_name,dbCAN_family, RFAM_accession,
                                PFAM_accession,EC, COG,COG_cat,
                                KEGG,Description,best_tax_level,
                                CAT_assigned_taxonomy,MMseq2_assigned_taxonomy
    ) %>% rename("best_tax_level"="Eggnog_best_tax_level") %>% partition(cluster) %>% collect()
  
  #saveRDS(big_tbl, file=paste(directorio_db,"big_tbl_red.rds", sep="/"))

  btmp <- big_tbl %>% dplyr::filter(MMseq2_assigned_taxonomy != "unclassified") %>% dplyr::select(GeneID,MMseq2_assigned_taxonomy)
  write_parquet(btmp,paste(directorio_db, paste0("big_tbl","MMseq2_assigned_taxonomy_S",".parquet"), sep="/") )
  rm(btmp)
  btmp <- big_tbl %>% dplyr::filter(CAT_assigned_taxonomy != "unclassified") %>% dplyr::select(GeneID,CAT_assigned_taxonomy)
  write_parquet(btmp,paste(directorio_db, paste0("big_tbl","CAT_assigned_taxonomy_S",".parquet"), sep="/") )
  rm(btmp)
  btmp <-big_tbl %>% dplyr::filter(!is.na(RFAM_accession)) %>% dplyr::filter(RFAM_accession != "-") %>% dplyr::select(GeneID,RFAM_accession)
  write_parquet(btmp,paste(directorio_db, paste0("big_tbl","RFAM_accession_S",".parquet"), sep="/") )
  rm(btmp)
  btmp <-big_tbl %>% dplyr::filter(!is.na(PFAM_accession) ) %>% dplyr::filter(PFAM_accession != "-") %>% dplyr::select(GeneID,PFAM_accession)
  write_parquet(btmp,paste(directorio_db, paste0("big_tbl","PFAM_accession_S",".parquet"), sep="/") )
  rm(btmp)
  btmp <-big_tbl %>% dplyr::filter(!is.na(KEGG) ) %>% dplyr::filter(KEGG != "-") %>% dplyr::select(GeneID,KEGG)
  write_parquet(btmp,paste(directorio_db, paste0("big_tbl","KEGG_S",".parquet"), sep="/") )
  rm(btmp)
  btmp <-big_tbl %>% dplyr::filter(!is.na(COG) ) %>% dplyr::filter( COG != "-") %>% dplyr::select(GeneID,COG)
  write_parquet(btmp,paste(directorio_db, paste0("big_tbl","COG_S",".parquet"), sep="/") )
  rm(btmp)
  btmp <-big_tbl %>% dplyr::filter(!is.na(dbCAN_family) ) %>% dplyr::filter(dbCAN_family != "-") %>% dplyr::select(GeneID,dbCAN_family)
  write_parquet(btmp,paste(directorio_db, paste0("big_tbl","dbCAN_family_S",".parquet"), sep="/") )
  rm(btmp)
  btmp <-big_tbl %>% dplyr::filter(!is.na(Preferred_name) ) %>% dplyr::filter(Preferred_name != "-" ) %>% dplyr::select(GeneID,Preferred_name)
  write_parquet(btmp,paste(directorio_db, paste0("big_tbl","Preferred_name_S",".parquet"), sep="/") )
  rm(btmp)
  btmp <-big_tbl %>% dplyr::filter(!is.na(Eggnog_best_tax_level) ) %>% dplyr::filter(Eggnog_best_tax_level != "-" ) %>% dplyr::select(GeneID,Eggnog_best_tax_level)
  write_parquet(btmp,paste(directorio_db, paste0("big_tbl","Eggnog_best_tax_level_S",".parquet"), sep="/") )
  rm(btmp)
  
  write_parquet(big_tbl, paste(directorio_db,"big_tbl_red.parquet", sep="/"))
  
}

if ( !file.exists(paste(path_to_db, ".ntf", sep=""))) {
  #Before creating the database, remove no relevant information from the headers of fasta files 
  #with the command sed -i s/".#.*"//g <path_to_db>
  makeblastdb(path_to_db, dbtype = "nucl",args = "-parse_seqids" )
}

if ( !file.exists(paste(path_to_db_p, ".ptf", sep=""))) {
  #Before creating the database, remove no relevant information from the headers of fasta files 
  #with the command sed -i s/".#.*"//g <path_to_db>
  makeblastdb(path_to_db_p, dbtype = "prot",args = "-parse_seqids")
}


#get the coordinates
if (!file.exists(paste(directorio_db, "ref_coords.rds", sep="/"))) {
  ref_coords<-as.data.frame(t(read.table(path_to_metadata)))%>%
    rownames_to_column(.,var="ref_id")%>%
    select(ref_id, Lon, Lat,Sal,Depth,Temp ) %>%
    mutate(Lon=as.numeric(Lon), Lat=as.numeric(Lat)) %>% partition(cluster) %>% collect()
  
  saveRDS(ref_coords, file=paste(directorio_db, "ref_coords.rds", sep="/"))
  #write_parquet(ref_coords, paste(directorio_db, "ref_coords.parquet", sep="/"))
}


if (file.exists(path_to_clusters) & !file.exists(paste(directorio_db, "cluster_dfco.parquet", sep="/") )) {
  cluster_df<-read_tsv(path_to_clusters, col_names = c("Rep", "genes_in_cluster"), num_threads = C_pus, show_col_types = FALSE)
  
  #saveRDS(cluster_df, file=paste(directorio_db, "cluster_df.rds", sep="/"))
  cluster_dfco <- cluster_df[grep("CO::", cluster_df[[1]]),]
  cluster_dfin <- cluster_df[-grep("CO::", cluster_df[[1]]),]
  
  write_parquet(cluster_dfin, paste(directorio_db, "cluster_dfin.parquet", sep="/"))
  write_parquet(cluster_dfco, paste(directorio_db, "cluster_dfco.parquet", sep="/"))
  
}

