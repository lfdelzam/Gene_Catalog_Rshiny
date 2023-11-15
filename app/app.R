list_of_packages <- c("tidyverse", "data.table", "rBLAST", "leaflet", "sp", "magrittr", "KEGGREST", "shiny", "shinythemes", "DT","shinyjs","vembedr", "multidplyr","shinybusy")

for (s in list_of_packages) { suppressPackageStartupMessages(library(s, character.only = TRUE)) }

source("app_functions.R")


shinyApp(
  ui=fluidPage( theme = shinytheme("spacelab"),
                useShinyjs(),
                titlePanel(div(
                  img(src='SciLifeLab_Logotype_Green_POS.png', align='left', height='30px'), img(src='Logos.png', align='right', height='40px'),
                  #br(),
                  h1(strong("BAltic Gene Set - BAGS.v1.1",
                            style='color:steelblue;font-family:Avenir;font-size:28px;padding:10px'),
                     # class="btn-primary",
                     align="center"),
                  hr()
                )
                ),
                sidebarPanel(width = 3,
                             ##INPUT
                             fluidRow(style='font-family:Avenir;front-size:24px',
                                      p("Search using filters on functional annotation and/or taxonomy:"),
                                      actionButton("searchFT", "Function/Taxonomy", icon=icon("magnifying-glass"),class="btn btn-primary"),hr(),
                                      p("Search using Basic Local Alignment Search Tool (BLAST):"),
                                      textAreaInput('query', paste('Enter FASTA sequence(s) - up to', Max_num_query), value = "", placeholder = "", width = "600px", height="150px"),
                                      tags$head(tags$style(".btn-file {
                                                color: #fff;
                                                background: #286090;
                                                border-color: #122b40;
                                                }"
                                      )
                                      ),
                                      fileInput("file_query","or, upload file", buttonLabel = "Browse"),
                                      div(sliderInput("N_hits",label="Max number of target sequences",value=10,min=1,max=100),
                                          sliderInput("minPerIden",label="Percentage Identity",value=95,min=1,max=100),
                                          sliderInput("minPeralg",label="Percentage query coverage",value=95,min=1,max=100), style="front-size:20px"),
                                      div(style="display:inline-block;front-size:20px",
                                          selectInput("B_type", "Program:", choices=c("blastn","blastp"), width="100px")),
                                      div(style="display:inline-block;front-size:20px",
                                          selectInput("eval", "e-value:",choices=c(1,0.01,0.001,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10), width="100px")),
                                      br(),
                                      actionButton("run_blast", "BLAST", icon=icon("magnifying-glass"),class="btn btn-primary"),
                                      actionButton("refresh", "reset", icon=icon("redo"),class="btn btn-primary")
                             ), #fluidrow
                ), #sidebarPanel
                ##OUTPUT
                div(id = "main",style='font-family:Avenir',
                    mainPanel(width = 9,
                              tabsetPanel(id="results_tabs",
                                          type="pills",
                                          tabPanel("Blast Results", br(),uiOutput("blast_results_selection"),
                                                   em("Select hit(s) and click on GO! to obtain functional annotation, taxonomy affiliation and geolocalisation"),hr(),
                                                   DT::dataTableOutput("blastResults"),
                                                   checkboxInput("dt_sel", "select/deselect all"),
                                                   actionButton("go_anotax", "GO!", icon=icon("magnifying-glass-location"),class="btn btn-primary"),
                                                   downloadButton("downloadblast", "Download",class="btn btn-primary"),
                                                   hr(),
                                                   icon= icon("folder-open")
                                                   
                                          ),
                                          tabPanel("Filter", icon=icon("filter"), br(),
                                                   p(strong("Select one or several filters, one item per filter."),em("Regular expresion can be used when Exact matching is set to FALSE.")),
                                                   br(),br(),
                                                   fixedRow(column(width = 4, uiOutput("assigned_taxonomy")),column(width = 4,uiOutput("PFAM_accession"))),
                                                   fixedRow(column(width = 4, uiOutput("KEGG")),column(width = 4,uiOutput("COG"))),
                                                   fixedRow(column(width = 4, uiOutput("dbCAN_family")),column(width = 4,uiOutput("Rfam_accession"))),
                                                   fixedRow(column(width = 4, uiOutput("Others")),column(width = 4,uiOutput("And_or"))),
                                                   fixedRow(column(width = 4, uiOutput("ignore_case")),column(width = 4, uiOutput("exact_matching"))),
                                                   actionButton("searchFT_2", "search", icon=icon("magnifying-glass"),class="btn btn-primary"),hr(),
                                                   em("Select hit(s) and click on GO! to obtain further KEGG functional annotation, and geolocalisation"),
                                                   br(),br(),
                                                   DT::dataTableOutput("AnnotaxResults"),
                                                   checkboxInput("dt_sel_all", "select/deselect all"),
                                                   actionButton("go_anotax_all", "GO!", icon=icon("magnifying-glass-location"),class="btn btn-primary"),
                                                   hr()
                                          ),
                                          tabPanel("Annot. & Tax", icon= icon("tag"),
                                                   br(),
                                                   strong("Functional annotation and taxonomy affiliation of selected hit(s) :"), br(),DT::dataTableOutput("Anotax_hits"),
                                                   downloadButton("downloadData", "Download",class="btn btn-primary"),
                                                   hr(),
                                                   em("Select KEGG annotation to obtain KEGG BRITTE and PATHWAY information"),
                                                   hr(),
                                                   strong("KEGG BRITTE on selected reference:"),br(),DT::dataTableOutput("Kegg_BRITTE"),
                                                   downloadButton("downloadbritte", "Download",class="btn btn-primary"),
                                                   hr(),
                                                   strong("KEGG PATHWAY on selected reference:"),br(),DT::dataTableOutput("Kegg_PATH"),
                                                   downloadButton("downloadpath", "Download",class="btn btn-primary"),
                                                   hr()
                                          ),
                                          tabPanel("Map", uiOutput("Map_ggplot"), icon= icon("location-dot"), hr()) #map
                              ) #ends tabsetPanel
                    ) #ens mainPanel
                ) %>% shinyjs::hidden(),
                
                div(id="footer", style="font-family:Avenir",
                    div(align="left", style= "font-size:16px",
                        p(
                          "The BAltic Gene Set gene catalogue encompasses 66,530,673 genes.",br(),
                          "The 66 million genes are based on metagenomic data from Alneberg at al. (2020) from 124 seawater samples",br(),
                          "that span the salinity and oxygen gradients of the Baltic Sea and capture seasonal dynamics at two locations.",br(),
                          "To obtain the gene catalogue, we used a mix-assembly approach described in Delgado et al. (2022).",br(),
                          "The gene catalogue has been functionally and taxonomically annotated, using the",
                          a(href="https://github/lfdelzam/mix_assembly_pipeline", " Mix-assembly Gene Catalog pipeline." ),br(),br(),
                          "The BAltic Gene Set gene catalogue can be downloaded ",a(href="https://figshare.scilifelab.se/articles/dataset/BAGS_v1_BAltic_Gene_Set_gene_catalogue/16677252", "here")
                        ),
                        br(),
                        embed_url("https://www.youtube.com/watch?v=PNqEjTsNxIs&pp=ygUrZXZhbHVhdGluZyBtZXRhZ2Vub21pYyBhc3NlbWJseSBhcHByb2FjaGVzIA%3D%3D"),
                        br(),
                        p(
                          "When using the BAGS gene catalogue, please cite:",br(),br(),
                          a(style= "font-size:14px",href= 'https://doi.org/10.1186/s40168-022-01259-2', "1. Delgado LF, Andersson AF. Evaluating metagenomic assembly approaches for biome-specific gene catalogues. Microbiome 10, 72 (2022)", br()),
                          a(style= "font-size:14px",href= 'https://doi.org/10.1038/s42003-020-0856-x', "2. Alneberg J, et. al. Ecosystem-wide metagenomic binning enables prediction of ecological niches from genomes. Commun Biol 3, 119 (2020)",br()),
                        ),
                        
                        
                        hr(),
                        span(#style= "font-size:16px",
                          "Taxonomic composition of the BAGS.", a(href="rep_genes_taxonomy_krona.html", "Click to open the chart" )
                        ),
                        hr()
                    ), #div
                    
                    fluidRow(style="font-size:12px;color:#aeabab",
                             column(2,"FUNDING",br(),"Swedish Biodiversity Data Infrastructure",br(),"Swedish Research Council"),
                             column(2,"PUBLISHER", br(),"Swedish Biodiversity Data Infrastructure (SBDI)"),
                             column(2,"CONTACT EMAIL", br(), "anders.andersson@scilifelab.se"),
                             column(2,"© 2023 authored by Luis F. Delgado, Anders Andersson")
                    ),
                    p(style="font-size:14px;color:#aeabab",align="center",
                      "Rshiny app developped by Luis F. Delgado, Simon Kebede Merid, Gilbert Osena, Marco Vicari, Samah Abousharieha"
                    ),
                ) #div(id="footer"
                
  ), #close fluidPage
  
  server=function(input,output, session) {
    observeEvent(input$refresh, {
      session$reload()
    },
    priority=1000)
    
    event_action <- reactive({
      list(input$run_blast, input$searchFT)
    })
    
    observeEvent(ignoreInit = TRUE, event_action() ,{
      shinyjs::toggle("main")
      shinyjs::toggle("footer")
    }, once = F)
    
    observeEvent(input$searchFT, {
      updateTabsetPanel(session, "results_tabs",selected = "Filter")
    })
    
    observeEvent(input$run_blast, {
      updateTabsetPanel(session, "results_tabs", selected = "Blast Results")
    },priority=4)
    
    
    blastresults <- eventReactive(input$run_blast, {
      text<-NULL
      if (grepl(">",isolate(input$query))){
        texti=isolate(input$query)
        dir.create(directorio_qr)
        if (input$B_type == "blastn") {
          path_to_query_n <- paste(directorio_qr, paste(session$token,"queryp.fna", sep="_"), sep="/")
          fwrite(as.list(texti), file=path_to_query_n, quote = F)
          text=isolate(readDNAStringSet(path_to_query_n))
        } else if (input$B_type == "blastp") {
          path_to_query_p <- paste(directorio_qr, paste(session$token,"queryp.faa", sep="_"), sep="/")
          fwrite(as.list(texti), file=path_to_query_p, quote = F)
          text=isolate(readAAStringSet(path_to_query_p))
        } 
      } else if (!is.null(isolate(input$file_query$datapath))) {
        if (input$B_type == "blastn") {
          text=readDNAStringSet(isolate(input$file_query$datapath))
        } else if (input$B_type == "blastp") {
          text=readAAStringSet(isolate(input$file_query$datapath))
        }
      }
      if (is.null(text)) { stop("Please provide query sequences in fasta format")}
      if (length(text) > Max_num_query) { stop(paste("Please provide a maximum of", Max_num_query, "queries", sep=" "))}
      showModal( modalDialog( "Running blast", add_busy_spinner(spin = "fingerprint") )   )
      
      LISTA=Hit_blast(input$B_type, text, C_pus, input$N_hits, input$eval, input$minPerIden, input$minPeralg)
      removeModal()
      return(LISTA)
    },  ignoreNULL= T, ignoreInit=F)
    
    output$assigned_taxonomy <- renderUI({
      textInput("selected_tax",label = "Taxonomy", value=NULL)
    })
    output$dbCAN_family <- renderUI({
      textInput("selected_dbCAN_family",label = "dbCAN family", value=NULL)
    })
    output$PFAM_accession <- renderUI({
      textInput("selected_PFAM_accession",label = "PFAM accession", value=NULL)
    })
    output$Rfam_accession <- renderUI({
      textInput("selected_Rfam_accession",label = "Rfam accession", value=NULL)
    })
    output$COG <- renderUI({
      textInput("selected_COG",label = "COG", value=NULL)
    })
    output$KEGG <- renderUI({
      textInput("selected_KEGG",label = "KEGG", value=NULL)
    })
    output$Others <- renderUI({
      textInput("selected_Others",label = "Preferred name", value=NULL)
    })
    output$And_or<- renderUI({
      selectInput("selected_And_or", "Logical operator", choices=c("and","or"), width="200px")
    })
    output$ignore_case<- renderUI({
      selectInput("selected_ignore_case", "Ignore case", choices=c("TRUE","FALSE"), width="200px")
    })
    output$exact_matching<- renderUI({
      selectInput("selected_exact_matching", "Exact matching", choices=c("FALSE","TRUE"), width="200px")
    })
    
    output$blast_results_selection <- renderUI({
      if(!is.null(blastresults())){
        bs <- blastresults()
        selectInput("selected_table",label = "Select Query",choices = names(bs),selected=1,multiple = F)
      } else {stop("Please provide query sequences in fasta format")}
    })
    
    output$blastResults <- renderDataTable({
      if (!is.null(blastresults()) && !is.null(input$selected_table)){
        blastresults()[[input$selected_table]]
      }
    }, selection="multiple", future = F, server=T)
    
    
    Anotaxhits <- eventReactive({ input$go_anotax && !is.null(input$selected_table) && !is.null(input$blastResults_rows_selected)}, {
      s = input$blastResults_rows_selected
      if (!is.null(s) ) {
        hit=as.character(blastresults()[[input$selected_table]][s, "SubjectID"])
        showModal(modalDialog( "Retrieving information", add_busy_spinner(spin = "fingerprint") ) )
        thehits<-get_annotation(hit)
        removeModal()
        return(thehits)
      } else { return(NULL)}
    },  ignoreNULL= T, ignoreInit=F)
    
    Filter_table <- eventReactive({ input$searchFT_2}, {
      showModal(modalDialog( "Retrieving information", add_busy_spinner(spin = "fingerprint")) )
      
      FThits <- NULL
      if ((input$selected_And_or == "and") && (input$selected_tax != "" |input$selected_PFAM_accession != "" |input$selected_KEGG != "" |input$selected_COG != "" | input$selected_dbCAN_family != ""| input$selected_Rfam_accession != "" | input$selected_Others != "" ) ) {
        
        FThits <- big_tbl %>% select(GeneID,
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
                                     assigned_taxonomy
        )
        
        if (input$selected_tax != "") {
          try(FThitsTa <- FThits %>% filter(grepl(input$selected_tax,assigned_taxonomy, ignore.case = input$selected_ignore_case, fixed=input$selected_exact_matching)) %>% partition(cluster) %>% collect() , silent = T)
          if (exists("FThitsTa")) { FThits <- FThitsTa } else {FThits <- NULL}
        }
        if (input$selected_PFAM_accession != "" && !is.null(FThits)) {
          
          try(FThitsTb <- FThits %>% filter(grepl(input$selected_PFAM_accession,PFAM_accession, ignore.case = input$selected_ignore_case, fixed=input$selected_exact_matching)) %>% partition(cluster) %>% collect(), silent = T)
          if (exists("FThitsTb")) { FThits <- FThitsTb } else {FThits <- NULL}
        }
        if (input$selected_KEGG != "" && !is.null(FThits)) {
          try(FThitsTc <- FThits %>% filter(grepl(input$selected_KEGG,KEGG), ignore.case = input$selected_ignore_case, fixed=input$selected_exact_matching) %>% partition(cluster) %>% collect(), silent = T)
          if (exists("FThitsTc")) { FThits <- FThitsTc } else {FThits <- NULL}
        }
        if (input$selected_COG != "" && !is.null(FThits)) {
          try(FThitsTd <- FThits %>% filter(grepl(input$selected_COG,COG), ignore.case = input$selected_ignore_case, fixed=input$selected_exact_matching) %>% partition(cluster) %>% collect(), silent = T)
          if (exists("FThitsTd")) { FThits <- FThitsTd } else {FThits <- NULL}
        }
        if (input$selected_dbCAN_family != "" && !is.null(FThits)) {
          try(FThitsTe <- FThits %>% filter(grepl(input$selected_dbCAN_family,dbCAN_family, ignore.case = input$selected_ignore_case, fixed=input$selected_exact_matching)) %>% partition(cluster) %>% collect(), silent = T)
          if (exists("FThitsTe")) { FThits <- FThitsTe } else {FThits <- NULL}
        }
        if (input$selected_Rfam_accession != "" && !is.null(FThits)) {
          try(FThitsTf <- FThits %>% filter(grepl(input$selected_Rfam_accession,RFAM_accession, ignore.case = input$selected_ignore_case, fixed=input$selected_exact_matching)) %>% partition(cluster) %>% collect(), silent = T)
          if (exists("FThitsTf")) { FThits <- FThitsTf } else {FThits <- NULL}
        }
        if (input$selected_Others != "" && !is.null(FThits)) {
          try(FThitsTg <- FThits %>% filter(grepl(input$selected_Others,Preferred_name, ignore.case = input$selected_ignore_case, fixed=input$selected_exact_matching)) %>% partition(cluster) %>% collect(), silent = T)
          if (exists("FThitsTg")) { FThits <- FThitsTg
          } else {FThits <- NULL}
        }
      } else if ((input$selected_And_or == "or") && (input$selected_tax != "" |input$selected_PFAM_accession != "" |input$selected_KEGG != "" |input$selected_COG != "" | input$selected_dbCAN_family != ""| input$selected_Rfam_accession != "" | input$selected_Others != "") ) {
        big_tbls = big_tbl%>% select(GeneID,
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
                                     assigned_taxonomy
        )
        
        list_df = vector("list", 7)
        if (input$selected_tax != "") {
          try(list_df1 <- big_tbls %>% filter(grepl(input$selected_tax,assigned_taxonomy, ignore.case = input$selected_ignore_case, fixed=input$selected_exact_matching)) %>% partition(cluster) %>% collect(), silent =T)
          if (exists("list_df1")) { list_df[[1]] <- list_df1 }
        }
        if (input$selected_PFAM_accession != "") {
          try(list_df2 <- big_tbls %>% filter(grepl(input$selected_PFAM_accession,PFAM_accession, ignore.case = input$selected_ignore_case, fixed=input$selected_exact_matching)) %>% partition(cluster) %>% collect(), silent =T)
          if (exists("list_df2")) { list_df[[2]] <- list_df2 }
        }
        if (input$selected_KEGG != "") {
          try(list_df3 <- big_tbls %>% filter(grepl(input$selected_KEGG,KEGG, ignore.case = input$selected_ignore_case, fixed=input$selected_exact_matching)) %>% partition(cluster) %>% collect(), silent =T)
          if (exists("list_df3")) { list_df[[3]] <- list_df3 }
        }
        if (input$selected_COG != "") {
          try(list_df4 <- big_tbls %>% filter(grepl(input$selected_COG,COG, ignore.case = input$selected_ignore_case, fixed=input$selected_exact_matching)) %>% partition(cluster) %>% collect(), silent =T)
          if (exists("list_df4")) { list_df[[4]] <- list_df4}
        }
        if (input$selected_dbCAN_family != "") {
          try(list_df5<- big_tbls %>% filter(grepl(input$selected_dbCAN_family,dbCAN_family, ignore.case = input$selected_ignore_case, fixed=input$selected_exact_matching)) %>% partition(cluster) %>% collect(), silent =T)
          if (exists("list_df5")) { list_df[[5]] <- list_df5 }
        }
        if (input$selected_Rfam_accession != "") {
          try(list_df6 <- big_tbls %>% filter(grepl(input$selected_Rfam_accession,RFAM_accession, ignore.case = input$selected_ignore_case, fixed=input$selected_exact_matching)) %>% partition(cluster) %>% collect(), silent =T)
          if (exists("list_df6")) { list_df[[6]] <- list_df6 }
        }
        if (input$selected_Others != "") {
          try(list_df7 <- big_tbls %>% filter(grepl(input$selected_Others,Preferred_name, ignore.case = input$selected_ignore_case, fixed=input$selected_exact_matching)) %>% partition(cluster) %>% collect(), silent =T)
          if (exists("list_df7")) { list_df[[7]] <- list_df7 }
        }
        
        list_df = list_df[lapply(list_df,length)>0]
        if (length(list_df) > 1) {
          FThits<- list_df %>% purrr::reduce(full_join, by='GeneID')
          names(FThits) <-names(big_tbls)
        } else if (length(list_df) == 1) {
          FThits<- list_df[[1]]
          names(FThits) <-names(big_tbls)
        }
      }
      
      removeModal()
      return(FThits)
    },  ignoreNULL= T, ignoreInit=F)
    
    output$AnnotaxResults <- renderDataTable({
      if (!is.null(Filter_table()) ) {
        return(Filter_table())
      } else {stop("No hits found")}
    }, selection="multiple", future = F, server=T)
    
    
    Anotaxhits_all <- eventReactive({ input$go_anotax_all && !is.null(input$AnnotaxResults_rows_selected)}, {
      s = input$AnnotaxResults_rows_selected
      if (!is.null(s) ) {
        hit=Filter_table()$GeneID[s]
        showModal(modalDialog( "Retrieving information", add_busy_spinner(spin = "fingerprint") ) )
        
        thehits<-get_annotation(hit)
        removeModal()
        return(thehits)
      } else { return(NULL)}
    },  ignoreNULL= T, ignoreInit=F)
    
    output$Anotax_hits <- renderDataTable({
      if (!is.null(Anotaxhits())){
        return(Anotaxhits())
      } else if (!is.null(Anotaxhits_all())){
        return(Anotaxhits_all())
      } else { return(NULL)}
    }, selection="single", future = F)
    
    
    dt_proxy <- DT::dataTableProxy("blastResults")
    observeEvent(input$dt_sel, {
      
      if (isTRUE(input$dt_sel)) {
        DT::selectRows(dt_proxy, input$blastResults_rows_all)
      } else {
        DT::selectRows(dt_proxy, NULL)
      }
    },ignoreNULL = TRUE,
    ignoreInit = T)
    
    
    observeEvent(input$go_anotax , {
      updateTabsetPanel(session, "results_tabs", selected ="Annot. & Tax") #c("Annot. & Tax","Map")
    },priority=4)
    
    observeEvent(input$go_anotax_all , {
      updateTabsetPanel(session, "results_tabs", selected ="Annot. & Tax")
    },priority=4)
    
    dt_proxy_all <- DT::dataTableProxy("AnnotaxResults")
    observeEvent(input$dt_sel_all, {
      
      if (isTRUE(input$dt_sel_all)) {
        DT::selectRows(dt_proxy_all, input$AnnotaxResults_rows_all)
      } else {
        DT::selectRows(dt_proxy_all, NULL)
      }
    },ignoreNULL = TRUE,
    ignoreInit = T)
    
    observeEvent(input$searchFT_2, {
      updateCheckboxGroupInput(session, "dt_sel_all", selected = character(1))
    },ignoreNULL = F,
    ignoreInit = T)
    
    
    observeEvent(input$selected_table, {
      updateCheckboxGroupInput(session, "dt_sel", selected = character(1))
    },ignoreNULL = F,
    ignoreInit = T)
    
    output$downloadData <- downloadHandler(
      filename = function() {
        texto=""
        if (!is.null(Anotaxhits())) {
          
          RENAME=unlist(strsplit(input$selected_table," "))[1]
          RENAME=gsub("::","_",RENAME)
          texto=paste(RENAME, "_annot_and_tax_output.tsv", sep = "")
        } else if (!is.null(Anotaxhits_all())) {
          texto=paste("Filter","annot_and_tax_output.tsv", sep = "_")
        }
        return(texto)
      },
      content = function(file) {
        if (!is.null(Anotaxhits())) {
          fwrite(Anotaxhits(), file, sep="\t", row.names = FALSE)
        } else if (!is.null(Anotaxhits_all())) {
          fwrite(Anotaxhits_all(), file, sep="\t", row.names = FALSE)
        }
      }
    )
    
    
    output$downloadblast <- downloadHandler(
      filename = function() {
        RENAME=unlist(strsplit(input$selected_table," "))[1]
        RENAME=gsub("::","_",RENAME)

        return(paste(RENAME, "_blast_output.tsv", sep = ""))
      },
      content = function(file) {
        fwrite(blastresults()[[input$selected_table]], file, sep="\t", row.names = FALSE)
      }
    )
    
    
    Kegg <- eventReactive(input$Anotax_hits_rows_selected , {
      tmp <- NULL
      s = input$Anotax_hits_rows_selected
      if (length(s)) {
        if (!is.null(Anotaxhits())) {
          keggR=as.character(Anotaxhits()[s, "KEGG"])
        } else if (!is.null(Anotaxhits_all())) {  #new
          keggR=as.character(Anotaxhits_all()[s, "KEGG"])  #new
        }
        if (!is.na(keggR)) {
          showModal(modalDialog( "Retrieving information", add_busy_spinner(spin = "fingerprint") ) )
          tmp <- expand_kegg(keggR)
          removeModal()
        }
        
      }
      
      return(tmp)
    },  ignoreNULL= T, ignoreInit = T )
    
    
    output$Kegg_BRITTE <- renderDataTable({
      outp<-NULL
      if (!is.null(Kegg())) {
        Kegg_BRITTE = Kegg()$BRITE
        if (!is.null(Kegg_BRITTE)){
          outp<-Kegg_BRITTE
        }
      }
      return(outp)
    }, selection="single", future = F)
    
    
    output$downloadbritte <- downloadHandler(
      filename = function() {
        texto=""
        s = input$Anotax_hits_rows_selected
        if (length(s)) {
          if (!is.null(Anotaxhits())) {
            
            REFNAME=gsub("::","_",as.character(Anotaxhits()[s, "GeneID"]))
            QNAME=unlist(strsplit(input$selected_table," "))[1]
            QNAME=gsub("::","_",QNAME)
            
            texto=paste("query",QNAME,"ref",REFNAME, "KEGG_BRITTE_output.tsv", sep = "_")
          } else if (!is.null(Anotaxhits_all())) {
            
            REFNAME=gsub("::","_",as.character(Anotaxhits_all()[s, "GeneID"]))
            texto=paste("ref",REFNAME, "KEGG_BRITTE_output.tsv", sep = "_")
            
          }
        }
        return(texto)
      },
      content = function(file) {
        fwrite(Kegg()$BRITE, file, sep="\t", row.names = FALSE)
      }
    )
    
    output$Kegg_PATH <- renderDataTable({
      outp <- NULL
      if (!is.null(Kegg())) {
        Kegg_PATH = Kegg()$PATHWAY
        if (!is.null(Kegg_PATH)){
          outp<-Kegg_PATH
        }
      }
      return(outp)
    }, selection="single", future = F)
    
    output$downloadpath <- downloadHandler(
      filename = function() {
        texto=""
        s = input$Anotax_hits_rows_selected
        if (length(s)) {
          
          if (!is.null(Anotaxhits())) {
            
            REFNAME=gsub("::","_",as.character(Anotaxhits()[s, "GeneID"]))
            QNAME=unlist(strsplit(input$selected_table," "))[1]
            QNAME=gsub("::","_",QNAME)
            texto=paste("query",QNAME,"ref",REFNAME, "KEGG_PATHWAY_output.tsv", sep = "_")            
          } else if (!is.null(Anotaxhits_all())) {
            
            REFNAME=gsub("::","_",as.character(Anotaxhits_all()[s, "GeneID"]))
            texto=paste("ref",REFNAME, "KEGG_BRITTE_output.tsv", sep = "_")
            
          }
        }
        return(texto)
      },
      content = function(file) {
        fwrite(Kegg()$PATHWAY, file, sep="\t", row.names = FALSE)
      }
    )
    
    output$Map_ggplot = renderUI({
      showModal(modalDialog( "Retrieving information", add_busy_spinner(spin = "fingerprint") ) )
      themap<-NULL
      s = input$blastResults_rows_selected
      if (length(s) ) {
        hit=as.character(blastresults()[[input$selected_table]][s, "SubjectID"])
        themap=get_localisation(hit)
        if (!is.null(themap)) {
          output$leafletMap <- renderLeaflet({themap})
          removeModal()
          return(leafletOutput('leafletMap', width = "100%", height = 700))
        } else {themap <- "Sorry, we don't have a localisation for the selected hit(s)"}
      }
      A = input$AnnotaxResults_rows_selected
      if (length(A) ) {
        hit=Filter_table()$GeneID[A]
        themap=get_localisation(hit)
        if (!is.null(themap)) {
          output$leafletMap <- renderLeaflet({themap})
          removeModal()
          return(leafletOutput('leafletMap', width = "100%", height = 700))
        } else {themap <- "Sorry, we don't have a localisation for the selected hit(s)"}
      }
      
      removeModal()
      return(themap)
    })
    
  })
