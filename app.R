library(shiny)

# ---------------------------------------------------------------------------
# LOOKUP TABLES
# ---------------------------------------------------------------------------

# Genotype mode: two alleles picked independently -> paste them -> look up
abo_genotype_table <- c(
  "AA" = "Blood Type A", "AO" = "Blood Type A", "OA" = "Blood Type A",
  "BB" = "Blood Type B", "BO" = "Blood Type B", "OB" = "Blood Type B",
  "AB" = "Blood Type AB", "BA" = "Blood Type AB",
  "OO" = "Blood Type O"
)

rh_genotype_table <- c(
  "DD" = "Rh Positive (+)", "Dd" = "Rh Positive (+)", "dD" = "Rh Positive (+)",
  "dd" = "Rh Negative (-)"
)

# Serology mode: reactions with antisera (TRUE/FALSE) -> paste them -> look up
abo_serology_table <- c(
  "TRUE.TRUE"   = "Blood Type AB",
  "TRUE.FALSE"  = "Blood Type A",
  "FALSE.TRUE"  = "Blood Type B",
  "FALSE.FALSE" = "Blood Type O"
)

rh_serology_table <- c(
  "TRUE"  = "Rh Positive (+)",
  "FALSE" = "Rh Negative (-)"
)

# ---------------------------------------------------------------------------
# HELPER: cross two parents' alleles -> % probability of each phenotype
# ---------------------------------------------------------------------------

cross_probabilities <- function(mom_alleles, dad_alleles, lookup_table) {
  combos <- expand.grid(a1 = mom_alleles, a2 = dad_alleles, stringsAsFactors = FALSE)
  keys <- paste0(combos$a1, combos$a2)
  phenotypes <- lookup_table[keys]
  probs <- round(prop.table(table(phenotypes)) * 100)
  paste0(names(probs), ": ", probs, "%", collapse = "   |   ")
}

# ---------------------------------------------------------------------------
# HELPER: build a colored blood-type badge (final visualization)
# ---------------------------------------------------------------------------

abo_color_table <- c(
  "A" = "#e74c3c", "B" = "#2980b9", "AB" = "#8e44ad", "O" = "#27ae60"
)

build_badge <- function(abo_phenotype, rh_phenotype = NULL) {
  letter <- sub("Blood Type ", "", abo_phenotype)
  bg_color <- abo_color_table[[letter]]
  symbol <- ""
  if (!is.null(rh_phenotype)) {
    symbol <- ifelse(grepl("Positive", rh_phenotype), "+", "-")
  }
  div(
    style = paste0(
      "display:inline-block; padding:20px 35px; border-radius:12px;",
      "background-color:", bg_color, "; color:white;",
      "font-size:36px; font-weight:bold; text-align:center;",
      "box-shadow:0 4px 8px rgba(0,0,0,0.2);"
    ),
    paste0(letter, symbol)
  )
}

# ---------------------------------------------------------------------------
# UI  (base shiny only)
# ---------------------------------------------------------------------------

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body { background-color: #f4f6f7; font-family: 'Segoe UI', sans-serif; }
      .navbar, .well { border-radius: 8px; }
      .well { background-color: #ffffff; box-shadow: 0 2px 6px rgba(0,0,0,0.08); border: none; }
      .nav-tabs > li > a { color: #c0392b; font-weight: 600; }
      .nav-tabs > li.active > a { color: #ffffff !important; background-color: #c0392b !important; border-color: #c0392b !important; }
      .btn, .form-control, .selectize-input { border-radius: 6px; }
      h4 { color: #c0392b; font-weight: 700; }
      titlePanel h1, .title-panel h1 { color: #c0392b; }
    "))
  ),

  titlePanel(
    div(style = "color:#c0392b; font-weight:800;", "Blood Group Predictor")
  ),

  tabsetPanel(

    # ---- TAB 1: Genotype-based prediction ---------------------------------
    tabPanel(
      "Genotype Predictor",
      sidebarLayout(
        sidebarPanel(
          h4("ABO Alleles"),
          selectInput("g_allele1", "Allele 1", choices = c("A", "B", "O")),
          selectInput("g_allele2", "Allele 2", choices = c("A", "B", "O")),

          checkboxInput("g_include_rh", "Also predict Rh factor?", value = FALSE),

          conditionalPanel(
            condition = "input.g_include_rh == true",
            h4("Rh Alleles"),
            selectInput("g_rh1", "Allele 1", choices = c("D", "d")),
            selectInput("g_rh2", "Allele 2", choices = c("D", "d"))
          )
        ),
        mainPanel(
          h4("Result"),
          h5("ABO Type:"),
          textOutput("g_abo_result"),
          br(),
          conditionalPanel(
            condition = "input.g_include_rh == true",
            h5("Rh Type:"),
            textOutput("g_rh_result")
          ),
          br(),
          uiOutput("g_badge")
        )
      )
    ),

    # ---- TAB 2: Serology-based prediction ----------------------------------
    tabPanel(
      "Serology Predictor",
      sidebarLayout(
        sidebarPanel(
          h4("ABO Reactions"),
          checkboxInput("s_anti_a", "Agglutinates with Anti-A?", value = FALSE),
          checkboxInput("s_anti_b", "Agglutinates with Anti-B?", value = FALSE),

          checkboxInput("s_include_rh", "Also test Rh factor?", value = FALSE),

          conditionalPanel(
            condition = "input.s_include_rh == true",
            h4("Rh Reaction"),
            checkboxInput("s_anti_d", "Agglutinates with Anti-D?", value = FALSE)
          )
        ),
        mainPanel(
          h4("Result"),
          h5("ABO Type:"),
          textOutput("s_abo_result"),
          br(),
          conditionalPanel(
            condition = "input.s_include_rh == true",
            h5("Rh Type:"),
            textOutput("s_rh_result")
          ),
          br(),
          uiOutput("s_badge")
        )
      )
    ),

    # ---- TAB 3: Parental cross -> offspring probability --------------------
    tabPanel(
      "Parental Cross Predictor",
      sidebarLayout(
        sidebarPanel(
          h4("Mother's ABO Alleles"),
          selectInput("p_mom1", "Allele 1", choices = c("A", "B", "O")),
          selectInput("p_mom2", "Allele 2", choices = c("A", "B", "O")),
          h4("Father's ABO Alleles"),
          selectInput("p_dad1", "Allele 1", choices = c("A", "B", "O")),
          selectInput("p_dad2", "Allele 2", choices = c("A", "B", "O")),

          checkboxInput("p_include_rh", "Also predict Rh factor?", value = FALSE),

          conditionalPanel(
            condition = "input.p_include_rh == true",
            h4("Mother's Rh Alleles"),
            selectInput("p_mom_rh1", "Allele 1", choices = c("D", "d")),
            selectInput("p_mom_rh2", "Allele 2", choices = c("D", "d")),
            h4("Father's Rh Alleles"),
            selectInput("p_dad_rh1", "Allele 1", choices = c("D", "d")),
            selectInput("p_dad_rh2", "Allele 2", choices = c("D", "d"))
          )
        ),
        mainPanel(
          h4("Offspring Probability"),
          h5("ABO Type Odds:"),
          textOutput("p_abo_result"),
          br(),
          conditionalPanel(
            condition = "input.p_include_rh == true",
            h5("Rh Type Odds:"),
            textOutput("p_rh_result")
          )
        )
      )
    ),

    # ---- TAB 4: About / Disclaimer -----------------------------------------
    tabPanel(
      "About",
      wellPanel(
        h4("About This App"),
        p("This app demonstrates ABO and Rh blood group inheritance using ",
          "simplified Mendelian genetics. It is built for educational purposes only."),
        h5("Limitations to keep in mind:"),
        tags$ul(
          tags$li("Rh factor here is modeled as a single gene (D/d). Real-world Rh ",
                  "biology involves multiple genes (RHD, RHCE) and weak/partial D variants."),
          tags$li("ABO genotype combinations shown assume standard codominant/recessive ",
                  "inheritance and don't account for rare variants (e.g. cis-AB, Bombay phenotype)."),
          tags$li("This tool is not a diagnostic or clinical device and should not be used ",
                  "for medical, paternity, or legal decisions.")
        )
      )
    )
  )
)

# ---------------------------------------------------------------------------
# SERVER
# ---------------------------------------------------------------------------

server <- function(input, output, session) {

  # ---- Genotype tab ----
  output$g_abo_result <- renderText({
    key <- paste0(input$g_allele1, input$g_allele2)
    abo_genotype_table[[key]]
  })

  output$g_rh_result <- renderText({
    req(input$g_include_rh)
    key <- paste0(input$g_rh1, input$g_rh2)
    rh_genotype_table[[key]]
  })

  output$g_badge <- renderUI({
    abo_key <- paste0(input$g_allele1, input$g_allele2)
    abo_pheno <- abo_genotype_table[[abo_key]]
    rh_pheno <- NULL
    if (isTRUE(input$g_include_rh)) {
      rh_key <- paste0(input$g_rh1, input$g_rh2)
      rh_pheno <- rh_genotype_table[[rh_key]]
    }
    build_badge(abo_pheno, rh_pheno)
  })

  # ---- Serology tab ----
  output$s_abo_result <- renderText({
    key <- paste(input$s_anti_a, input$s_anti_b, sep = ".")
    abo_serology_table[[key]]
  })

  output$s_rh_result <- renderText({
    req(input$s_include_rh)
    key <- paste(input$s_anti_d)
    rh_serology_table[[key]]
  })

  output$s_badge <- renderUI({
    abo_key <- paste(input$s_anti_a, input$s_anti_b, sep = ".")
    abo_pheno <- abo_serology_table[[abo_key]]
    rh_pheno <- NULL
    if (isTRUE(input$s_include_rh)) {
      rh_key <- paste(input$s_anti_d)
      rh_pheno <- rh_serology_table[[rh_key]]
    }
    build_badge(abo_pheno, rh_pheno)
  })

  # ---- Parental Cross tab ----
  output$p_abo_result <- renderText({
    mom_alleles <- c(input$p_mom1, input$p_mom2)
    dad_alleles <- c(input$p_dad1, input$p_dad2)
    cross_probabilities(mom_alleles, dad_alleles, abo_genotype_table)
  })

  output$p_rh_result <- renderText({
    req(input$p_include_rh)
    mom_alleles <- c(input$p_mom_rh1, input$p_mom_rh2)
    dad_alleles <- c(input$p_dad_rh1, input$p_dad_rh2)
    cross_probabilities(mom_alleles, dad_alleles, rh_genotype_table)
  })
}

shinyApp(ui, server)
