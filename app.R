library(shiny)
library(dplyr)
library(rcdk)
library(randomForest)

# --------- Load pre-trained model ----------
model_path <- "phenol_pyridine_model.rds"
if (!file.exists(model_path)) stop("phenol_pyridine_model.rds not found in app folder!")
trained_model <- readRDS(model_path)

p_external <- 0.35
threshold <- 0.50

# --------- Safe RDKit descriptor function ----------
safe_rdkit_descriptors <- function(smiles, prefix) {
  tryCatch({
    mol <- parse.smiles(smiles)[[1]]
    if (is.null(mol)) return(NULL)
    desc_names <- get.desc.names(type = "all")
    desc_values <- eval.desc(mol, desc_names)
    names(desc_values) <- paste0(prefix, names(desc_values))
    as.list(desc_values)
  }, error = function(e) {
    message("Error parsing SMILES: ", smiles)
    return(NULL)
  })
}

# ---------- PubChem image helper ----------
smiles_to_url <- function(smiles) {
  encoded <- URLencode(smiles, reserved = TRUE)
  paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/",
    encoded,
    "/PNG?image_size=300x300"
  )
}

# --------- UI ----------
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body { background-color: #f5f5f5; }
      .centered { display:flex; flex-direction:column; align-items:center; margin-top:20px; }
      input[type=text] { font-size:16px; padding:8px; }
      #predict_btn { font-size:18px; padding:10px 25px; border-radius:12px; }
      .output_box { border-radius:12px; padding:20px; color:white; text-align:center; font-weight:bold; font-size:32px; min-width:160px; }
      #prob_box { font-size:16px; margin-top:10px; text-align:center; }
      #mix_box { font-size:14px; margin-top:8px; font-style:italic; text-align:center; }
      .mol_img { border:2px solid #2c2f36; border-radius:12px; padding:10px; background:#fff; width:300px; height:300px; }
      .input_col { padding:10px; display:flex; flex-direction:column; }
      .mol_row { display:flex; gap:50px; flex-wrap:wrap; justify-content:center; }
      .reaction_row { display:flex; align-items:center; gap:40px; margin-top:20px; }
      h1 { font-size:28px; margin-bottom:20px; }
    "))
  ),
  
  h1("Machine Learning Prediction for Photochemical Phenol Heteroarylation"),
  
  div(class="centered",
      div(class="mol_row",
          div(class="input_col",
              strong("Phenol SMILES:"),
              textInput("phenol", NULL, "", width="300px"),
              uiOutput("phenol_img")
          ),
          div(class="input_col",
              strong("Pyridine SMILES (must contain Cl, Br, or I):"),
              textInput("pyridine", NULL, "", width="300px"),
              uiOutput("pyridine_img")
          )
      ),
      div(class="reaction_row",
          actionButton("predict_btn", "Predict"),
          div(
            uiOutput("reaction_box"),
            uiOutput("prob_box"),
            uiOutput("mix_box")
          )
      ),
      uiOutput("validation_box")
  ),
  
  div(
    style="
      margin: 40px auto 20px auto;
      max-width: 900px;
      border: 1px solid #ccc;
      border-radius: 8px;
      padding: 12px;
      background-color: #fefefe;
      font-size: 12px;
      line-height: 1.4;
      color: #333;",
    HTML(
      "Credit: Matthew C. Carson and Kalyana B. Duggal for web design.<br>
       Model obtained from Carson, M. C.; Wu, A.; Duggal, K. B.; Rotella, M. E.; Kozlowski, M. C.,
       <i>Machine Learning-Guided Photocatalytic Cross-Coupling of Phenols and Heteroaryl Halides</i>.
       ChemRxiv. 2025, doi:10.26434/chemrxiv-2025-5vcsz.<br>
       This content is a preprint and has not been peer-reviewed."
    )
  )
)

# --------- Server ----------
server <- function(input, output, session) {
  
  # Render molecule images
  observe({
    if (nchar(input$phenol) > 0) {
      output$phenol_img <- renderUI({
        tags$img(src = smiles_to_url(input$phenol), class = "mol_img")
      })
    }
    if (nchar(input$pyridine) > 0) {
      output$pyridine_img <- renderUI({
        tags$img(src = smiles_to_url(input$pyridine), class = "mol_img")
      })
    }
  })
  
  # Prediction
  prediction <- eventReactive(input$predict_btn, {
    
    if (input$phenol == "" || input$pyridine == "")
      return(list(error=TRUE, msg="Both SMILES must be entered."))
    
    if (!grepl("Cl|Br|I", input$pyridine))
      return(list(error=TRUE, msg="Pyridine must contain Cl, Br, or I."))
    
    phenol_input <- safe_rdkit_descriptors(input$phenol, "Phenol_")
    pyridine_input <- safe_rdkit_descriptors(input$pyridine, "Pyridine_")
    
    if (is.null(phenol_input) || is.null(pyridine_input))
      return(list(error=TRUE, msg="Cannot process this molecule. Try a simpler SMILES."))
    
    X_input <- bind_cols(
      as.data.frame(list(phenol_input)),
      as.data.frame(list(pyridine_input))
    )
    
    # Ensure all model columns exist
    for (c in setdiff(trained_model$X_train_names, colnames(X_input)))
      X_input[[c]] <- 0
    
    X_input <- X_input[, trained_model$X_train_names]
    X_input[is.na(X_input)] <- 0
    
    # ML prediction
    p_raw <- predict(trained_model$rf, X_input, type="prob")[,2]
    p_raw <- pmin(pmax(p_raw, 1e-6), 1 - 1e-6)
    
    odds <- p_raw/(1-p_raw) *
      (p_external/(1-p_external)) /
      (trained_model$p_train/(1-trained_model$p_train))
    
    p_corr <- odds/(1+odds)
    
    list(probability = p_corr, error = FALSE)
  })
  
  # Render outputs
  output$reaction_box <- renderUI({
    pred <- prediction()
    if (is.null(pred) || pred$error) return(NULL)
    label <- ifelse(pred$probability >= threshold, "C–C", "S<sub>N</sub>Ar")
    color <- ifelse(pred$probability >= threshold, "#1f3a93", "#800000")
    div(class="output_box", style=paste0("background:", color), HTML(label))
  })
  
  output$prob_box <- renderUI({
    pred <- prediction()
    if (is.null(pred) || pred$error) return(NULL)
    div(id="prob_box", paste0("C–C probability: ", round(pred$probability*100,1), "%"))
  })
  
  output$mix_box <- renderUI({
    pred <- prediction()
    if (is.null(pred) || pred$error) return(NULL)
    if (pred$probability >= 0.25 && pred$probability <= 0.55)
      HTML("A mix of C–C and S<sub>N</sub>Ar may be observed.")
  })
  
  output$validation_box <- renderUI({
    pred <- prediction()
    if (!is.null(pred$error) && pred$error)
      div(style="color:red; font-weight:bold;", pred$msg)
  })
}

# ---------- Run App ----------
shinyApp(ui, server)



