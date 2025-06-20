library(Biostrings)
setwd("~/Documents/SignalNoiseBiasPaper/BirdAnalysisRedo/Fasta")


extract_model_params <- function(iqtree_file, output_file) {
  lines <- readLines(iqtree_file)

  # --- Rate parameters ---
  rate_labels <- c("A-C", "A-G", "A-T", "C-G", "C-T", "G-T")
  rate_map <- c(
    "A-C" = "rACp", "A-G" = "rAGp", "A-T" = "rATp",
    "C-G" = "rCGp", "C-T" = "rTCp", "G-T" = "rTGp"
  )

  rate_entries <- unlist(lapply(rate_labels, function(label) {
    match_line <- grep(paste0("^\\s*", label, ":"), lines, value = TRUE)
    if (length(match_line) > 0) {
      val <- as.numeric(trimws(strsplit(match_line, ":")[[1]][2]))
      paste0(rate_map[[label]], "=", format(val, digits=6), ";")
    } else {
      NULL
    }
  }))

  # --- State frequencies ---
  freq_line_idx <- grep("State frequencies:", lines)
  freq_entries <- c()

  if (length(freq_line_idx) > 0) {
    if (grepl("equal frequencies", lines[freq_line_idx])) {
      freqs <- rep(0.25, 4)
      names(freqs) <- c("piAp", "piCp", "piGp", "piTp")
    } else {
      freq_lines <- lines[(freq_line_idx + 1):(freq_line_idx + 5)]
      freqs <- unlist(lapply(freq_lines, function(line) {
        matches <- regmatches(line, regexec("pi\\((.)\\)\\s*=\\s*([0-9.]+)", line))[[1]]
        if (length(matches) == 3) {
          base <- matches[2]
          val <- as.numeric(matches[3])
          setNames(val, paste0("pi", base, "p"))
        } else {
          NULL
        }
      }))
    }
    freq_entries <- paste0(names(freqs), "=", format(freqs, digits=6), ";")
  }

  # --- Write to file ---
  writeLines(c(rate_entries, freq_entries), con = output_file)
  message("✅ Written model parameters to: ", output_file)
}


calculate_base_frequencies <- function(fasta_file, taxon1, taxon2, taxon3, taxon4, output_file) {
  # Read the alignment
  aln <- readDNAStringSet(fasta_file)
  aln_names <- names(aln)

  # Helper: match a group of candidate names to the first matching alignment name
  match_taxon <- function(taxon_group) {
    for (candidate in taxon_group) {
      match_idx <- grep(candidate, aln_names, ignore.case = TRUE)
      if (length(match_idx) > 0) {
        return(match_idx[1])
      }
    }
    warning(sprintf("⚠️ None of the candidates in %s were found in alignment.", paste(taxon_group, collapse = ", ")))
    return(NA)
  }

  # Match each taxon group
  matched_indices <- c(
    match_taxon(taxon1),
    match_taxon(taxon2),
    match_taxon(taxon3),
    match_taxon(taxon4)
  )

  if (any(is.na(matched_indices))) {
    stop("❌ One or more taxon groups could not be matched. Check the FASTA and taxon lists.")
  }

  # Function to compute base frequencies
  compute_freqs <- function(seq) {
    seq <- toupper(as.character(seq))
    seq <- gsub("-", "", seq)
    bases <- table(strsplit(seq, "")[[1]])
    total <- sum(bases[c("A", "C", "G", "T")], na.rm = TRUE)
    freqs <- sapply(c("T", "C", "A", "G"), function(b) {
      if (!is.na(bases[b])) bases[b] / total else 0
    })
    names(freqs) <- c("T", "C", "A", "G")
    return(freqs)
  }

  output_lines <- c()
  for (i in seq_along(matched_indices)) {
    idx <- matched_indices[i]
    seq <- aln[[idx]]
    freqs <- compute_freqs(seq)

    output_lines <- c(output_lines,
                      sprintf("piT%d=%.9f;", i, freqs["T"]),
                      sprintf("piC%d=%.9f;", i, freqs["C"]),
                      sprintf("piA%d=%.9f;", i, freqs["A"]),
                      sprintf("piG%d=%.9f;", i, freqs["G"]))
  }

  write(output_lines, file = output_file, append = TRUE, sep = "\n")
  message("✅ Base frequencies for selected taxa written to: ", output_file)
}

# Get all model files (edit the pattern if needed)

model_files <- list.files(
  pattern = "\\.iqtree$",
  full.names = TRUE
)

# Loop through each model file
for (mdl_file in model_files) {
  # Strip extension to get the base name
  mdl_base <- tools::file_path_sans_ext(basename(mdl_file))
  
  # Construct full paths
  outputmodel_file <- file.path(
    "~/Documents/SignalNoiseBiasPaper/BirdAnalysisRedo/Fasta",
    paste0(mdl_base, ".params.txt")
  )
  fasta_file <- file.path(
    "~/Documents/SignalNoiseBiasPaper/BirdAnalysisRedo/Fasta",mdl_base
  )

  cat("🔁 Processing:", mdl_base, "\n")

  # Run the model extraction
  extract_model_params(paste0(mdl_base,".iqtree"), outputmodel_file)

  # Then calculate base frequencies for the matching alignment
  # Try catch to not stop the loop
  tryCatch({
    calculate_base_frequencies(
      fasta_file = fasta_file,
      taxon1 = c("I0369_848_Opisthocomidae_Opisthocomus_hoazin"),
      taxon2 = c("I0202_1066_Cathartidae_Cathartes_burrovianus","I0276_5411_Cathartidae_Vultur_gryphus"),
      taxon3 = c("I0361_9545_Sagittariidae_Sagittarius_serpentarius" ,"I0450_81668_Pandionidae_Pandion_haliaetus","I0385_9309_Accipitridae_Elanus_leucurus","I0268_3790_Accipitridae_Accipiter_superciliosus","I0268_3790_Accipitridae_Accipiter_superciliosus","I0185_86_Accipitridae_Buteo_jamaicensis"),
      taxon4 = c("I0208_1208_Strigidae_Strix_varia", "I0452_88784_Tytonidae_Tyto_alba", "I0413_67124_Coliidae_Colius_indicus","I0410_52764_Coliidae_Colius_colius","I0367_449184_Leptosomidae_Leptosomus_discolor",
"I0190_295_Trogonidae_Apaloderma_aequatoriale","I0187_115_Trogonidae_Trogon_elegans","I0199_955_Meropidae_Merops_muelleri",
"I0364_384759_Brachypteraciidae_Atelornis_pittoides","I0284_6031_Coraciidae_Coracias_cyanogaster","I0234_1994_Todidae_Todus_mexicanus",
"I0186_111_Momotidae_Momotus_lessonii","I0272_5303_Alcedinidae_Chloroceryle_inda","I0194_495_Alcedinidae_Alcedo_quadribrachys",
"I0273_5315_Galbulidae_Galbula_dea","I0274_5339_Bucconidae_Chelidoptera_tenebrosa","I0267_3690_Bucconidae_Bucco_capensis",
"I0191_344_Indicatoridae_Indicator_exilis","I0213_1334_Picidae_Jynx_torquilla","I0259_3344_Picidae_Picus_canus",
"I0397_17790_Ramphastidae_Megalaima_chrysopogon","I0195_522_Ramphastidae_Buccanodon_duchaillui","I0266_3617_Ramphastidae_Capito_niger",
"I0286_6061_Ramphastidae_Ramphastos_ambiguus","I2828_3250_Bucerotidae_Bucorvus_leadbeateri","I0198_937_Bucerotidae_Tockus_camurus",
"I0392_15430_Phoeniculidae_Phoeniculus_purpureus","I0233_1920_Upupidae_Upupa_epops","I0208_1208_Strigidae_Strix_varia",
"I0452_88784_Tytonidae_Tyto_alba"),
      output_file = outputmodel_file
    )
  },
  error = function(e) {
    message("⚠️ Skipping ", mdl_base, " due to missing taxon: ", e$message)
  })
}

