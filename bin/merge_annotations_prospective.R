#!/usr/bin/env Rscript
# Version 2: Handles flexible AF position in VCF FORMAT field
# Supports both GT:GQ:DP:AF:... (AF at position 4) and GT:GQ:DP:AD:AF:... (AF at position 5)

merge_variant_caller_output <- function(Merged_file, Pileup_file, Somatic_file, output_file, occgenes, filtered_output_file) {

  library(dplyr)
  library(tidyr)

  # Helper function to check if file exists and has content
  check_file <- function(file_path, file_name) {
    if (!file.exists(file_path)) {
      print(paste("Warning:", file_name, "file does not exist:", file_path))
      return(FALSE)
    }

    file_size <- file.size(file_path)
    if (file_size == 0) {
      print(paste("Warning:", file_name, "file is empty:", file_path))
      return(FALSE)
    }

    # Try to read first few lines to check if file has content
    tryCatch({
      test_read <- readLines(file_path, n = 2)
      if (length(test_read) < 2) {
        print(paste("Warning:", file_name, "file has no data rows:", file_path))
        return(FALSE)
      }
      return(TRUE)
    }, error = function(e) {
      print(paste("Warning: Error reading", file_name, "file:", e$message))
      return(FALSE)
    })
  }

  # Helper function to dynamically parse FORMAT field and extract AF
  # Handles both scenarios: AF at position 4 or 5
  parse_format_field <- function(df, format_column, caller_prefix = "") {
    # Detect FORMAT structure from first non-NA value
    sample_value <- df[[format_column]][!is.na(df[[format_column]])][1]

    if (is.na(sample_value) || sample_value == "") {
      print(paste("Warning: Cannot detect FORMAT structure for", caller_prefix))
      return(df)
    }

    # Split to count fields
    fields <- strsplit(sample_value, ":")[[1]]
    num_fields <- length(fields)

    print(paste(caller_prefix, "FORMAT detected with", num_fields, "fields:", paste(fields[1:min(5, num_fields)], collapse=":")))

    # Determine parsing strategy based on number of fields
    if (num_fields >= 9) {
      # ClairS-TO format: GT:GQ:DP:AF:AD:AU:CU:GU:TU (9 fields)
      print(paste(caller_prefix, "- Using ClairS-TO format (AF at position 4)"))
      df <- df %>%
        separate(!!sym(format_column),
                 c("GT", paste0(caller_prefix, "_GQ"), paste0(caller_prefix, "_Depth"),
                   paste0(caller_prefix, "_AF"), paste0(caller_prefix, "_AD"),
                   "AU", "CU", "GU", "TU"),
                 sep = ":", fill = "right", extra = "drop")
    } else if (num_fields == 5) {
      # Could be GT:GQ:DP:AF:AD or GT:GQ:DP:AD:AF
      # Check if 4th field looks like AF (decimal < 1) or AD (integer/comma-separated)
      fourth_field <- fields[4]

      # Try to determine if 4th field is AF or AD
      # AF is typically 0.xxxx (decimal < 1)
      # AD is typically integer or comma-separated integers like "6,1"
      is_af_fourth <- FALSE
      if (grepl("^0\\.", fourth_field)) {
        is_af_fourth <- TRUE
      } else {
        # Try numeric conversion safely
        fourth_num <- suppressWarnings(as.numeric(fourth_field))
        if (!is.na(fourth_num) && fourth_num < 1 && fourth_num > 0) {
          is_af_fourth <- TRUE
        }
      }

      if (is_af_fourth) {
        print(paste(caller_prefix, "- Using format GT:GQ:DP:AF:AD (AF at position 4)"))
        df <- df %>%
          separate(!!sym(format_column),
                   c("GT", paste0(caller_prefix, "_GQ"), "Depth", "AF", "AD"),
                   sep = ":", fill = "right", extra = "drop")
      } else {
        print(paste(caller_prefix, "- Using format GT:GQ:DP:AD:AF (AF at position 5)"))
        df <- df %>%
          separate(!!sym(format_column),
                   c("GT", paste0(caller_prefix, "_GQ"), "Depth", "AD", "AF"),
                   sep = ":", fill = "right", extra = "drop")
      }
    } else if (num_fields == 4) {
      # GT:GQ:DP:AF format
      print(paste(caller_prefix, "- Using format GT:GQ:DP:AF (AF at position 4)"))
      df <- df %>%
        separate(!!sym(format_column),
                 c("GT", paste0(caller_prefix, "_GQ"), "Depth", "AF"),
                 sep = ":", fill = "right", extra = "drop")
    } else {
      print(paste("Warning:", caller_prefix, "- Unknown FORMAT structure with", num_fields, "fields"))
      # Default to simple 4-field format
      df <- df %>%
        separate(!!sym(format_column),
                 c("GT", paste0(caller_prefix, "_GQ"), "Depth", "AF"),
                 sep = ":", fill = "right", extra = "drop")
    }

    return(df)
  }

  # Helper function to create empty data frame with expected columns
  create_empty_df <- function(caller_type) {
    if (caller_type == "Merged") {
      return(data.frame(
        Chr = character(),
        Start = character(),
        End = character(),
        Ref = character(),
        Alt = character(),
        Func.refGene = character(),
        Gene.refGene = character(),
        GeneDetail.refGene = character(),
        ExonicFunc.refGene = character(),
        AAChange.refGene = character(),
        Otherinfo1 = character(),
        Otherinfo2 = character(),
        Otherinfo3 = character(),
        Otherinfo4 = character(),
        Otherinfo5 = character(),
        Otherinfo6 = character(),
        Otherinfo7 = character(),
        Otherinfo8 = character(),
        Otherinfo9 = character(),
        Otherinfo10 = character(),
        Otherinfo11 = character(),
        Otherinfo12 = character(),
        Otherinfo13 = character(),
        Otherinfo14 = character(),
        Otherinfo15 = character(),
        Otherinfo16 = character(),
        Otherinfo17 = character(),
        Otherinfo18 = character(),
        Otherinfo19 = character(),
        GT = character(),
        Merged_GQ = character(),
        Depth = character(),
        AD = character(),
        AF = character(),
        callerM = character(),
        stringsAsFactors = FALSE
      ))
    } else if (caller_type == "Pileup") {
      return(data.frame(
        Chr = character(),
        Start = character(),
        End = character(),
        Ref = character(),
        Alt = character(),
        Func.refGene = character(),
        Gene.refGene = character(),
        GeneDetail.refGene = character(),
        ExonicFunc.refGene = character(),
        AAChange.refGene = character(),
        Otherinfo1 = character(),
        Otherinfo2 = character(),
        Otherinfo3 = character(),
        Otherinfo4 = character(),
        Otherinfo5 = character(),
        Otherinfo6 = character(),
        Otherinfo7 = character(),
        Otherinfo8 = character(),
        Otherinfo9 = character(),
        Otherinfo10 = character(),
        Otherinfo11 = character(),
        Otherinfo12 = character(),
        Otherinfo13 = character(),
        Otherinfo14 = character(),
        Otherinfo15 = character(),
        Otherinfo16 = character(),
        Otherinfo17 = character(),
        Otherinfo18 = character(),
        Otherinfo19 = character(),
        GT = character(),
        Pileup_GQ = character(),
        Depth = character(),
        AD = character(),
        AF = character(),
        callerP = character(),
        stringsAsFactors = FALSE
      ))
    } else if (caller_type == "Somatic") {
      return(data.frame(
        Chr = character(),
        Start = character(),
        End = character(),
        Ref = character(),
        Alt = character(),
        Func.refGene = character(),
        Gene.refGene = character(),
        GeneDetail.refGene = character(),
        ExonicFunc.refGene = character(),
        AAChange.refGene = character(),
        Otherinfo1 = character(),
        Otherinfo2 = character(),
        Otherinfo3 = character(),
        Otherinfo4 = character(),
        Otherinfo5 = character(),
        Otherinfo6 = character(),
        Otherinfo7 = character(),
        Otherinfo8 = character(),
        Otherinfo9 = character(),
        Otherinfo10 = character(),
        Otherinfo11 = character(),
        Otherinfo12 = character(),
        Otherinfo13 = character(),
        Otherinfo14 = character(),
        Otherinfo15 = character(),
        Otherinfo16 = character(),
        Otherinfo17 = character(),
        Otherinfo18 = character(),
        Otherinfo19 = character(),
        GT = character(),
        ClairS_GQ = character(),
        ClairS_Depth = character(),
        ClairS_AF = character(),
        ClairS_AD = character(),
        AU = character(),
        CU = character(),
        GU = character(),
        TU = character(),
        callerS = character(),
        stringsAsFactors = FALSE
      ))
    }
  }

  # Check and load Merged file (Clair3)
  has_merged <- check_file(Merged_file, "Merged")
  if (has_merged) {
    tryCatch({
      Merged <- read.delim(Merged_file, header = TRUE, colClasses = c("character"))
      if (nrow(Merged) == 0) {
        print("Warning: Merged file has no data rows")
        has_merged <- FALSE
        Merged <- create_empty_df("Merged")
      } else {
        Merged <- Merged %>% select(7,1:6,8:16,19)
        # Use dynamic parsing for Merged (Clair3)
        Merged <- parse_format_field(Merged, "Otherinfo13", "Merged")
        Merged$callerM <- "Merged"
      }
    }, error = function(e) {
      print(paste("Error processing Merged file:", e$message))
      has_merged <- FALSE
      Merged <- create_empty_df("Merged")
    })
  } else {
    Merged <- create_empty_df("Merged")
  }

  # Check and load Pileup file (Clair3 pileup mode)
  has_pileup <- check_file(Pileup_file, "Pileup")
  if (has_pileup) {
    tryCatch({
      Pileup <- read.delim(Pileup_file, header = TRUE, colClasses = c("character"))
      if (nrow(Pileup) == 0) {
        print("Warning: Pileup file has no data rows")
        has_pileup <- FALSE
        Pileup <- create_empty_df("Pileup")
      } else {
        Pileup <- Pileup %>% select(7,1:6,8:16,19)
        # Use dynamic parsing for Pileup
        Pileup <- parse_format_field(Pileup, "Otherinfo13", "Pileup")
        Pileup$callerP <- "Pileup"
      }
    }, error = function(e) {
      print(paste("Error processing Pileup file:", e$message))
      has_pileup <- FALSE
      Pileup <- create_empty_df("Pileup")
    })
  } else {
    Pileup <- create_empty_df("Pileup")
  }

  # Check and load Somatic file (ClairS-TO)
  has_somatic <- check_file(Somatic_file, "Somatic")
  if (has_somatic) {
    tryCatch({
      Somatic <- read.delim(Somatic_file, header = TRUE, colClasses = c("character"))
      if (nrow(Somatic) == 0) {
        print("Warning: Somatic file has no data rows")
        has_somatic <- FALSE
        Somatic <- create_empty_df("Somatic")
      } else {
        Somatic <- Somatic %>% select(7,1:6,9:16,18)
        # Use dynamic parsing for Somatic (ClairS-TO)
        Somatic <- parse_format_field(Somatic, "Otherinfo10", "ClairS")

        # Select only the required columns (remove AU, CU, GU, TU if present)
        Somatic <- Somatic %>% select(-any_of(c("AU", "CU", "GU", "TU")))
        Somatic$callerS <- "ClairS_TO"
      }
    }, error = function(e) {
      print(paste("Error processing Somatic file:", e$message))
      has_somatic <- FALSE
      Somatic <- create_empty_df("Somatic")
    })
  } else {
    Somatic <- create_empty_df("Somatic")
  }

  # Load OCC genes
  tryCatch({
    occgenes <- readLines(occgenes)
    occgenes <- occgenes[nchar(trimws(occgenes)) > 0]
    print(paste("Loaded", length(occgenes), "OCC genes"))
  }, error = function(e) {
    print(paste("Error loading OCC genes:", e$message))
    occgenes <- character(0)
  })

  # Check if we have any data to process
  if (!has_merged && !has_pileup && !has_somatic) {
    print("Warning: No valid data found in any input files. Creating empty output files.")

    # Create empty output files
    empty_df <- data.frame(
      cosmic100_ID = character(),
      Chr = character(),
      Start = character(),
      End = character(),
      Ref = character(),
      Alt = character(),
      Func.refGene = character(),
      Gene.refGene = character(),
      ExonicFunc.refGene = character(),
      AAChange.refGene = character(),
      Variant_caller = character(),
      GQ = character(),
      Depth = character(),
      stringsAsFactors = FALSE
    )

    write.table(empty_df, file = output_file, row.names = FALSE, sep = "\t", quote = FALSE, na = "0")
    write.table(empty_df, file = filtered_output_file, row.names = FALSE, sep = "\t", quote = FALSE, na = "0")

    print(paste("Empty output files created:", output_file, "and", filtered_output_file))
    return()
  }

  # Join on variant position only (Chr, Start, End, Ref, Alt).
  # Joining on annotation columns (Gene, AAChange, etc.) causes split entries
  # when VEP picks different transcripts for the same position across callers.
  # For ANNOVAR the result is identical since annotations are always consistent.
  join_keys <- c("Chr", "Start", "End", "Ref", "Alt")

  coalesce_join <- function(df, suffix1 = ".x", suffix2 = ".y") {
    dup_base <- sub(paste0("\\", suffix1, "$"), "",
                    grep(paste0("\\", suffix1, "$"), names(df), value = TRUE))
    for (col in dup_base) {
      c1 <- paste0(col, suffix1); c2 <- paste0(col, suffix2)
      df[[col]] <- coalesce(df[[c1]], df[[c2]])
      df[[c1]] <- NULL; df[[c2]] <- NULL
    }
    df
  }

  # Merge files - handle empty data frames
  if (has_pileup && has_merged) {
    All_calls <- full_join(Pileup, Merged, by = join_keys, suffix = c(".x", ".y")) %>%
      coalesce_join() %>%
      unite(Variant_caller, callerP, callerM, sep = ", ", na.rm = TRUE)
  } else if (has_pileup) {
    All_calls <- Pileup %>%
      mutate(Variant_caller = callerP)
  } else if (has_merged) {
    All_calls <- Merged %>%
      mutate(Variant_caller = callerM)
  } else {
    All_calls <- data.frame()
  }

  # Only join with Somatic data if it exists and has data
  if (has_somatic && nrow(Somatic) > 0) {
    if (nrow(All_calls) > 0) {
      All_calls <- All_calls %>%
        full_join(Somatic, by = join_keys, suffix = c(".x", ".y")) %>%
        coalesce_join() %>%
        unite(Variant_caller, Variant_caller, callerS, sep = ", ", na.rm = TRUE)
    } else {
      All_calls <- Somatic %>%
        mutate(Variant_caller = callerS)
    }
  }

  # Handle case where no data was merged
  if (nrow(All_calls) == 0) {
    print("Warning: No data after merging. Creating empty output files.")
    empty_df <- data.frame(
      cosmic100_ID = character(),
      Chr = character(),
      Start = character(),
      End = character(),
      Ref = character(),
      Alt = character(),
      Func.refGene = character(),
      Gene.refGene = character(),
      ExonicFunc.refGene = character(),
      AAChange.refGene = character(),
      Variant_caller = character(),
      GQ = character(),
      Depth = character(),
      stringsAsFactors = FALSE
    )

    write.table(empty_df, file = output_file, row.names = FALSE, sep = "\t", quote = FALSE, na = "0")
    write.table(empty_df, file = filtered_output_file, row.names = FALSE, sep = "\t", quote = FALSE, na = "0")

    print(paste("Empty output files created:", output_file, "and", filtered_output_file))
    return()
  }

  # Unite GQ values based on availability of data
  if (has_somatic && nrow(Somatic) > 0) {
    All_calls <- All_calls %>%
      unite(GQ, any_of(c("Pileup_GQ", "Merged_GQ")), sep = ",", na.rm = TRUE, remove = FALSE) %>%
      unite(GQ, GQ, any_of("ClairS_GQ"), sep = ",", na.rm = TRUE)
  } else {
    # If no somatic data, do not unite ClairS_GQ
    All_calls <- All_calls %>%
      unite(GQ, any_of(c("Pileup_GQ", "Merged_GQ")), sep = ",", na.rm = TRUE)
  }

  # Handle COSMIC100 column if it exists
  if ("COSMIC100" %in% colnames(All_calls)) {
    All_calls <- All_calls %>%
      separate(COSMIC100, c("cosmic100_ID", "y"), sep = ";", remove = FALSE, fill = "right") %>%
      mutate(cosmic100_ID = gsub("ID=", "", cosmic100_ID))
  } else {
    print("COSMIC100 column is missing, skipping separation and mutation steps for this column.")
    All_calls$cosmic100_ID <- ""
  }

  # Continue with the rest of the pipeline - select relevant columns
  # Handle ClairS-TO only variants: coalesce ClairS-specific columns with general columns
  if ("ClairS_AF" %in% colnames(All_calls) && !"AF" %in% colnames(All_calls)) {
    All_calls$AF <- All_calls$ClairS_AF
  } else if ("ClairS_AF" %in% colnames(All_calls) && "AF" %in% colnames(All_calls)) {
    All_calls <- All_calls %>% mutate(AF = coalesce(AF, ClairS_AF))
  }

  if ("ClairS_Depth" %in% colnames(All_calls) && !"Depth" %in% colnames(All_calls)) {
    All_calls$Depth <- All_calls$ClairS_Depth
  } else if ("ClairS_Depth" %in% colnames(All_calls) && "Depth" %in% colnames(All_calls)) {
    All_calls <- All_calls %>% mutate(Depth = coalesce(Depth, ClairS_Depth))
  }

  if ("ClairS_AD" %in% colnames(All_calls) && !"AD" %in% colnames(All_calls)) {
    All_calls$AD <- All_calls$ClairS_AD
  } else if ("ClairS_AD" %in% colnames(All_calls) && "AD" %in% colnames(All_calls)) {
    All_calls <- All_calls %>% mutate(AD = coalesce(AD, ClairS_AD))
  }

  # Select only user-specified columns in the requested order
  required_columns <- c("ClairS_Depth", "Gene.refGene", "Chr", "Start", "End",
                       "Ref", "Alt", "Func.refGene", "ExonicFunc.refGene",
                       "AAChange.refGene", "CLNSIG", "cosmic100_ID", "GQ",
                       "Depth", "AD", "GT", "AF", "Variant_caller")

  # Select only columns that exist in the dataframe, maintain order
  All_calls <- All_calls %>% select(any_of(required_columns))

  # Rename cosmic100_ID to COSMIC100 if it exists
  if ("cosmic100_ID" %in% colnames(All_calls)) {
    All_calls <- All_calls %>% rename(COSMIC100 = cosmic100_ID)
  }

  # Filter by OCC genes if we have any
  if (length(occgenes) > 0) {
    All_calls <- All_calls %>% filter(Gene.refGene %in% occgenes)
    print(paste("Filtered to", length(unique(All_calls$Gene.refGene)), "OCC genes"))
  } else {
    print("No OCC genes loaded, skipping gene filtering")
  }

  # Create filtered version for additional output
  All_calls_filtered <- All_calls %>%
    # Filter 1: Remove "Pileup only" from Variant_caller
    filter(!grepl("^Pileup$", Variant_caller)) %>%
    # Filter 2: Remove rows with depth below 10 (only if Depth column exists and has numeric values)
    filter(is.na(Depth) | Depth == "" | as.numeric(Depth) >= 10) %>%
    # Filter 3: Remove specific TERT variant row
    filter(!(Gene.refGene == "TERT" & Chr == "chr5" & Start == "1295957" & End == "1295957"))

  # Write original output
  All_calls[] <- lapply(All_calls, function(x) {
    if (is.factor(x)) {
      x <- as.character(x)
    }
    # Replace empty strings with "0"
    x[x == ""] <- "0"
    return(x)
  })

  # Write the table, replacing any NA values with "0"
  write.table(All_calls,
              file = output_file,
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = "0")

  # Write filtered output
  All_calls_filtered[] <- lapply(All_calls_filtered, function(x) {
    if (is.factor(x)) {
      x <- as.character(x)
    }
    # Replace empty strings with "0"
    x[x == ""] <- "0"
    return(x)
  })

  write.table(All_calls_filtered,
              file = filtered_output_file,
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = "0")

  print(paste("Original output written to:", output_file))
  print(paste("Filtered output written to:", filtered_output_file))
  print(paste("Original rows:", nrow(All_calls)))
  print(paste("Filtered rows:", nrow(All_calls_filtered)))
}

# Parsing arguments passed from the Nextflow process
args <- commandArgs(trailingOnly = TRUE)
Merged_file <- args[1]
Pileup_file <- args[2]
Somatic_file <- args[3]
output_file <- args[4]
occgenes <- args[5]
filtered_output_file <- args[6]

# Call the function with provided arguments
merge_variant_caller_output(Merged_file, Pileup_file, Somatic_file, output_file, occgenes, filtered_output_file)
