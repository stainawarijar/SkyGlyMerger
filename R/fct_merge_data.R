# This file contains functions used in the module `mod_merge_data.R`


# TODO: documentation
read_skyline_csv <- function(file_path) {
  if (grepl(";", readLines(file_path, n = 1))) {
    read.csv(file_path, header = TRUE, sep = ";")
  }
  else {
    read.csv(file_path, header = TRUE, sep = ",")
  }
}


# TODO: documentation
rename_skyline_isomers <- function(data_raw) {
  # Look for isomers in the glycan compositions
  data <- data_raw |>
    dplyr::group_by(
      `Protein.Name`, `Peptide`, `Precursor`, `Precursor.Charge`
    ) |>
    dplyr::mutate(n = dplyr::n(), .after = `Precursor.Charge`) |>
    dplyr::ungroup()

  # n == 1 implies unique glycan composition
  data_unique <- data |>
    dplyr::filter(n == 1)

  # n > 1 implies presence of isomers
  data_isomers <- data |>
    dplyr::filter(n > 1) |>
    dplyr::group_by(
      `Protein.Name`, `Peptide`, `Precursor`, `Precursor.Charge`
    ) |>
    dplyr::mutate(glycan_unique = make.unique(Peptide)) |>
    dplyr::ungroup() |>
    # Instead of ".1", ".2", etc at the end of duplicates, add "_a","_b", to
    # the ends of all isomers, including the first one.
    dplyr::mutate(
      Peptide = dplyr::case_when(
        endsWith(glycan_unique, ".1") ~ paste0(Peptide, "_b"),
        endsWith(glycan_unique, ".2") ~ paste0(Peptide, "_c"),
        endsWith(glycan_unique, ".3") ~ paste0(Peptide, "_d"),
        endsWith(glycan_unique, ".4") ~ paste0(Peptide, "_e"),
        endsWith(glycan_unique, ".5") ~ paste0(Peptide, "_f"),
        endsWith(glycan_unique, ".6") ~ paste0(Peptide, "_g"),
        endsWith(glycan_unique, ".7") ~ paste0(Peptide, "_h"),
        endsWith(glycan_unique, ".8") ~ paste0(Peptide, "_i"),
        endsWith(glycan_unique, ".9") ~ paste0(Peptide, "_j"),
        .default = paste0(Peptide, "_a")
      )
    ) |>
    dplyr::select(-glycan_unique)

  # Combine the data again
  data_renamed <- dplyr::bind_rows(data_unique, data_isomers) |>
    dplyr::select(-n)

  return(data_renamed)
}



# TODO: documentation
transform_skyline_data <- function(data_renamed) {
  # Assume variable columns for samples start after `Precursor.Charge`
  # Everything after that point belongs to one sample-variable combination such
  # as "L20252001084c Best Retention Time" or "L20252001084c Min Start Time".
  # These column names will be split into:
  #   sample   = L20252001084c
  #   variable = Best.Retention.Time / Min.Start.Time / ...
  start_idx = which(colnames(data_renamed) == "Precursor.Charge") + 1
  variable_cols = colnames(data_renamed)[
    start_idx:length(colnames(data_renamed))
  ]

  data_long <- data_renamed |>
    dplyr::mutate(
      # Convert Skyline measurement columns to numeric.
      # Notes:
      # - "#N/A" should become NA
      # - some Skyline exports prefix scientific notation with "*"
      #   (for example "*2.4246E+7"), which must be removed before conversion
      dplyr::across(
        .cols = all_of(variable_cols),
        .fns = ~ .x |>
          # Force to character
          as.character() |>
          # Turn "#N/A" into `NA`
          dplyr::na_if("#N/A") |>
          # Skyline can prefix scientific notation with "*"
          stringr::str_remove("^\\*") |>
          # Back to numeric
          as.numeric()
      )
    ) |>
    tidyr::pivot_longer(tidyr::all_of(variable_cols), names_to = "sample_variable") |>
    tidyr::extract(
      # Split the original Skyline column name into sample ID and measurement type
      col   = sample_variable,
      into  = c("sample", "variable"),
      regex = paste0(
        "^(.+)\\.(",
        "Best\\.Retention\\.Time|",
        "Total\\.Area\\.MS1|",
        "Isotope\\.Dot\\.Product|",
        "Average\\.Mass\\.Error\\.PPM|",
        "Normalized\\.Area|",
        "Min\\.Start\\.Time|",
        "Max\\.End\\.Time",
        ")$"
      ),
      remove = TRUE
    )

  # Reshape back to wide format, but now with one row per analyte-sample
  # combination and one column per measurement variable.
  data_wide <- data_long |>
    tidyr::pivot_wider(names_from = "variable", values_from = "value")

  return(data_wide)
}


# TODO: documentation
load_skyline_data <- function(file_path) {
  data_raw <- read_skyline_csv(file_path)
  data_renamed <- rename_skyline_isomers(data_raw)
  data_transformed <- transform_skyline_data(data_renamed)
  return(data_transformed)
}


# TODO: documentation
# `files`: named list containing original filenames (names) and those stored
# in memory (values)
load_glycounter_data <- function(files) {
  # Read and process all OxoSignal files
  purrr::imap_dfr(files, function(file, name) {
    read.delim(
      file = file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ) |>
      dplyr::select(
        # Required output variables
        ScanNumber,
        RetentionTime,
        PrecursorMZ,
        LikelyGlycoSpectrum,
        DissociationType,
        # Fragment columns are named like "204.0867, HexNAc"
        tidyr::matches("^\\d+\\.\\d+,\\s")
      ) |>
      dplyr::mutate(
        # Recover the sample name from the original filename
        sample = sub("_?GlyCounter.*$", "", basename(name)),
        # Convert "True"/"False" text values to logical values
        LikelyGlycoSpectrum = tolower(LikelyGlycoSpectrum) == "true",
        # Sum all fragment-ion intensities within each scan.
        # The individual fragment columns are retained as well.
        fragment_sum = rowSums(
          dplyr::pick(tidyr::matches("^\\d+\\.\\d+,\\s")),
          na.rm = TRUE
        )
      ) |>
      dplyr::relocate(sample, .before = ScanNumber) |>
      dplyr::relocate(fragment_sum, .after = DissociationType)
  })
}

