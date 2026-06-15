# This file contains functions used in the module `mod_merge_data.R`
# It also calls functions in `fct_isotopes.R`


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
rename_skyline_isomers <- function(
    data_raw
    # TODO Adjustable column names for Protein Name, Charge etc
  ) {
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
# TODO: Adjustable column names for Precursor.Charge etc...
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
    tidyr::pivot_longer(
      tidyr::all_of(variable_cols), names_to = "sample_variable"
    ) |>
    tidyr::extract(
      # Split original Skyline column name into sample ID and measurement type
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
calculate_skyline_isotopic_patterns <- function(
    skyline_data,
    molecular_formula_col = "Molecule.Formula",
    charge_col = "Precursor.Charge",
    charge_carrier = "H"  # Fix for now
  ) {

  # Extract unique composition + charge combinations
  compositions <- skyline_data |>
    dplyr::select(dplyr::all_of(c(molecular_formula_col, charge_col))) |>
    dplyr::distinct() |>
    # Ensure charge is integer.
    dplyr::mutate(dplyr::across(dplyr::all_of(charge_col), as.integer))

  # Get nominal isotopic pattern(M, M+1, M+2, ...) for each ion
  patterns <- list()
  for (i in seq_len(nrow(compositions))) {
    formula <- compositions[[molecular_formula_col]][[i]]
    charge <- compositions[[charge_col]][[i]]
    charge_name <- as.character(charge)

    fine_structure_pattern <- calculate_ion_fine_structure(
      formula = formula, charge = charge, carrier = charge_carrier
    )

    nominal_pattern <- collapse_to_nominal_pattern(fine_structure_pattern)

    patterns[[formula]][[charge_name]] <- nominal_pattern
  }

  return(patterns)
}



# TODO: Documentation
extract_top_isotopic_mz <- function(
    isotopic_patterns,
    n_peaks = 3,
    molecular_formula_col = "Molecule.Formula",
    charge_col = "Precursor.Charge"
  ) {
  purrr::imap_dfr(
    isotopic_patterns, function(charge_list, formula) {
      purrr::imap_dfr(
        charge_list, function(peak_list, charge) {

          charge_int <- as.integer(charge)

          peak_df <- purrr::map_dfr(peak_list, tibble::as_tibble)

          peak_df |>
            dplyr::arrange(dplyr::desc(prob)) |>
            dplyr::slice_head(n = n_peaks) |>
            dplyr::mutate(
              rank = dplyr::row_number(),
              !!molecular_formula_col := formula,
              !!charge_col := charge_int,
              # For negative charges, m/z is still reported as positive value
              mz = abs(mass / charge_int)
            ) |>
            dplyr::select(
              dplyr::all_of(c(molecular_formula_col, charge_col)),
              rank,
              mz
            )
        }
      )
    }
  ) |>
    tidyr::pivot_wider(
      names_from = rank, values_from = mz,
      names_prefix = "mz_prob"
    )
}



# TODO: Documentation
extract_isotopic_mz_candidates <- function(
    isotopic_patterns,
    n_peaks = 3,
    min_relative_prob = NULL,
    molecular_formula_col = "Molecule.Formula",
    charge_col = "Precursor.Charge"
  ) {

  candidates <- purrr::imap_dfr(
    isotopic_patterns, function(charge_list, formula) {
      purrr::imap_dfr(
        charge_list, function(peak_list, charge) {

          charge_int <- as.integer(charge)

          peak_df <- purrr::map_dfr(peak_list, tibble::as_tibble)

          peak_df |>
            dplyr::arrange(dplyr::desc(prob)) |>
            dplyr::mutate(
              isotope_rank = dplyr::row_number(),
              !!molecular_formula_col := formula,
              !!charge_col := charge_int,
              isotope_mz = abs(mass / charge_int),
              isotope_prob = prob,
              isotope_prob_relative = prob_relative
            ) |>
            dplyr::select(
              dplyr::all_of(c(molecular_formula_col, charge_col)),
              isotope_rank,
              isotope_group,
              extra_neutrons,
              isotope_mz,
              isotope_prob,
              isotope_prob_relative,
              n_fine_structure_peaks
            )
        }
      )
    }
  )

  if (!is.null(min_relative_prob)) {
    candidates <- candidates |>
      dplyr::filter(isotope_prob_relative >= min_relative_prob)
  }

  slice <- candidates |>
    dplyr::group_by(
      dplyr::across(dplyr::all_of(c(molecular_formula_col, charge_col)))
    ) |>
    dplyr::slice_head(n = n_peaks) |>
    dplyr::ungroup()

  return(slice)
}



# TODO: Documentation
expand_skyline_isotope_candidates <- function(
    skyline_prepped,
    isotope_mz_candidates,
    molecular_formula_col = "Molecule.Formula",
    charge_col = "Precursor.Charge"
) {
  skyline_prepped |>
    dplyr::left_join(
      isotope_mz_candidates,
      by = c(molecular_formula_col, charge_col),
      relationship = "many-to-many"
    ) |>
    dplyr::mutate(
      isotope_mz_min = isotope_mz * (1 - ppm_tolerance * 1e-6),
      isotope_mz_max = isotope_mz * (1 + ppm_tolerance * 1e-6)
    )
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


# TODO: documentation
extract_fragment_cols <- function(glycounter_data) {
  names <- colnames(glycounter_data)
  fragment_cols <- names[stringr::str_detect(names, "^\\d+\\.\\d+,\\s")]
  return(fragment_cols)
}


# TODO: documentation
prepare_skyline_data <- function(
    skyline_data,
    ppm_tolerance = 10,
    charge_col = "Precursor.Charge"
  ) {
  skyline_data |>
    dplyr::mutate(
      # Add a stable row identifier so GlyCounter summaries can be calculated
      # per Skyline row and safely joined back afterwards.
      skyline_row_id = dplyr::row_number(),
      # Add ppm tolerance, and charge as integer
      ppm_tolerance = ppm_tolerance,
      `Precursor.Charge` = as.integer(.data[["Precursor.Charge"]])  # TODO: adjustable colname
    )
}



# TODO: documentation
extract_glycounter_candidates <- function(
    skyline_isotope_candidates,
    glycounter_data
) {
  skyline_isotope_candidates |>
    dplyr::select(
      skyline_row_id,
      sample,
      `Best.Retention.Time`,
      Min.Start.Time,
      Max.End.Time,
      isotope_rank,
      isotope_group,
      extra_neutrons,
      isotope_mz,
      isotope_prob,
      isotope_prob_relative,
      isotope_mz_min,
      isotope_mz_max
    ) |>
    dplyr::left_join(
      glycounter_data,
      by = "sample",
      relationship = "many-to-many"
    ) |>
    dplyr::filter(
      dplyr::between(PrecursorMZ, isotope_mz_min, isotope_mz_max),
      dplyr::between(RetentionTime, Min.Start.Time, Max.End.Time)
    ) |>
    dplyr::mutate(
      glycounter_mz_error_ppm = (
        1e6 * abs(PrecursorMZ - isotope_mz) / isotope_mz
      )
    ) |>
    dplyr::group_by(skyline_row_id, sample, ScanNumber) |>
    dplyr::slice_min(
      order_by = glycounter_mz_error_ppm,
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::ungroup()
}



# TODO: documentation
summarize_glycounter_data <- function(
    glycounter_candidates,
    fragment_cols,
    relative_abundances = TRUE
) {
  # Aggregate the scan-level GlyCounter data to one summary per Skyline row.
  # This is where repeated fragmentation events are collapsed into one total
  # signal profile for each sample/analyte combination.
  #
  # Note: GlyCounter precursor m/z values may match M, M+1, M+2, etc. rather
  # than only the monoisotopic precursor m/z. Therefore, isotope-match metadata
  # is retained here for traceability.
  summary <- glycounter_candidates |>
    dplyr::group_by(skyline_row_id) |>
    dplyr::summarize(
      # Number of GlyCounter scans contributing to this Skyline row.
      glycounter_scan_count = dplyr::n(),

      # Keep contributing scan numbers for traceability.
      glycounter_scan_numbers = stringr::str_c(
        sort(unique(ScanNumber)), collapse = ";"
      ),

      # DissociationType can vary across contributing scans.
      glycounter_dissociation_types = stringr::str_c(
        sort(unique(DissociationType)), collapse = ";"
      ),

      # Keep track of which isotope-envelope peaks were matched.
      # These can be M, M+1, M+2, etc.
      glycounter_matched_isotope_groups = stringr::str_c(
        sort(unique(isotope_group)), collapse = ";"
      ),

      # Keep track of the probability ranks of the matched isotope peaks.
      # Rank 1 is the most probable isotope peak for that molecular formula
      # and charge state, not necessarily the monoisotopic peak.
      glycounter_matched_isotope_ranks = stringr::str_c(
        sort(unique(isotope_rank)), collapse = ";"
      ),

      # Count how many contributing scans were flagged as likely glyco spectra.
      glycounter_likely_glyco_count = sum(LikelyGlycoSpectrum, na.rm = TRUE),

      # Summaries of matched retention times.
      glycounter_retention_time_mean = mean(RetentionTime, na.rm = TRUE),
      glycounter_retention_time_min = min(RetentionTime, na.rm = TRUE),
      glycounter_retention_time_max = max(RetentionTime, na.rm = TRUE),

      # Summaries of observed GlyCounter precursor m/z values.
      glycounter_precursor_mz_mean = mean(PrecursorMZ, na.rm = TRUE),

      # Summaries of the theoretical isotope m/z values that were matched.
      # These are candidate isotope-envelope m/z values, not necessarily
      # monoisotopic m/z values.
      glycounter_matched_isotope_mz_mean = mean(isotope_mz, na.rm = TRUE),

      # Precursor m/z error is calculated relative to the matched isotope m/z,
      # not relative to the Skyline monoisotopic m/z.
      glycounter_mz_error_ppm_mean = mean(glycounter_mz_error_ppm, na.rm = TRUE),
      glycounter_mz_error_ppm_min = min(glycounter_mz_error_ppm, na.rm = TRUE),

      # Sum total fragment signal across all matched scans.
      fragment_sum = sum(fragment_sum, na.rm = TRUE),

      # Sum every individual fragment-ion column across matched scans so each
      # Skyline row ends up with one total intensity per fragment ion.
      dplyr::across(dplyr::all_of(fragment_cols), ~ sum(.x, na.rm = TRUE)),

      .groups = "drop"
    )

  # Optional: report relative abundances for fragments.
  # This converts each fragment-ion intensity to a percentage of the total
  # fragment signal for that Skyline row.
  if (relative_abundances) {
    summary <- summary |>
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(fragment_cols),
          ~ dplyr::if_else(fragment_sum > 0, .x / fragment_sum * 100, NA_real_)
        )
      )
  }

  return(summary)
}



# TODO: documentation
merge_data <- function(
    skyline_prepped,
    glycounter_summary
  ) {
  # Add the GlyCounter summaries back to the full Skyline table.
  # Because this is a left join starting from Skyline, the number of rows
  # remains identical to the original Skyline analyte-sample table.
  skyline_prepped |>
    dplyr::left_join(glycounter_summary, by = "skyline_row_id") |>
    # Remove helper columns only needed for matching.
    dplyr::select(-skyline_row_id, -ppm_tolerance)
}

