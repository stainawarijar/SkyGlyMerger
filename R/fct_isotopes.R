# This file contains functions that are used for calculating isotopic fine
# structure based on a molecular formula (elemental composition), and for
# collapsing the fine structure into M, M+1, M+2, etc. peaks.

# It uses constants defined in `constants.R`
# (ISOTOPES, PROTON_MASS, ELECTRON_MASS, EXTRA_NEUTRON_LOOKUP)



# Returns full name of an element based on symbol.
symbol_to_element <- function(symbol) {

  symbol_element_map <- list(
    C = "carbon",
    H = "hydrogen",
    O = "oxygen",
    N = "nitrogen",
    S = "sulfur",
    Na = "sodium",
    K = "potassium",
    Fe = "iron"
  )

  return(symbol_element_map[[symbol]])
}



# Converts a molecular formula in character format to a dictionary.
# Example: "C6H12O6" is converted to `list("C" = 6, "H" = 12, "O" = 6)`
parse_molecular_formula <- function(formula) {
  composition <- list()

  matches <- gregexpr("([A-Z][a-z]*)([0-9]*)", formula, perl = TRUE)
  matched_text <- regmatches(formula, matches)[[1]]

  for (match_text in matched_text) {
    symbol <- sub("^([A-Z][a-z]*)([0-9]*)$", "\\1", match_text, perl = TRUE)
    count_text <- sub("^([A-Z][a-z]*)([0-9]*)$", "\\2", match_text, perl = TRUE)

    if (count_text == "") {
      count <- 1
    }
    else {
      count <- as.integer(count_text)
    }

    if (is.null(composition[[symbol]])) {
      composition[[symbol]] <- count
    }
    else {
      composition[[symbol]] <- composition[[symbol]] + count
    }
  }

  return(composition)
}



# Generate all ways to distribute n identical atoms over k isotopes.
# Example: n = 2, k = 3 gives combinations:
# (2, 0, 0), (1, 1, 0), (1, 0, 1), ..., (0, 0, 2)
# This example applies to 2 oxygen atoms, which has 3 stable isotopes.
isotope_count_combis <- function(n, k) {
  if (k == 1) {
    # Base case: only one isotope.
    return(list(c(n)))
  }

  # Store all combinations here.
  combinations <- list()

  # Try every possible count for the first isotope.
  for (first_count in 0:n) {
    # Distribute remaining atoms over the remaining isotopes.
    remaining_combinations <- isotope_count_combis(
      n = n - first_count,
      k = k - 1
    )

    # Add first_count in front of each remaining combination.
    for (remaining_counts in remaining_combinations) {
      new_combination <- c(first_count, remaining_counts)
      combinations[[length(combinations) + 1]] <- new_combination
    }
  }

  return(combinations)
}



# Calculate the multinomial probability for isotope counts.
# For example, for carbon with counts (5, 1):
# probability = 6! / (5! 1!) * P(C12)^5 * P(C13)^1
#
# This is done in log-space, because direct factorials and products can become
# numerically unstable for larger molecules.
multinomial_prob <- function(counts, probs) {
  total <- sum(counts)

  # The gamma function generalizes factorials.
  # `lgamma(n + 1)` is therefore log(n!)
  log_prob <- lgamma(total + 1)

  # Subtract the log-factorials of the isotope counts.
  log_prob <- log_prob - sum(lgamma(counts + 1))

  # Add the probability terms.
  # Only include counts > 0 to avoid problems such as 0 * log(0).
  positive <- counts > 0
  log_prob <- log_prob + sum(counts[positive] * log(probs[positive]))

  return(exp(log_prob))
}



# Calculate the isotopic fine structure for a single element.
element_fine_structure <- function(symbol, atom_count, min_prob) {

  element <- symbol_to_element(symbol)
  isotope_data <- ISOTOPES[[element]]

  isotope_labels <- names(isotope_data)
  isotope_masses <- purrr::map_dbl(
    isotope_labels, function(label) isotope_data[[label]]$mass
  )
  isotope_probs <- purrr::map_dbl(
    isotope_labels, function(label) isotope_data[[label]]$abundance
  )

  pattern = list()

  for (counts in isotope_count_combis(atom_count, length(isotope_labels))) {
    prob = multinomial_prob(counts, isotope_probs)

    # Skip masses with very low probabilities, to reduce computation.
    if (prob <= min_prob) {
      next
    }

    mass <- sum(counts * isotope_masses)

    isotope_counts <- counts[counts > 0]
    names(isotope_counts) <- isotope_labels[counts > 0]

    pattern[[length(pattern) + 1]] <- list(
      mass = mass, prob = prob, isotope_counts = isotope_counts
    )
  }

  return(pattern)
}



# Combine two isotope patterns.
# Every peak from pattern_a is combined with every peak from pattern_b.
# The masses add and the probabilities multiply.
convolve_patterns <- function(pattern_a, pattern_b, min_prob) {
  combined_pattern <- list()
  peak_index <- 1

  for (peak_a in pattern_a) {
    for (peak_b in pattern_b) {

      # When two independent isotope patterns are combined,
      # their probabilities are multiplied.
      prob <- peak_a$prob * peak_b$prob

      # Skip very small peaks.
      if (prob < min_prob) {
        next
      }

      # Start with the isotope counts from peak_a.
      isotope_counts <- peak_a$isotope_counts

      # Add the isotope counts from peak_b.
      for (isotope_label in names(peak_b$isotope_counts)) {
        count <- peak_b$isotope_counts[[isotope_label]]

        if (is.null(isotope_counts[[isotope_label]])) {
          isotope_counts[[isotope_label]] <- count
        }
        else {
          isotope_counts[[isotope_label]] <- (
            isotope_counts[[isotope_label]] + count
          )
        }
      }

      # Store the combined peak.
      combined_pattern[[peak_index]] <- list(
        mass = peak_a$mass + peak_b$mass,
        prob = prob,
        isotope_counts = isotope_counts
      )

      peak_index <- peak_index + 1
    }
  }

  return(combined_pattern)
}



calculate_fine_structure <- function(
    formula,
    charge,
    min_prob=1e-12
) {

  # Determine elemental composition
  composition <- parse_molecular_formula(formula)

  # Start with one empty mass. Each element is convolved into this pattern.
  molecular_pattern = list(
    list(mass = 0.0, prob = 1.0, isotope_counts = list())
  )

  # Loop over elements
  element_symbols <- names(composition)
  for (symbol in element_symbols) {
    atom_count <- composition[[symbol]]

    element_pattern <- element_fine_structure(symbol, atom_count, min_prob)

    molecular_pattern <- convolve_patterns(
      molecular_pattern, element_pattern, min_prob
    )
  }

  # Correct for electron masses based on charge state.
  mass_correction <- as.integer(charge) * ELECTRON_MASS

  for (i in seq_along(molecular_pattern)) {
    molecular_pattern[[i]]$mass <- molecular_pattern[[i]]$mass - mass_correction
  }

  # Sort final fine-structure peaks by mass.
  masses <- vapply(
    molecular_pattern,
    function(peak) peak$mass,
    numeric(1)
  )

  molecular_pattern <- molecular_pattern[order(masses)]

  # Add probabilities relative to most probable peak.
  probs <- vapply(
    molecular_pattern,
    function(peak) peak$prob,
    numeric(1)
  )

  max_prob <- max(probs)

  for (i in seq_along(molecular_pattern)) {
    molecular_pattern[[i]]$relative_prob <- (
      molecular_pattern[[i]]$prob / max_prob
    )
  }

  return(molecular_pattern)
}



# Collapse fine structure to nominal pattern.
collapse_to_nominal_pattern <- function(fine_structure_pattern) {

  grouped_pattern <- list()

  for (peak in fine_structure_pattern) {
    # Determine how many extra neutrons this fine-structure peak has.
    extra_neutrons <- 0

    for (isotope_label in names(peak$isotope_counts)) {
      isotope_count <- peak$isotope_counts[[isotope_label]]
      isotope_extra_neutrons <- EXTRA_NEUTRON_LOOKUP[[isotope_label]]

      extra_neutrons <- extra_neutrons +
        isotope_count * isotope_extra_neutrons
    }

    group_name <- paste0("M+", extra_neutrons)

    if (extra_neutrons == 0) {
      group_name <- "M"
    }

    # If this M+n group does not exist yet, create it.
    if (is.null(grouped_pattern[[group_name]])) {
      grouped_pattern[[group_name]] <- list(
        isotope_group = group_name,
        extra_neutrons = extra_neutrons,
        prob = 0,
        weighted_mass_sum = 0,
        n_fine_structure_peaks = 0
      )
    }

    # Sum probabilities.
    grouped_pattern[[group_name]]$prob <- (
      grouped_pattern[[group_name]]$prob + peak$prob
    )

    # Store sum(mass * probability), so we can calculate the
    # probability-weighted average mass later.
    grouped_pattern[[group_name]]$weighted_mass_sum <- (
      grouped_pattern[[group_name]]$weighted_mass_sum + peak$mass * peak$prob
    )

    grouped_pattern[[group_name]]$n_fine_structure_peaks <- (
      grouped_pattern[[group_name]]$n_fine_structure_peaks + 1
    )
  }

  # Convert weighted mass sums into weighted average masses.
  for (group_name in names(grouped_pattern)) {
    grouped_pattern[[group_name]]$mass <- (
      grouped_pattern[[group_name]]$weighted_mass_sum /
        grouped_pattern[[group_name]]$prob
    )

    grouped_pattern[[group_name]]$weighted_mass_sum <- NULL
  }

  # Sort by M, M+1, M+2, ...
  extra_neutrons <- vapply(
    grouped_pattern,
    function(group) group$extra_neutrons,
    numeric(1)
  )

  grouped_pattern <- grouped_pattern[order(extra_neutrons)]

  # Add probabilities relative to most probable peak.
  probs <- vapply(
    grouped_pattern,
    function(group) group$prob,
    numeric(1)
  )

  max_prob <- max(probs)

  for (i in seq_along(grouped_pattern)) {
    grouped_pattern[[i]]$prob_relative <- grouped_pattern[[i]]$prob / max_prob
  }

  return(grouped_pattern)
}



# Calculate the fine structure for an ion.
# This function adds charge carrier to the molecular formula,
# so the molecular formula is assumed to be for the neutral molecule.
calculate_ion_fine_structure <- function(
    formula,  # Molecular formula (elemental composition)
    charge,  # Integer, may be positive or negative
    carrier = "H"  # Fix to hydrogen for now
) {
  # Add charge carrier to the molecular formula
  charge_composition <- paste0(carrier, as.integer(charge))
  formula_incl_charge <- paste0(formula, charge_composition)

  # Get fine structure pattern for formula including charge carrier.
  ion_fine_pattern <- calculate_fine_structure(
    formula_incl_charge, charge = charge
  )

  return(ion_fine_pattern)
}


