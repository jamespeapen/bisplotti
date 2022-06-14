#' Create metagene intervals
#'
#' "metagene" intervals associated with a given interval (input, in bed format) is
#' the set of intervals positioned with respect to the input interval. This allows
#' for metagene/waterfall plots to be made.
#'
#' @param bedfile         Input bed file
#' @param middle          Use middle point of each probe window
#' @param outer           Use outer point of each probe window (w.r.t. to the target region)
#' @param collapse        Collapse internal interval to middle - all probe windows will now be relative to this location
#' @param flankbygene     Allow the size of steps to vary according to the region of interest length
#' @param flanktoneighbor Set step size to X bases (default: 100)
#' @param flankstep       Set step size to X bases (default: 100)
#' @param flanknumber     Size of step is dependent on the closest region of interest to current region
#' @param expansion       Number of bases to expand window, expands window in both directions (default: 0)
#' @param numinternal     Number of points to sample in the region of interest, \code{middle} ignores this (default: 30)
#' @param fold            Use the same index for intervals on both sides of the target region (usually used when strand is irrelevant)
#' @param strand          The column which contains strand information in input file - if NULL then ignore strand (default: NULL),
#' @param ignoreend       Ignore the end of the input interval
#'
#' @return
#'
#' @examples
wzmetagene <- function(bedfile,
                       middle = F,
                       outer = F,
                       collapse = F,
                       flankbygene = F,
                       flanktoneighbor = F,
                       flankstep = 100,
                       flanknumber = 30,
                       expansion = 0,
                       numinternal = 30,
                       fold = F,
                       strand = NULL,
                       ignoreend = F) {

  opts <- list(
    middle = middle,
    outer = outer,
    collapse = collapse,
    flankbygene = flankbygene,
    flanktoneighbor = flanktoneighbor,
    flankstep = flankstep,
    flanknumber = flanknumber,
    expansion = expansion,
    numinternal = numinternal,
    fold = fold,
    strand = strand,
    ignoreend = ignoreend
  )

  if (opts$collapse) {
    if (opts$outer) {
      opts$numinternal <- 1
    } else {
      opts$numinternal <- 0
    }
  }

  r0 <- NULL
  r <- NULL
  r2 <- NULL

  header <- c("chr", "start", "end", "strand")
  input_bed <- read.table(bedfile, header = F)
  if (ncol(input_bed) == 3) {
    colnames(input_bed) <- header[1:3]
    opts$strand <- NULL
  } else {
      colnames(input_bed) <- header
  }

  # preserve original start and end
  input_bed$og_start <- input_bed$start
  input_bed$og_end <- input_bed$end

  if (opts$collapse) {
    input_bed$mid <- (input_bed$start + input_bed$end + 1) / 2
    input_bed$start <- input_bed$mid
    input_bed$end <- input_bed$mid
  }

  input_bed$step1 <- flankstep
  input_bed$step2 <- flankstep

  for (i in 1:nrow(input_bed)) {
    cur_record <- input_bed[i, ]
    prev_record <- input_bed[i - 1, ]
    next_record <- input_bed[i + 1, ]

    process_record(cur_record, process_record, next_record, opts)
  }
}

#' Parse target region, creating windows based on the flanking inputs given on the command line
#'
#' @param cur_record      Record to process
#' @param previous_record Previous record in file
#' @param next_record     Next record in file
#' @param opts            List of options for wzmetagene
#'
#' @return NULL
process_record <- function(cur_record, previous_record, next_record, opts) {

  if (opts$flankbygene) {
    cur_record$step1 <- (cur_record$end - cur_record$start + 1) / opts$numinternal
    cur_record$step2 <- cur_record$step1
  } else if (opts$flanktoneighbor) {
    if (is.null(previous_record) || cur_record$start < cur_record$end) {
      cur_record$step1 <- -1
    } else {
      cur_record$step1 <- (cur_record$start - cur_record$end + 1) / opts$flanknumber
      if (opts$fold) {
        cur_record$step1 <- cur_record$step1 / 2
      }
    }

    if (is.null(next_record)) {
      cur_record$step2 <- -1
    } else {
      cur_record$step2 <- (next_record$start - cur_record$end + 1) / opts$flanknumber

      if (opts$fold) {
        cur_record$step2 <- cur_record$step2 / 2
      }
    }
  }

  if (is.null(opts$strand)) {
    strand <- "+"
  } else {
    strand <- cur_record$strand
  }

  if (opts$ignoreend) {
    opts$numinternal <- 1
    if (strand == "+") {
      cur_record$end <- cur_record$start + 1
    } else {
      cur_record$start <- cur_record$end - 1
    }
  }

  if (opts$fold) {
    sample_backward(cur_record, function(i) { -i - 1 }, opts)
    sample_internal(cur_record, function(i) { min(i, opts$numinternal - 1 - i) }, opts)
    sample_forward(cur_record, function(i) { -i - 1 } , opts)
  }

  if (strand == "+") {
    sample_backward(cur_record, function(i) { -i }, opts) # removed -1 since 1 base
    sample_internal(cur_record, function(i) { i - 1 }, opts) # -1 to offset 1 base and get a 0 start
    sample_forward(cur_record, function(i) { opts$numinternal + i }, opts)
  } else {
    sample_forward(cur_record, function(i) { -i - 1 }, opts)
    sample_internal(cur_record, function(i) { opts$numinternal - i}, opts)
    sample_backward(cur_record, function(i) { opts$numinternal + i - 1 }, opts)
  }
}

#' Move towards the start of the chromosome (relative to the reference) from target region
#'
#' @param cur_record the record being sampled
#' @param index_func function for adjusting the index value for the bin
#' @param opts       List of options for wzmetagene
#'
#' @return NULL
sample_backward <- function(cur_record, index_func, opts) {

  if (cur_record$step1 == 0 || opts$flanknumber == 0) {
    return()
  }

  temp_window_end <- cur_record$start

  for (i in 1:opts$flanknumber) {
    temp_window_beg <- temp_window_end - cur_record$step1

    if (opts$outer) {
      window_beg <- as.integer(temp_window_beg)
      window_end <- window_beg + 1
    } else if (opts$middle) {
        window_mid <- as.integer((temp_window_beg + temp_window_end) / 2)
        window_beg <- window_mid - 1
        window_end <- window_mid
    } else {
        window_beg <- as.integer(temp_window_beg)
        window_end <- as.integer(temp_window_end)
    }

    index <- index_func(i)
    if (index < 0) {
      if (opts$flankbygene) {
        reg <- paste0("(", -i - 1, ")-(", -i, ")")
      } else if (opts$flanktoneighbor) {
        if (opts$fold) {
          reg <- paste0(
            "(",
            as.integer(-i - 1 / opts$flanknumber / 2 * 100),
            ")-(",
            as.integer(-i / opts$flanknumber / 2 * 100),
            ")"
          )
        } else {
          reg <- paste0(
            "(",
            as.integer(-i - 1 / opts$flanknumber * 100),
            ")-(",
            as.integer(-i / opts$flanknumber * 100),
            ")"
          )
        }
      } else {
        reg <- paste0(
            "(",
            as.integer(window_beg - cur_record$start),
            ")-(",
            as.integer(window_end - cur_record$start),
            ")"
        )
      }
    } else {
      if (opts$flankbygene) {
        reg <- paste0(i, "-", i + 1)
      } else if (opts$flanktoneighbor) {
        if (opts$fold) {
          reg <- paste0(
            as.integer(i / opts$flanknumber / 2 * 100),
            "-",
            as.integer(i + 1 / opts$flanknumber / 2 * 100)
          )
        } else {
          reg <- paste0(
            as.integer(i / opts$flanknumber * 100),
            "-",
            as.integer(i + 1 / opts$flanknumber * 100)
          )
        }
      } else {
        reg <- paste0(
          as.integer(cur_record$start - window_end),
          "-",
          as.integer(cur_record$start - window_beg)
        )
      }
    }

    if (index >= 0) {
      area <- 1
    } else {
      area <- -1
    }

    if (opts$expansion > 0) {
      window_beg <- max(window_beg - opts$expansion, 0)
      window_end <- window_end + opts$expansion
    }

    if (window_beg > 0 && window_end > window_beg) {
      cat(paste(
        cur_record$chr,
        as.integer(window_beg),
        as.integer(window_end),
        index,
        reg,
        area,
        cur_record$chr,
        cur_record$og_start,
        cur_record$og_end,
        cur_record$strand,
        sep = "\t"
      ), sep="\n")
    }

    temp_window_end <- temp_window_beg
  }
}

#' Move towards the end of the chromosome (relative to the reference) from target region
#'
#' @param cur_record the record being sampled
#' @param index_func function for adjusting the index value for the bin
#' @param opts       List of options for wzmetagene
#'
#' @return NULL
sample_forward <- function(cur_record, index_func, opts) {

  if (cur_record$step2 < 0 || opts$flanknumber == 0) {
    return ()
  }

  temp_window_beg <- cur_record$end

  for (i in 1:opts$flanknumber - 1) {
    temp_window_end <- temp_window_beg + cur_record$step2

    if (opts$outer) {
      window_end <- as.integer(temp_window_end)
      window_beg <- window_end - 1
    } else if (opts$middle) {
      window_mid <- as.integer((temp_window_beg + temp_window_end) / 2)
      window_beg <- window_mid - 1
      window_end <- window_mid
    } else {
      window_beg <- as.integer(temp_window_beg)
      window_end <- as.integer(temp_window_end)
    }

    index <- index_func(i)

    if (index < 0) {
      if (opts$flankbygene) {
        reg <- paste0("(", -i - 1, ")-(", -i, ")")
      } else if (opts$flanktoneighbor) {
        if (opts$fold) {
          reg <- paste0(
            "(",
            as.integer(-i - 1 / opts$flanknumber / 2 * 100),
            ")-(",
            as.integer(-i / opts$flanknumber / 2 * 100),
            ")"
          )
        } else {
          reg <- paste0(
            "(",
            as.integer(-i - 1 / opts$flanknumber * 100),
            ")-(",
            as.integer(-i / opts$flanknumber * 100),
            ")"
          )
        }
      } else {
        reg <- paste0(
          "(",
          as.integer(cur_record$end - window_end),
          ")-(",
          as.integer(cur_record$end - window_beg),
          ")"
        )
      }
    } else {
      if (opts$flankbygene) {
        reg <- paste0(i, "-", i + 1)
      } else if (opts$flanktoneighbor) {
        reg <- paste0(
          as.integer(i / opts$flanknumber * 100),
          as.integer((i + 1) / opts$flanknumber * 100)
        )
      } else {
        reg <- paste0(
            as.integer(window_beg - cur_record$end),
            "-",
            as.integer(window_end - cur_record$end)
        )
      }
    }

    if (index >= 0) {
      area <- 1
    } else {
      area <- -1
    }

    if (opts$expansion > 0) {
      window_beg <- max(window_beg - opts$expansion, 0)
      window_end <- window_end + opts$expansion
    }

    if (window_beg > 0 && window_end > window_beg) {
      cat(paste(
        cur_record$chr,
        as.integer(window_beg),
        as.integer(window_end),
        index,
        reg,
        area,
        cur_record$chr,
        cur_record$og_start,
        cur_record$og_end,
        cur_record$strand,
        sep = "\t"
      ), sep="\n")
    }

    temp_window_beg <- temp_window_end
  }
}

#' Move within the target region
#'
#' @param cur_record the record being sampled
#' @param index_func function for adjusting the index value for the bin
#' @param opts       List of options for wzmetagene
#'
#' @return NULL
sample_internal <- function(cur_record, index_func, opts) {

  if (opts$outer) {

    sentinels <- seq(cur_record$start + 1, cur_record$end, length.out = opts$numinternal)

    for (i in 1:length(sentinels)) {
      window_end <- as.integer(sentinels[i])
      window_beg <- window_end - 1

      if (opts$expansion) {
        window_beg <- max(window_beg - opts$expansion, 0)
        window_end <- window_end + opts$expansion
      }

      index <- index_func(i)

      if (window_beg > 0 && window_end > window_beg) {
        cat(paste(
          cur_record$chr,
          as.integer(window_beg),
          as.integer(window_end),
          index,
          paste0(
            as.integer(index / opts$numinternal * 100),
            "-",
            as.integer((index + 1) / opts$numinternal * 100),
            "%"),
          0,
          cur_record$chr,
          cur_record$og_start,
          cur_record$og_end,
          cur_record$strand,
          sep = "\t"
        ), sep="\n")
      }
    }
  } else {

    sentinels <- seq(cur_record$start, cur_record$end, length.out = opts$numinternal + 1)
    if (length(sentinels) < 2) return()

      for (i in 1:(length(sentinels) - 1)) {
        window_beg <- as.integer(sentinels[i])
        window_end <- as.integer(sentinels[i + 1])

        if (opts$middle) {
          window_mid <- as.integer((sentinels[i] + sentinels[i + 1]) / 2)
          window_beg <- window_mid - 1
          window_end <- window_mid
        }

        if (opts$expansion > 0) {
          window_beg <- max(window_beg - opts$expansion, 0)
          window_end <- window_end + opts$expansion
        }

        index <- index_func(i)

        if (window_beg > 0 && window_end > window_beg) {
          cat(paste(
            cur_record$chr,
            as.integer(window_beg),
            as.integer(window_end),
            index,
            paste0(
              as.integer(index / opts$numinternal * 100),
              "-",
              as.integer((index + 1) / opts$numinternal * 100),
              "%"),
            0,
            cur_record$chr,
            cur_record$og_start,
            cur_record$og_end,
            cur_record$strand,
            sep = "\t"
          ), sep = "\n")
        }
      }
  }
}
