#' Compile tables with min, mean and max GC contents per group of baits
#' 
#' @param x The output of main_function
#' @param combine Logical: Should results be displayed per chromosome?
#' 
#' @keywords internal
#' 
gc_table <- function(x, combine = FALSE) {
  if (!is.list(x) || is.null(x$baits))
    stop("could not recognise x as the output of main_function\n", call. = FALSE)

  baits <- x$baits

  if (combine)
    input <- list(data.table::rbindlist(baits))
  else
    input <- baits
  
	output <- lapply(input, function(i) {
    aux <- as.data.frame(with(i, table(bait_type)))
    colnames(aux)[2] <- "n"
    aux$min_GC <-  with(i, aggregate(pGC, list(bait_type), min))$x
    aux$mean_GC <- with(i, aggregate(pGC, list(bait_type), mean))$x
    aux$max_GC <-  with(i, aggregate(pGC, list(bait_type), max))$x
    return(aux)
  })

  if (combine)
    return(as.data.frame(output))
  else
    return(output)
}

#' Compile tables with coverage values per chromosome
#' 
#' @param x The output of main_function
#' @param combine Logical: Should results be displayed per type of bait?
#' 
#' @keywords internal
#' 
coverage <- function(x, combined = TRUE) {
  if (!is.list(x) || is.null(x$baits))
    stop("could not recognise x as the output of main_function\n", call. = FALSE)

  chr.lengths <- x$input.summary$chr.lengths
  baits <- x$baits
  exclusions <- x$input.summary$exclusions

  output <- lapply(names(baits), function(i) {
    input <- baits[[i]]
    
    if (combined)
      input$bait_type <- "All"

    input <- split(input, input$bait_type)

    total_length <- chr.lengths$size[chr.lengths$name == i]
    if (!is.null(exclusions)) {
      if (any(exclusions$chr == i)) {
        aux <- exclusions[exclusions$chr == i, ]
        aux <- apply(aux, 1, function(j) j["start"]:j["stop"])
        if (is.list(aux))
          excluded_bps <- length(unique(unlist(aux)))
        else
          excluded_bps <- length(unique(as.vector(aux)))
        valid_length <- total_length - excluded_bps
      } else {
        valid_length <- total_length
      }
    } else {
      valid_length <- total_length
    }

    capture <- lapply(names(input), function(x) {
      bp_covered <- length(unique(as.vector(apply(input[[x]], 1, function(j) j["bait_start"]:j["bait_stop"]))))
      output <- data.frame(bait_type = x, 
                           bp = bp_covered, 
                           valid_coverage = bp_covered/valid_length,
                           total_coverage = bp_covered/total_length)
      return(output)
    })
    return(data.table::rbindlist(capture))
  })
  names(output) <- names(baits)
  return(output)
}

#' Print graphics with coverage per sequence
#' 
#' @param x The output of main_function.
#' @param sequence the sequence to be printed.
#' 
#' @keywords internal
#' 
print_coverage <- function(x, seq.name) {
  if (!is.list(x) || is.null(x$baits))
    stop("could not recognise x as the output of main_function\n", call. = FALSE)

  if (is.na(match(seq.name, names(x$baits))))
    stop("Could not find sequence '", seq.name, "' in the input data\n", call. = FALSE)

# extract inputs  
  # bait size
  size <- x$input.summary$size

  # length for chosen sequence
  chr.length <- x$input.summary$chr.lengths$size[x$input.summary$chr.lengths == seq.name]

  # exclusions for chosen sequence
  if (!is.null(x$input.summary$exclusions))
    exclusions <- x$input.summary$exclusions[x$input.summary$exclusions$chr == seq.name, ]
  else
    exclusions <- NULL

  # targets for chosen sequence
  if (!is.null(x$input.summary$targets))
    targets <- x$input.summary$targets[x$input.summary$targets$chr == seq.name, ]
  else
    targets <- NULL
  
  # bait positions
  baitpoints <- x$baits[[seq.name]][, c("bait_no", "bait_type", "bait_start", "bait_stop")]
# -----

  plotdata <- reshape2::melt(baitpoints, id.vars = c("bait_type", "bait_no"))
  plotdata$Colour <- "covered zones"

  baitpoints$X <- apply(baitpoints[, c("bait_start", "bait_stop")], 1, mean)
  baitpoints$Colour <- "covered zones"

  # Add solid line for whole sequence at the top of the graphic
  aux <- data.frame(bait_type = rep(seq.name, 2), 
    bait_no = rep(-999, 2), 
    variable = c("bait_start", "bait_stop"), 
    value = c(1, chr.length),
    Colour = rep("chromosome", 2))
  plotdata <- rbind(plotdata, aux)
  # ---

  # add exclusions, if relevant, on top of the whole sequence
  if (!is.null(exclusions) && nrow(exclusions) > 0) {
    colnames(exclusions) <- c("bait_type", "bait_start", "bait_stop")
    exclusions$bait_no <- c(-998:(-999 + nrow(exclusions)))
    exclusions <- reshape2::melt(exclusions, id.vars = c("bait_type", "bait_no"))
    exclusions$Colour <- "excluded zones"
    plotdata <- rbind(plotdata, exclusions)
  }

  # final preparations
  plotdata$bait_type <- factor(plotdata$bait_type, levels = c("random", "target", "region", seq.name))
  used.levels <- levels(droplevels(plotdata$bait_type))
  baitpoints$bait_type <- factor(baitpoints$bait_type, levels = c("random", "target", "region", seq.name))
  plotdata$Colour <- factor(plotdata$Colour, levels = c("chromosome", "excluded zones", "covered zones"))

  p <- ggplot2::ggplot()
  if (size / chr.length > 0.003) {
    p <- p + ggplot2::geom_line(data = plotdata, ggplot2::aes(x = value, y = bait_type, group = bait_no, colour = Colour), size = 2)
  } else {
    p <- p + ggplot2::geom_line(data = plotdata[plotdata$bait_no < 0, ], ggplot2::aes(x = value, y = bait_type, group = bait_no, colour = Colour), size = 2)
    p <- p + ggplot2::geom_point(data = baitpoints, ggplot2::aes(x = X, y = bait_type, colour = Colour), shape = "|", size = 1)
  }
  p <- p + ggplot2::scale_colour_manual(values = c("chromosome" = "#56B4E9", "excluded zones" = "grey", "covered zones" = "#009E73"))
  p <- p + ggplot2::scale_y_discrete(limits = used.levels)
  p <- p + ggplot2::labs(x = "bp", y = "")
  if (!is.null(targets) && nrow(targets) > 0) {
    targets$bait_type <- "target"
    p <- p + ggplot2::geom_point(data = targets, ggplot2::aes(x = target, y = bait_type), colour = "#E69F00", shape = "I", size = 5)
  }
  p <- p + ggplot2::theme(legend.position = "none")

  return(p)
}


print_report <- function(gc.table) {
	 inst.ver <- utils::packageVersion("supeRbaits")
  inst.ver.short <- substr(inst.ver, start = 1, stop = nchar(as.character(inst.ver)) - 5) 
  if (file.exists(reportname <- "Report/supeRbaits_report.Rmd")) {
    continue <- TRUE
    index <- 1
    while (continue) {
      if(file.exists(reportname <- paste0("Report/supeRbaits_report.", index, ".Rmd"))) {
        index <- index + 1
      } else {
        continue <- FALSE
      }
    }
    message("M: A superBaits report is already present in the current directory.\n   Saving new report as 'supeRbaits_report.", index, ".html'.")
    rm(continue,index)
  } else {
    appendTo("Screen", "M: Saving supeRbaits report as 'supeRbaits_report.html'.")
  }
  report <- readr::read_file("temp_log.txt")
  sink(reportname)
  cat(paste0(
'---
title: "Summary of generated baits"
author: "supeRbaits R package (', inst.ver.short, ')"
output: 
  html_document:
    includes:
      after_body: toc_menu.html
---

### GC content

', paste(knitr::kable(gc.table, row.names = FALSE), collapse = "\n"),' 

}

'))
 sink()
}