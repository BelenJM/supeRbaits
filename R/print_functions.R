#' Compile tables with min, mean and max GC contents per group of baits
#' 
#' @param baits A list of baits per chromosome
#' @param combine Logical: Should results be displayed per chromosome?
#' 
#' @keywords internal
#' 
gc_table <- function(baits, combine = FALSE) {
  if (combine) {
    if (!is.list(baits))
      stop("'combine' is set to TRUE but 'baits' is not a list.\n", call. = FALSE)
    else
      input <- data.table::rbindlist(baits)
  } else {
    input <- baits
  }
  
  if(is.data.frame(input))
    input <- list(input)

	output <- lapply(input, function(i) {
    aux <- as.data.frame(with(i, table(bait_type)))
    colnames(aux)[2] <- "n"
    aux$min_GC <-  with(i, aggregate(pGC, list(bait_type), min))$x
    aux$mean_GC <- with(i, aggregate(pGC, list(bait_type), mean))$x
    aux$max_GC <-  with(i, aggregate(pGC, list(bait_type), max))$x
    return(aux)
  })

  if (length(output) == 1)
    return(as.data.frame(output))
  else
    return(output)
}

#' Compile tables with coverage values per chromosome
#' 
#' @param chr.lengths The lengths of the chromosomes
#' @param baits A list of baits per chromosome
#' @param exclusions The exclusions table
#' @param combine Logical: Should results be displayed per type of bait?
#' 
#' @keywords internal
#' 
coverage <- function(chr.lengths, baits, exclusions = NULL, combined = TRUE) {
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
      output <- data.frame(bait_type = x, bp = bp_covered, Valid_coverage = bp_covered/valid_length, Total_coverage = bp_covered/total_length)
      return(output)
    })
    return(data.table::rbindlist(capture))
  })
  names(output) <- names(baits)
  return(output)
}

#' Print graphics with coverage per chromosome
#' 
#' @param chr.lengths The lengths of the chromosomes
#' @param baits A list of baits per chromosome
#' @param exclusions The exclusions table
#' @param targets The targets table
#' 
#' @keywords internal
#' 
print_coverage <- function(chr.lengths, baits, size, exclusions = NULL, targets = NULL) {
  capture <- lapply(names(baits), function(i) {
    baitpoints <- baits[[i]][, c("bait_no", "bait_type", "bait_start", "bait_stop")]
    plotdata <- reshape2::melt(baitpoints, id.vars = c("bait_type", "bait_no"))
    plotdata$Colour <- "covered zones"

    baitpoints$X <- apply(baitpoints[, c("bait_start", "bait_stop")], 1, mean)
    baitpoints$Colour <- "covered zones"

    aux <- data.frame(bait_type = rep(i, 2), 
      bait_no = rep(-999, 2), 
      variable = c("bait_start", "bait_stop"), 
      value = c(1, chr.lengths$size[chr.lengths$name == i]),
      Colour = rep("chromosome", 2))

    plotdata <- rbind(plotdata, aux)

    if (!is.null(exclusions)) {
      if (any(exclusions$chr == i)) {
        aux <- exclusions[exclusions$chr == i, ]
        colnames(aux) <- c("bait_type", "bait_start", "bait_stop")
        aux$bait_no <- c(-998:(-999 + nrow(aux)))
        aux <- reshape2::melt(aux, id.vars = c("bait_type", "bait_no"))
        aux$Colour <- "excluded zones"
      }
      plotdata <- rbind(plotdata, aux)
    }

    plotdata$bait_type <- factor(plotdata$bait_type, levels = c("random", "target", "region", i))
    used.levels <- levels(droplevels(plotdata$bait_type))
    baitpoints$bait_type <- factor(baitpoints$bait_type, levels = c("random", "target", "region", i))
    plotdata$Colour <- factor(plotdata$Colour, levels = c("chromosome", "excluded zones", "covered zones"))

    p <- ggplot2::ggplot()
    if (size / chr.lengths$size[chr.lengths$name == i] > 0.003) {
      p <- p + ggplot2::geom_line(data = plotdata, ggplot2::aes(x = value, y = bait_type, group = bait_no, colour = Colour), size = 2)
    } else {
      p <- p + ggplot2::geom_line(data = plotdata[plotdata$bait_no < 0, ], ggplot2::aes(x = value, y = bait_type, group = bait_no, colour = Colour), size = 2)
      p <- p + ggplot2::geom_point(data = baitpoints, ggplot2::aes(x = X, y = bait_type, colour = Colour), shape = "I", size = 2)
    }
    p <- p + ggplot2::scale_colour_manual(values = c("chromosome" = "#56B4E9", "excluded zones" = "grey", "covered zones" = "#009E73"))
    p <- p + ggplot2::scale_y_discrete(limits = used.levels)
    p <- p + ggplot2::labs(x = "bp", y = "")
    if (any(targets$chr == i)) {
      aux <- targets[targets$chr == i, ]
      aux$bait_type <- "target"
      p <- p + ggplot2::geom_point(data = aux, ggplot2::aes(x = target, y = bait_type), colour = "#E69F00", shape = "I", size = 5)
    }
    p <- p + ggplot2::theme(legend.position = "none")
    p
    ggplot2::ggsave(paste0("test_", i, ".png"), width = 10, height = 1.2)
    return(p)
  })
  names(capture) <- names(baits)
  return(capture)
}
