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
    aux <- as.data.frame(with(i, table(Type)))
    colnames(aux)[2] <- "n"
    aux$min_GC <-  with(i, aggregate(pGC, list(Type), min))$x
    aux$mean_GC <- with(i, aggregate(pGC, list(Type), mean))$x
    aux$max_GC <-  with(i, aggregate(pGC, list(Type), max))$x
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
      input$Type <- "All"

    input <- split(input, input$Type)

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
      bp_covered <- length(unique(as.vector(apply(input[[x]], 1, function(j) j["Start_bp"]:j["End_bp"]))))
      output <- data.frame(Type = x, bp = bp_covered, Valid_coverage = bp_covered/valid_length, Total_coverage = bp_covered/total_length)
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
print_coverage <- function(chr.lengths, baits, exclusions = NULL, targets = NULL) {
  capture <- lapply(names(baits), function(i) {
    aux <- baits[[i]][, c("Bait_no", "Type", "Start_bp", "End_bp")]
    plotdata <- reshape2::melt(aux, id.vars = c("Type", "Bait_no"))
    plotdata$Colour <- "covered zones"

    aux <- data.frame(Type = rep(i, 2), 
      Bait_no = rep(-999, 2), 
      variable = c("Start_bp", "End_bp"), 
      value = c(1, chr.lengths$size[chr.lengths$name == i]),
      Colour = rep("chromosome", 2))

    plotdata <- rbind(plotdata, aux)

    if (!is.null(exclusions)) {
      if (any(exclusions$chr == i)) {
        aux <- exclusions[exclusions$chr == i, ]
        colnames(aux) <- c("Type", "Start_bp", "End_bp")
        aux$Bait_no <- c(-998:(-999 + nrow(aux)))
        aux <- reshape2::melt(aux, id.vars = c("Type", "Bait_no"))
        aux$Colour <- "excluded zones"
      }
      plotdata <- rbind(plotdata, aux)
    }

    plotdata$Type <- factor(plotdata$Type, levels = c("random", "target", "region", i))
    plotdata$Colour <- factor(plotdata$Colour, levels = c("chromosome", "excluded zones", "covered zones"))

    p <- ggplot2::ggplot()
    p <- p + ggplot2::geom_line(data = plotdata, ggplot2::aes(x = value, y = Type, group = Bait_no, colour = Colour), size = 2)
    p <- p + ggplot2::scale_colour_manual(values = c("chromosome" = "#56B4E9", "excluded zones" = "grey", "covered zones" = "#009E73"))
    p <- p + ggplot2::labs(x = "bp", y = "")
    if (any(targets$chr == i)) {
      aux <- targets[targets$chr == i, ]
      aux$Type <- "target"
      p <- p + ggplot2::geom_point(data = aux, ggplot2::aes(x = target, y = Type), colour = "#E69F00", shape = "I", size = 5)
    }
    p <- p + ggplot2::theme(legend.position = "none")
    ggplot2::ggsave(paste0("test_", i, ".png"), width = 10, height = 1.2)
    return(p)
  })
  names(capture) <- names(baits)
  return(capture)
}
