specificity_plotv2 <- function(mk, subclass = NULL,
                             group = NULL,
                             type = 1,
                             use_filter = FALSE,
                             nrank = 8,
                             nsubclass = NULL,
                             expfilter = NULL,
                             scheme = NULL,
                             add_labels = NULL,
                             label_pos = "right",
                             axis_extend = 0.4,
                             nudge_x = NULL, nudge_y = NULL,
                             ...) {
  if (!inherits(mk, "cellMarkers")) stop("not a 'cellMarkers' class object")
  if (is.null(subclass) & is.null(group))
    stop("Either subclass or group must be specified")
  
  if (!is.null(subclass)) {
    if (is.numeric(subclass)) subclass <- colnames(mk$genemeans)[subclass]
    if (!subclass %in% colnames(mk$genemeans))
      stop("subclass ", subclass, " not found")
    genemeans <- if (use_filter) mk$genemeans_filtered else mk$genemeans
    if (is.null(nsubclass)) nsubclass <- mk$nsubclass[subclass]
    if (is.null(nsubclass)) nsubclass <- 5
    #labs <- rownames(mk$best_angle[[subclass]][1L:nsubclass, ])
    subc <- subclass
  } else {
    if (is.numeric(group)) group <- colnames(mk$groupmeans)[group]
    if (!group %in% colnames(mk$groupmeans))
      stop("group ", group, " not found")
    genemeans <- if (use_filter) mk$groupmeans_filtered else mk$groupmeans
    if (is.null(nsubclass)) nsubclass <- 5
    #labs <- rownames(mk$group_angle[[group]][1L:nsubclass, ])
    subc <- group
  }
  
  vecLength <- sqrt(rowSums(genemeans^2))
  genemeans_scaled <- genemeans / vecLength
  genemeans_angle <- acos(genemeans_scaled)
  gene_rank <- apply(-genemeans, 1, rank)[subc, ]
  nrank <- pmin(ncol(genemeans), nrank)
  gene_rank[gene_rank > nrank] <- nrank
  gene_rank <- factor(floor(gene_rank), levels = 1:nrank)
  if (type == 1) {
    if (is.null(expfilter)) expfilter <- mk$opt$expfilter
    low <- genemeans[, subc] < expfilter #& !rownames(genemeans) %in% labs
    gene_rank[low] <- nrank
    levels(gene_rank)[nrank] <- paste0(nrank, "+/low")
  } else {
    levels(gene_rank)[nrank] <- paste0(nrank, "+")
  }
  
  df <- data.frame(angle = genemeans_angle[, subc],
                   angle.deg = genemeans_angle[, subc] * 180/pi,
                   mean = genemeans[, subc],
                   rank = gene_rank)
  df$x <- vecLength * sin(df$angle)
  df$y <- vecLength * cos(df$angle)
  df <- df[vecLength != 0, ]
  labs <- unique(c(#labs, 
                   add_labels))
  df$label <- ""
  df$label[match(labs, rownames(df))] <- labs
  df <- df[rev(order(df$rank)), ]
  
  if (is.null(scheme)) {
    scheme <- c(hue_pal(h = c(0, 270), c = 120)(nrank -1),
                adjustcolor("grey", 0.5))
    scheme[1] <- "red"
  }
  
  if (type == 2) {
    # use actual angle; radius is vecLength
    xlim <- xr <- range(df$x, na.rm = TRUE)
    yr <- range(df$y, na.rm = TRUE)
    if (label_pos == "left") {
      xlim[1] <- xlim[1] - diff(xr) * axis_extend
      if (is.null(nudge_x)) nudge_x <- -diff(xr) * 0.2
    } else {
      if (is.null(nudge_x)) nudge_x <- diff(xr) * 0.5
    }
    if (is.null(nudge_y)) nudge_y <- 0
    
    ggplot(df, aes(x = .data$x, y = .data$y, color = .data$rank,
                   label = .data$label)) +
      (if (label_pos == "left") geom_vline(xintercept = 0, lty = 2)) +
      geom_point(show.legend = TRUE) +
      scale_color_manual(values = scheme, drop = FALSE) +
      (if (label_pos == "left") {
        geom_label_repel(size = 3, color = "black",
                         nudge_x = nudge_x, nudge_y = nudge_y,
                         hjust = 1, label.size = NA, direction = "y", ...)
      } else {
        geom_text_repel(size = 3, color = "black",
                        nudge_x = nudge_x, nudge_y = nudge_y,
                        hjust = 0, direction = "y", na.rm = TRUE, ...)
      }) +
      xlim(xlim) + ylim(yr) +
      xlab("Non-specific gene expression") +
      ylab(paste(subc, "mean expression")) +
      theme_classic() +
      theme(axis.text = element_text(colour = "black"))
  } else {
    # angle on x, mean exp on y
    xr <- range(df$angle.deg, na.rm = TRUE)
    yr <- range(df$mean, na.rm = TRUE)
    if (is.null(nudge_x)) nudge_x <- 0
    if (is.null(nudge_y)) nudge_y <- 0.1
    
    ggplot(df, aes(x = .data$angle.deg, y = .data$mean, color = .data$rank,
                   label = .data$label)) +
      rasterize(geom_point(show.legend = TRUE), scale = 0.6, dpi = 600) +
      scale_color_manual(values = scheme, drop = FALSE) +
      geom_text_repel(size = 3, color = "black",
                      nudge_x = nudge_x, nudge_y = nudge_y, ...) +
      xlab("Vector angle") +
      #ylab(paste(subc, "mean expression")) +
      ylab("mean expression") +
      ggtitle(subc) +
      theme_classic() +
      theme(axis.text = element_text(colour = "black"),
            axis.ticks = element_line(color = "black"),
            plot.title = element_text(size = 10)) +
      guides(color=guide_legend(nrow=1)) #for legend = "bottom"
  }
}