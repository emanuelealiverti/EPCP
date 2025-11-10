# --- Required Libraries ---
require(ggplot2)
require(latex2exp) # For math labels in facets
require(grid)
# --- Helper Function for Facet Labels ---
appender <- function(string) {
  TeX(sprintf("$\\beta_%s$", string))
}

#' Plot Univariate Marginals for Gaussian Approximations and MCMC
#'
#' This function creates a faceted plot comparing the univariate marginal 
#' densities of one or more Gaussian distributions (as lines) against an 
#' optional MCMC sample (as a filled polygon).
#'
#' @param m_list A list of mean vectors.
#' @param S_list A list of covariance matrices. Must be the same length 
#'   as m_list.
#' @param mcmc_sam (Optional) A matrix or data frame of MCMC samples.
#' @param methods_labels (Optional) A character vector of names for 
#'   each Gaussian approximation.
#' @param mcmc_label (Optional) A character label for the MCMC sample. 
#'   Defaults to "MCMC".
#' @param param_names (Optional) A character vector of names for the 
#'   parameters (e.g., c("Intercept", "Slope")).
#' @param n_points (Optional) The number of points to use for 
#'   evaluating the density curves. Defaults to 201.
#'
#' @return A ggplot object.
#'
plot_marginals <- function(m_list, S_list,
                                     mcmc_sam = NULL,
                                     methods_labels = NULL,
                                     mcmc_label = "MCMC",
                                     param_names = NULL,
                                     n_points = 201) {

  # --- 1. Input Validation ---
  if (length(m_list) != length(S_list)) {
    stop("m_list and S_list must have the same length.")
  }
  n_gaussians <- length(m_list)
  if (n_gaussians == 0) {
    stop("At least one mean vector (m_list) and covariance matrix (S_list) must be provided.")
  }
  
  p <- length(m_list[[1]]) # Number of parameters

  if (is.null(methods_labels)) {
    methods_labels <- paste("Gaussian", 1:n_gaussians)
  }
  if (length(methods_labels) != n_gaussians) {
    stop("Length of methods_labels must match length of m_list.")
  }
  
  # Flag for whether MCMC data is present
  has_mcmc <- !is.null(mcmc_sam)

  # --- 2. Build the Plotting Data Frame ---
  plot_data_list <- list()

  for (id in 1:p) {
    
    # --- 2a. Determine x-axis range ---
    all_ranges <- list()
    for (k in 1:n_gaussians) {
      mean_k <- m_list[[k]][id]
      sd_k <- sqrt(S_list[[k]][id, id])
      all_ranges[[k]] <- c(mean_k - 4 * sd_k, mean_k + 4 * sd_k)
    }
    
    if (has_mcmc) {
      all_ranges[["mcmc"]] <- range(mcmc_sam[, id], na.rm = TRUE)
    }
    
    final_range <- range(unlist(all_ranges), na.rm = TRUE)
    x_grid <- seq(final_range[1], final_range[2], length.out = n_points)

    # --- 2b. Calculate Gaussian densities ---
    for (k in 1:n_gaussians) {
      y_values <- dnorm(x_grid, 
                        mean = m_list[[k]][id], 
                        sd = sqrt(S_list[[k]][id, id]))
      
      plot_data_list[[length(plot_data_list) + 1]] <- data.frame(
        x = x_grid,
        y = y_values,
        parameter_id = id,
        method = methods_labels[k],
        type = "Gaussian"
      )
    }

    # --- 2c. Calculate MCMC density (if provided) ---
    if (has_mcmc) {
      kde <- density(mcmc_sam[, id], 
                     n = n_points, 
                     from = final_range[1], 
                     to = final_range[2],
                     na.rm = TRUE)
      
      plot_data_list[[length(plot_data_list) + 1]] <- data.frame(
        x = kde$x,
        y = kde$y,
        parameter_id = id,
        method = mcmc_label,
        type = "MCMC"
      )
    }
  }

  plot_data <- do.call(rbind, plot_data_list)

  # --- 3. Create the Plot ---

  # --- 3a. Set up facet labels ---
  if (is.null(param_names)) {
    plot_data$parameter_label <- factor(plot_data$parameter_id)
    facet_labeller <- as_labeller(appender, default = label_parsed)
  } else {
    plot_data$parameter_label <- factor(plot_data$parameter_id, 
                                        levels = 1:p, 
                                        labels = param_names)
    facet_labeller <- "label_value"
  }

  # --- 3b. Set up scales for a combined legend ---
  
  # Palette from your original code
  palette <- c("#ef476f", "#f4a261", "#457b9d", "grey85", "black")
  
  # Define labels
  all_method_labels <- c(methods_labels, if(has_mcmc) mcmc_label else NULL)

  # 1. Colors (for lines)
  gauss_colors <- palette[1:n_gaussians]
  # Ensure we have enough colors
  if (n_gaussians > 5) { 
    gauss_colors <- scales::hue_pal()(n_gaussians)
  }
  color_vals <- c(gauss_colors, if(has_mcmc) NA else NULL)
  
  # 2. Linetypes (for lines)
  line_types <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
  # Repeat if necessary
  gauss_linetypes <- rep(line_types, length.out = n_gaussians)
  linetype_vals <- c(gauss_linetypes, if(has_mcmc) NA else NULL)

  # 3. Fills (for ribbon)
  mcmc_fill <- palette[4] # "grey85"
  fill_vals <- c(rep(NA, n_gaussians), if(has_mcmc) mcmc_fill else NULL)
  
  # Set names for all scales
  names(color_vals) <- all_method_labels
  names(linetype_vals) <- all_method_labels
  names(fill_vals) <- all_method_labels
  
  
  # --- 3c. Generate ggplot object ---
  p_out <- ggplot(plot_data, aes(x = x, group = method)) +
    
    # Add MCMC density polygon (if present)
    {if (has_mcmc) 
      geom_ribbon(data = subset(plot_data, type == "MCMC"),
                  aes(ymin = 0, ymax = y, fill = method),
                  alpha = 0.8)
    } +
    
    # Add Gaussian density lines
    geom_line(data = subset(plot_data, type == "Gaussian"),
              aes(y = y, color = method, linetype = method),
              lwd = 1,key_glyph = draw_key_colored_box_line) +
    
    # Facets and theme
    facet_wrap(~ parameter_label, labeller = facet_labeller, scales = "free") +
    theme_bw(base_size = 14) +
    xlab("") +
    ylab("") +
    scale_y_continuous(expand = c(0.01, 0.01)) +
    
    # Add combined scales
    scale_color_manual(name = "Method", 
                       values = color_vals) +
    scale_linetype_manual(name = "Method", 
                          values = linetype_vals) +
    {if (has_mcmc) 
    scale_fill_manual(name = "Method", 
                      values = fill_vals)
    } +
    
    theme(legend.position = "bottom")

  return(p_out)
}


# Helper
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Custom legend key: colored box background + white line using mapped linetype
draw_key_colored_box_line <- function(data, params, size) {
  #lwd <- (data$linewidth %||% data$size %||% 1.1)
  lwd <- 2
  grobTree(
    rectGrob(gp = gpar(fill = data$colour %||% "grey60", col = NA)),
    linesGrob(
      x = unit(c(0.12, 0.88), "npc"),
      y = unit(c(0.5, 0.5), "npc"),
      gp = gpar(col = "white", lty = data$linetype %||% 1, lwd = lwd)
    )
  )
}

