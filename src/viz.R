#' Visualization functions
#' 

plot_pairs <- function(guess_points, df_results, log_lik_diff = 10) {
  cols <- names(guess_points)
  cols_to_plot <- cols[!(cols %in% c("s_0", "e_0", "a_0", "i_0", "h_0", "r_0", "ra_0", "ri_0"))]
  res_to_plot <- (
    df_results
    %>% filter(log_lik>max(log_lik) - log_lik_diff)
    %>% bind_rows(guess_points)
    %>% mutate(
      type = if_else(is.na(log_lik), "guess", "result")
    )
    %>% arrange(type)
    %>% select(log_lik, all_of(cols_to_plot), type)
  )
  pairs(
    res_to_plot[1: dim(res_to_plot)[2] - 1],
    col=ifelse(res_to_plot$type == "guess", grey(0.5), "red"),
    cex.axis = 1.3
  )
}


plot_data_and_sample_sims <- function(data, target, sims) {
    min_year <- min(year(data$date))
    max_year <- max(year(data$date))
    years <- seq(min_year, max_year, by = 1)
    num_years <- length(years)
    base_plot <- (
      ggplot()
      + theme(
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.ticks.x = element_line(color = c(rep("black", num_years), rep(NA, num_years)))
      )
    )
    p1 <- (
      base_plot
      + geom_line(
        data = data,
        mapping = aes(x = date, y = .data[[target]]), 
        color = "firebrick", 
        linetype = "solid",
        linewidth = 1
      )
      + scale_x_date(
        breaks = c(as.Date(paste0(years, "-01-01")), as.Date(paste0(years, "-07-01"))),
        labels = c(rep(c(""), num_years), rep(years)),
      )
      + labs(
        title = "True suspected cholera cases",
        x = "Date",
        y = "Counts"
      )
    )
    p2 <- (
      base_plot
      + geom_line(
        data = sims %>% filter(.id %in% sample.int(Nsim, 4)),
        mapping = aes(x = date, y = y),
        color = "cornflowerblue",
        linetype = "solid",
        linewidth = 1
      )
      + facet_wrap(~.id, ncol = 2)
      + scale_x_date(
        breaks = c(as.Date(paste0(years, "-01-01")), as.Date(paste0(years, "-07-01"))),
        labels = c(rep(c(""), num_years), rep(years)),
      )
      + labs(
        title = "Simulated reported cases - best Log Likelihood",
        x = "Date",
        y = "Counts"
      )
    )
    (
      p1/p2
      + plot_layout(heights = c(1,2))
    )
}

plot_data_and_summary_sims <- function(data, target, sims) {
  min_year <- min(year(data$date))
  max_year <- max(year(data$date))
  years <- seq(min_year, max_year, by = 1)
  num_years <- length(years)
  (
    ggplot()
    + theme(
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
      panel.grid = element_blank(),
      plot.background = element_blank(),
      plot.title = element_text(size = 16, hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      axis.ticks.x = element_line(color = c(rep("black", num_years), rep(NA, num_years))),
      legend.background = element_blank(),
      legend.text = element_text(size = 13),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing.y = unit(-15, "pt"),
      legend.position = "inside",
      legend.position.inside = c(0.75, 0.88)
    )
    + geom_line(
      data = data,
      mapping = aes(x = date, y = .data[[target]], color = "data"), 
      linetype = "dashed"
    )
    + geom_line(
      data = sims,
      mapping = aes(x = date, y = q50, color = "median", ), 
      linetype = "solid",
      linewidth = 1
    )
    + geom_ribbon(
      data = sims,
      aes(x = date, ymin = q07, ymax = q93, fill = "band"),
      alpha=0.3
    )
    + scale_color_manual(
      name = "",
      values = c("data" = "black", "median" = "cornflowerblue"),
      labels = c("data" = "Data", "median" = "Median of simulated trajectories")
    )
    + scale_fill_manual(
      name = "",
      values = c("band" = "cornflowerblue"),
      labels = c("band" = "86% coverage")
    )
    + scale_x_date(
      breaks = c(as.Date(paste0(years, "-01-01")), as.Date(paste0(years, "-07-01"))),
      labels = c(rep(c(""), num_years), rep(years)),
    )
    + guides(
      color = guide_legend(order = 1),
      fill = guide_legend(order = 2)
    )
    + labs(
      x = "Date",
      y = "Cases"
    )
  )
}

plot_compartments_sims <- function(data, sims, compartments, n = 20) {
  compartment_map <- c(
    "Population" = "N",
    "Susceptible" = "S",
    "Exposed" = "E",
    "Asymptomatic" = "A",
    "Symptomatic" = "I",
    "Hospitalized" = "H",
    "Recovered Asymptomatic" = "RA",
    "Recovered Infectious" = "RI",
    "Vaccinated" = "VE",
    "Vaccinated Immune" = "V",
    "Cases" = "C"
  )
  sims <- (
    sims
    %>% select(date, .id, all_of(compartments))
    %>% rename(all_of(compartment_map[compartment_map %in% compartments]))
    %>% pivot_longer(
      cols = -c(date, .id),
      names_to = "compartment",
      values_to = "counts"
    )
    %>% group_by(date, compartment)
    %>% summarise(
      q07 = quantile(counts, 0.07, na.rm = TRUE),
      q50 = quantile(counts, 0.5, na.rm = TRUE),
      mean = mean(counts, na.rm = TRUE),
      q93 = quantile(counts, 0.93, na.rm = TRUE),
      .groups = "drop"
    )
    %>% ungroup()
  )
  min_year <- min(year(data$date))
  max_year <- max(year(data$date))
  years <- seq(min_year, max_year, by = 1)
  num_years <- length(years)
  (
    ggplot()
    + theme(
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
      panel.grid = element_blank(),
      plot.background = element_blank(),
      plot.title = element_text(size = 16, hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      axis.ticks.x = element_line(color = c(rep("black", num_years), rep(NA, num_years))),
      legend.background = element_blank(),
      legend.text = element_text(size = 13),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing.y = unit(-15, "pt"),
      legend.position = "inside",
      legend.position.inside = c(0.4, 0.93),
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold")
    )
    + geom_line(
      data = sims,
      mapping = aes(x = date, y = q50, color = "median"), 
      linetype = "solid",
      linewidth = 1
    )
    + geom_ribbon(
      data = sims,
      aes(x = date, ymin = q07, ymax = q93, fill = "band"),
      alpha=0.3
    )
    + scale_color_manual(
      name = "",
      values = c("median" = "forestgreen"),
      labels = c("median" = "Median of simulated trajectories")
    )
    + scale_fill_manual(
      name = "",
      values = c("band" = "forestgreen"),
      labels = c("band" = "86% coverage")
    )
    + facet_wrap(
      ~compartment, 
      scales = "free_y", 
      ncol = 1
    )
    + scale_x_date(
      breaks = c(as.Date(paste0(years, "-01-01")), as.Date(paste0(years, "-07-01"))),
      labels = c(rep(c(""), num_years), rep(years)),
    )
    + guides(
      color = guide_legend(order = 1),
      fill = guide_legend(order = 2)
    )
    + labs(
      x = "Date",
      y = "Counts"
    )
  )
}
