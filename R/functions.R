#' create_table_descriptive_stats
#'
#' @param data
#'
#' @returns “A data.frame/tibble.”
#'  summarised mean + SD rounded metabolite data values

create_table_descriptive_stats <-
  function(data) {
    data |>
      dplyr::group_by(metabolite) |>
      dplyr::summarise(dplyr::across(value, list(mean = mean, sd = sd))) |>
      dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), \(x) round(x, digits = 2))) |>
      dplyr::mutate(MeanSD = glue::glue("{value_mean} ({value_sd})")) |>
      dplyr::select(Metabolite = metabolite, "Mean SD" = "MeanSD")
  }

#' create_plot_distributions
#'
#' @param data
#'
#' @returns create plot distributions

create_plot_distributions <- function(data) {
  data |>
    ggplot2::ggplot(
      ggplot2::aes(x = value)
    ) +
    ggplot2::geom_histogram() +
    ggplot2::facet_wrap(ggplot2::vars(metabolite), scales = "free") +
    ggplot2::theme_minimal()
}
