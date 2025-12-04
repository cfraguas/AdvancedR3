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
      dplyr::summarise(dplyr::across(value, list(mean = mean, median = median, sd = sd, iqr = IQR))) |>
      dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), \(x) round(x, digits = 2))) |>
      dplyr::mutate(MeanSD = glue::glue("{value_mean} ({value_sd})")) |>
      dplyr::mutate(MedianIQR = glue::glue("{value_median} ({value_iqr})")) |>
      dplyr::select(Metabolite = metabolite, "Mean SD" = MeanSD, "Median IQR" = MedianIQR)
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

#' cleaning_data
#'
#' @param data
#'
#' @returns clean lipidomics data
#'
clean <- function(data) {
  data |>
    dplyr::group_by(dplyr::pick(-value)) |>
    dplyr::summarise(value = mean(value), .groups = "keep") |>
    dplyr::ungroup()
}

#' preprocessing
#'
#' @param data
#'
#' @returns preprocessing
preprocessing <- function(data) {
  data |>
    dplyr::mutate(
      class = as.factor(class),
      value = scale(value)
    )
}

#' Fitting model
#'
#' @param data
#' @param model
#'
#' @returns fitting model

fit_model <- function(data, model) {
  stats::glm(model, data = data, family = binomial) |>
    broom::tidy(exponentiate = TRUE) |>
    dplyr::mutate(
      metabolite = unique(data$metabolite),
      model = format(model),
      .before = tidyselect::everything()
    )
}


#' fit all models
#'
#' @param data
#'
#' @returns fit any model 1 or 2
fit_all_models <- function(data){
  list(
    class ~ value,
    class ~ value + gender + age
  ) |>
    purrr::map(\(model) fit_model(data, model = model)) |>
    purrr::list_rbind()
}

#' split groups and fit models and bind again
#'
#' @param data
#'
#' @returns final models
create_model_results <- function(data) {
  data |>
    dplyr::group_split(metabolite)|>
    purrr::map(preprocessing) |>
    purrr::map(fit_all_models) |>
    purrr::list_rbind()
}
