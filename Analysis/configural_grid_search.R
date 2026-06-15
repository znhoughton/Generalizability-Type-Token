## configural_grid_search.R
## Full configural cue saliency grid search.
## Run from the Analysis/ directory:
##   Rscript configural_grid_search.R

N_WORKERS = 16
N_SIMS    = 50

library(tidyverse)
library(ndl)
library(arm)
library(furrr)
library(future)
library(progressr)
source("RescorlaLogistic.R")

# ── Load and prepare condition data ──────────────────────────────────────────

load_cond = function(path) {
  d = read_csv(path, show_col_types = FALSE) %>%
    mutate(stem = orthoCoding(stem),
           type = case_when(
             suffix %in% c('dan', 'sil')   ~ 'big_pl',
             suffix == 'none'               ~ 'big_sg',
             suffix %in% c('nem', 'shoon') ~ 'dim_sg'
           ),
           key = paste0(stem, '_', type))
  d %>%
    group_by(key, suffix) %>%
    summarize(frequency = table(key), .groups = 'drop') %>%
    setNames(c('Cues', 'Outcomes', 'Frequency')) %>%
    mutate(Frequency = as.numeric(Frequency) * 5) %>%
    uncount(Frequency) %>%
    mutate(Frequency = 1)
}

add_configural = function(d) {
  d %>% mutate(Cues = case_when(
    str_ends(Cues, "big_pl") ~ paste0(Cues, "_big.pl"),
    str_ends(Cues, "big_sg") ~ paste0(Cues, "_big.sg"),
    str_ends(Cues, "dim_sg") ~ paste0(Cues, "_dim.sg"),
    TRUE ~ Cues
  ))
}

cat("Loading condition data...\n")
data_cond1_configural = load_cond('../Conditions/Condition1.csv') %>% add_configural()
data_cond2_configural = load_cond('../Conditions/Condition2.csv') %>% add_configural()
data_cond3_configural = load_cond('../Conditions/Condition3.csv') %>% add_configural()
data_cond4_configural = load_cond('../Conditions/Condition4.csv') %>% add_configural()
data_cond5_configural = load_cond('../Conditions/Condition5.csv') %>% add_configural()
data_cond6_configural = load_cond('../Conditions/Condition6.csv') %>% add_configural()

human_data = read_csv('../Data/fixed_effects_df_all.csv', show_col_types = FALSE) %>%
  filter(Task == 'prod')

# Pre-load production condition to avoid repeated file reads inside the function
data_prod_base = read_csv('../Conditions/production_condition.csv', show_col_types = FALSE) %>%
  mutate(
    stem = orthoCoding(prod_label),
    type = case_when(
      str_detect(file, "dim_sg") ~ "dim_sg_dim.sg",
      str_detect(file, "_2")     ~ "big_pl_big.pl",
      str_detect(file, "dim_pl") ~ "dim_pl",
      TRUE                       ~ NA_character_
    ),
    row_id = row_number()
  )

# ── Model function ────────────────────────────────────────────────────────────

run_rw_production_conf = function(
  n_sims = 5,
  Alpha, SemAlpha, ConfAlpha, Beta,
  configural_list,
  prod_base = data_prod_base,
  hum = human_data
) {
  run_rw = function(config_data) {
    map_dfr(1:n_sims, function(i) {
      Rescorla_modified_alpha_for_configural_cues(
        cuesOutcomes = config_data[sample(nrow(config_data)), ],
        logistic  = TRUE,
        Alpha     = Alpha,
        SemAlpha  = SemAlpha,
        ConfAlpha = ConfAlpha,
        Beta      = Beta
      ) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("cue") %>%
        mutate(sim = i)
    })
  }

  rw_results = map(configural_list, run_rw)

  data_prod = map_dfr(1:n_sims, ~ mutate(prod_base, sim = .x))

  process_one_condition = function(prod_data, rw_data, cond_num) {
    prod_cond = prod_data %>%
      mutate(key = paste0(stem, "_", type)) %>%
      group_by(sim) %>%
      mutate(row_id = row_number())

    res = prod_cond %>%
      separate_rows(key, sep = "_") %>%
      left_join(rw_data, by = c("key" = "cue", "sim")) %>%
      group_by(row_id, sim) %>%
      summarize(
        dan_activations   = sum(dan,   na.rm = TRUE),
        sil_activations   = sum(sil,   na.rm = TRUE),
        nem_activations   = sum(nem,   na.rm = TRUE),
        shoon_activations = sum(shoon, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      left_join(prod_data, by = c("row_id", "sim")) %>%
      mutate(
        condition = cond_num,
        freq_minus_infreq =
          if (cond_num %in% c(1, 2, 3))
            dan_activations - nem_activations
          else
            nem_activations - dan_activations
      )
    res
  }

  prod_results = map2_dfr(seq_along(rw_results), rw_results,
                          ~ process_one_condition(data_prod, .y, .x))

  prod_results = prod_results %>%
    group_by(prod_label, sim, condition) %>%
    mutate(dv = case_when(
      type == "dim_pl" ~ freq_minus_infreq,
      condition %in% c(1, 2, 3) ~
        dan_activations[type == "big_pl_big.pl"][1] -
        nem_activations[type == "dim_sg_dim.sg"][1],
      condition %in% c(4, 5, 6) ~
        nem_activations[type == "dim_sg_dim.sg"][1] -
        dan_activations[type == "big_pl_big.pl"][1]
    ))

  out = prod_results %>%
    mutate(
      meaning = case_when(
        type %in% c("dim_sg_dim.sg", "big_pl_big.pl") ~ "original",
        type == "dim_pl"                               ~ "novel",
        TRUE                                           ~ NA_character_
      ),
      condition = case_when(
        condition %in% c(1, 4) ~ 1,
        condition %in% c(2, 5) ~ 2,
        condition %in% c(3, 6) ~ 3
      )
    ) %>%
    group_by(meaning, prod_condition, condition, sim) %>%
    summarize(mean_freq_diff = mean(freq_minus_infreq), .groups = "keep") %>%
    group_by(meaning, prod_condition, condition) %>%
    summarize(
      mean_freq_diff_v2 = mean(mean_freq_diff),
      ci2.5  = quantile(mean_freq_diff, 0.025, na.rm = TRUE),
      ci97.5 = quantile(mean_freq_diff, 0.975, na.rm = TRUE),
      .groups = "keep"
    )

  hum2 = hum %>%
    mutate(
      Meaning   = tolower(Meaning),
      Condition = case_when(Condition == 'cond1' ~ 1,
                            Condition == 'cond2' ~ 2,
                            Condition == 'cond3' ~ 3),
      Stem      = tolower(Stem)
    ) %>%
    dplyr::select(Meaning, Condition, Stem, Estimate)

  out %>%
    left_join(hum2, by = join_by(meaning       == Meaning,
                                 condition      == Condition,
                                 prod_condition == Stem)) %>%
    mutate(mse = (mean_freq_diff_v2 - Estimate)^2)
}

# ── Grid search ───────────────────────────────────────────────────────────────

param_grid_conf = tidyr::crossing(
  ConfAlpha = seq(0, 1, by = 0.1),
  SemAlpha  = seq(0, 1, by = 0.1),
  Beta      = seq(0, 1, by = 0.1),
  Alpha     = seq(0, 1, by = 0.1)
)

cat(sprintf("Running grid search: %d combinations, %d sims, %d workers\n",
            nrow(param_grid_conf), N_SIMS, N_WORKERS))

plan(multisession, workers = N_WORKERS)
handlers(global = TRUE)
handlers("txtprogressbar")

results_all_conf = param_grid_conf %>%
  mutate(
    model_results = future_pmap(
      list(ConfAlpha, SemAlpha, Beta, Alpha),
      ~ run_rw_production_conf(
          n_sims    = N_SIMS,
          Alpha     = ..4,
          SemAlpha  = ..2,
          ConfAlpha = ..1,
          Beta      = ..3,
          configural_list = list(
            data_cond1_configural, data_cond2_configural, data_cond3_configural,
            data_cond4_configural, data_cond5_configural, data_cond6_configural
          )
      ),
      .options = furrr_options(packages = c("tidyverse", "ndl", "arm"),
                               seed = 964),
      .progress = TRUE
    )
  ) %>%
  unnest(model_results)

plan(sequential)

results_all_conf = results_all_conf %>%
  group_by(ConfAlpha, SemAlpha, Beta, Alpha) %>%
  mutate(mean_mse = mean(mse, na.rm = TRUE))

# ── Score and rank parameter sets ─────────────────────────────────────────────

pattern_scores_conf = results_all_conf %>%
  group_by(Alpha, Beta, SemAlpha, ConfAlpha) %>%
  summarize(
    spearman = cor(Estimate, mean_freq_diff_v2, method = "spearman"),
    kendall  = cor(Estimate, mean_freq_diff_v2, method = "kendall"),
    mean_mse = mean(mse, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(spearman), desc(kendall))

cat(sprintf("\nBest parameters: Alpha=%.2f  Beta=%.2f  SemAlpha=%.2f  ConfAlpha=%.2f\n",
            pattern_scores_conf$Alpha[1], pattern_scores_conf$Beta[1],
            pattern_scores_conf$SemAlpha[1], pattern_scores_conf$ConfAlpha[1]))
cat(sprintf("  Spearman=%.4f  Kendall=%.4f  MSE=%.4f\n",
            pattern_scores_conf$spearman[1], pattern_scores_conf$kendall[1],
            pattern_scores_conf$mean_mse[1]))

# ── Save outputs ──────────────────────────────────────────────────────────────

cat("\nSaving results...\n")
write_csv(results_all_conf,    '../Data/grid_search_results_conf_alpha.csv')
write_csv(pattern_scores_conf, '../Data/grid_search_scores_conf_alpha.csv')

cat("Done. Outputs saved to ../Data/\n")
