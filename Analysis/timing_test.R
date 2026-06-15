## timing_test.R
## Estimates how long the full configural cue grid search will take.
## Run on the target machine and check the final printed estimate.
## Change N_WORKERS to match the number of cores available on that machine.

N_WORKERS = 20   # <-- set this to the number of available cores
N_SAMPLE  = 40   # combinations to time (keep at 40 for a reliable estimate)
N_SIMS    = 50   # must match the real grid search

# Run from the Analysis/ directory:
#   Rscript timing_test.R
# Or set the path explicitly:
# setwd("path/to/Generalizability-Type-Token/Analysis")

library(tidyverse)
library(ndl)
library(arm)
library(furrr)
library(future)
source("RescorlaLogistic.R")

# ── Replicate data setup from grid_search.Rmd ────────────────────────────────

load_cond = function(path) {
  d = read_csv(path, show_col_types = FALSE) %>%
    mutate(stem = orthoCoding(stem),
           type = case_when(
             suffix %in% c('dan', 'sil') ~ 'big_pl',
             suffix == 'none'             ~ 'big_sg',
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

data_cond1_configural = load_cond('../Conditions/Condition1.csv') %>% add_configural()
data_cond2_configural = load_cond('../Conditions/Condition2.csv') %>% add_configural()
data_cond3_configural = load_cond('../Conditions/Condition3.csv') %>% add_configural()
data_cond4_configural = load_cond('../Conditions/Condition4.csv') %>% add_configural()
data_cond5_configural = load_cond('../Conditions/Condition5.csv') %>% add_configural()
data_cond6_configural = load_cond('../Conditions/Condition6.csv') %>% add_configural()

human_data = read_csv('../Data/fixed_effects_df_all.csv', show_col_types = FALSE) %>%
  filter(Task == 'prod')

# ── run_rw_production_conf (exact copy from grid_search.Rmd) ─────────────────

run_rw_production_conf = function(
  n_sims = 5,
  Alpha, SemAlpha, ConfAlpha, Beta,
  configural_list,
  prod_file = '../Conditions/production_condition.csv',
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

  data_prod = read_csv(prod_file, show_col_types = FALSE) %>%
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

  data_prod = map_dfr(1:n_sims, ~ mutate(data_prod, sim = .x))

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
    left_join(hum2, by = join_by(meaning    == Meaning,
                                 condition  == Condition,
                                 prod_condition == Stem)) %>%
    mutate(mse = (mean_freq_diff_v2 - Estimate)^2)
}

# ── Timing run ───────────────────────────────────────────────────────────────

N_TOTAL = 21^4  # 194,481 combinations in the full grid

set.seed(964)
full_grid = tidyr::crossing(
  ConfAlpha = seq(0, 1, by = 0.05),
  SemAlpha  = seq(0, 1, by = 0.05),
  Beta      = seq(0, 1, by = 0.05),
  Alpha     = seq(0, 1, by = 0.05)
)
sample_grid = full_grid[sample(nrow(full_grid), N_SAMPLE), ]

plan(multisession, workers = N_WORKERS)

cat(sprintf("Timing %d combinations across %d workers (n_sims = %d)...\n",
            N_SAMPLE, N_WORKERS, N_SIMS))
t_start = proc.time()

sample_grid %>%
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
          ),
          hum = human_data
      ),
      .options = furrr_options(packages = c("tidyverse", "ndl", "arm"),
                               seed = 964)
    )
  )

t_elapsed = (proc.time() - t_start)[["elapsed"]]
plan(sequential)

secs_per_combo = t_elapsed / N_SAMPLE * N_WORKERS
estimated_secs = secs_per_combo * N_TOTAL / N_WORKERS
estimated_hrs  = estimated_secs / 3600

cat(sprintf("\nSample wall-clock time    : %.1f s for %d combos (%d workers)\n",
            t_elapsed, N_SAMPLE, N_WORKERS))
cat(sprintf("Time per combo (1 worker) : %.2f s\n", secs_per_combo))
cat(sprintf("Estimated total time      : %.1f s  (~%.1f hours)\n",
            estimated_secs, estimated_hrs))
