# install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
# install.packages(c("spdep", "dplyr", "tidyr", "geodata", "loo"))
library(INLA)
library(spdep)
library(dplyr)
library(tidyr)
library(fmesher)
library(geodata)
library(loo)
geodata_path("/home/tvd/data/rgeodata")
dk <- gadm("DNK", level = 0)
dkbdry <- inla.sp2segment(sf::st_as_sf(dk))

data <- read.csv("/home/tvd/K/FasciolaDK/predictors.csv")
data$øko = as.numeric(data$øko == "true")
data$chr <- dense_rank(data$CHR_ID)
data$prevalence <- data$positive / data$count
data$Intercept <- 1 # to avoid any shenanigans with intercepts

coords <- as.matrix(data[, c("LON", "LAT")])
mesh <- fm_mesh_2d_inla(
  loc = coords,
  boundary = dkbdry,
  max.edge = c(0.5, 2.0),
  cutoff = 0.2 # coarser mesh results in better WAIC/DIC scores
)

### Prepare data into INLA data structures
startyr = min(data$yr)
nyr = max(data$yr) - startyr + 1
spde <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = c(0.5, 0.1), # push towards bigger range = smoothness in space
                            prior.sigma = c(1, 0.5)) # allow fairly large std
# SPDE with space and time component
spde_index <- inla.spde.make.index("space", spde$n.spde, n.group = nyr)
A <- inla.spde.make.A(mesh = mesh, loc = coords, group = data$yr - startyr + 1, n.group = nyr)

stack <- inla.stack(
  data = list(y = data$positive, n = data$count),
  A = list(A, 1),
  effects = list(
    list(space = spde_index$space, year = spde_index$space.group),
    data
  ),
  tag = "fit"
)

# Convenience function to call inla with some settings and the data set
runinla <- function(x) inla(
  x,
  data = inla.stack.data(stack),
  family = "binomial",
  Ntrials = n,
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE),
  control.inla = list(strategy = "simplified.laplace"),
)

####### Test many models with many possible variable combinations
# normal variables - include or not
random_vars <- c(
  øko = "øko",
  chr = 'f(chr, model = "iid")',
  slagtid = 'f(slagteri_id, model = "iid")',
  time = "year"
)

vars <- c("tavg", "ppt", "soil", "def")
seasons <- c("spring_year1", "summer_year1", "autumn_year1", "winter_year2", "spring_year2", "summer_year2")
env_vars <- c(
  as.vector(outer(vars, seasons, function(x,y) paste(x,y,sep="_"))),
  "ollerenshaw_year1", "ollerenshaw_year2", "1"
)
env_vars <- setdiff(env_vars, c("soil_winter_year2", "def_winter_year2")) # these have basically no pattern at all

# spacetime - multiple options
vars_spacetime <- c(
  none = '1',
  time = 'f(year, model = "ar1")',
  space = 'f(space, model = spde)',
  space_and_time = 'f(space, model = spde) + f(year, model = "ar1")',
  spacetime = 'f(space, model = spde, group = year, control.group = list(
    model = "ar1", 
    hyper = list(rho = list(
        prior = "pc.cor1", 
        param = c(0.0, 0.8)
    ))
  ))' # prior shrinks towards bigger Rho = smoother over time
)

# Generate all combinations of binary variables
binary_combos <- rbind(
  # øko,  chr,   slagtid, time
  c(TRUE, FALSE, TRUE, FALSE, "none"), 
  c(TRUE, FALSE, TRUE, FALSE, "space"),
  c(TRUE, TRUE, TRUE, FALSE, "space"),
  c(TRUE, TRUE, TRUE, TRUE, "space"),
  c(TRUE, TRUE, TRUE, TRUE, "spacetime"))

# Combine binary and multi-option combinations
colnames(binary_combos) <- c(names(random_vars), "spacetime")
backbones <- apply(binary_combos, 1, function(row) {
  terms <- c()
  for (var in names(random_vars)) {
    if (row[[var]]) {
      terms <- c(terms, random_vars[[var]])
    }
  }
  terms <- c(terms, vars_spacetime[[row[["spacetime"]]]])
  as.formula(paste("y ~ ", paste(terms, collapse = " + "), "-1 + Intercept"))
})
binary_combos <- as.data.frame(binary_combos, stringsAsFactors = FALSE)

model_list <- merge(binary_combos, env_vars) # cross join; unnamed env_vars column is auto-named "y"

# update backbone with formula
model_list$formula <- lapply(seq_len(nrow(model_list)), function(i) {
  v <- model_list$y[i]
  bb <- backbones[[((i-1) %% 5)+1]]
  if (v == "1") return(bb)
  update(bb, as.formula(paste("y ~ . +", v)))
})

# Actually run all models
model_results <- lapply(model_list$formula, runinla)
# Get estimates for goodness of fit
model_list$mlik <- sapply(model_results, function(x) x$mlik[2])
model_list$waic <- sapply(model_results, function(x) x$waic$waic)
model_list$dic <- sapply(model_results, function(x) x$dic$dic)
model_list$dic_sat <- sapply(model_results, function(x) x$dic$dic.sat)

# calculate WAIC and CPO differences and SE relative to equivalent model w/o environmental var
model_list$row_id <- seq_len(nrow(model_list))
model_list <- model_list %>%
  group_by(øko, chr, slagtid, time, spacetime) %>%
  mutate(
    backbone_row_id = row_id[y == "1"],
    waic_diff = waic - waic[y == "1"],
    dic_diff = dic - dic[y == "1"]
  ) %>%
  ungroup() %>%
  as.data.frame()

model_list$waic_diff_se <- mapply(function(i, j) {
  pointwise_diff <- model_results[[i]]$waic$local.waic - model_results[[j]]$waic$local.waic
  sqrt(length(pointwise_diff) * var(pointwise_diff))
}, model_list$row_id, model_list$backbone_row_id)


cpo_diffs <- mapply(function(i, j) {
  res_i <- model_results[[i]]
  res_j <- model_results[[j]]
  bad <- res_i$cpo$failure > 0 | res_j$cpo$failure > 0
  d <- log(res_i$cpo$cpo) - log(res_j$cpo$cpo)
  bad <- bad | !is.finite(d)
  d[!bad]
  c(cpo_diff = sum(d), cpo_diff_se = sqrt(length(d) * var(d)), bad_idx = list(which(bad)))
}, model_list$row_id, model_list$backbone_row_id)

model_list$cpo_diff <- unlist(cpo_diffs["cpo_diff", ])
model_list$cpo_diff_se <- unlist(cpo_diffs["cpo_diff_se", ])
#model_list$cpo_bad <- cpo_diffs["bad_idx", ]

#### Write CSV file with full posteriors for every one of the 1080 models
# convert to named vector with names "columnname_rowname"
mat_to_vect <- function(m) setNames(as.vector(m),
              paste0(rep(colnames(m), each = nrow(m)),
                     "_",
                     rep(rownames(m), times = ncol(m))))

# get posterior variance from precision
parse_precision <- function(res, name){
  m <- res$internal.marginals.hyperpar[[paste0("Log precision for ", name)]]
  m.var <- inla.tmarginal(function(x) 1/exp(x), m)
  post_var <- inla.zmarginal(m.var, silent = TRUE)
  post_var <- unlist(post_var[c("mean", "quant0.025", "quant0.975")])
  names(post_var) <- paste0(gsub("_", "", name), "_", names(post_var))
  post_var
}

# Build an n_obs x n_models matrix of LPD's from INLA's CPO, for loo weights
cpo_lpd_matrix <- function(results_list) {
  failure_flags <- sapply(results_list, function(x) x$cpo$failure > 0)  # n_obs x K
  failure_counts <- colSums(failure_flags, na.rm = TRUE)
  if (any(failure_counts > 0)) {
    warning(
      sum(failure_counts > 0), " of ", length(results_list), " model(s) have a nonzero ",
      "cpo$failure count (", paste(failure_counts[failure_counts > 0], collapse = ", "),
      "); excluding the flagged observations from the weight computation for every model."
    )
  }
  lpd <- sapply(results_list, function(x) log(x$cpo$cpo))
  bad_obs <- which(rowSums(!is.finite(lpd)) > 0 | rowSums(failure_flags, na.rm = TRUE) > 0)
  if (length(bad_obs) > 0) {
    warning(
      length(bad_obs), " observation(s) have a non-finite or failure-flagged CPO in at ",
      "least one model; dropping them from the weight computation."
    )
    lpd <- lpd[-bad_obs, , drop = FALSE]
  }
  lpd
}

posteriors_summarized <- mapply(function(res, m){
  # Fixed posteriors
  posterior_fixed <- apply(res$summary.fixed, 1, function(x) c(x["mean"], x["0.025quant"], x["0.975quant"]))
  colnames(posterior_fixed) <- replace(colnames(posterior_fixed), colnames(posterior_fixed) == m$y, "envvar")
  rownames(posterior_fixed) <- c("mean", "quant0.025", "quant0.975")
  posterior_fixed_v <- mat_to_vect(posterior_fixed)

  # i.i.d. posteriors
  chr_var <- if (m$chr) parse_precision(res, "chr") else c()
  slagtid_var <- if (m$slagtid) parse_precision(res, "slagteri_id") else c()
  c(posterior_fixed_v, chr_var, slagtid_var)
},
  model_results,
  asplit(model_list, 1)
)

posteriors_fixed_df <- data.frame(bind_rows(posteriors_summarized))
model_results_full <- cbind(
  model_list[, -which(names(model_list) == "formula")],
  posteriors_fixed_df
)
write.csv(
  model_results_full, 
  paste0("model_runs_posteriors", Sys.Date(), ".csv"), row.names = FALSE, na = ""
)

##################
# Posterior draws for environmental predictors for each model
n_draws <- 1000
draw_rows <- which(model_list$y != "1")
posterior_draws <- bind_rows(lapply(draw_rows, function(i) {
  v <- model_list$y[i]
  samples <- inla.posterior.sample(n_draws, model_results[[i]], selection = setNames(list(1), v))
  values <- sapply(samples, function(s) s$latent[paste0(v, ":1"), ])
  data.frame(
    øko = model_list$øko[i], chr = model_list$chr[i], slagtid = model_list$slagtid[i], time = model_list$time[i],
    spacetime = model_list$spacetime[i], y = v, draw = seq_len(n_draws), value = values
  )
})
)

write.csv(posterior_draws, paste0("posterior_draws_envvars_", Sys.Date(), ".csv"), row.names = FALSE)

###############
# Two-variable models
top_envvars <- model_list %>%
  filter(y != "1") %>%
  group_by(y) %>%
  summarise(waic_diff_mean = mean(waic_diff)) %>%
  arrange(waic_diff_mean) %>%
  slice_head(n = 10) %>%
  pull(y)

backbone3 <- binary_combos[3, ]
backbone3_f <- backbones[[3]]
backbone3_row <- which(
  model_list$øko == backbone3$øko & model_list$chr == backbone3$chr &
    model_list$slagtid == backbone3$slagtid & model_list$time == backbone3$time &
    model_list$spacetime == backbone3$spacetime & model_list$y == "1"
)
backbone3_result <- model_results[[backbone3_row]]

var_pairs <- as.data.frame(t(combn(top_envvars, 2)), stringsAsFactors = FALSE)
colnames(var_pairs) <- c("envvar1", "envvar2")
two_var_list <- var_pairs
two_var_list$formula <- lapply(seq_len(nrow(two_var_list)), function(i) {
  update(backbone3_f, 
    as.formula(paste(". ~ ", two_var_list$envvar1[i], two_var_list$envvar2[i], sep = "+")))
})

two_var_results <- lapply(two_var_list$formula, runinla)

two_var_list$mlik <- sapply(two_var_results, function(x) x$mlik[2])
two_var_list$waic <- sapply(two_var_results, function(x) x$waic$waic)
two_var_list$dic <- sapply(two_var_results, function(x) x$dic$dic)
two_var_list$dic_sat <- sapply(two_var_results, function(x) x$dic$dic.sat)

# waic/dic/cpo differences relative to the backbone-3 null model (no env vars)
compute_diff_stats <- function(res, ref) {
  pointwise_diff <- res$waic$local.waic - ref$waic$local.waic
  waic_diff_se <- sqrt(length(pointwise_diff) * var(pointwise_diff))

  bad <- res$cpo$failure > 0 | ref$cpo$failure > 0
  d <- log(res$cpo$cpo) - log(ref$cpo$cpo)
  bad <- bad | !is.finite(d)
  d <- d[!bad]

  c(
    waic_diff = res$waic$waic - ref$waic$waic,
    dic_diff = res$dic$dic - ref$dic$dic,
    waic_diff_se = waic_diff_se,
    cpo_diff = sum(d),
    cpo_diff_se = sqrt(length(d) * var(d))
  )
}
diff_stats <- t(sapply(two_var_results, compute_diff_stats, ref = backbone3_result))
two_var_list <- cbind(two_var_list, diff_stats)

# posterior summaries (mean + 95% CI) for both env var coefficients, plus
# the backbone-3 iid random effects, same layout as the single-var export
two_var_posteriors <- mapply(function(res, row) {
  posterior_fixed <- apply(res$summary.fixed, 1, function(x) c(x["mean"], x["0.025quant"], x["0.975quant"]))
  colnames(posterior_fixed) <- replace(colnames(posterior_fixed), colnames(posterior_fixed) == row[["envvar1"]], "envvar1")
  colnames(posterior_fixed) <- replace(colnames(posterior_fixed), colnames(posterior_fixed) == row[["envvar2"]], "envvar2")
  rownames(posterior_fixed) <- c("mean", "quant0.025", "quant0.975")
  posterior_fixed_v <- mat_to_vect(posterior_fixed)

  chr_var <- parse_precision(res, "chr")
  slagtid_var <- parse_precision(res, "slagteri_id")
  c(posterior_fixed_v, chr_var, slagtid_var)
},
  two_var_results,
  asplit(two_var_list, 1),
  SIMPLIFY = FALSE
)

two_var_results_full <- cbind(
  two_var_list[, -which(names(two_var_list) == "formula")],
  data.frame(bind_rows(two_var_posteriors))
)
write.csv(
  two_var_results_full,
  paste0("model_runs_posteriors_twovars_", Sys.Date(), ".csv"), row.names = FALSE, na = ""
)

##################
# Posterior draws for the two env var coefficients in each two-variable model
two_var_draws <- bind_rows(lapply(seq_len(nrow(two_var_list)), function(i) {
  v1 <- two_var_list$envvar1[i]
  v2 <- two_var_list$envvar2[i]
  samples <- inla.posterior.sample(
    n_draws, two_var_results[[i]],
    selection = setNames(list(1, 1), c(v1, v2))
  )
  value1 <- sapply(samples, function(s) s$latent[paste0(v1, ":1"), ])
  value2 <- sapply(samples, function(s) s$latent[paste0(v2, ":1"), ])
  data.frame(
    envvar1 = v1, envvar2 = v2, draw = seq_len(n_draws),
    value1 = value1, value2 = value2
  )
}))

write.csv(two_var_draws, paste0("posterior_draws_twovars_", Sys.Date(), ".csv"), row.names = FALSE)