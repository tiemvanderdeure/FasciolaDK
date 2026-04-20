library(INLA)
library(spdep)
library(dplyr)

library(geodata)
geodata_path("/home/tvd/data/rgeodata")
dk <- gadm("DNK", level = 0)
dkbdry <- inla.sp2segment(sf::st_as_sf(dk))

data <- read.csv("/home/tvd/K/FasciolaDK/predictors.csv")
data$øko = as.numeric(data$øko == "true")
data$besid <- dense_rank(data$BES_ID)
data$logsize <- log(data$bes_size)
data$prevalence <- data$positive / data$count
data$Intercept <- 1 # to avoid any shenanigans with intercepts

coords <- as.matrix(data[, c("LON", "LAT")])
mesh <- inla.mesh.2d(
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
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE),
  control.inla = list(strategy = "adaptive"), # for speed up
) 

####### Test many models with many possible variable combinations
# normal variables - include or not
random_vars <- c(
  øko = "øko", 
  besid = 'f(besid, model = "iid", hyper = list(prec = list(prior = "pc.prec", param = c(1.0, 0.10))))',
  slagtid = 'f(slagteri_id, model = "iid")'
)

vars <- c("tavg", "ppt", "soil", "def")
seasons <- c("spring_lag1", "summer_lag1", "autumn_lag1", "winter", "spring", "summer")
env_vars <- c(
  as.vector(outer(vars, seasons, function(x,y) paste(x,y,sep="_"))), 
  "ollerenshaw", "ollerenshaw_lag1", "1"
)

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
        param = c(0.8, 0.8)
    ))
  ))' # prior shrinks towards bigger Rho = smoother over time
)

# Generate all combinations of binary variables
binary_combos <- expand.grid(rep(list(c(FALSE, TRUE)), length(random_vars)))

# Combine binary and multi-option combinations
model_list <- data.frame(merge(binary_combos, names(vars_spacetime)))
colnames(model_list) <- c(names(random_vars), "spacetime")

model_list <- merge(model_list, env_vars)

# Generate formulas
model_list$formula <- apply(model_list, 1, function(row) {
  terms <- c()
    # Add binary variables if included
  for (var in names(random_vars)) {
    if (row[[var]]) {
      terms <- c(terms, random_vars[[var]])
    }
  }
  # Add the selected spacetime term
  terms <- c(terms, row[["y"]], vars_spacetime[[row[["spacetime"]]]])
  # Build formula
  as.formula(paste("y ~ ", paste(terms, collapse = " + "), "-1 + Intercept"))
})

# Actually run all models
model_results <- sapply(model_list$formula, runinla)
# Get estimates for goodness of fit
model_list$mlik <- sapply(model_results, function(x) x$mlik[2])
model_list$waic <- sapply(model_results, function(x) x$waic$waic)
model_list$dic <- sapply(model_results, function(x) x$dic$dic)
model_list$dic_sat <- sapply(model_results, function(x) x$dic$dic.sat)

write.csv(
  model_list[,-which(names(model_list) == "formula")], 
  paste0("model_runs_", Sys.Date(), ".csv"), row.names = FALSE
)

#### Write CSV file with posteriors
# convert to named vector with names "columnname_rowname"
mat_to_vect <- function(m) setNames(as.vector(m),
              paste0(rep(colnames(m), each = nrow(m)),
                     "_",
                     rep(rownames(m), times = ncol(m))))

# get posterior variance from precision
parse_precision <- function(res, name){
  m <- res$internal.marginals.hyperpar[[paste0("Log precision for ", name)]]
  m.var <- inla.tmarginal(function(x) 1/exp(x), m)
  post_var <- post_var <- inla.zmarginal(m.var, silent =TRUE)
  post_var <- unlist(post_var[c("mean", "quant0.025", "quant0.975")])
  names(post_var) <- paste0(gsub("_", "", name), "_", names(post_var))
  post_var
}

posteriors_summarized <- mapply(function(res, m){
  # Fixed posteriors
  posterior_fixed <- apply(res$summary.fixed, 1, function(x) c(x["mean"], x["0.025quant"], x["0.975quant"]))
  colnames(posterior_fixed) <- replace(colnames(posterior_fixed), colnames(posterior_fixed)==m$y, "envvar")
  rownames(posterior_fixed) <- c("mean", "quant0.025", "quant0.975")
  posterior_fixed_v <- mat_to_vect(posterior_fixed)

  # i.i.d. posteriors
  besid_var <- if(m$besid) parse_precision(res, "besid") else c()
  slagtid_var <- if(m$slagtid) parse_precision(res, "slagteri_id") else c()
  vars <- c(posterior_fixed_v, besid_var, slagtid_var)
},
  model_results, 
  asplit(model_list, 1)
)

posteriors_fixed_df <- data.frame(bind_rows(posteriors_summarized))

model_results_full <- cbind(
  model_list[,-which(names(model_list) == "formula")], 
  posteriors_fixed_df
)

write.csv(
  model_results_full, 
  paste0("model_runs_posteriors", Sys.Date(), ".csv"), row.names = FALSE, na = ""
)

##################
# Explore combinations of variables

# Identified by earlier exploratory analysis as best candidates 
candidate_vars <- c(
  "ppt_summer_lag1", "tavg_summer_lag1", "soil_summer_lag1", 
  "soil_autumn_lag1", "def_spring_lag1") 
  
# all combos of size = 2
candidate_combos <- apply(combn(candidate_vars, 2), 2, paste, collapse = " + ")

forms <- sapply(candidate_combos, function(x) as.formula(
  paste("y ~ ", 
  vars_spacetime["spacetime"], "+", 
  paste(random_vars, collapse = "+"), 
  "-1 + Intercept + ", x)))

results <- lapply(forms, runinla)
result_summary <- lapply(results, function(res) {
  # fixed effects: pick mean, 0.025 and 0.975 quantiles for each fixed-effect row
  fixed_names <- rownames(res$summary.fixed)
  posterior_fixed <- sapply(fixed_names, function(v) {
    c(res$summary.fixed[v, "mean"],
      res$summary.fixed[v, "0.025quant"],
      res$summary.fixed[v, "0.975quant"])
  }, simplify = "array")

  # make rownames consistent with earlier code
  rownames(posterior_fixed) <- c("mean", "quant0.025", "quant0.975")
  colnames(posterior_fixed) <- c("øko", "Intercept", "envvar1", "envvar2")
  posterior_fixed_v <- mat_to_vect(posterior_fixed)

  # i.i.d. random effects (all your candidate models include these)
  besid_var <- parse_precision(res, "besid")
  slagtid_var <- parse_precision(res, "slagteri_id")

  # combine into a single named numeric vector
  c(posterior_fixed_v, besid_var, slagtid_var)
}
)
results_df <- bind_rows(result_summary)
results_df$formula <- candidate_combos


results_df$mlik <- sapply(results, function(x) x$mlik[2])
results_df$waic <- sapply(results, function(x) x$waic$waic)
results_df$dic <- sapply(results, function(x) x$dic$dic)
results_df$dic_sat <- sapply(results, function(x) x$dic$dic.sat)

which(results_df$mlik == max(results_df$mlik))
which(results_df$waic == min(results_df$waic))

write.csv(results_df, paste0("candidate_models_two_vars_", Sys.Date(), ".csv"), row.names = FALSE)



#############
## Write for spatio-temporal effect
spacetime_formula_base <- as.formula(paste(
  "y ~ ", paste(random_vars, collapse = " + "), "-1 + Intercept + ", vars_spacetime["spacetime"]))

optimal_model <- update(spacetime_formula_base, y ~ . + ppt_summer_lag1 + tavg_summer_lag1)

# Run inla with the model with besid, øko and spatiotemporal effect
#spacetime_formula <- last(model_list$formula)
res <- runinla(optimal_model)

# export for estimates of abattoir random effect
write.csv(res$summary.random$slagteri_id, "slagteri_posterior.csv", row.names = FALSE)

#### Export posteroir spatiotemporal component
library(fmesher)
library(terra)
library(ncdf4)

nsamples <- 500
samples <- inla.posterior.sample(nsamples, res)

space_indices <- grep("space", rownames(samples[[1]]$latent))
space_samples <- lapply(samples, function(s) s$latent[space_indices])

dkrast <- rast("data/terraclimate_dk.nc", lyrs = 1)
dkpoints <- crds(dkrast, na.rm = TRUE)
ev <- fm_evaluator(mesh, loc = dkpoints)

rast_template <- rast(dkrast, nlyrs = 14)
namask <- !is.na(dkrast)

for (i in 1:nsamples){
  print(i)
  mmat <- matrix(space_samples[[i]], nrow = mesh$n)
  vals <- fm_evaluate(ev, mmat)
  rast_template[namask] <- c(vals)
  writeRaster(rast_template, paste0("/home/tvd/K/fasciolaDK/posterior_spatiotemporal", i, ".tif"), overwrite = TRUE)
}

## Write estimates for the effects of the two environmental variables in each sample
vars_to_export <- c("øko", "Intercept", "ppt_summer_lag1", "tavg_summer_lag1")
effects <- lapply(samples, function(s) s$latent[paste0(vars_to_export, ":1"),])
effects_df <- data.frame(do.call(rbind, effects))
colnames(effects_df) <- vars_to_export
write.csv(effects_df, paste0("posterior_effects_envvars_samples", Sys.Date(), ".csv"), row.names = FALSE)

####### Effect of year
# Run models with year as linear predictor, with and without environmental variables
formula_year <- update(spacetime_formula_base, y ~ . + year)
optimal_formula_year <- update(optimal_model, y ~ . + year)

res_optimal_year <- runinla(optimal_formula_year)
res_year <- runinla(formula_year)

# Export posteriors
year_posterior <- rbind(res_year$summary.fixed["year",], res_optimal_year$summary.fixed["year",])[c("mean", "0.025quant", "0.975quant")]
rownames(year_posterior) <- c("year_simple", "year_with_environment")

write.csv(year_posterior, paste0("posteriors_year", Sys.Date(), ".csv"))
