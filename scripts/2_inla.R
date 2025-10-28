library(INLA)
library(geodata)
library(spdep)
library(dplyr)
# 
geodata_path("/home/tvd/data/rgeodata")

dk <- gadm("DNK", level = 0)
dkbdry <- inla.sp2segment(sf::st_as_sf(dk))

data <- read.csv("/home/tvd/K/predictors.csv")
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
  besid = 'f(besid, model = "iid")',
  slagtid = 'f(slagteri_id, model = "iid")'
  #tempanomaly = "tavg_winter + tavg_spring + tavg_summer",
  #precanomaly = "ppt_winter + ppt_spring + ppt_summer"
)

vars <- c("tavg", "ppt", "soil", "def")
seasons <- c("winter", "spring", "summer", "autumn")
env_vars <- as.vector(outer(vars, seasons, function(x,y) paste(x,y,sep="_")))

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
  for (var in names(vars)) {
    if (row[[var]]) {
      terms <- c(terms, vars[[var]])
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
  head(model_list[,-which(names(model_list) == "formula")]), 
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


#############
## Write for spatio-temporal effect
spacetime_formula <- y ~ øko + f(besid, model = "iid") + f(slagteri_id, model = "iid") + 
    1 + f(space, model = spde, group = year, control.group = list(model = "ar1", 
    hyper = list(rho = list(prior = "pc.cor1", param = c(0.8, 
        0.8))))) - 1 + Intercept


# Run inla with the model with besid, øko and spatiotemporal effect
#spacetime_formula <- last(model_list$formula)
res <- runinla(spacetime_formula)

write.csv(res$summary.random$slagteri_id, "slagteri_posterior.csv", row.names = FALSE)

samples <- inla.posterior.sample(500, res)
tail(samples[[1]]$latent)

predict(res, fm_pixels(mesh))

#### Attempts to properly export posteroir spaciotemporal component
library(fmesher)
library(terra)
library(ncdf4)

nsamples <- 100


sapply(samples, function(s) s$hyperpar["Precision for besid"])

samples <- inla.posterior.sample(nsamples, res)
space_indices <- grep("space", rownames(samples[[1]]$latent))
space_samples <- lapply(samples, function(s) s$latent[space_indices])
intercepts <- sapply(samples, function(s) s$latent[nrow(s$latent)])

dkrast <- rasterize(dk, rast(dk, res = 0.05), touches = TRUE)
dkpoints <- crds(dkrast, na.rm = TRUE)
ev <- fm_evaluator(mesh, loc = dkpoints)
ev2 <- fm_evaluator(mesh, loc = crds(vect(data[1:1, c("LON", "LAT")], geom = c("LON", "LAT"))))

rast_template <- rast(dkrast, nlyrs = 14)
namask <- !is.na(dkrast)

for (i in 1:nsamples){
  mmat <- matrix(space_samples[[i]], nrow = mesh$n)
  vals <- fm_evaluate(ev, mmat) + intercepts[[i]]
  rast_template[namask] <- c(vals)
  writeRaster(rast_template, paste0("/home/tvd/K/posterior_spatiotemporal", i, ".tif"), overwrite = TRUE)
}


###### Playground
# No clear effect of year
r5 <- runinla(y ~ 1 + f(space, model = spde) - 1 + Intercept)

r2 <- runinla(y ~ øko)
samples2 <- inla.posterior.sample(2, r2)


r3 <- runinla(y ~ 0 + øko)
samples2 <- inla.posterior.sample(2, r4)

# Including slagteri_id does improve the fit!
res <- runinla(y ~ øko + yr + f(space, model = spde) + f(year, model = "ar1") + f(slagteri_id, model = "iid"))

sample = samples[[1]]
predictor_indices <- grep("^Predictor", rownames(sample$latent))
apredictor_indices <- grep("APredictor", rownames(sample$latent))
space_indices <- grep("space", rownames(sample$latent))
slagteri_indices <- grep("slagteri_id", rownames(sample$latent))
bes_indices <- grep("besid", rownames(sample$latent))

r2$summary.fitted.values

head(sample$latent[predictor_indices], 20)
head(sample$latent[apredictor_indices], 20)

pred <- sample$latent[apredictor_indices]
space_posterior <- sample$latent[space_indices]
slagteri_posterior <- sample$latent[slagteri_indices]
bes_posterior <- sample$latent[bes_indices]

intercept <- sample$latent[length(sample$latent)]
øko <- sample$latent[length(sample$latent)-1]

data$mean_pred <- r3$summary.fitted.values[idx$data, "mean"]
res$summary.random$slagteri_id[, "mean"]

res$summary.fitted.values[idx$data[9228],]
res$summary.random$besid[data$besid[9228],]
res$summary.random$slagteri_id[data$slagteri_id[9228]+1,]

data$mean_pred <- pred
data$besid_pred <- bes_posterior[data$besid]
data$slagteri_pred <- slagteri_posterior[data$slagteri_id+1]
data$spat_pred <- (A %*% space_posterior)[,1]

data$spat_pred + data$slagteri_pred + data$besid_pred - data$mean_pred

p.pred<-exp(post.mean.pred.logit)/(1 + exp(post.mean.pred.logit))

rownames(res$summary.random$space)
nrow(res$summary.random$space)



head(data)

predictor_indices[idx$data]
