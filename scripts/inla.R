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
vars <- c(
  øko = "øko", 
  besid = 'f(besid, model = "iid")',
  slagtid = 'f(slagteri_id, model = "iid")'
  #tempanomaly = "tavg_winter + tavg_spring + tavg_summer",
  #precanomaly = "ppt_winter + ppt_spring + ppt_summer"
)
env_vars <- c(
  "tavg_winter", "tavg_spring", "tavg_summer",
  "ppt_winter", "ppt_spring", "ppt_summer",
  "1" # 1 = no environmental variables
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
binary_combos <- expand.grid(rep(list(c(FALSE, TRUE)), length(vars)))

# Combine binary and multi-option combinations
model_list <- data.frame(merge(binary_combos, names(vars_spacetime)))
colnames(model_list) <- c(names(vars), "spacetime")

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

# get posteriors of hyperparameters
m <- res$internal.marginals.hyperpar$`Log precision for besid`
m.var <- inla.tmarginal(function(x) 1/exp(x), m)
inla.zmarginal(m.var) 


m <- res$internal.marginals.hyperpar$`Log precision for slagteri_id`
m.var <- inla.tmarginal(function(x) 1/exp(x), m)
inla.zmarginal(m.var) 

model_results_order <- model_results[order(model_list$mlik, decreasing = TRUE)]
summary(model_results_order[[6]])

mean(abs(data$tavg_winter))

r <- runinla(y ~ 1 + f(besid, model = "iid") + øko + f(space, model = spde, 
    group = year, control.group = list(model = "ar1", hyper = list(rho = list(prior = "pc.cor1", 
        param = c(0.8, 0.8))))))

write.csv(model_list[,-5], paste0("model_runs_", Sys.Date(), ".csv"), row.names = FALSE)

# Run inla with the model with besid, øko and spatiotemporal effect
spacetime_formula <- last(model_list$formula)
res <- runinla(spacetime_formula)

samples <- inla.posterior.sample(10, res)
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

