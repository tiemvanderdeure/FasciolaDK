library(INLA)
library(geodata)
library(spdep)
library(dplyr)
# 
geodata_path("/home/tvd/data/rgeodata")

dk <- gadm("DNK", level = 0)
dkbdry <- inla.sp2segment(sf::st_as_sf(dk))

data <- read.csv("/home/tvd/K/predictors.csv")
data$øko = data$øko == "true"
data$besid <- dense_rank(data$BES_ID)
data$logsize <- log(data$bes_size)
data$prevalence <- data$positive / data$count

#data <- data[data$besid <= 200,]

coords <- as.matrix(data[, c("LON", "LAT")])
mesh <- inla.mesh.2d(
  loc = coords, 
  boundary = dkbdry, 
  max.edge = c(0.3, 1.0), 
  cutoff = 0.15
)

### Prepare data into INLA data structures
startyr = min(data$yr)
nyr = max(data$yr) - startyr + 1
spde <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = c(0.1, 0.2),
                            prior.sigma = c(1, 0.2))
# SPDE with space and time component
spde_index <- inla.spde.make.index("space", spde$n.spde, n.group = nyr)
A <- inla.spde.make.A(mesh = mesh, loc = coords, group = data$yr - startyr + 1, n.group = nyr)

stack <- inla.stack(
  data = list(y = data$positive, n = data$count),
  A = list(A, 1),
  effects = list(
    list(space = spde_index$space, year = spde_index$space.group),
    data
  )
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
  tempanomaly = "tavg_winter + tavg_spring + tavg_summer",
  climateanomaly = "tavg_winter + tavg_spring + tavg_summer + ppt_winter + ppt_spring + ppt_summer"
)
# spacetime - multiple options
vars_spacetime <- c(
  none = '1',
  time = 'f(year, model = "ar1")',
  space = 'f(space, model = spde)',
  space_and_time = 'f(space, model = spde) + f(year, model = "ar1")',
  spacetime = 'f(space, model = spde, group = year, control.group = list(model = "ar1"))'
)

# Generate all combinations of binary variables
binary_combos <- expand.grid(rep(list(c(FALSE, TRUE)), length(vars)))

# Combine binary and multi-option combinations
model_list <- data.frame(merge(binary_combos, names(vars_spacetime)))
colnames(model_list) <- c(names(vars), "spacetime")

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
  terms <- c("1", terms, vars_spacetime[[row[["spacetime"]]]])
  # Build formula
  as.formula(paste("y ~ ", paste(terms, collapse = " + ")))
})

model_results <- sapply(model_list$formula, runinla)
model_list$mlik <- sapply(model_results, function(x) x$mlik[2])
model_list$waic <- sapply(model_results, function(x) x$waic$waic)
model_list$dic <- sapply(model_results, function(x) x$dic$dic)
model_list$dic_sat <- sapply(model_results, function(x) x$dic$dic.sat)

arrange(model_list, desc(mlik))

write.csv(model_list[,-6], "model_runs_3107.csv", row.names = FALSE)


