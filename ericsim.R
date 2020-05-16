library(DSsim)


## ----region, warning=FALSE, message=FALSE, fig.width=4, fig.cap="The study region."--------
# Create a polgon
poly1 <- data.frame(x = c(0,0,2800,2800,0), y = c(0,500,500,0,0))

# Create an empty list
coords <- list()
# Store the polygon inside a list in the first element of the coords list referring to strata 1.
coords[[1]] <- list(poly1)

# Create the survey region
region <- make.region(region.name = "study area", 
                      units = "m",
                      coords = coords)
# The plot function allows plotting in km or m.
#plot(region, plot.units = "m")



## ----density, warning=FALSE, message=FALSE, fig.width=4, fig.cap="The density surface."----
# Create the density surface
density <- make.density(region.obj = region, 
                        x.space = 50, 
                        y.space = 200, 
                        constant = 1)

## ----popdesc1, warning=FALSE, message=FALSE------------------------------------------------
# # Create the population description, with a population size N = 200
# pop.desc <- make.population.description(region.obj = region, 
#                                             density.obj = density, 
#                                             N = 200)


## ----popdesc2, warning=FALSE, message=FALSE------------------------------------------------
# Create the covariate list
covariate.list <- list()
# The population will be 50% males and 50% females
covariate.list$sex <- list(data.frame(level = c("female", "male"), prob = c(0.3,0.7)))

# Create the population description, with a population size N = 200
pop.desc.cov <- make.population.description(region.obj = region, 
                                            density.obj = density, 
                                            covariates = covariate.list, 
                                            N = 180)


## ----detect1, warning=FALSE, message=FALSE, fig.width=4, fig.cap="The detection functions for males and females."----
# Make a simple half normal detection function with a scale parameter of 200
# detect.hn <- make.detectability(key.function = "hn",
#                                  scale.param = 20, 
#                                  truncation = 100)
# # We can now visualise these detection functions
# plot(detect.hn, pop.desc)


## ----detect2, warning=FALSE, message=FALSE, fig.width=4, fig.cap="The detection functions for males and females."----
# Create the covariate parameter list
cov.params <- list()
# Note the covariate parameters are supplied on the log scale
cov.params$sex = data.frame(level = c("female", "male"), 
                            param = c(0, 1))

detect.cov <- make.detectability(key.function = "hn" ,
                                 scale.param = 12,
                                 cov.param = cov.params, 
                                 truncation = 80)

# This setup gives a scale parameter of around 120 for the females and 540 for 
# the males. We can calculate the sigma for the males using the formula above:
# exp(log(scale.param) + sex.male)
exp(log(12) + .1)   
# We can now visualise these detection functions
plot(detect.cov, pop.desc.cov)


## ----design, warning=FALSE, message=FALSE--------------------------------------------------

# Define the design
design <- make.design(transect.type = "line",
                      design.details = c("parallel", "systematic"),
                      region.obj = region,
                      spacing = 220)


## ----transects, warning=FALSE, message=FALSE, echo = TRUE, fig.width=4, fig.cap="Example survey transects."----
transects <- generate.transects(design, region = region)
plot(region, plot.units = "m")
plot(transects, col = 4, lwd = 2)


## ----analyses1-----------------------------------------------------------------------------
ddf.analyses <- make.ddf.analysis.list(dsmodel = list(~cds(key = "hn", formula = ~1),
                                                      ~cds(key = "hr", formula = ~1)), 
                                       method = "ds",
                                       criteria = "AIC",
                                       truncation = 80)


## ----analyses2-----------------------------------------------------------------------------
ddf.analyses.cov <- make.ddf.analysis.list(dsmodel = list(~mcds(key = "hn", formula = ~sex)), 
                                           method = "ds",
                                           truncation = 80)


## ----set.seed, echo = FALSE, eval = TRUE---------------------------------------------------
set.seed(474)


# ## ----checksim, fig.height=5.5, fig.width=7.2, fig.cap="Example survey. Top left - the density suface with an example population. Top right - an example set of transects. Bottom left - the detections from the transects. Bottom right - A histogram of the distances from these observations to the transect it was detected."----
# sim <- make.simulation(reps = 99, 
#                        region.obj = region,
#                        design.obj = design,
#                        population.description.obj = pop.desc,
#                        detectability.obj = detect.hn,
#                        ddf.analyses.list = ddf.analyses)
# # Produce simulation setup plots
# check.sim.setup(sim)


## ----checksim2, fig.height=5.5, fig.width=7.2, fig.cap="Example survey."-------------------
sim.cov <- make.simulation(reps = 99, 
                       region.obj = region,
                       design.obj = design,
                       population.description.obj = pop.desc.cov,
                       detectability.obj = detect.cov,
                       ddf.analyses.list = ddf.analyses)
# Produce simulation setup plots
check.sim.setup(sim.cov)


fredcov <- create.survey.results(sim.cov)
billcov <- get.distance.data(fredcov)
library(Distance)
billcov$Region.Label <- "OneStratum"
billcov$Area <- region@area
billcov$Effort <- unname(region@box[4])
names(billcov)[2] <- "Sample.Label"
newbill <- billcov[billcov$object %% 1 == 0, ]
notrun <- ds(newbill)
notrun$dht$individuals$N
plot(notrun, main="No truncation")
trun60 <- ds(newbill, truncation = 60)
trun60$dht$individuals$N
plot(trun60, main="Trunc 60")
covtrunc <- ds(newbill, truncation=60, formula=~sex)
covtrunc$dht$individuals$N
plot(covtrunc, main="Covar truncation")
covtrunc.hr <- ds(newbill, truncation=60, key="hr")
covtrunc.hr$dht$individuals$N
covtruncsex.hr <- ds(newbill, truncation=60, key="hr", formula=~sex)
covtruncsex.hr$dht$individuals$N
AIC(trun60, covtrunc, covtrunc.hr, covtruncsex.hr)
trun60$dht$individuals$summary
un <- ds(newbill, key="uni", adj="cos", truncation=60)
un$dht$individuals$N
summarize_ds_models(un, trun60, covtrunc, covtrunc.hr, covtruncsex.hr, output="plain")

colorder <- c("Region.Label", "Area", "Sample.Label", "Effort", "distance", "sex")
newbill$distance <- round(newbill$distance)
newbill$Area <- newbill$Area/1e6  # m^2 to km^2
newbill$Effort <- newbill$Effort/1000  # m to km
dataset <- newbill[order(newbill$Sample.Label) , colorder]
write.csv(dataset, file="mystery.csv", row.names = FALSE, quote=FALSE)
