library(DSsim)
## ----region, warning=FALSE, message=FALSE, fig.width=4, fig.cap="The study region."--------
# Create a polgon
poly1 <- data.frame(x = c(0,0,2800,2800,0), y = c(0,500,500,0,0))
coords <- list()
coords[[1]] <- list(poly1)
# Create the survey region
region <- make.region(region.name = "study area", 
                      units = "m",
                      coords = coords)
## ----density
# Create the density surface
density <- make.density(region.obj = region, 
                        x.space = 50, 
                        y.space = 200, 
                        constant = 1)
## ----popdesc2, warning=FALSE, message=FALSE------------------------------------------------
# Create the covariate list
covariate.list <- list()
# The population will be 70% males and 30% females
covariate.list$sex <- list(data.frame(level = c("female", "male"), prob = c(0.3,0.7)))
# Create the population description, with a population size N = 200
pop.desc.cov <- make.population.description(region.obj = region, 
                                            density.obj = density, 
                                            covariates = covariate.list, 
                                            N = 180)
## ----detect2
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
## ----transects
transects <- generate.transects(design, region = region)
plot(region, plot.units = "m")
plot(transects, col = 4, lwd = 2)
## ----analyses2-----------------------------------------------------------------------------
ddf.analyses.cov <- make.ddf.analysis.list(dsmodel = list(~mcds(key = "hn", formula = ~sex)), 
                                           method = "ds", truncation = 80)
## ----set.seed, echo = FALSE, eval = TRUE---------------------------------------------------
set.seed(474)
## ----checksim2
sim.cov <- make.simulation(reps = 99, 
                       region.obj = region,
                       design.obj = design,
                       population.description.obj = pop.desc.cov,
                       detectability.obj = detect.cov,
                       ddf.analyses.list = ddf.analyses)
# Produce simulation setup plots
check.sim.setup(sim.cov)
#  keep one realisation, munge it for formatting, write to file
mysurvey <- create.survey.results(sim.cov)
mydata <- get.distance.data(mysurvey)
library(Distance)
mydata$Region.Label <- "OneStratum"
mydata$Area <- region@area
mydata$Effort <- unname(region@box[4])
names(mydata)[2] <- "Sample.Label"
mydata.nodoubles <- mydata[mydata$object %% 1 == 0, ]
notrun <- ds(mydata.nodoubles)
plot(notrun, main="No truncation")
trun60 <- ds(mydata.nodoubles, truncation = 60)
plot(trun60, main="Trunc 60")
covtrunc <- ds(mydata.nodoubles, truncation=60, formula=~sex)
plot(covtrunc, main="Covar truncation")
covtrunc.hr <- ds(mydata.nodoubles, truncation=60, key="hr")
covtruncsex.hr <- ds(mydata.nodoubles, truncation=60, key="hr", formula=~sex)
AIC(trun60, covtrunc, covtrunc.hr, covtruncsex.hr)
un <- ds(mydata.nodoubles, key="uni", adj="cos", truncation=60)
summarize_ds_models(un, trun60, covtrunc, covtrunc.hr, covtruncsex.hr, output="plain")
colorder <- c("Region.Label", "Area", "Sample.Label", "Effort", "distance", "sex")
mydata.nodoubles$distance <- round(mydata.nodoubles$distance)
mydata.nodoubles$Area <- mydata.nodoubles$Area/1e6  # m^2 to km^2
mydata.nodoubles$Effort <- mydata.nodoubles$Effort/1000  # m to km
dataset <- mydata.nodoubles[order(mydata.nodoubles$Sample.Label) , colorder]
write.csv(dataset, file="mystery.csv", row.names = FALSE, quote=FALSE)
