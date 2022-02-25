#### PRACTICAL 999 ####
#######################

#Distance package
library(Distance)
library(knitr)
conversion.factor <- convert_units("meter", "kilometer", "square kilometer")

# Load data
data <- read.csv("mystery.csv")

# Check that it has been imported correctly
head(data)
 
# How many observations (note: detections on all lines)
nrow(data)


# Histogram with 8 bins
hist(data$distance, xlab="Distance (m)",
     main="Perpendicular distances duck nests")


# Truncate at 60metres
data.uf.t60m <- ds(data, key="unif", adjustment="cos",truncation=60,
              convert.units=conversion.factor)

data.hr.t60m <- ds(data, key="hr", adjustment=NULL,truncation=60, 
              convert.units=conversion.factor)

data.hn.t60m <- ds(data=data, key="hn", adjustment=NULL,truncation=60, 
                 convert.units=conversion.factor)


knitr::kable(summarize_ds_models(data.hn.t60m, data.uf.t60m, data.hr.t60m, output="plain"), 
             caption="Model results for data set.", digits=3)

par(mfrow=c(2,2))
plot(data.uf.t60m, nc=8, main="Uniform, cosine adjustement")
plot(data.hr.t60m, nc=8, main="Hazard-rate, no adjustement")
plot(data.hn.t60m, nc=8, main="Half-normal, no adjustement")


####Best selected model is with keyfunction uniform and cosine adjustment. But could be also Half normal with no adj. 
gof_ds(data.uf.t60m)

## ---------------------------------------------------------------------------

###Plot for covariates
boxplot(data$distance~data$sex, xlab="Sex", ylab="Distance (m)")
boxplot(data$Sample.Label~data$sex, xlab="Sex", ylab="Distance (m)")

# Adjusting the raw data
# Convert sex to a factor
data$sex <- factor(data$sex)
# Set the reference level 
data$sex <- relevel(data$sex, ref="male")

# Fit model selected by Marques et al (2007)
data.hn1<- ds(data, transect="line", key="hn",convert.units = conversion.factor,
              truncation=60)

data.hn.sex<- ds(data, transect="line", key="hn", formula=~factor(sex), convert.units = conversion.factor,
                      truncation=60)

# Compare models (using pretty printing)
knitr::kable(summarize_ds_models(data.hn1, data.hn.sex),
           digits=3) 


summary(data.hn.sex)

###best model selected is data.hn.sex

###performing bootstrap for a better estimation of the variance

est.boot <- bootdht(model=data.hn.sex, flatfile=data,
                    summary_fun=bootdht_Nhat_summarize, 
                    convert.units=conversion.factor, nboot=50)
summary(est.boot)
