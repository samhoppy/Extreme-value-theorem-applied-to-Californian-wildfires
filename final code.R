#
#Code explanation 
#
library(ismev)
library(extRemes)
library(evd)

x <- seq(-3, 5, by = 0.01)

#simulate distributions
pdf_gumbel <- dgev(x, loc = 0, scale = 1,shape=0)
pdf_weibull <- dgev(x ,loc = 0 ,shape = -0.4, scale = 1)
pdf_frechet <- dgev(x, loc = 0 ,shape = 0.4, scale = 1)

df <- data.frame(x = x, 
                 Gumbel = pdf_gumbel, 
                 Weibull = pdf_weibull, 
                 Frechet = pdf_frechet)
#install.packages("reshape2")
library(reshape2)
df_melted <- melt(df, id.vars = "x", variable.name = "Distribution", value.name = "Density")

ggplot(df_melted, aes(x = x, y = Density, color = Distribution)) +
  annotate("text", x = 4.5, y = 0.3, label = expression(mu == 0)) +
  annotate("text", x = 4.5, y = 0.275, label = expression(sigma == 1)) +
  geom_line() +
  labs(title = "Probability Density Functions Comparison",
       x = "x",
       y = "Density") +
  theme_bw()
#GEV dists comparisons
sim<-rgev(1000,loc = 0,scale = 1,shape = 0.4)
sim2<-rnorm(10000,0,1)
sim_fit<-fevd(sim)
plot(sim_fit)
threshrange.plot(sim2,c(0,3),
                 type = "GP",nint = 100)
sim_fit2<-fevd(sim2,type = "GP",threshold = 2)
plot(sim_fit2)

#
#Data preprep
#
mywd <- paste0(Sys.getenv('USERPROFILE'), "\\OneDrive - University of Plymouth\\Year 3\\Extreme Value Theory\\Code")      



setwd(mywd)

getwd()
#install.packages("lubridate")
library(lubridate)
citation("lubridate")
### Read in the data
california.wildfires<- read.csv("California_Fire_Perimeters_(all) (1).csv", 
                                header = TRUE, na.strings = "NA")


# Load the dplyr library
#install.packages("dplyr")
library(dplyr)
citation("dplyr")

names(california.wildfires)


# Convert 'ALARM_DATE' to a date type
california.wildfires$ALARM_DATE <- as.Date(california.wildfires$ALARM_DATE,
                                           format = "%Y/%m/%d")

# Calculate the median date for each year
median_dates <- california.wildfires %>%
  group_by(YEAR_) %>%
  summarize(median_date = median(ALARM_DATE, na.rm = TRUE))
median_dates
# Merge the median dates with the original data set
california.wildfires <- california.wildfires %>%
  left_join(median_dates, by = c("YEAR_" = "YEAR_")) %>%
  mutate(ALARM_DATE = ifelse(is.na(ALARM_DATE), median_date, ALARM_DATE))%>% 
  select(-median_date)#gets red of the median variable that was introduced 
#with join
?is.na
#back to date 
california.wildfires$ALARM_DATE <- as.Date(california.wildfires$ALARM_DATE,
                                           format = "%Y/%m/%d")
#there was a date that was entered as 0219 in the original data set that mi 
#assuming is meant to be 2019.
min(california.wildfires$ALARM_DATE)#

#
# Replace "0219" with "2019" in the ALARM_DATE column
california.wildfires$ALARM_DATE <- as.character(california.wildfires$ALARM_DATE)
california.wildfires$ALARM_DATE <- sub("219-05-29", "2019-05-29",
                                       california.wildfires$ALARM_DATE)
?sub

# Convert the ALARM_DATE column to the Date data type
california.wildfires$ALARM_DATE <- as.Date(california.wildfires$ALARM_DATE)

#
california.wildfires <- california.wildfires %>%
  arrange(ALARM_DATE)


california.wildfires.2<-california.wildfires[c("YEAR_","ALARM_DATE","GIS_ACRES")]

cal.f.2000<- california.wildfires.2%>%
  filter(YEAR_<2001)


#Check reliability of early data
# Load the extRemes package
library(extRemes)
c1<-california.wildfires.2%>%
  filter(YEAR_<1973)
c2<-california.wildfires.2%>%
  filter(1972<YEAR_,YEAR_<1995)

threshrange.plot(c1$GIS_ACRES,c(10000,70000),
                 type = "GP",nint = 200)

threshrange.plot(c2$GIS_ACRES,c(10000,50000),
                 type = "GP",nint = 200)

fevd(c1$GIS_ACRES,c1,type = "GP", threshold = 40000)
fevd(c2$GIS_ACRES,c2,type = "GP", threshold = 40000)

# Fit generalized Pareto distributions to the datasets
fit1 <- fevd(c1$GIS_ACRES, data = c1, type = "GP", threshold = 40000)
plot(fit1)
fit2 <- fevd(c2$GIS_ACRES, data = c2, type = "GP", threshold = 40000)
plot(fit2)


sumf.11<- summary(fit1)
se.11<- sumf.11$se.theta
CI.f11<- confindence(fit1$results$par,se.11)
sumf.12<- summary(fit2)
se.12<- sumf.12$se.theta
CI.f12<- confindence(fit2$results$par,se.12)


combined_data.comp <- rbind(CI.f11, CI.f12)

# Create xtable object
xtable(combined_data.comp)

#
#block maxima
#

library(extRemes)
library(ismev)
library(evd)
library(lubridate)
library(ggplot2)


###-2000
##Year

#create a new data set that only includes the largest wild fires per year up to 
#2000
bm<-blockmaxxer(california.wildfires.2,which = 3,
                blocks = california.wildfires.2$YEAR_)%>%
  filter(YEAR_<2001)
#check there are no missing values from the new data 
t<-is.na(california.wildfires.2$GIS_ACRES)
if (any(t)) {
  cat("There are TRUE values in the list.")
} else {
  cat("There are no TRUE values in the list.")
}
length(california.wildfires.2$YEAR_)
length(california.wildfires.2$GIS_ACRES)
?blockmaxxer
#create a new data set that is only up to the year 2000
cal.f.2000<- california.wildfires.2%>%
  filter(YEAR_<2001)
# plot of the wildfire data being blocked into years, with the max value of each 
#block 
unique_years <- unique(format(cal.f.2000$ALARM_DATE, "%Y"))
base<- ggplot(cal.f.2000,aes(x= ALARM_DATE,y= GIS_ACRES))+
  geom_line()
full<- base+
  annotate("point",x= bm$ALARM_DATE,y=bm$GIS_ACRES,color = "red", size = 3)+
  geom_vline(xintercept = as.Date(unique_years, format = "%Y"), color = "blue",
             linetype = "dashed")+
  labs(title = "Wildfire data blocked into years")

full
#test if there are many non extreme values 
un_exreme<- bm%>%
  filter(bm$GIS_ACRES<5000)
un_exreme
#fit the GEV model 
f<-fevd(bm$GIS_ACRES,bm)
#Confidence Interval=MLE±z×Standard Error
#Here, "MLE" is the maximum likelihood estimate of the parameter, "Standard Error"
#is the standard error of the estimate, and z is the critical value from the 
#standard normal distribution corresponding to the desired confidence level.
#For a 95% confidence interval, z is typically around 1.96.
confindence <- function(x,se){
  lower<-rep(NA,length(x))
  upper<-rep(NA,length(x))
  for(i in 1:length(x)){
    lower[i]<- x[i]-(1.96*se[i])
    upper[i]<- x[i]+ (1.96*se[i])
  }
  df<-data.frame(Lower=lower,Estimate=x,Upper=upper)
  print(df)
}
sumf<- summary(f)
se.f<- sumf$se.theta
CI.f<- confindence(f$results$par,se.f)
f
library(xtable)
xtable(CI.f)


mean(bm$GIS_ACRES)
median(bm$GIS_ACRES)
quantile(bm$GIS_ACRES,c(0.25,0.5,0.75))
#plot the diagnostic graphs 
plot(f)
#plots show that the gev model is actualy quite a good fit with a blocksize of a year 
#find the 2,20 and 100 year rturn levels with a confindence interval 
block.return.levelyear.2000<-return.level(f,return.period = c(2,50),do.ci=T)
block.return.levelyear.2000.df<-data.frame(
  lower_CI=block.return.levelyear.2000[,1],
  Estimate =block.return.levelyear.2000[,2],
  upper_CI=block.return.levelyear.2000[,3])

xtable(block.return.levelyear.2000.df)

ggplot(cal.f.2000, aes(x = ALARM_DATE, y = GIS_ACRES)) +
  geom_point()+
  geom_hline(yintercept = block.return.levelyear.2000[2,2],color="red")+
  geom_ribbon(aes(ymin=block.return.levelyear.2000[2,1],
                  ymax=block.return.levelyear.2000[2,3]),
              fill = "pink", alpha = 0.5)+
  geom_hline(yintercept =block.return.levelyear.2000[1,2],color="green")+
  geom_ribbon(aes(ymin=block.return.levelyear.2000[1,1],
                  ymax=block.return.levelyear.2000[1,3]),
              fill = "lightgreen", alpha = 0.5)+
  labs(title = "The 2 and 50-year Return- levels and their CI's aginst the data  ",
       x="Alarm Date",y="Wildfire Perimeter(Acres)")+
  theme_bw()






#
#plot pdf of GEV dist for year blocks as its best fit
#

hist(bm$GIS_ACRES, freq = F, main = "Histogram with PDF of Fitted GEV Distribution",
     ylim = c(0,max(dgev(bm$GIS_ACRES, loc = f$results$par[1], 
                         scale = f$results$par[2], 
                         shape = f$results$par[3]))),
     xlab = "Wildfire perimeter size (acres)")
curve(dgev(x, loc = f$results$par[1], 
           scale = f$results$par[2], 
           shape = f$results$par[3]), add = TRUE, col = "blue",)



#
#Gumbel
#


f.G<-fevd(bm$GIS_ACRES,bm,type = "Gumbel")
sumf.g<-summary(f.G)
#negative log likelihood

neg_log_likelihood <- sumf.g$nllh

# Extract the number of parameters in the model
num_parameters <- length(f.G$results$par)


#
#Perform a deviance test to test the nested gumbel distribtution agaisnt the GEV
#
#set geres of freedom = difference in the number of parameters
#
df <- 1
# Confidence level
confidence_level <- 0.95
# Critical value
critical_value <- qchisq(confidence_level, df)

#
#the deviance test
#
Deviance.t<- function(L,S){
  -2*(L-S)
}
devicance.test.yearly.bloking<- Deviance.t(sumf$nllh,sumf.g$nllh)
?logLik
f
f.G
print(c(devicance.test.yearly.bloking,critical_value))

# parmeter estimates fo gumbel dist with yealy blocking 
sumf.g
se.f.g<- sumf.g$se.theta
CI.f.g<- confindence(f.G$results$par,se.f.g)
xtable(CI.f.g)





#
#2 years
#
#

#blovk the data into two years
even.years<- seq(1950,2000,by = 2)
blocked.by.2.years<-rep(NA,length(cal.f.2000$YEAR_))
for(i in 1:length(cal.f.2000$YEAR_)){
  ifelse(cal.f.2000$YEAR_[i] %in% even.years,
         blocked.by.2.years[i]<-cal.f.2000$YEAR_[i],
         blocked.by.2.years[i]<-cal.f.2000$YEAR_[i]-1)
}
blocked.by.2.years


df.2.y.to.be.blocked<- data.frame(block=blocked.by.2.years,
                                  acres=cal.f.2000$GIS_ACRES,
                                  alarm_date=cal.f.2000$ALARM_DATE)
blocked.by.2.years<- blockmaxxer(df.2.y.to.be.blocked,which = 2,
                                 blocks = df.2.y.to.be.blocked$block)
length(blocked.by.2.years$acres)


#plot od bi yealry blocked data

#
#Data split into blocks of 2 years
#
even.years <- as.character(even.years)

# Create date strings
date_strings <- paste(even.years, "0101", sep = "")

# Convert to Date
dates <- as.Date(date_strings, format = "%Y%m%d")

ggplot(cal.f.2000,aes(x= ALARM_DATE,y= GIS_ACRES))+
  geom_line()+
  annotate("point",x= blocked.by.2.years$alarm_date,
           y=blocked.by.2.years$acres,color = "red", size = 3)+
  geom_vline(xintercept = as.Date(dates), color = "blue",
             linetype = "dashed")+
  labs(title = "Wildfire data blocked by every 2 years",
       x="Alarm Date",y="Wildfire Perimeter(Acres)")+
  theme_bw()


#Model fitting
f.2.years<-fevd(blocked.by.2.years$acres,blocked.by.2.years)
f.2.years
plot(f.2.years)###############################################################################
r.2.b<-return.level(f.2.years,return.period = c(2,25,50,72),do.ci=T)
sumf.2<-summary(f.2.years)
se.f.2.years<- sumf.2$se.theta
CI.f.2<- confindence(f.2.years$results$par,se.f.2.years)
CI.f.2
xtable(CI.f.2)
#histogram 
hist(blocked.by.2.years$acres, freq = F, main = "Histogram with PDF of Fitted GEV Distribution",
     xlab = "Wildfire perimeter size (acres)")
curve(dgev(x, loc = f.2.years.gumbel$results$par[1], 
           scale = f.2.years.gumbel$results$par[2],
           shape = 0), 
      add = TRUE, col = "blue",)

#gumbel 
f.2.years.gumbel<-fevd(blocked.by.2.years$acres,blocked.by.2.years,type = "Gumbel")
sumf.2.g<- summary(f.2.years.gumbel)
se.f.2.g<- sumf.2.g$se.theta
CI.f.2.g<- confindence(f.2.years.gumbel$results$par,se.f.2.g)
xtable(CI.f.2.g)
plot(f.2.years.gumbel)
return.level(f.2.years.gumbel,return.period = c(2,5,50),do.ci=T)




#Deviance tests
sumf.2$nllh
devicance.test.bi.yearly.bloking<- Deviance.t(sumf.2$nllh,sumf.2.g$nllh)
print(c(devicance.test.bi.yearly.bloking,critical_value))

#
#
#
#data per-process again 

mywd <- paste0(Sys.getenv('USERPROFILE'), "\\OneDrive - University of Plymouth\\Year 3\\Extreme Value Theory\\Code")      



setwd(mywd)

getwd()
#install.packages("lubridate")
library(lubridate)
citation("lubridate")
### Read in the data
california.wildfires<- read.csv("California_Fire_Perimeters_(all) (1).csv", 
                                header = TRUE, na.strings = "NA")


# Load the dplyr library
#install.packages("dplyr")
library(dplyr)
citation("dplyr")

names(california.wildfires)


# Convert 'ALARM_DATE' to a date type
california.wildfires$ALARM_DATE <- as.Date(california.wildfires$ALARM_DATE,
                                           format = "%Y/%m/%d")

# Calculate the median date for each year
median_dates <- california.wildfires %>%
  group_by(YEAR_) %>%
  summarize(median_date = median(ALARM_DATE, na.rm = TRUE))
median_dates
# Merge the median dates with the original data set
california.wildfires <- california.wildfires %>%
  left_join(median_dates, by = c("YEAR_" = "YEAR_")) %>%
  mutate(ALARM_DATE = ifelse(is.na(ALARM_DATE), median_date, ALARM_DATE))%>% 
  select(-median_date)#gets red of the median variable that was introduced 
#with join
?is.na
#back to date 
california.wildfires$ALARM_DATE <- as.Date(california.wildfires$ALARM_DATE,
                                           format = "%Y/%m/%d")
#there was a date that was entered as 0219 in the original data set that mi 
#assuming is meant to be 2019.
min(california.wildfires$ALARM_DATE)#

#
# Replace "0219" with "2019" in the ALARM_DATE column
california.wildfires$ALARM_DATE <- as.character(california.wildfires$ALARM_DATE)
california.wildfires$ALARM_DATE <- sub("219-05-29", "2019-05-29",
                                       california.wildfires$ALARM_DATE)
?sub

# Convert the ALARM_DATE column to the Date data type
california.wildfires$ALARM_DATE <- as.Date(california.wildfires$ALARM_DATE)

#
california.wildfires <- california.wildfires %>%
  arrange(ALARM_DATE)


california.wildfires.2<-california.wildfires[c("YEAR_","ALARM_DATE","GIS_ACRES")]

cal.f.2000<- california.wildfires.2%>%
  filter(YEAR_<2001)
#instal and load packages 
#install.packages("lubridate")
library(extRemes)
library(ismev)
library(evd)
library(lubridate)
library(ggplot2)


#
#valie over threshold approach 
#
#### Value over thresholds 
#2000
?threshrange.plot()

#find the 0.95 quantile to use as an estimate of the threshold 
quantile(cal.f.2000$GIS_ACRES,0.99)
#this did not help 
#Create parameter estimate plots for a range of given thresholds to aid in 
#threshold seletion 
threshrange.plot(cal.f.2000$GIS_ACRES,c(10000,60000),
                 type = "GP",nint = 200)

# Plot a graph illustrating the chosen threshold 
ggplot(cal.f.2000, aes(x = ALARM_DATE, y = GIS_ACRES)) +
  geom_point() +
  geom_hline(yintercept = 40000, color = "red", linetype = "solid")+
  labs(x = "Alarm Date", y = "GIS_ACRES", 
       title = "Visualisation of values over the threshold")+
  theme_bw()

#see if parameters stay relativly unchanged 
threshrange.plot(cal.f.2000$GIS_ACRES,c(39800,40200),
                 type = "GP",nint = 80)

#
#plots the diagnostic graphs for a range of values to help in threshold
#selection 
#
for(i in 1:13){
  plot(fevd(cal.f.2000$GIS_ACRES,cal.f.2000,
            threshold = 15000+(i*5000),type = "GP",))
}



#
#check to see if there is a big change in the shape and scale parameters around 
#possible threshold choices 
#

fevd(cal.f.2000$GIS_ACRES,cal.f.2000,
     threshold = 39990,type = "GP",)$results$par
fevd(cal.f.2000$GIS_ACRES,cal.f.2000,
     threshold = 39995,type = "GP",)$results$par
fevd(cal.f.2000$GIS_ACRES,cal.f.2000,
     threshold = 40000,type = "GP",)$results$par
fevd(cal.f.2000$GIS_ACRES,cal.f.2000,
     threshold = 40005,type = "GP",)$results$par






length(seq(min(cal.f.2000$YEAR_),max(cal.f.2000$YEAR_)))

confindence <- function(x,se){
  lower<-rep(NA,length(x))
  upper<-rep(NA,length(x))
  for(i in 1:length(x)){
    lower[i]<- x[i]-(1.96*se[i])
    upper[i]<- x[i]+ (1.96*se[i])
  }
  df<-data.frame(Lower=lower,Estimate=x,Upper=upper)
  print(df)
}
#Fit the general Pareto distribution 
f.gp.2000<-fevd(cal.f.2000$GIS_ACRES,cal.f.2000,
                threshold = 40000,type = "GP"
)
?fevd
f.gp.2000$results$par[2]
period_span <- 51 
time_units <- "years"
f.gp.2000
sumf.gp.2000<- summary(f.gp.2000)
se.f.gp.2000<- sumf.gp.2000$se.theta
CI.f.gp.2000<- confindence(f.gp.2000$results$par,se.f.gp.2000)
# Plot the diagnostics 
plot(f.gp.2000)
plot(f)
plot(fevd(cal.f.2000$GIS_ACRES,cal.f.2000,
          threshold = 50000,type = "GP",))
#find the return levels 

N<-51
u<-40000
n_y<-sum(cal.f.2000$GIS_ACRES > u)/N

return.level.threshold.2000<-return.level(f.gp.2000,
                                          return.period = c(2*n_y,
                                                            50*n_y,
                                                            72*n_y,
                                                            21*n_y),do.ci=T)
return.level(f.gp.2000,
             return.period =c(2,50,72,21))
#
return.level.threshold.2000[2,2]
?return.level

VOT.return.levelyear.2000.df<-data.frame(
  Period = c("2-year-return-level","50-year-return-level"),
  lower_CI=c(return.level.threshold.2000[1,1],
             return.level.threshold.2000[2,1]),
  Estimate =c(return.level.threshold.2000[1,2],
              return.level.threshold.2000[2,2]),
  upper_CI=c(return.level.threshold.2000[1,3],
             return.level.threshold.2000[2,3]))

max(cal.f.2000$YEAR_)

ggplot(cal.f.2000,aes(x= ALARM_DATE,y=GIS_ACRES))+
  geom_point()+
  geom_hline(yintercept =return.level.threshold.2000[2,2] ,color="red")+
  geom_ribbon(aes(ymin=return.level.threshold.2000[2,1],
                  ymax=return.level.threshold.2000[2,3]), 
              fill = "pink", alpha = 0.5)+
  geom_hline(yintercept =return.level.threshold.2000[1,2] ,color="green")+
  geom_ribbon(aes(ymin=return.level.threshold.2000[1,1],
                  ymax=return.level.threshold.2000[1,3]), 
              fill = "lightgreen", alpha = 0.5)+
  labs(x="Alarm date", y="GIS acres", 
       title = "50 and 2 year return levels with CIs plotted againts the wildfire data up to 2001")+
  theme_bw()
#
#
#PDF of distribution 
#
#find data over the threshold 
thing <- cal.f.2000%>%
  filter(GIS_ACRES>20000)
#use this to create a list of 1000 cvalue from the min ro maxe of the values 
#over threshold 


par(mfrow = c(1, 1))
hist(thing$GIS_ACRES, freq = FALSE, main = "Histogram with Fitted GPD Distribution"
     ,breaks = 20
)
curve(dgpd(x, scale = f.gp.2000$results$par[1], 
           shape = f.gp.2000$results$par[2]), add = TRUE, col = "blue")

?return.level
#
#



####2022




#choose a thrershold 
threshrange.plot(california.wildfires.2$GIS_ACRES,c(10000,120000),
                 type = "GP",nint = 200)


ggplot(california.wildfires.2, aes(x = ALARM_DATE, y = GIS_ACRES)) +
  geom_point() +
  geom_hline(yintercept = 47500, color = "red", linetype = "solid")+
  labs(x = "Alarm Date", y = "GIS_ACRES", 
       title = "Visualisation of values over the threshold")



f.gp.2022<-fevd(california.wildfires.2$GIS_ACRES,california.wildfires.2,
                threshold = 40000,type = "GP")
f.gp.2022

fevd(california.wildfires.2$GIS_ACRES,california.wildfires.2,
     threshold = 20010,type = "GP")$results$par
fevd(california.wildfires.2$GIS_ACRES,california.wildfires.2,
     threshold = 20005,type = "GP")$results$par
fevd(california.wildfires.2$GIS_ACRES,california.wildfires.2,
     threshold = 20000,type = "GP")$results$par
fevd(california.wildfires.2$GIS_ACRES,california.wildfires.2,
     threshold = 19995,type = "GP")$results$par
fevd(california.wildfires.2$GIS_ACRES,california.wildfires.2,
     threshold = 19910,type = "GP")$results$par

plot(f.gp.2022)
#still not an awful fir an gives all 34 diagnostic plots 
N<-72
u<-40000
n_y_2<-sum(california.wildfires.2$GIS_ACRES > u)/N


r.level.full<-return.level(f.gp.2022,return.period = c(2*n_y_2,72*n_y_2
                                                       ,100*n_y_2,
                                                       21*n_y_2),do.ci=T)

#
#Value over threshold approach for only 2001-2022
#

#
#using own threshold
#

cal.f.2000.plus<-california.wildfires.2%>%
  filter(YEAR_>2000)

threshrange.plot(cal.f.2000.plus$GIS_ACRES,c(60000,110000),
                 type = "GP",nint = 200)


#
#plots the diagnostic graphs for a range of values to help in threshold
#selection 
#
for (i in 1:15) {
  tryCatch({
    # Attempt to execute the code within this block
    plot(fevd(cal.f.2000.plus$GIS_ACRES, cal.f.2000.plus, threshold = 50000 + (i * 5000), type = "GP"))
  }, error = function(e) {
    # If an error occurs, execute the following block
    cat("Error occurred in iteration", i, ":", conditionMessage(e), "\n")
    # Print an error message with the iteration number and the error message
  })
}
#

fevd(cal.f.2000.plus$GIS_ACRES,cal.f.2000.plus,
     threshold = 80900,type = "GP")$results$par
fevd(cal.f.2000.plus$GIS_ACRES,cal.f.2000.plus,
     threshold = 80950,type = "GP")$results$par
fevd(cal.f.2000.plus$GIS_ACRES,cal.f.2000.plus,
     threshold = 81000,type = "GP")$results$par
fevd(cal.f.2000.plus$GIS_ACRES,cal.f.2000.plus,
     threshold = 81050,type = "GP")$results$par

#
#these show that 70000 would be a suitable threshold choice for this part of 
#the data set  
#

fit.2001.2022.2<-fevd(cal.f.2000.plus$GIS_ACRES,cal.f.2000.plus,
                      threshold = 90000,type = "GP")

plot(fit.2001.2022.2)

N<-22
u<-90000
n_y_3<-sum(cal.f.2000.plus$GIS_ACRES > u)/N

return.level.last.part<-return.level(fit.2001.2022.2,return.period = c(2*n_y_3,
                                                                       21*n_y_3,
                                                                       77*n_y_3)
                                     ,do.ci=T)



#
#Plot the PDF 
#

#1950-2001



par(mfrow = c(2, 2))
hist(thing$GIS_ACRES, freq = FALSE, main = "1950 - 2001"
     ,breaks =20,xlab = "Wildfire size(Acres)")
curve(dgpd(x, scale = f.gp.2000$results$par[1], 
           shape = f.gp.2000$results$par[2]), add = TRUE, col = "blue")
hist(t.3$GIS_ACRES, freq = FALSE, main = "1950-2022"
     ,breaks = 20,xlab = "Wildfire size(Acres)",ylim = c(0,0.000025)
)
curve(dgpd(x, scale = f.gp.2022$results$par[1], 
           shape = f.gp.2022$results$par[2]), add = TRUE, col = "blue")
hist(thing.2$GIS_ACRES, freq = FALSE, main = "2001-2022"
     ,breaks = 20,xlab = "Wildfire size(Acres)"
)
curve(dgpd(x, scale = fit.2001.2022.2$results$par[1], 
           shape = fit.2001.2022.2$results$par[2]), add = TRUE, col = "blue")


# Extract parameters for the 1950-2022 data
params_1950_2022 <- c(scale = f.gp.2022$results$par[1], shape = f.gp.2022$results$par[2])

# Extract parameters for the 2001-2022 data
params_2001_2022 <- c(scale = fit.2001.2022.2$results$par[1], shape = fit.2001.2022.2$results$par[2])

# Extract parameters for the 1950-2001 data
params_1950_2001 <- c(scale = f.gp.2000$results$par[1], shape = f.gp.2000$results$par[2])

# Create a data frame to store the parameters
params_table <- data.frame(
  Time_Period = c("1950-2022", "2001-2022","1950-2001"),
  Scale_Parameter = c(params_1950_2022[1], params_2001_2022[1],params_1950_2001[1]),
  Shape_Parameter = c(params_1950_2022[2], params_2001_2022[2],params_1950_2001[2])
)




r.level.full
return.level.last.part
return.level.threshold.2000




#
#comparison of 21-year return level
#

ggplot(cal.f.2000.plus, aes(x = ALARM_DATE, y = GIS_ACRES)) +
  geom_point() +
  geom_hline(aes(yintercept = r.level.full[4, 2], color = "Full dataset 50-year RL"), show.legend = TRUE) +
  geom_ribbon(aes(ymin = r.level.full[4, 1], ymax = r.level.full[4, 3]), 
              fill = "pink", alpha = 0.5) +
  geom_hline(aes(yintercept = return.level.last.part[2, 2], color = "Last Part 50-year RL"), show.legend = TRUE) +
  geom_ribbon(aes(ymin = return.level.last.part[2, 1], ymax = return.level.last.part[2, 3]), 
              fill = "lightgreen", alpha = 0.3) +
  geom_hline(aes(yintercept = return.level.threshold.2000[4, 2], color = "Threshold 2000 50-year RL"), show.legend = TRUE) +
  geom_ribbon(aes(ymin = return.level.threshold.2000[4, 1], ymax = return.level.threshold.2000[4, 3]), 
              fill = "lightblue", alpha = 0.5) +
  labs(x = "Alarm date", y = "Wildfire size(acres)", 
       title = "Return level comparison graph") +
  scale_color_manual(name = "Return Level", 
                     values = c("red", "green", "blue"),
                     labels = c("1950-2022", "2001-2022", "1950-2001")) +
  theme_bw() 

max(cal.f.2000.plus$YEAR_)-min(cal.f.2000.plus$YEAR_)

max(cal.f.2000.plus$GIS_ACRES)
max(cal.f.2000$GIS_ACRES)

#
#Dependent stuff


#
#
#
uncluttered_return_level <- function(nexceed, index, nwildfires, threshold, scale, shape, period) {
  proportion <- nexceed / nwildfires
  #formula is n year level so has to by timesed by numbe rof observations per year
  #i.e N-year return level = u+(sigma/xi)(((number of cluster/number of observations)*
  #number of observations per yearI extremal index) -1)
  return(threshold + ((scale / shape) * ((((proportion * index*(nwildfires/50) )*period)^shape) - 1)))
}

head(california.wildfires.2)
# Calculate the total number of events
total_events <- nrow(california.wildfires.2)
# Calculate the total number of unique days
total_unique_days <- n_distinct(california.wildfires.2$ALARM_DATE)
# Calculate the average number of events per day
average_events_per_day <- total_events / total_unique_days


ext.index<-rep(0, 30)
theshold.list<-rep(0, 30)
cluster.cutoff<-rep(0, 30)
mean.100.year.rl<-rep(0, 30)
number.of.clusters<-rep(0, 30)
scale.estimate<-rep(0, 30)
shape.estimate<-rep(0, 30)


cluster_ranges<-seq(3,18,by = 3)
u<-seq(39000,41000,by = 500)
for (i in 1:6) {
  cluster.range<-cluster_ranges[i]*average_events_per_day
  #((i-1)*5)+j
  for (j in 1:5) {
    t<-u[j]
    e<-extremalindex(california.wildfires.2,threshold = t
                     ,which.cols =3,method = "runs"
                     ,run.length = cluster.range)
    ext.index[((i-1)*5)+j]<- as.numeric(e[1])
    declusted.i<-decluster(california.wildfires.2,threshold = t
                           ,which.cols =3,method = "runs"
                           ,r = cluster.range)
    fit.dec.i<- fevd(declusted.i,type = "GP",threshold = t)
    theshold.list[((i-1)*5)+j]<-fit.dec.i$threshold
    scale.estimate[((i-1)*5)+j] <-fit.dec.i$results$par[1]
    shape.estimate[((i-1)*5)+j] <-fit.dec.i$results$par[2]
    number.of.clusters[((i-1)*5)+j]<-as.numeric(e[2])
    ave.y<-as.numeric(e[2])/72
    rl<-return.level(fit.dec.i,ave.y*100)
    mean.100.year.rl[((i-1)*5)+j]<-rl
    cluster.cutoff[((i-1)*5)+j]<-cluster_ranges[i]
  }
}


declustering_data <- data.frame(
  days = cluster.cutoff,
  theshold = theshold.list,
  theta = ext.index,
  mean.100.year.rl = mean.100.year.rl,
  number.of.clusters = number.of.clusters,
  scale = scale.estimate,
  shape = shape.estimate
)


# Calculate coefficient of variation (CV)
cv.rl <- as.numeric(tapply(declustering_data$mean.100.year.rl, declustering_data$days, 
                           function(x) sd(x) / mean(x)))
cv.theta <- as.numeric(tapply(declustering_data$theta, declustering_data$days, 
                              function(x) sd(x) / mean(x)))
cv.nc <- as.numeric(tapply(declustering_data$number.of.clusters, declustering_data$days, 
                           function(x) sd(x) / mean(x)))
cv.scale <- as.numeric(tapply(declustering_data$scale, declustering_data$days, 
                              function(x) sd(x) / mean(x)))
cv.shape <- as.numeric(tapply(declustering_data$shape, declustering_data$days, 
                              function(x) sd(x) / mean(x)))
cv_data<-data.frame(
  #days =cluster_ranges ,
  theta = cv.theta,
  mean.100.year.rl = cv.rl,
  number.of.clusters = cv.nc,
  scale = cv.scale,
  shape = cv.shape
)
average<-rowMeans(cv_data)

cv_data.2<-data.frame(
  days =cluster_ranges ,
  theta = cv.theta,
  mean.100.year.rl = cv.rl,
  number.of.clusters = cv.nc,
  scale = cv.scale,
  shape = cv.shape,
  ave= average
)


# Find the day choice with the smallest CV
optimal_day <- c(3*which.min(cv.rl),3*which.min(cv.theta)
                 ,3*which.min(cv.nc),3*which.min(cv.scale),3*which.min(cv.shape))

# Print the result
print("The choice of days with the smallest average coefficient of variation (CV) is: 6 days")
# Check the lengths of the vectors



declusted.finale<-decluster(california.wildfires.2,threshold = 40000
                            ,which.cols =3,method = "runs"
                            ,r = 6*average_events_per_day)
fit.dec.finale<- fevd(declusted.finale,type = "GP",threshold = 40000)
sumf.gp.dep<- summary(fit.dec.finale)
se.f.gp.dep<- sumf.gp.dep$se.theta
CI.f.gp.dep<- confindence(fit.dec.finale$results$par,se.f.gp.dep)
xtable(CI.f.gp.dep)
VOT.return.levelyear.2000.df<-data.frame(
  Period = c("2-year-return-level","50-year-return-level"),
  lower_CI=c(return.level.threshold.2000[1,1],
             return.level.threshold.2000[2,1]),
  Estimate =c(return.level.threshold.2000[1,2],
              return.level.threshold.2000[2,2]),
  upper_CI=c(return.level.threshold.2000[1,3],
             return.level.threshold.2000[2,3]))
plot(fit.dec.finale)
ave.y<-95/72
return.level(fit.dec.finale,return.period = c(ave.y*72*0.6842105
                                              ,ave.y*100*0.6842105),do.ci=T)
uncluttered_return_level(nexceed =  95,
                         index =  0.6842105,
                         nwildfires =  length(cal.f.2000$YEAR_),
                         threshold =  fit.dec.finale$threshold,
                         scale =  fit.dec.finale$results$par[1],
                         shape = fit.dec.finale$results$par[2],
                         period =  100*(95/72))


#install.packages("xtable")
library(xtable)
latex_table <- xtable(declustering_data)
print(latex_table)
latex_table.2<-xtable(cv_data.2,digits = 4)
print(latex_table.2)

#
#nonstionarity 
#
?fevd
ggplot(california.wildfires.2, aes(x = ALARM_DATE, y = GIS_ACRES)) +
  geom_point()+
  labs(x = "Alarm Date", y = "Wildfire Perimeter(Acres) ", 
       title = "Wildire size Over Time")+
  theme_bw()





california.wildfires.2$numeric_time.2 <- as.numeric(difftime(california.wildfires.2$ALARM_DATE, min(california.wildfires.2$ALARM_DATE), units = "days"))
head(california.wildfires.2)

extremes <- california.wildfires.2%>%
  filter(GIS_ACRES>40000)

ggplot(extremes, aes(x = ALARM_DATE, y = GIS_ACRES)) +
  geom_point()+
  labs(x = "Alarm Date", y = "Wildfire Perimeter(Acres) ", 
       title = "Wildire size Over Time")+
  theme_bw()
?decluster
aa<-decluster(california.wildfires.2,threshold = 40000
              ,which.cols =3,method = "runs"
              ,r = 6*average_events_per_day)
# Use fevd function with numeric_time in scale.fun
non.stationary.f.gp.2022 <- fevd(aa,data = california.wildfires.2,
                                 type = "GP",
                                 threshold = 40000,
                                 scale.fun =  ~ numeric_time.2
                                 +I(numeric_time.2^2)
                                 ,use.phi = T,
                                 threshold.fun =~ exp((numeric_time.2-19000)/3650)
                                 
)



non.stationary.f.gp.2022



# Extract the diagonal elements of the inverted Hessian matrix
se.f.gp.ns <- sqrt(diag(non.stationary.f.gp.2022$results$hessian))
CI.f.gp.ns<- confindence(non.stationary.f.gp.2022$results$par,se.f.gp.ns)
#View(confindence)
xtable(CI.f.gp.ns)
parameter_test <- data.frame( Estimate = estimates, 'Std. Error' = se.f.gp.ns, 'Z-value' = z_values, 'P-value' = p_values)
print(parameter_test)
# Assuming non.stationary.f.gp.2022$threshold is your threshold vector
threshold <- non.stationary.f.gp.2022$threshold

# Create a plot
par(mfrow = c(1, 1))
#plot(threshold, type = "l", xlab = "Index", ylab = "Threshold Value", main = "Threshold Values")
plot(california.wildfires.2$numeric_time.2, threshold, type = "l", 
     xlab = "Numeric Time", ylab = "Threshold", main = "Threshold vs Numeric Time")
# Add grid lines
grid()

# Add a horizontal line for better visualization
abline(h = 40000, col = "red", lty = 2)

non.stationary.f.gp.2022
summary(non.stationary.f.gp.2022)
summary(f.gp.2022)
plot(non.stationary.f.gp.2022)

par(mfrow = c(1, 2))
# Plotting trace plot
plot(non.stationary.f.gp.2022, type = "qq",main = "")

# Plotting qq plot
plot(non.stationary.f.gp.2022, type = "qq2",main = "")
par(mfrow = c(1, 1))
# Plotting return level plot
plot(non.stationary.f.gp.2022, type = "trace")

fit.dec.finale

?decluster

threshold.data <- data.frame(ALARM_DATE = california.wildfires.2$ALARM_DATE,
                             threshold = non.stationary.f.gp.2022$threshold)
ggplot(california.wildfires.2, aes(x = ALARM_DATE, y = GIS_ACRES)) +
  geom_point() +
  geom_line(data = threshold.data,aes(y=threshold),colour="red")+
  labs(x = "Alarm Date", y = "Wildfire Perimeter (Acres)", 
       title = "Wildfire Size Over Time") +
  theme_bw()
?fevd
head(california.wildfires.2)
head(threshold.data)
exceeds <- 0

exccedances<-california.wildfires$GIS_ACRES - threshold.data$threshold
for (i in 1:length(exccedances)) {
  if (exccedances[i] >0) {
    exceeds <- exceeds + 1
  }
}

exceeds


max(exccedances)




r.l.non.stationary.vot.1<-return.level(non.stationary.f.gp.2022,return.period =0.6842105*50*(exceeds/72))
california.wildfires.2$r.l.1<- r.l.non.stationary.vot.1
ggplot(california.wildfires.2, aes(x = ALARM_DATE, y = GIS_ACRES)) +
  geom_point() +
  geom_line(aes(x= ALARM_DATE,y=r.l.1,color="non-stationary"),
            show.legend = T
  )+
  geom_hline(aes(yintercept = return.level.threshold.2000[2,2],color="stationary")
             ,show.legend = T)+
  geom_ribbon(aes(ymin=return.level.threshold.2000[2,1],
                  ymax=return.level.threshold.2000[2,3]),
              fill = "lightgreen", alpha = 0.5)+
  #geom_line(data = threshold.data,aes(y=threshold),colour="blue")+
  labs(x = "Alarm Date", y = "Wildfire Perimeter (Acres)", 
       title = "Return level comparison ") +
  scale_color_manual(name = "50-year Return Level", 
                     values = c("red", "green"),
                     labels = c("non-stationary", "stationary")) +
  theme_bw()
min(california.wildfires.2$ALARM_DATE)
head(r.l.non.stationary.vot.1,1)
max(california.wildfires.2$ALARM_DATE)
tail(r.l.non.stationary.vot.1,1)


cal.plus<-california.wildfires.2%>%
  filter(YEAR_>2001)
extremalindex(cal.plus,threshold = 90000
              ,which.cols =3,method = "runs"
              ,run.length =  6*average_events_per_day)
bb<-decluster(cal.plus,threshold = 90000
              ,which.cols =3,method = "runs"
              ,r = 6*average_events_per_day)
cal.plus$numeric_time.2 <- as.numeric(difftime(cal.plus$ALARM_DATE, min(cal.plus$ALARM_DATE), units = "days"))
non.stationary.f.gp.2012_2022 <- fevd(bb,data = cal.plus,
                                      type = "GP",
                                      threshold = 90000,
                                      scale.fun =  ~ numeric_time.2
                                      #+I(numeric_time.2^2)
                                      #+I(numeric_time.2^3)
                                      #,threshold.fun =~ I(numeric_time.2/10000)
)
non.stationary.f.gp.2012_2022
plot(non.stationary.f.gp.2012_2022)
par(mfrow = c(1, 2))
plot(non.stationary.f.gp.2012_2022,type="qq",main="")
plot(non.stationary.f.gp.2012_2022,type="qq2")


exceeds <- 0

exccedances<-cal.plus$GIS_ACRES - 90000
for (i in 1:length(exccedances)) {
  if (exccedances[i] >0) {
    exceeds <- exceeds + 1
  }
}
r.l.non.stationary.vot.2<-return.level(non.stationary.f.gp.2012_2022,return.period = 0.6818182*20*(exceeds/20))
r.l.non.stationary.vot.3<-return.level(non.stationary.f.gp.2012_2022,return.period = 0.6818182*50*(exceeds/20))
r.l.non.stationary.vot.4<-return.level(non.stationary.f.gp.2012_2022,return.period = 0.6818182*5*(exceeds/20))
cal.plus$r.l.1<- r.l.non.stationary.vot.2
cal.plus$r.l.2<- r.l.non.stationary.vot.3
cal.plus$r.l.3<- r.l.non.stationary.vot.4
ggplot(cal.plus, aes(x = ALARM_DATE, y = GIS_ACRES)) +
  geom_point() +
  geom_line(aes(x = ALARM_DATE, y = r.l.1, colour = "20 yrear return level")) +
  geom_line(aes(x = ALARM_DATE, y = r.l.2, colour = "50 yrear return level")) +
  geom_line(aes(x = ALARM_DATE, y = r.l.3, colour = "5 yrear return level")) +
  labs(title = "Retun levels over time",
       x = "Date",
       y = "Wildfire Perimeter (Acres)",
       colour = "Return level") +
  theme_bw()
non.stationary.f.gp.2012_2022
max(cal.plus$numeric_time.2)
date.2075 <-as.numeric(difftime(seq(min(cal.plus$ALARM_DATE), max(cal.plus$ALARM_DATE)+365*53,by="day") , min(cal.plus$ALARM_DATE), units = "days"))
model.rl<-lm(r.l.2~numeric_time.2,data = cal.plus)
model.rl2<-lm(r.l.1~numeric_time.2,data = cal.plus)
model.rl3<-lm(r.l.3~numeric_time.2,data = cal.plus)
summary(model.rl)
# Create a data frame for prediction
newdata <- data.frame(numeric_time.2 = date.2075)
newdata$date<-as.Date(as.numeric(min(cal.plus$ALARM_DATE))+date.2075)
# Predict using the model
prediction <- predict(model.rl, newdata = newdata)
newdata$predict<-prediction
newdata$pre2<-predict(model.rl2, newdata = newdata)
newdata$pre3<-predict(model.rl3, newdata = newdata)
# View the prediction
print(prediction)

tail(newdata)

CI.f.gp.ns.2<- confindence(non.stationary.f.gp.2012_2022$results$par,
                           sqrt(diag(non.stationary.f.gp.2012_2022$results$hessian)))
xtable(CI.f.gp.ns.2,digits = 4)

ggplot() +
  geom_line(data = newdata ,aes(x = date, y = predict,colour = "50-year")) +
  geom_line(data = newdata ,aes(x = date, y = pre2,colour = "20-year")) +
  geom_line(data = newdata ,aes(x = date, y = pre3,colour = "5-year")) +
  labs(title = "Retun levels over time",
       x = "Date",
       y = "Return level(Acres)",
       colour = "Return level") +
  theme_bw()







