rm(list=ls())

library(unmarked)
library(ggplot2)
library(grid)
library(RColorBrewer)
library(ggpubr)
library(gtable)
library(AICcmodavg)
library(MuMIn)

##################################################
## Detection Histories - 2015 and 2016
##################################################

## BROWN leech
brown15 <- read.csv("data/Occupancy/brown.y.2015.csv")
brown16 <- read.csv("data/Occupancy/brown.y.2016.csv")
tiger15 <- read.csv("data/Occupancy/tiger.y.2015.csv")
tiger16 <- read.csv("data/Occupancy/tiger.y.2016.csv")

brown15.y <- brown15[,2:5]
brown16.y <- brown16[,2:5]
tiger15.y <- tiger15[,2:5]
tiger16.y <- tiger16[,2:5]

# Remove the columns with all NAs
brown16.y2 <- brown16.y[!is.na(brown16.y$y.1) | !is.na(brown16.y$y.2) | !is.na(brown16.y$y.3) | !is.na(brown16.y$y.4), ]
tiger16.y2 <- tiger16.y[!is.na(tiger16.y$y.1) | !is.na(tiger16.y$y.2) | !is.na(tiger16.y$y.3) | !is.na(tiger16.y$y.4), ]

##################################################
## Covariates 
##################################################

# Extract the covariates 
# Site
cov <- read.csv("data/EnvData/RosieLiDARSites_scale_50.csv")
height <- cov[,c('canopy_height_mean')]
moran <-  cov[,c('canopy_height_moran')]
pai <- cov[,c('pai_mean')]

# Survey
date15 <- brown15[,8:11]
date16 <- brown16[,8:11]
effort15 <- brown15[,12:15]
effort16 <- brown16[,12:15]

# for 2015
height15 <- height[-c(170:179)]
moran15 <- moran[-c(170:179)]
pai15 <- pai[-c(170:179)]
date15 <- date15[-c(170:179),]
effort15 <- effort15[-c(170:179),]

# match up for 2016 sites
height16 <- height[-c(1:4,17:20, 33:36, 49:52, 65:72, 81:84, 170:179)] 
moran16 <- moran[-c(1:4,17:20, 33:36, 49:52, 65:72, 81:84, 170:179)] 
pai16 <- pai[-c(1:4,17:20, 33:36, 49:52, 65:72, 81:84, 170:179)]
date16 <- date16[-c(1:4,17:20, 33:36, 49:52, 65:72, 81:84, 170:179),]
effort16 <- effort16[-c(1:4,17:20, 33:36, 49:52, 65:72, 81:84, 170:179),]

##################################################
## Make the unmarked data frame
##################################################

brown15.umf <- unmarkedFrameOccu(y = brown15.y, siteCovs = data.frame(height = height15, moran = moran15), obsCovs = list(date = date15, effort = effort15))
tiger15.umf <- unmarkedFrameOccu(y = tiger15.y, siteCovs = data.frame(height = height15, moran = moran15), obsCovs = list(date = date15, effort = effort15))
brown16.umf <- unmarkedFrameOccu(y = brown16.y2, siteCovs = data.frame(height = height16, moran = moran16), obsCovs = list(date = date16, effort = effort16))
tiger16.umf <- unmarkedFrameOccu(y = tiger16.y2, siteCovs = data.frame(height = height16, moran = moran16), obsCovs = list(date = date16, effort = effort16))

summary(brown15.umf)
summary(tiger15.umf)
summary(brown16.umf)
summary(tiger16.umf)

# Scale site covariates
sc1 <- scale(siteCovs(brown15.umf))
siteCovs(brown15.umf) <- sc1

sc1 <- scale(siteCovs(tiger15.umf))
siteCovs(tiger15.umf) <- sc1

sc2 <- scale(siteCovs(brown16.umf))
siteCovs(brown16.umf) <- sc2

sc2 <- scale(siteCovs(tiger16.umf))
siteCovs(tiger16.umf) <- sc2

# Scale obs covariates
sc3 <- scale(obsCovs(brown15.umf))
obsCovs(brown15.umf) <- sc3

sc3 <- scale(obsCovs(tiger15.umf))
obsCovs(tiger15.umf) <- sc3

sc4 <- scale(obsCovs(brown16.umf))
obsCovs(brown16.umf) <- sc4

sc4 <- scale(obsCovs(tiger16.umf))
obsCovs(tiger16.umf) <- sc4

#################################################################
# Model fitting
######################################################################
## Model1 - NULL MODEL - constant occupancy and detection
######################################################################

(fm1.brown15 <- occu(~1 ~1, brown15.umf))
(fm1.tiger15 <- occu(~1 ~1, tiger15.umf))

(fm1.brown16 <- occu(~1 ~1, brown16.umf))
(fm1.tiger16 <- occu(~1 ~1, tiger16.umf))

# For standard error
# Extract occupancy estimates 
psi.brown15 <- backTransform(fm1.brown15, type="state")
psi.brown16 <- backTransform(fm1.brown16, type="state")

psi.tiger15 <- backTransform(fm1.tiger15, type="state")
psi.tiger16 <- backTransform(fm1.tiger16, type="state")

# Extract detection estimates 
p.brown15 <- backTransform(fm1.brown15, type="det")
p.brown16 <- backTransform(fm1.brown16, type="det")

p.tiger15 <- backTransform(fm1.tiger15, type="det")
p.tiger16 <- backTransform(fm1.tiger16, type="det")

# Bind in a dataframe
df <- as.data.frame(rbind(B15 = psi.brown15@estimate, B16 = psi.brown16@estimate, T15 = psi.tiger15@estimate, T16 = psi.tiger16@estimate))
df <- as.data.frame(cbind(psi = df, SE.psi = c(0.0356, 0.0505, 0.0321, 0.0518)))
colnames(df) <- c("estimate", "SE")
df$year <- c("Dry", "Wet","Dry", "Wet")
df$species <- c("Brown", "Brown", "Tiger", "Tiger")
df$type <- c(rep("psi",4))
df

# Bind all into a dataframe
df2 <- as.data.frame(rbind(B15 = p.brown15@estimate, B16 = p.brown16@estimate, T15 = p.tiger15@estimate, T16 = p.tiger16@estimate))
df2 <- as.data.frame(cbind(p = df2, SE.p = c(0.0251, 0.0338, 0.0244, 0.0348)))
colnames(df2) <- c("estimate", "SE")
df2$year <- c("Dry", "Wet","Dry", "Wet")
df2$species <- c("Brown", "Brown", "Tiger", "Tiger")
df2$type <- c(rep("p",4))
df2


null.plot.psi <- ggplot(df, aes(species, estimate, colour = species)) + 
                geom_errorbar(mapping = aes(x = species, ymin = estimate - SE, ymax = estimate + SE), width = 0.1, size = 1.5)
null.plot.psi <- null.plot.psi + geom_point(size = 3.5, shape = 21, fill = "white") +
                facet_wrap(.~year, nrow =1) + theme_bw()
null.plot.psi <- null.plot.psi + ylab("Occupancy probability") + xlab("") + 
                labs(tag = "A)") + ylim(0,1) +
                scale_colour_manual(values = c("saddlebrown", "forestgreen", "saddlebrown", "forestgreen"))
null.plot.psi <- null.plot.psi + theme(text = element_text(size = 25),
                               strip.background = element_blank(), strip.text = element_blank(),
                               axis.title.y = element_text(margin = margin(r = 20)),
                               legend.title = element_blank(),
                               legend.text = element_text(margin = margin(c(5, 5, 5, 0))),
                               axis.text.x = element_text(margin = margin(c(5, 0, 0, 0))),
                               axis.text.y = element_text(margin = margin(c(0, 5, 0, 0))),
                               plot.title = element_text(margin = margin(c(0, 0, 5, 0))), legend.position = "none")
null.plot.psi
ggsave("figures/Fig.2.SeasonalOccupancy.png", height = 12.8, width = 10.9)

null.plot.p <- ggplot(df2, aes(species, estimate, colour = species)) + 
  geom_errorbar(aes(ymin = estimate - SE, ymax = estimate + SE), width = 0.1, size = 1.5)
null.plot.p <- null.plot.p + geom_point(size = 3.5, shape = 21, fill = "white") +
  facet_wrap(.~year, nrow =1) + theme_bw()
null.plot.p <- null.plot.p + ylim(0,1) + ylab("Detection probability") + xlab("") + 
  labs(tag = "B)")+
  scale_colour_manual(values = c("saddlebrown", "forestgreen", "saddlebrown", "forestgreen"))
null.plot.p <- null.plot.p + theme(text = element_text(size = 25),
                                       strip.background = element_blank(), strip.text = element_blank(),
                                       axis.title.y = element_text(margin = margin(r = 20)),
                                       legend.title = element_blank(),
                                       legend.text = element_text(margin = margin(c(5, 5, 5, 0))),
                                       axis.text.x = element_text(margin = margin(c(5, 0, 0, 0))),
                                       axis.text.y = element_text(margin = margin(c(0, 5, 0, 0))),
                                       plot.title = element_text(margin = margin(c(0, 0, 5, 0))), legend.position = "none")
null.plot.p
ggsave("figures/Fig.2.SeasonalDetection.png", height = 12.8, width = 10.9)

ggarrange(null.plot.psi, null.plot.p, nrow = 2, legend = "none")
ggsave("figures/Fig.2.Combi.png", height = 12.8, width = 10.9)
ggsave("figures/Fig.2.Combi.jpeg", height = 12.8, width = 10.9)
ggsave("figures/PaperFigures/Fig.2.Combi.pdf", height = 12.8, width = 10.9)

##################################################
## Model2 - global model
##################################################

(fm2.brown15 <- occu(~date+effort+height+moran ~height*moran, brown15.umf))
(fm2.tiger15 <- occu(~date+effort+height+moran ~height*moran, tiger15.umf))
(fm2.brown16 <- occu(~date+effort+height+moran ~height*moran, brown16.umf))
(fm2.tiger16 <- occu(~date+effort+height+moran ~height*moran, tiger16.umf))

# (brown15.GOF <- mb.gof.test(fm2.brown15, nsim = 100)) # c-hat = 3.16
# (tiger15.GOF <- mb.gof.test(fm2.tiger15, nsim = 100)) # c-hat = 2.08
# (brown16.GOF <- mb.gof.test(fm2.brown16, nsim = 100)) # c-hat = 2.84
# (tiger16.GOF <- mb.gof.test(fm2.tiger16, nsim = 100)) # c-hat = 2.26

ranef(fm2.brown15)
LRT(fm1.brown15,fm2.brown15)
LRT(fm1.tiger15,fm2.tiger15)
LRT(fm1.brown16,fm2.brown16)
LRT(fm1.tiger16,fm2.tiger16)

##################################################
## 
##################################################

(fm3.brown15 <- occu(~date+effort+height+moran ~1, brown15.umf))
(fm3.tiger15 <- occu(~date+effort+height+moran ~1, tiger15.umf))
(fm3.brown16 <- occu(~date+effort+height+moran ~1, brown16.umf))
(fm3.tiger16 <- occu(~date+effort+height+moran ~1, tiger16.umf))

########################################################################
# Detection model selection - observation process
########################################################################

(fm4.brown15 <- occu(~effort+date+height ~1, brown15.umf))
(fm4.tiger15 <- occu(~effort+date+height ~1, tiger15.umf))
(fm4.brown16 <- occu(~effort+date+height ~1, brown16.umf))
(fm4.tiger16 <- occu(~effort+date+height ~1, tiger16.umf))

(fm5.brown15 <- occu(~effort+date+moran ~1, brown15.umf))
(fm5.tiger15 <- occu(~effort+date+moran ~1, tiger15.umf))
(fm5.brown16 <- occu(~effort+date+moran ~1, brown16.umf))
(fm5.tiger16 <- occu(~effort+date+moran ~1, tiger16.umf))

(fm6.brown15 <- occu(~effort+height+moran ~1, brown15.umf))
(fm6.tiger15 <- occu(~effort+height+moran ~1, tiger15.umf))
(fm6.brown16 <- occu(~effort+height+moran ~1, brown16.umf))
(fm6.tiger16 <- occu(~effort+height+moran ~1, tiger16.umf))

(fm7.brown15 <- occu(~date+height+moran ~1, brown15.umf))
(fm7.tiger15 <- occu(~date+height+moran ~1, tiger15.umf))
(fm7.brown16 <- occu(~date+height+moran ~1, brown16.umf))
(fm7.tiger16 <- occu(~date+height+moran ~1, tiger16.umf))

(fm8.brown15 <- occu(~effort+date ~1, brown15.umf))
(fm8.tiger15 <- occu(~effort+date ~1, tiger15.umf))
(fm8.brown16 <- occu(~effort+date ~1, brown16.umf))
(fm8.tiger16 <- occu(~effort+date ~1, tiger16.umf))

(fm9.brown15 <- occu(~effort+height ~1, brown15.umf))
(fm9.tiger15 <- occu(~effort+height ~1, tiger15.umf))
(fm9.brown16 <- occu(~effort+height ~1, brown16.umf))
(fm9.tiger16 <- occu(~effort+height ~1, tiger16.umf))

(fm10.brown15 <- occu(~effort+moran ~1, brown15.umf))
(fm10.tiger15 <- occu(~effort+moran ~1, tiger15.umf))
(fm10.brown16 <- occu(~effort+moran ~1, brown16.umf))
(fm10.tiger16 <- occu(~effort+moran ~1, tiger16.umf))

(fm11.brown15 <- occu(~date+height ~1, brown15.umf))
(fm11.tiger15 <- occu(~date+height ~1, tiger15.umf))
(fm11.brown16 <- occu(~date+height ~1, brown16.umf))
(fm11.tiger16 <- occu(~date+height ~1, tiger16.umf))

(fm12.brown15 <- occu(~date+moran ~1, brown15.umf))
(fm12.tiger15 <- occu(~date+moran ~1, tiger15.umf))
(fm12.brown16 <- occu(~date+moran ~1, brown16.umf))
(fm12.tiger16 <- occu(~date+moran ~1, tiger16.umf))

(fm13.brown15 <- occu(~height+moran ~1, brown15.umf))
(fm13.tiger15 <- occu(~height+moran ~1, tiger15.umf))
(fm13.brown16 <- occu(~height+moran ~1, brown16.umf))
(fm13.tiger16 <- occu(~height+moran ~1, tiger16.umf))

########################################################################
# Detection fitlist
########################################################################
fms.brown15 <- fitList("m1" = fm1.brown15,
                       "m2" = fm2.brown15,
                       "m3" = fm3.brown15,
                       "m4" = fm4.brown15,
                       "m5" = fm5.brown15,
                       "m6" = fm6.brown15,
                       "m7" = fm7.brown15,
                       "m8" = fm8.brown15,
                       "m9" = fm9.brown15,
                       "m10" = fm10.brown15,
                       "m11" = fm11.brown15,
                       "m12" = fm12.brown15,
                       "m13" = fm13.brown15)

fms.tiger15 <- fitList("m1" = fm1.tiger15,
                       "m2" = fm2.tiger15,
                       "m3" = fm3.tiger15,
                       "m4" = fm4.tiger15,
                       "m5" = fm5.tiger15,
                       "m6" = fm6.tiger15,
                       "m7" = fm7.tiger15,
                       "m8" = fm8.tiger15,
                       "m9" = fm9.tiger15,
                       "m10" = fm10.tiger15,
                       "m11" = fm11.tiger15,
                       "m12" = fm12.tiger15,
                       "m13" = fm13.tiger15)

fms.brown16 <- fitList("m1" = fm1.brown16,
                       "m2" = fm2.brown16,
                       "m3" = fm3.brown16,
                       "m4" = fm4.brown16,
                       "m5" = fm5.brown16,
                       "m6" = fm6.brown16,
                       "m7" = fm7.brown16,
                       "m8" = fm8.brown16,
                       "m9" = fm9.brown16,
                       "m10" = fm10.brown16,
                       "m11" = fm11.brown16,
                       "m12" = fm12.brown16,
                       "m13" = fm13.brown16)

fms.tiger16 <- fitList("m1" = fm1.tiger16,
                       "m2" = fm2.tiger16,
                       "m3" = fm3.tiger16,
                       "m4" = fm4.tiger16,
                       "m5" = fm5.tiger16,
                       "m6" = fm6.tiger16,
                       "m7" = fm7.tiger16,
                       "m8" = fm8.tiger16,
                       "m9" = fm9.tiger16,
                       "m10" = fm10.tiger16,
                       "m11" = fm11.tiger16,
                       "m12" = fm12.tiger16,
                       "m13" = fm13.tiger16)

(ms.brown15 <- modSel(fms.brown15)) #m9  -detection structure
(ms.tiger15 <- modSel(fms.tiger15)) #m11 - detection structure
(ms.brown16 <- modSel(fms.brown16)) #m9 - detection structure
(ms.tiger16 <- modSel(fms.tiger16)) #m6 - detection structure

##################################################
## Occupancy structure
##################################################

(fm14.brown15 <- occu(~effort+height ~height*moran, brown15.umf))
(fm14.tiger15 <- occu(~date+height ~height*moran, tiger15.umf))
(fm14.brown16 <- occu(~effort+height ~height*moran, brown16.umf))
(fm14.tiger16 <- occu(~effort+height+moran ~height*moran, tiger16.umf))

(fm15.brown15 <- occu(~effort+height ~height+moran, brown15.umf))
(fm15.tiger15 <- occu(~date+height ~height+moran, tiger15.umf))
(fm15.brown16 <- occu(~effort+height ~height+moran, brown16.umf))
(fm15.tiger16 <- occu(~effort+height+moran ~height+moran, tiger16.umf))

(fm16.brown15 <- occu(~effort+height ~height, brown15.umf))
(fm16.tiger15 <- occu(~date+height ~height, tiger15.umf))
(fm16.brown16 <- occu(~effort+height ~height, brown16.umf))
(fm16.tiger16 <- occu(~effort+height+moran ~height, tiger16.umf))

(fm17.brown15 <- occu(~effort+height ~moran, brown15.umf))
(fm17.tiger15 <- occu(~date+height ~moran, tiger15.umf))
(fm17.brown16 <- occu(~effort+height ~moran, brown16.umf))
(fm17.tiger16 <- occu(~effort+height+moran ~moran, tiger16.umf))

##################################################
## Model selection, prediction and averaging
##################################################

######################################################
# Model selection - USING STEPWISE - DETECTION FIRST
######################################################

# Put the fitted models in a "fitList"
fms.brown15.all <- fitList("m9"  = fm9.brown15,
                           "m14"  = fm14.brown15,
                           "m15"  = fm15.brown15,
                           "m16"  = fm16.brown15,
                           "m17"  = fm17.brown15)

fms.tiger15.all <- fitList("m11"  = fm11.tiger15,
                           "m14"  = fm14.tiger15,
                           "m15"  = fm15.tiger15,
                           "m16"  = fm16.tiger15,
                           "m17"  = fm17.tiger15)

fms.brown16.all <- fitList("m9"  = fm9.brown16,
                           "m14"  = fm14.brown16,
                           "m15"  = fm15.brown16,
                           "m16"  = fm16.brown16,
                           "m17"  = fm17.brown16)

fms.tiger16.all <- fitList("m6"  = fm6.tiger16,
                           "m14"  = fm14.tiger16,
                           "m15"  = fm15.tiger16,
                           "m16"  = fm16.tiger16,
                           "m17"  = fm17.tiger16)

# Rank them by AIC and choose 95% confidence set
(ms.brown15.all <- modSel(fms.brown15.all)) #m17, m15, m14
(ms.tiger15.all <- modSel(fms.tiger15.all)) #m16, m17, m15
(ms.brown16.all <- modSel(fms.brown16.all)) #m16, m15, m14
(ms.tiger16.all <- modSel(fms.tiger16.all)) #m17, m16, m14

# Using qAICc for small and overdispersed 
QAIC.B15 <- QAIC(fm9.brown15, fm14.brown15, fm15.brown15, fm16.brown15, fm17.brown15, chat = 3.16)
QAIC.B15[order(QAIC.B15$QAIC),]

QAIC.T15 <- QAIC(fm9.tiger15, fm14.tiger15, fm15.tiger15, fm16.tiger15, fm17.tiger15, chat = 2.08)
QAIC.T15[order(QAIC.T15$QAIC),]

QAIC.B16 <- QAIC(fm9.brown16, fm14.brown16, fm15.brown16, fm16.brown16, fm17.brown16, chat = 2.84)
QAIC.B16[order(QAIC.B16$QAIC),]

QAIC.T16 <- QAIC(fm9.tiger16, fm14.tiger16, fm15.tiger16, fm16.tiger16, fm17.tiger16, chat = 2.26)
QAIC.T16[order(QAIC.T16$QAIC),]

# Average B15
aveB15 <- model.avg(list("m17"  = fm17.brown15, "m16"  = fm16.brown15, "m9"  = fm9.brown15), fit = TRUE, 
                    rank = "QAIC", rank.args = list(chat = 3.16))
# For plotting 
fmsAve.B15 <- fitList("m17"  = fm17.brown15, "m16"  = fm16.brown15, "m9"  = fm9.brown15)

# Average T15
aveT15 <- model.avg(list("m16"  = fm16.tiger15, "m17"  = fm17.tiger15, "m15"  = fm15.tiger15), fit = TRUE,
                    rank = "QAIC", rank.args = list(chat = 2.08))
# For plotting 
fmsAve.T15 <- fitList("m16"  = fm16.tiger15, "m17"  = fm17.tiger15, "m15"  = fm15.tiger15)
summary(aveT15)

# Average B16
aveB16 <- model.avg(list("m16"  = fm16.brown16, "m17"  = fm17.brown16, "m9"  = fm9.brown16), fit = TRUE,
                    rank = "QAIC", rank.args = list(chat = 2.84))
# For plotting 
fmsAve.B16 <- fitList("m16"  = fm16.brown16, "m17"  = fm17.brown16, "m9"  = fm9.brown16)
summary(aveB16)

# Average T16
aveT16 <- model.avg(list("m17"  = fm17.tiger16, "m9"  = fm9.tiger16), subset = cumsum(weight) <= .95, fit = TRUE, rank = "QAIC", rank.args = list(chat = 2.26))

# For plotting 
fmsAve.T16 <- fitList("m17"  = fm17.tiger16, "m9"  = fm9.tiger16)
summary(aveT16)

########################################################################
## Brown dry prediction plots
########################################################################

# predict new data using the averaged model - select values based on min and max of the umf
newDataBrownDry1 <- data.frame(height = seq(-2.23265, 2.80800, by = 0.01), moran = 0)
pred.occ.heightB15 <- predict(fmsAve.B15, type = "state", newdata = newDataBrownDry1, appendData = TRUE)

newDataBrownDry2 <- data.frame(height = 0, moran = seq(-2.62339, 2.51316, by = 0.01))
pred.occ.moranB15 <- predict(fmsAve.B15, type = "state", newdata = newDataBrownDry2, appendData = TRUE)

########################################################################

newDataBrownDry3 <- data.frame(height = seq(-2.23265, 2.80800, by = 0.01), date = 0, effort = 0, moran = 0)
pred.det.heightB15 <- predict(fmsAve.B15, type = "det", newdata = newDataBrownDry3, appendData = TRUE)

newDataBrownDry4 <- data.frame(height = 0, date = seq(-1.6578, 1.6194, by = 0.01), effort = 0, moran = 0)
pred.det.dateB15 <- predict(fmsAve.B15, type = "det", newdata = newDataBrownDry4, appendData = TRUE)

newDataBrownDry5 <- data.frame(height = 0, date = 0, effort = seq(-1.80297, 2.56551, by = 0.01), moran = 0)
pred.det.effB15 <- predict(fmsAve.B15, type = "det", newdata = newDataBrownDry5, appendData = TRUE)

newDataBrownDry6 <- data.frame(height = 0, date = 0, effort = 0, moran = seq(-2.62339, 2.51316, by = 0.01))
pred.det.moranB15 <- predict(fmsAve.B15, type = "det", newdata = newDataBrownDry6, appendData = TRUE)

########################################################################
## Tiger dry prediction plots
########################################################################

# predict new data using the averaged model - select values based on min and max of the umf
newDataTigerDry1 <- data.frame(height = seq(-2.23265, 2.80800, by = 0.01), moran = 0)
pred.occ.heightT15 <- predict(fmsAve.T15, type = "state", newdata = newDataTigerDry1, appendData = TRUE)

newDataTigerDry2 <- data.frame(height = 0, moran = seq(-2.62339, 2.51316, by = 0.01))
pred.occ.moranT15 <- predict(fmsAve.T15, type = "state", newdata = newDataTigerDry2, appendData = TRUE)

########################################################################

newDataTigerDry3 <- data.frame(height = seq(-2.23265, 2.80800, by = 0.01), date = 0, effort = 0, moran = 0)
pred.det.heightT15 <- predict(fmsAve.T15, type = "det", newdata = newDataTigerDry3, appendData = TRUE)

newDataTigerDry4 <- data.frame(height = 0, date = seq(-1.6578, 1.6194, by = 0.01), effort = 0, moran = 0)
pred.det.dateT15 <- predict(fmsAve.T15, type = "det", newdata = newDataTigerDry4, appendData = TRUE)

newDataTigerDry5 <- data.frame(height = 0, date = 0, effort = seq(-1.80297, 2.56551, by = 0.01), moran = 0)
pred.det.effT15 <- predict(fmsAve.T15, type = "det", newdata = newDataTigerDry5, appendData = TRUE)

newDataTigerDry6 <- data.frame(height = 0, date = 0, effort = 0, moran = seq(-2.62339, 2.51316, by = 0.01))
pred.det.moranT15 <- predict(fmsAve.T15, type = "det", newdata = newDataTigerDry6, appendData = TRUE)

########################################################################
## Brown wet prediction plots
########################################################################

# predict new data using the averaged model - select values based on min and max of the umf
newDataBrownWet1 <- data.frame(height = seq(-2.23265, 2.80800, by = 0.01), moran = 0)
pred.occ.heightB16 <- predict(fmsAve.B16, type = "state", newdata = newDataBrownWet1, appendData = TRUE)

newDataBrownWet2 <- data.frame(height = 0, moran = seq(-2.62339, 2.51316, by = 0.01))
pred.occ.moranB16 <- predict(fmsAve.B16, type = "state", newdata = newDataBrownWet2, appendData = TRUE)

########################################################################

newDataBrownWet3 <- data.frame(height = seq(-2.23265, 2.80800, by = 0.01), date = 0, effort = 0, moran = 0)
pred.det.heightB16 <- predict(fmsAve.B16, type = "det", newdata = newDataBrownWet3, appendData = TRUE)

newDataBrownWet4 <- data.frame(height = 0, date = seq(-1.6578, 1.6194, by = 0.01), effort = 0, moran = 0)
pred.det.dateB16 <- predict(fmsAve.B16, type = "det", newdata = newDataBrownWet4, appendData = TRUE)

newDataBrownWet5 <- data.frame(height = 0, date = 0, effort = seq(-1.80297, 2.56551, by = 0.01), moran = 0)
pred.det.effB16 <- predict(fmsAve.B16, type = "det", newdata = newDataBrownWet5, appendData = TRUE)

newDataBrownWet6 <- data.frame(height = 0, date = 0, effort = 0, moran = seq(-2.62339, 2.51316, by = 0.01))
pred.det.moranB16 <- predict(fmsAve.B16, type = "det", newdata = newDataBrownWet6, appendData = TRUE)

########################################################################
## Tiger wet prediction plots
########################################################################

# predict new data using the averaged model - select values based on min and max of the umf
newDataTigerWet1 <- data.frame(height = seq(-2.23265, 2.80800, by = 0.01), moran = 0)
pred.occ.heightT16 <- predict(fmsAve.T16, type = "state", newdata = newDataTigerWet1, appendData = TRUE)

newDataTigerWet2 <- data.frame(height = 0, moran = seq(-2.62339, 2.51316, by = 0.01))
pred.occ.moranT16 <- predict(fmsAve.T16, type = "state", newdata = newDataTigerWet2, appendData = TRUE)

########################################################################

newDataTigerWet3 <- data.frame(height = seq(-2.23265, 2.80800, by = 0.01), date = 0, effort = 0, moran = 0)
pred.det.heightT16 <- predict(fmsAve.T16, type = "det", newdata = newDataTigerWet3, appendData = TRUE)

newDataTigerWet5 <- data.frame(height = 0, date = 0, effort = seq(-1.80297, 2.56551, by = 0.01), moran = 0)
pred.det.effT16 <- predict(fmsAve.T16, type = "det", newdata = newDataTigerWet5, appendData = TRUE)

newDataTigerWet6 <- data.frame(height = 0, date = 0, effort = 0, moran = seq(-2.62339, 2.51316, by = 0.01))
pred.det.moranT16 <- predict(fmsAve.T16, type = "det", newdata = newDataTigerWet6, appendData = TRUE)

#################################
# Plots for the paper

######################################################
# Calculate the means and sd of the covariates - to transform from scaled data
######################################################

height.mean <- mean(cov$canopy_height_mean)
height.sd <- sd(cov$canopy_height_mean)
height.z <- (cov$canopy_height_mean-height.mean)/height.sd

xticksH <- -2:3
xlabsH <- xticksH*height.sd + height.mean
xlabsH <- round(xlabsH, 0)

moran.mean <- mean(cov$canopy_height_moran)
moran.sd <- sd(cov$canopy_height_moran)
moran.z <- (cov$canopy_height_moran-moran.mean)/moran.sd

xticksM <- -2:2
xlabsM <- xticksM*moran.sd + moran.mean
xlabsM <- round(xlabsM, 1)

moran.mean <- mean(cov$canopy_height_moran)
moran.sd <- sd(cov$canopy_height_moran)
moran.z <- (cov$canopy_height_moran-moran.mean)/moran.sd

xticksM <- -2:2
xlabsM <- xticksM*moran.sd + moran.mean
xlabsM <- round(xlabsM, 1)

sc6 <- as.numeric(rbind(brown15[,12], brown15[,13], brown15[,14], brown15[,15]))
effort.mean <- mean(sc6)
effort.sd <- sd(sc6)
effort.z <- (sc6 - effort.mean) / effort.sd

xticksE <- -2:2
xlabsE <- xticksE*effort.sd + effort.mean
xlabsE <- round(xlabsE, 0)

######################################################
# OCCUPANCY PLOTS
######################################################

# plots 
B1 <- ggplot(pred.occ.heightB15) + geom_line(aes(height,Predicted), size = 1, colour = "saddlebrown") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, x = height), alpha = 0.2, fill = "saddlebrown") +
  ylim(0,1) + scale_x_continuous(breaks = xticksH, labels = xlabsH) + labs(tag = "A)") + 
  theme_bw() + xlab("") + ylab("Occupancy probability") + theme(axis.title.y = element_text(margin = margin(c(0, 20, 0, 0))),
                                                                text = element_text(size = 18),
                                                                plot.tag.position = c(0.038,0.99))

B2 <- ggplot(pred.occ.moranB15) + geom_line(aes(moran,Predicted), size = 1, colour = "saddlebrown") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, x = moran), alpha = 0.2, fill = "saddlebrown") +
  ylim(0,1) + scale_x_continuous(breaks = xticksM, labels = xlabsM) + labs(tag = "B)") +
  theme_bw() + xlab("") + ylab("") + theme(axis.title.y = element_text(margin = margin(c(0, 20, 0, 0))),
                                          text = element_text(size = 18),
                                          plot.tag.position = c(0.038,0.99))

T1 <- ggplot(pred.occ.heightT15) + geom_line(aes(height,Predicted), size = 1, colour = "forestgreen") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, x = height), alpha = 0.2, fill = "forestgreen") +
  ylim(0,1) + scale_x_continuous(breaks = xticksH, labels = xlabsH) + labs(tag = "C)") +
  theme_bw() + xlab("") + ylab("Occupancy probability") + theme(axis.title.y = element_text(margin = margin(c(0, 20, 0, 0))),
                                                                text = element_text(size = 18),
                                                                plot.tag.position = c(0.038,0.99))

T2 <- ggplot(pred.occ.moranT15) + geom_line(aes(moran,Predicted), size = 1, colour = "forestgreen") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, x = moran), alpha = 0.2, fill = "forestgreen") +
  ylim(0,1) + scale_x_continuous(breaks = xticksM, labels = xlabsM) + labs(tag = "D)") +
  theme_bw() + xlab("") + ylab("") + theme(axis.title.y = element_text(margin = margin(c(0, 20, 0, 0))),
                                           text = element_text(size = 18),
                                           plot.tag.position = c(0.04,0.99))

B7 <- ggplot(pred.occ.heightB16) + geom_line(aes(height,Predicted), size = 1, colour = "saddlebrown") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, x = height), alpha = 0.2, fill = "saddlebrown") +
  ylim(0,1) + scale_x_continuous(breaks = xticksH, labels = xlabsH) + labs(tag = "E)") +
  theme_bw() + xlab("") + ylab("Occupancy probability") + theme(axis.title.y = element_text(margin = margin(c(0, 20, 0, 0))),
                                                                text = element_text(size = 18),
                                                                plot.tag.position = c(0.038,0.99))

B8 <- ggplot(pred.occ.moranB16) + geom_line(aes(moran,Predicted), size = 1, colour = "saddlebrown") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, x = moran), alpha = 0.2, fill = "saddlebrown") +
  ylim(0,1) + scale_x_continuous(breaks = xticksM, labels = xlabsM) + labs(tag = "F)") + 
  theme_bw() + xlab("") + ylab("") + theme(axis.title.y = element_text(margin = margin(c(0, 20, 0, 0))),
                                           text = element_text(size = 18),
                                           plot.tag.position = c(0.04,0.99))

T7 <- ggplot(pred.occ.heightT16) + geom_line(aes(height,Predicted), size = 1, colour = "forestgreen") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, x = height), alpha = 0.2, fill = "forestgreen") +
  ylim(0,1) + scale_x_continuous(breaks = xticksH, labels = xlabsH) + labs(tag = "G)") +
  theme_bw() + xlab("Canopy height") + ylab("Occupancy probability") + theme(axis.title.y = element_text(margin = margin(c(0, 20, 0, 0))),
                                                                text = element_text(size = 18),
                                                                plot.tag.position = c(0.038,0.99),
                                                                axis.title.x = element_text(margin = margin(c(5, 0, 0, 0))))

T8 <- ggplot(pred.occ.moranT16) + geom_line(aes(moran,Predicted), size = 1, colour = "forestgreen") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, x = moran), alpha = 0.2, fill = "forestgreen") +
  ylim(0,1) + scale_x_continuous(breaks = xticksM, labels = xlabsM) + labs(tag = "H)") +
  theme_bw() + xlab("Habitat heterogeneity") + ylab("") + theme(axis.title.y = element_text(margin = margin(c(0, 20, 0, 0))),
                                           text = element_text(size = 18),
                                           plot.tag.position = c(0.038,0.99),
                                           axis.title.x = element_text(margin = margin(c(5, 0, 0, 0))))

Occ <- ggarrange(B1, B2, T1, T2,
                 B7, B8, T7, T8, nrow = 4, ncol = 2)
Occ

ggsave("figures/ModelAveFigures/OccPlotsColour.png", height = 17.0, width = 11.5)
ggsave("figures/ModelAveFigures/OccPlotsColour.pdf", height = 17.0, width = 11.5)
ggsave("figures/ModelAveFigures/OccPlotsColour.tiff", height = 17.0, width = 11.5)

######################################################
# DETECTION PLOTS
######################################################

B3 <- ggplot(pred.det.heightB15) + geom_line(aes(height,Predicted), size = 1, colour = "saddlebrown") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, x = height), alpha = 0.2, fill = "saddlebrown") +
  ylim(0,1) + scale_x_continuous(breaks = xticksH, labels = xlabsH) + labs(tag = "A)") +
  theme_bw() + xlab("") + ylab("Detection probability") + theme(axis.title.y = element_text(margin = margin(c(0, 20, 0, 0))),
                                                                             text = element_text(size = 18),
                                                                             plot.tag.position = c(0.038,0.99),
                                                                             axis.title.x = element_text(margin = margin(c(5, 0, 0, 0))))

B6 <- ggplot(pred.det.moranB15) + geom_line(aes(moran,Predicted), size = 1, colour = "saddlebrown") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, x = moran), alpha = 0.2, fill = "saddlebrown") +
  ylim(0,1) + scale_x_continuous(breaks = xticksM, labels = xlabsM) + labs(tag = "B)") +
  theme_bw() + xlab("") + ylab("") + theme(axis.title.y = element_text(margin = margin(c(0, 20, 0, 0))),
                                           text = element_text(size = 18),
                                           plot.tag.position = c(0.038,0.99),
                                           axis.title.x = element_text(margin = margin(c(5, 0, 0, 0))))

B5 <- ggplot(pred.det.effB15) + geom_line(aes(effort,Predicted), size = 1, colour = "saddlebrown") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, x = effort), alpha = 0.2, fill = "saddlebrown") +
  ylim(0,1) + scale_x_continuous(breaks = xticksE, labels = xlabsE) + labs(tag = "C)") +
  theme_bw() + xlab("") + ylab("") + theme(axis.title.y = element_text(margin = margin(c(0, 20, 0, 0))),
                                           text = element_text(size = 18),
                                           plot.tag.position = c(0.038,0.99),
                                           axis.title.x = element_text(margin = margin(c(5, 0, 0, 0))))

T3 <- ggplot(pred.det.heightT15) + geom_line(aes(height,Predicted), size = 1, colour = "forestgreen") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, x = height), alpha = 0.2, fill = "forestgreen") +
  ylim(0,1) + scale_x_continuous(breaks = xticksH, labels = xlabsH) + labs(tag = "D)") +
  theme_bw() + xlab("") + ylab("Detection probability") + theme(axis.title.y = element_text(margin = margin(c(0, 20, 0, 0))),
                                                                             text = element_text(size = 18),
                                                                             plot.tag.position = c(0.038,0.99),
                                                                             axis.title.x = element_text(margin = margin(c(5, 0, 0, 0))))

T6 <- ggplot(pred.det.moranT15) + geom_line(aes(moran,Predicted), size = 1, colour = "forestgreen") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, x = moran), alpha = 0.2, fill = "forestgreen") +
  ylim(0,1) + scale_x_continuous(breaks = xticksM, labels = xlabsM) + labs(tag = "E)") +
  theme_bw() + xlab("") + ylab("") + theme(axis.title.y = element_text(margin = margin(c(0, 20, 0, 0))),
                                           text = element_text(size = 18),
                                           plot.tag.position = c(0.038,0.99),
                                           axis.title.x = element_text(margin = margin(c(5, 0, 0, 0))))

T5 <- ggplot(pred.det.effT15) + geom_line(aes(effort,Predicted), size = 1, colour = "forestgreen") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, x = effort), alpha = 0.2, fill = "forestgreen") +
  ylim(0,1) + scale_x_continuous(breaks = xticksE, labels = xlabsE) + labs(tag = "F)") +
  theme_bw() + xlab("") + ylab("") + theme(axis.title.y = element_text(margin = margin(c(0, 20, 0, 0))),
                                                                    text = element_text(size = 18),
                                                                    plot.tag.position = c(0.038,0.99),
                                                                    axis.title.x = element_text(margin = margin(c(5, 0, 0, 0))))

B9 <- ggplot(pred.det.heightB16) + geom_line(aes(height,Predicted), size = 1, colour = "saddlebrown") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, x = height), alpha = 0.2, fill = "saddlebrown") +
  ylim(0,1) + scale_x_continuous(breaks = xticksH, labels = xlabsH) + labs(tag = "G)") +
  theme_bw() + xlab("") + ylab("Detection probability") + theme(axis.title.y = element_text(margin = margin(c(0, 20, 0, 0))),
                                                                             text = element_text(size = 18),
                                                                             plot.tag.position = c(0.038,0.99),
                                                                             axis.title.x = element_text(margin = margin(c(5, 0, 0, 0))))
B12 <- ggplot(pred.det.moranB16) + geom_line(aes(moran,Predicted), size = 1, colour = "saddlebrown") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, x = moran), alpha = 0.2, fill = "saddlebrown") +
  ylim(0,1) + scale_x_continuous(breaks = xticksM, labels = xlabsM) + labs(tag = "H)") +
  theme_bw() + xlab("") + ylab("") + theme(axis.title.y = element_text(margin = margin(c(0, 20, 0, 0))),
                                           text = element_text(size = 18),
                                           plot.tag.position = c(0.038,0.99),
                                           axis.title.x = element_text(margin = margin(c(5, 0, 0, 0))))

B11 <- ggplot(pred.det.effB16) + geom_line(aes(effort,Predicted), size = 1, colour = "saddlebrown") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, x = effort), alpha = 0.2, fill = "saddlebrown") +
  ylim(0,1) + scale_x_continuous(breaks = xticksE, labels = xlabsE) + labs(tag = "I)") +
  theme_bw() + xlab("") + ylab("") + theme(axis.title.y = element_text(margin = margin(c(0, 20, 0, 0))),
                                           text = element_text(size = 18),
                                           plot.tag.position = c(0.038,0.99),
                                           axis.title.x = element_text(margin = margin(c(5, 0, 0, 0))))


T9 <- ggplot(pred.det.heightT16) + geom_line(aes(height,Predicted), size = 1, colour = "forestgreen") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, x = height), alpha = 0.2, fill = "forestgreen") +
  ylim(0,1) + scale_x_continuous(breaks = xticksH, labels = xlabsH) + labs(tag = "J)") +
  theme_bw() + xlab("Canopy height") + ylab("Detection probability") + theme(axis.title.y = element_text(margin = margin(c(0, 20, 0, 0))),
                                                                             text = element_text(size = 18),
                                                                             plot.tag.position = c(0.038,0.99),
                                                                             axis.title.x = element_text(margin = margin(c(5, 0, 0, 0))))

T12 <- ggplot(pred.det.moranT16) + geom_line(aes(moran,Predicted), size = 1, colour = "forestgreen") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, x = moran), alpha = 0.2, fill = "forestgreen") +
  ylim(0,1) + scale_x_continuous(breaks = xticksM, labels = xlabsM) + labs(tag = "K)") +
  theme_bw() + xlab("Habitat heterogeneity") + ylab("") + theme(axis.title.y = element_text(margin = margin(c(0, 20, 0, 0))),
                                                                text = element_text(size = 18),
                                                                plot.tag.position = c(0.038,0.99),
                                                                axis.title.x = element_text(margin = margin(c(5, 0, 0, 0))))

T11 <- ggplot(pred.det.effT16) + geom_line(aes(effort,Predicted), size = 1, colour = "forestgreen") + 
  geom_ribbon(aes(ymin = lower, ymax = upper, x = effort), alpha = 0.2, fill = "forestgreen") +
  ylim(0,1) + scale_x_continuous(breaks = xticksE, labels = xlabsE) + labs(tag = "L)") +
  theme_bw() + xlab("Sampling effort") + ylab("") + theme(axis.title.y = element_text(margin = margin(c(0, 20, 0, 0))),
                                                                               text = element_text(size = 18),
                                                                               plot.tag.position = c(0.038,0.99),
                                                                               axis.title.x = element_text(margin = margin(c(5, 0, 0, 0))))

Det2 <- ggarrange(B3, B6, B5, 
                  T3, T6, T5,
                  B9, B12, B11,
                  T9, T12, T11,
                  nrow = 4, ncol = 3)
Det2
ggsave("figures/ModelAveFigures/DetPlotsColour.png", height = 15.0, width = 11.5)
ggsave("figures/ModelAveFigures/DetPlotsColour.pdf", height = 15.0, width = 11.5)
ggsave("figures/ModelAveFigures/DetPlotsColour.tiff", height = 15.0, width = 11.5)
