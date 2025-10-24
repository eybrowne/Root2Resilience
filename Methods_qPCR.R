
library(readxl)
library(car)
library(agricolae)
library(ggplot2)
library(scales)

#load data with averages of technical replicates
Pilot <- read_excel("Pilot Trial1.xlsx", 
                    sheet = "Average values ")
View(Pilot)

#subset by Gene
d.data_ITS <- subset(Pilot,Gene == 'ITS')
d.data_AOA <- subset(Pilot,Gene == 'AOA')
d.data_16S <- subset(Pilot,Gene == '16S')
d.data_NOSZ2 <- subset(Pilot,Gene == 'NOSZ2')


# Perform Shapiro-Wilk test
shapiro_resultNOSZ2 <- shapiro.test(d.data_NOSZ2$Abundance)
print(shapiro_resultNOZS2)

# Print the result
print(shapiro_resultITS)
#data:  d.data_ITS$Abundance
#W = 0.77801, p-value = 9.899e-06
print(shapiro_resultAOA)
#data:  d.data_AOA$Abundance
#W = 0.85216, p-value = 0.00038
print(shapiro_result16S)
#data:  d.data_16S$Abundance
#W = 0.87525, p-value = 0.001084
print(shapiro_resultNOSZ2)
#data:  d.data_NOSZ2$Abundance
#W = 0.95501, p-value = 0.1997

###### ITS, AOA, 16S, all non-normal, I will run an ANOVA and check residuals 

#ANOVA
r.outITS <-aov(Abundance ~ Protocol, data = d.data_ITS)
summary(r.outITS)
HSD.test(r.outITS, trt = c("Protocol"), console = TRUE)   #Tukey test for multiple comparisons

r.out16S <-aov(Abundance ~ Protocol, data = d.data_16S)
summary(r.out16S)
HSD.test(r.out16S, trt = c("Protocol"), console = TRUE) 

r.outAOA <-aov(Abundance ~ Protocol, data = d.data_AOA)
summary(r.outAOA)
HSD.test(r.outAOA, trt = c("Protocol"), console = TRUE) 

r.outNOSZ2 <-aov(Abundance ~ Protocol, data = d.data_NOSZ2)
summary(r.outNOSZ2)
HSD.test(r.outNOSZ2, trt = c("Protocol"), console = TRUE) 


#Plotting residuals
par(mfrow=c(2,2))
plot(fitted(r.outITS), resid(r.outITS)/sd(resid(r.outITS))) ; abline(h=0)
lines( lowess(fitted(r.outITS), resid(r.outITS)/sd(resid(r.outITS))), lty=2, col=2 )
title("Tuckey-Anscombe")
#identify(fitted(r.outITS), resid(r.outITS)/sd(resid(r.outITS)))
qqPlot(resid(r.outITS))
title("QQ")
lever <- lm.influence(r.outITS)$hat
plot(lever, resid(r.outITS)) ; abline(h=0) ; abline(v=2*mean(lever), lty=2, col=2)
title("Leverage")
#identify(lever, resid(r.outITS))
hist(resid(r.outITS)/sd(resid(r.outITS)), freq=FALSE, main="", ylim=c(0,0.5))
curve(dnorm, add=TRUE, lty=2, col=2)
title("Resid/sd(Resid)")

###Repeated this for the ANOVA outputs of all genes

#COOKS TESTS to determine if outliers are significantly affecting the model
par(mfrow=c(2,2))
cooks <- plot(cooks.distance(r.outITS))
title("ITS")
cooks <- plot(cooks.distance(r.out16S))
title("16S")
cooks <- plot(cooks.distance(r.outAOA))
title("AOA")
cooks <- plot(cooks.distance(r.outNOSZ2))
title("NOSZ2")

#to find outlier in AOA
cooks.distance(r.outAOA)[which.max(cooks.distance(r.outAOA))]
#12 
#0.6936515

d.data_AOA_pruned <- subset(d.data_AOA, `Sample Name` != "22.EO")
r.outAOA_pruned <-aov(Abundance ~ Protocol, data = d.data_AOA_pruned)
summary(r.outAOA_pruned)
HSD.test(r.outAOA_pruned, trt = c("Protocol"), console = TRUE)

#Assign order 
Pilot <- Pilot[is.finite(Pilot$Abundance), ]
Pilot$Protocol <- factor(Pilot$Protocol, levels = c("Bulk-dip", "Bulk-vortex", "Rinse-dip", "Shake-dip","Shake-vortex", "Shake-cut")) 


########Plotting box plots####

ggplot(Pilot, aes(x = Protocol, y = Abundance, fill = Protocol)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  facet_wrap(~factor(Gene, levels = c("16S", "ITS", "AOA", "NOSZ2")), scales = "free_y") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15)),  # Adds 50% space above y-axis in each facet
    labels = scales::scientific            # Sets y-axis labels to scientific notation
  ) +
  scale_fill_manual(values = c("#D55E00","#0072B2","#D55E00","#E69F00","#56B4E9","#0072B2")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(
    title = "Pre-Trial qPCR",
    x = "Protocol",
    y = "gene copies per ng DNA",
    fill = "Protocol"
  )






