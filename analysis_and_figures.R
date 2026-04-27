############################################################################
# Characterizing populations of ash specialist insects in forests that
# have experienced over two decades of ash mortality by emerald ash borer
#
# Analyses for abundance and richness of ash specialists
#
# Anthony R Garro & Kayla I Perry
#
# 27 April 2026
############################################################################

# read dataset
df <- read.csv("ash_specialist_data.csv")
str(df)

# change categorical variables to factors
df$plot.id <- as.factor(df$plot.id)
df$ash.species <- as.factor(df$ash.species)
df$canopy.condition <- as.factor(df$canopy.condition)

str(df)
summary(df)

# load packages
library(lme4) #
library(lmerTest) #
library(blmeco)
library(emmeans) #
library(car) #
library(vegan) #
library(hillR) #
library(tidyverse)
library(performance)
library(glmmTMB)
library(DHARMa) #
library(bbmle)
library(fossil) #
library(RColorBrewer) #

colSums(df[c(11:21)])

## Response variables:
# Abundance of specialists (possion)
# Total abundance of insects (only the relevant families) (possion)
# Species richness of specialists (possion)

#________________________________________________________________________________
# Total specialists ----

df$specialist.abundance <- rowSums(df[,c(13,15:16,18)])

# look at the data
dotchart(df$specialist.abundance, group = df$ash.species)
dotchart(df$specialist.abundance, group = df$canopy.condition)
hist(df$specialist.abundance)
boxplot(df$specialist.abundance ~ df$ash.species)
boxplot(df$specialist.abundance ~ df$canopy.condition)
plot(df$specialist.abundance ~ df$dbh)

# run the mod1: dbh, species
total.specialists.mod1 <- glmer(specialist.abundance ~ dbh + ash.species + 
                              offset(log(number.branches.surveyed)) + (1|plot.id), 
                            family = poisson, data = df)

# check assumptions
testDispersion(total.specialists.mod1)
res.total.specialists.mod.1 <- simulateResiduals(total.specialists.mod1)
plot(res.total.specialists.mod.1)
plotResiduals(res.total.specialists.mod.1, form = df$dbh)
testCategorical(res.total.specialists.mod.1, catPred = df$ash.species)
testZeroInflation(res.total.specialists.mod.1)
testOutliers(res.total.specialists.mod.1)

# model 1 outputs
summary(total.specialists.mod1)
Anova(total.specialists.mod1)

##
# run the mod2: canopy condition
total.specialists.mod2 <- glmer(specialist.abundance ~ canopy.condition + 
                              offset(log(number.branches.surveyed)) + (1|plot.id), 
                            family = poisson, data = df)

# check assumptions
testDispersion(total.specialists.mod2)
res.total.specialists.mod.2 <- simulateResiduals(total.specialists.mod2)
plot(res.total.specialists.mod.2)
testCategorical(res.total.specialists.mod.2, catPred = df$canopy.condition)
testZeroInflation(res.total.specialists.mod.2)

# model 2 outputs
summary(total.specialists.mod2)
Anova(total.specialists.mod2)
emmeans(total.specialists.mod2, pairwise ~ canopy.condition)

boxplot(specialist.abundance ~ canopy.condition, data = df, 
        xlab = "Canopy Condition", ylab = "Specialist Abundance")
stripchart(specialist.abundance ~ canopy.condition, data = df, pch = 19, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)

#________________________________________________________________________________
# Total insects ----

df$relevantfam.abundance <- rowSums(df[,c(11:12,14,17)])

# look at the data
dotchart(df$relevantfam.abundance, group = df$ash.species)
dotchart(df$relevantfam.abundance, group = df$canopy.condition)
hist(df$relevantfam.abundance)
boxplot(df$relevantfam.abundance ~ df$ash.species)
boxplot(df$relevantfam.abundance ~ df$canopy.condition)
plot(df$relevantfam.abundance ~ df$dbh)

# run the mod1: dbh, species
total.insects.mod1 <- glmer(relevantfam.abundance ~ dbh + ash.species + 
                              offset(log(number.branches.surveyed)) + (1|plot.id), 
                            family = poisson, data = df)

# check assumptions
testDispersion(total.insects.mod1)
res.total.insects.mod.1 <- simulateResiduals(total.insects.mod1)
plot(res.total.insects.mod.1)
testCategorical(res.total.insects.mod.1, catPred = df$ash.species)
plotResiduals(res.total.insects.mod.1, form = df$dbh)
testZeroInflation(res.total.insects.mod.1)

# model 1 outputs
summary(total.insects.mod1)
Anova(total.insects.mod1)

##
# run the mod2: canopy condition
total.insects.mod2 <- glmer(relevantfam.abundance ~ canopy.condition + 
                              offset(log(number.branches.surveyed)) + (1|plot.id), 
                            family = poisson, data = df)

# check assumptions
testDispersion(total.insects.mod2)
res.total.insects.mod.2 <- simulateResiduals(total.insects.mod2)
plot(res.total.insects.mod.2)
testCategorical(res.total.insects.mod.2, catPred = df$canopy.condition)
testZeroInflation(res.total.insects.mod.2)

# model 1 outputs
summary(total.insects.mod1)
Anova(total.insects.mod1)

# model 2 outputs
summary(total.insects.mod2)
Anova(total.insects.mod2)

boxplot(relevantfam.abundance ~ canopy.condition, data = df, 
        xlab = "Canopy Condition", ylab = "Total Insects")
stripchart(relevantfam.abundance ~ canopy.condition, data = df, pch = 19, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)

#________________________________________________________________________________
# Richness of specialists ----

df$specialist.richness <- hill_taxa(df[,c(13,15:16,18)], q = 0, MARGIN = 1)

# look at the data
dotchart(df$specialist.richness, group = df$ash.species)
dotchart(df$specialist.richness, group = df$canopy.condition)
hist(df$specialist.richness)
boxplot(df$specialist.richness ~ df$ash.species)
boxplot(df$specialist.richness ~ df$canopy.condition)
plot(df$specialist.richness ~ df$dbh)

# run the mod1: dbh, species
rich.mod1 <- glmer(specialist.richness ~ dbh + ash.species + 
                                  offset(log(number.branches.surveyed)) + (1|plot.id), 
                                family = poisson, data = df)

# check assumptions
testDispersion(rich.mod1)
res.rich.mod.1 <- simulateResiduals(rich.mod1)
plot(res.rich.mod.1)
testCategorical(res.rich.mod.1, catPred = df$ash.species)
plotResiduals(res.rich.mod.1, form = df$dbh)
testZeroInflation(res.rich.mod.1)

# model 1 outputs
summary(rich.mod1)
Anova(rich.mod1)

##
# run the mod2: canopy condition
rich.mod2 <- glmer(specialist.richness ~ canopy.condition + 
                                  offset(log(number.branches.surveyed)) + (1|plot.id), 
                                family = poisson, data = df)

# check assumptions
testDispersion(rich.mod2)
res.rich.mod.2 <- simulateResiduals(rich.mod2)
plot(res.rich.mod.2)
testCategorical(res.rich.mod.2, catPred = df$ash.species)
testCategorical(res.rich.mod.2, catPred = df$canopy.condition)
testZeroInflation(res.rich.mod.2)

# model 2 outputs
summary(rich.mod2)
Anova(rich.mod2)
emmeans(rich.mod2, pairwise ~ canopy.condition)

boxplot(specialist.richness ~ canopy.condition, data = df, 
        xlab = "Canopy Condition", ylab = "Specialist Richness")
stripchart(specialist.richness ~ canopy.condition, data = df, pch = 19, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)

################################################################################
# rarefaction analysis
str(df)

green <- df[which(df$ash.species == "green"),]
black <- df[which(df$ash.species == "black"),]

sp.green <- specaccum(green[,c(13,15:16,18)], method = "rarefaction", permutations = 100, gamma = "jack2")
sp.black <- specaccum(black[,c(13,15:16,18)], method = "rarefaction", permutations = 100, gamma = "jack2")


# make the plot

plot(sp.black, pch = 19, col = "black", xvar = c("individuals"), lty = 4, lwd = 2,
     ylab = "Species Richness", xlab = "Number of Individuals", xlim = c(0, 40), ylim = c(0, 5))
plot(sp.green, add = TRUE, pch = 15, xvar = c("individuals"), lty = 1, lwd = 2, col = "#73D055FF")

legend("bottomright", legend = c("Green ash", "Black ash"),
       lty = c(1,2,3,4), cex = 1.5, bty = "n", lwd = 4,
       col = c("#73D055FF", "black"))


# specialists
jack1(green[,c(13,15:16,18)], taxa.row = FALSE, abund = TRUE)
jack2(green[,c(13,15:16,18)], taxa.row = FALSE, abund = TRUE)

jack1(black[,c(13,15:16,18)], taxa.row = FALSE, abund = TRUE)
jack2(black[,c(13,15:16,18)], taxa.row = FALSE, abund = TRUE)


brewer.pal(9, "Greens")

png("Figures/Fig_Rarefaction_Specialists.png", width = 1800, height = 1000, pointsize = 30)

par(mar=c(5,7,4,2))

plot(sp.black, pch = 19, col = "black", xvar = c("individuals"), lty = 4, lwd = 4, cex.lab = 1.5, cex.axis = 1.2,
     ylab = "Species Richness", xlab = "Number of Individuals", xlim = c(0, 40), ylim = c(0, 5))
plot(sp.green, add = TRUE, pch = 15, xvar = c("individuals"), lty = 1, lwd = 4, cex.lab = 1.5, cex.axis = 1.2, col = "#74C476")

legend("topright", legend = c("Green ash", "Black ash"),
       lty = c(1,2,3,4), cex = 1.5, bty = "n", lwd = 4,
       col = c("#74C476", "black"))

dev.off()
