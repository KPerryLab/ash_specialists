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
library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(vegan)
library(hillR)
library(DHARMa)
library(fossil)
library(RColorBrewer)

colSums(df[c(11:18)])

## Response variables:
# Abundance of specialists
# Total abundance of insects (only the relevant families)
# Species richness of specialists

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
plotQQunif(res.total.specialists.mod.1)
plotResiduals(res.total.specialists.mod.1, form = df$dbh)
testCategorical(res.total.specialists.mod.1, catPred = df$ash.species)
testZeroInflation(res.total.specialists.mod.1)
testOutliers(res.total.specialists.mod.1)

# model 1 outputs
summary(total.specialists.mod1)
Anova(total.specialists.mod1)

# random intercepts
rand_ints_total.specialists.mod1 <- ranef(total.specialists.mod1)$plot.id$`(Intercept)`
hist(rand_ints_total.specialists.mod1)
qqnorm(rand_ints_total.specialists.mod1)
qqline(rand_ints_total.specialists.mod1)


##
# run the mod2: canopy condition
total.specialists.mod2 <- glmer(specialist.abundance ~ canopy.condition + 
                              offset(log(number.branches.surveyed)) + (1|plot.id), 
                            family = poisson, data = df)

# check assumptions
testDispersion(total.specialists.mod2)
res.total.specialists.mod.2 <- simulateResiduals(total.specialists.mod2)
plotQQunif(res.total.specialists.mod.2)
testCategorical(res.total.specialists.mod.2, catPred = df$canopy.condition)
testZeroInflation(res.total.specialists.mod.2)

# model 2 outputs
summary(total.specialists.mod2)
Anova(total.specialists.mod2)
emmeans(total.specialists.mod2, pairwise ~ canopy.condition)

rand_ints_total.specialists.mod2 <- ranef(total.specialists.mod2)$plot.id$`(Intercept)`
hist(rand_ints_total.specialists.mod2)
qqnorm(rand_ints_total.specialists.mod2)
qqline(rand_ints_total.specialists.mod2)

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
plotQQunif(res.total.insects.mod.1)
testCategorical(res.total.insects.mod.1, catPred = df$ash.species)
plotResiduals(res.total.insects.mod.1, form = df$dbh)
testZeroInflation(res.total.insects.mod.1)

# model 1 outputs
summary(total.insects.mod1)
Anova(total.insects.mod1)

rand_ints_total.insects.mod1 <- ranef(total.insects.mod1)$plot.id$`(Intercept)`
hist(rand_ints_total.insects.mod1)
qqnorm(rand_ints_total.insects.mod1)
qqline(rand_ints_total.insects.mod1)

##
# run the mod2: canopy condition
total.insects.mod2 <- glmer(relevantfam.abundance ~ canopy.condition + 
                              offset(log(number.branches.surveyed)) + (1|plot.id), 
                            family = poisson, data = df)

# check assumptions
testDispersion(total.insects.mod2)
res.total.insects.mod.2 <- simulateResiduals(total.insects.mod2)
plotQQunif(res.total.insects.mod.2)
testCategorical(res.total.insects.mod.2, catPred = df$canopy.condition)
testZeroInflation(res.total.insects.mod.2)

# model 2 outputs
summary(total.insects.mod2)
Anova(total.insects.mod2)

rand_ints_total.insects.mod2 <- ranef(total.insects.mod2)$plot.id$`(Intercept)`
hist(rand_ints_total.insects.mod2)
qqnorm(rand_ints_total.insects.mod2)
qqline(rand_ints_total.insects.mod2)

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
plotQQunif(res.rich.mod.1)
testCategorical(res.rich.mod.1, catPred = df$ash.species)
plotResiduals(res.rich.mod.1, form = df$dbh)
testZeroInflation(res.rich.mod.1)

# model 1 outputs
summary(rich.mod1)
Anova(rich.mod1)

rand_ints_rich.mod1 <- ranef(rich.mod1)$plot.id$`(Intercept)`
hist(rand_ints_rich.mod1)
qqnorm(rand_ints_rich.mod1)
qqline(rand_ints_rich.mod1)

##
# run the mod2: canopy condition
rich.mod2 <- glmer(specialist.richness ~ canopy.condition + 
                                  offset(log(number.branches.surveyed)) + (1|plot.id), 
                                family = poisson, data = df)

# check assumptions
testDispersion(rich.mod2)
res.rich.mod.2 <- simulateResiduals(rich.mod2)
plotQQunif(res.rich.mod.2)
testCategorical(res.rich.mod.2, catPred = df$ash.species)
testCategorical(res.rich.mod.2, catPred = df$canopy.condition)
testZeroInflation(res.rich.mod.2)

# model 2 outputs
summary(rich.mod2)
Anova(rich.mod2)
emmeans(rich.mod2, pairwise ~ canopy.condition)

rand_ints_rich.mod2 <- ranef(rich.mod2)$plot.id$`(Intercept)`
hist(rand_ints_rich.mod2)
qqnorm(rand_ints_rich.mod2)
qqline(rand_ints_rich.mod2)

boxplot(specialist.richness ~ canopy.condition, data = df, 
        xlab = "Canopy Condition", ylab = "Specialist Richness")
stripchart(specialist.richness ~ canopy.condition, data = df, pch = 19, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)


png("Figures/Specialists_canopy_condition_panel.png", width = 2500, height = 2500, pointsize = 30)

par(mfrow=c(2,1))
par(bty = "n")
par(mar=c(5,7,0.5,2))
# c(bottom, left, top, right)

boxplot(specialist.abundance ~ canopy.condition, data = df,
        col = c("#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#238B45"),
        ylim = c(0,5), ylab = "Specialist Abundance", xlab = "", cex.lab = 2, cex.axis = 1.7, xaxt = "n")
axis(side = 1, labels = FALSE)
stripchart(specialist.abundance ~ canopy.condition, data = df, pch = 19, cex = 2.2, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)
text(0.45,4.7, "A", pos = 3, font = 2, cex = 2.4)
text(1,4.2, "a", pos = 3, font = 1, cex = 2)
text(5,4.2, "b", pos = 3, font = 1, cex = 2)
text(2,3.2, "ab", pos = 3, font = 1, cex = 2)
text(3,3.2, "ab", pos = 3, font = 1, cex = 2)
text(4,2.2, "ab", pos = 3, font = 1, cex = 2)

boxplot(specialist.richness ~ canopy.condition, data = df,
        col = c("#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#2171B5"),
        ylim = c(0,4), ylab = "Specialist Richness", xlab = "Canopy Condition", cex.lab = 2, cex.axis = 1.7)
stripchart(specialist.richness ~ canopy.condition, data = df, pch = 19, cex = 2.2, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)
text(0.45,3.75, "B", pos = 3, font = 2, cex = 2.4)
text(1,3.2, "a", pos = 3, font = 1, cex = 2)
text(4,1.2, "b", pos = 3, font = 1, cex = 2)
text(5,1.2, "b", pos = 3, font = 1, cex = 2)
text(2,2.2, "ab", pos = 3, font = 1, cex = 2)
text(3,2.2, "ab", pos = 3, font = 1, cex = 2)

dev.off()


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
