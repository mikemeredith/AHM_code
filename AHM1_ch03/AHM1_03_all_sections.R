#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
# Chapter 3 Linear models, generalised linear models (GLMs) and random effects models:
#    the components of hierarchical models
# =========================================================================

# 3.1 Introduction
# ================

# Define data
pop <- factor(c(rep("Navarra", 3), rep("Aragon", 3), rep("Catalonia", 3)),
    levels = c("Navarra", "Aragon", "Catalonia"))         # Population
wing <- c(10.5, 10.6, 11.0, 12.1, 11.7, 13.5, 11.4, 13.0, 12.9) # Wing span
body <- c(6.8, 8.3, 9.2, 6.9, 7.7, 8.9, 6.9, 8.2, 9.2) # Body length
sex <- factor(c("M","F","M","F","M","F","M","F","M"), levels = c("M", "F"))
mites <- c(0, 3, 2, 1, 0, 7, 0, 9, 6)      # Number of ectoparasites
color <- c(0.45, 0.47, 0.54, 0.42, 0.54, 0.46, 0.49, 0.42, 0.57) # Color intensity
damage <- c(0,2,0,0,4,2,1,0,1)                 # Number of wings damaged

cbind(pop, sex, wing, body, mites, color, damage) # Look at data

str(pop)

op <- par(mfrow = c(1, 3), cex = 1.2)
colorM <- c("red", "red", "blue", "green", "green")  # Pop color code males
colorF <- c("red", "blue", "blue", "green", "green") # Pop color code females
plot(body[sex == "M"], wing[sex == "M"], col =colorM, xlim = c(6.5, 9.5),
    ylim = c(10, 14), lwd = 2, frame.plot = FALSE, las = 1, pch = 17,
    xlab = "Body length", ylab = "Wing span")
points(body[sex == "F"], wing[sex == "F"], col =colorF, pch = 16)
text(6.8, 13.8, "A", cex = 1.5)
plot(body[sex == "M"], mites[sex == "M"], col = colorM, xlim = c(6.5, 9.5),
    ylim = c(0, 10), lwd = 2, frame.plot = FALSE, las = 1, pch = 17,
    xlab = "Body length", ylab = "Parasite load")
points(body[sex == "F"], mites[sex == "F"], col = colorF, pch = 16)
text(6.8, 9.7, "B", cex = 1.5)
plot(body[sex == "M"], damage[sex == "M"], col = colorM, xlim = c(6.5, 9.5),
    ylim = c(0, 4), lwd = 2, frame.plot = FALSE, las = 1, pch = 17,
    xlab = "Body length", ylab = "Damaged wings")
points(body[sex == "F"], damage[sex == "F"], col = colorF, pch = 16)
text(6.8, 3.9, "C", cex = 1.5)
par(op)


# 3.2 Linear models
# =================

# 3.2.1 Linear models with main effects of one factor and one continuous covariate
# --------------------------------------------------------------------------------
summary(fm1 <- lm(wing ~ pop + body))

summary(fm2 <- lm(wing ~ pop-1 + body))

cbind(model.matrix(~pop+body) %*% fm1$coef, predict(fm1)) # Compare two solutions

model.matrix(~ pop + body) # Effects parameterisation

model.matrix(~ pop-1 + body) # Means parameterization

op <- par(mfrow = c(1, 3), mar = c(5,4,2,2), cex = 1.2, cex.main = 1)
plot(body[sex == "M"], wing[sex == "M"], col = colorM, xlim = c(6.5, 9.5),
    ylim = c(10, 14), lwd = 2, frame.plot = FALSE, las = 1, pch = 17,
    xlab = "Body length", ylab = "Wing span")
points(body[sex == "F"], wing[sex == "F"], col = colorF, pch = 16)
abline(coef(fm2)[1], coef(fm2)[4], col = "red", lwd = 2)
abline(coef(fm2)[2], coef(fm2)[4], col = "blue", lwd = 2)
abline(coef(fm2)[3], coef(fm2)[4], col = "green", lwd = 2)
text(6.8, 14, "A", cex = 1.5)


# 3.2.2 Linear models with interaction between one factor and one continuous covariate
# ------------------------------------------------------------------------

model.matrix(~ pop*body)  # Effects parameterisation

model.matrix(~ pop*body-1-body)  # Means parameterisation

summary(fm3 <- lm(wing ~ pop*body-1-body))

# Plot
plot(body[sex == "M"], wing[sex == "M"], col = colorM, xlim = c(6.5, 9.5),
    ylim = c(10, 14), lwd = 2, frame.plot = FALSE, las = 1, pch = 17,
    xlab = "Body length", ylab = "")
points(body[sex == "F"], wing[sex == "F"], col = colorF, pch = 16)
abline(coef(fm3)[1], coef(fm3)[4], col = "red", lwd = 2)
abline(coef(fm3)[2], coef(fm3)[5], col = "blue", lwd = 2)
abline(coef(fm3)[3], coef(fm3)[6], col = "green", lwd = 2)
text(6.8, 14, "B", cex = 1.5)

# Create new design matrix
(DM0 <- model.matrix(~ pop*body-1-body)) # Original DM for means param
DM0[7:9,5] <- DM0[7:9,6]                 # Combine slopes for Ar and Cat
(DM1 <- DM0[, -6])                       # Delete former slope column for Cat

# Fit model with partial interaction
summary(fm4 <- lm(wing ~ DM1-1))

# Do significance test
anova(fm3, fm4)             # F test between two models

# Plot
plot(body[sex == "M"], wing[sex == "M"], col = colorM, xlim = c(6.5, 9.5),
    ylim = c(10, 14), lwd = 2, frame.plot = FALSE, las = 1, pch = 17,
    xlab = "Body length", ylab = "")
points(body[sex == "F"], wing[sex == "F"], col = colorF, pch = 16)
abline(coef(fm4)[1], coef(fm4)[4], col = "red", lwd = 2)
abline(coef(fm4)[2], coef(fm4)[5], col = "blue", lwd = 2)
abline(coef(fm4)[3], coef(fm4)[5], col = "green", lwd = 2)
text(6.8, 14, "C", cex = 1.5)
par(op)

# 3.2.3 Linear models with two factors
# ------------------------------------------------------------------------
model.matrix(~ pop+sex)  # Design matrix of main-effects 2-way ANOVA

# Fit linear model with that design matrix
summary(fm5 <- lm(wing ~ pop + sex))

model.matrix(~ pop+sex-1)  # Design matrix of the main-effects 2-way ANOVA

# Fit linear model with that design matrix
summary(fm6 <- lm(wing ~ pop + sex-1))

# Variant 1: Effects parameterisation (R default)
model.matrix(~ pop*sex)
#model.matrix(~ pop + sex + pop:sex)     # Same

# Variant 2: Means param. for pop, effects param. for sex
model.matrix(~ pop*sex-1)

# Variant 3 (output slightly trimmed): full means parameterisation
model.matrix(~ pop:sex-1)


# 3.2.4 Linear models with two continuous covariates and including polynomials
# ------------------------------------------------------------------------
model.matrix(~ body + color)  # main-effects of covariates

summary(fm7 <- lm(wing ~ body + color))  # Fit that model

model.matrix(~ body*color)  # Interaction between two covariates

summary(fm8 <- lm(wing ~ body*color))  # Fit that model

# Cubic polynomial of body in R
body2 <- body^2           # Squared body length
body3 <- body^3           # Cubed body length
model.matrix(~ body + body2 + body3)

summary(fm9 <- lm(wing ~ body + body2 + body3))  # Fit that model
# summary(fm9 <- lm(wing ~ body + I(body^2) + I(body^3))) # same



# 3.3 Generalised linear models (GLMs)
# ====================================

# 3.3.1 Poisson generalised linear model (GLM) for unbounded counts
# ------------------------------------------------------------------------
summary(fm10 <- glm(mites ~ pop-1 + body, family = poisson))


# 3.3.2 Offsets in the Poisson GLM
# ------------------------------------------------------------------------
summary(fm10a <- glm(mites ~ pop-1 + wing, offset = log(body), family = poisson))
# summary(fm10a <- glm(mites ~ offset(log(body)) + pop-1 + wing, family = poisson))     # same


# 3.3.3 Overdispersion and underdispersion (no code)
# 3.3.4 Zero-inflation (no code)

# 3.3.5 Bernoulli GLM: logistic regression for a binary response
# ------------------------------------------------------------------------
presence <- ifelse(mites > 0, 1, 0)  # convert abundance to presence/absence
summary(fm11 <- glm(presence ~ pop-1 + body, family = binomial))

# 3.3.6 Modeling a Poisson process from presence/absence data using a Bernoulli GLM with cloglog link
# ------------------------------------------------------------------------
summary(fm11a <- glm(presence ~ pop-1 + body, family = binomial(link = "cloglog")))

summary(fm10)

# 3.3.7 Binomial GLM: logistic regression for bounded counts
# -------------------------------------------------------------------------
summary(fm12 <- glm(cbind(damage, 4-damage) ~ pop + body -1, family = binomial))


# 3.3.8 The GLM as the quintessential statistical model (no code)

# 3.4 Random effects (mixed) models
# =================================

# 3.4.1 Random effects for a normal data distribution: normal-normal generalised linear mixed model (GLMM)
# ------------------------------------------------------------------------
# Plot data without distinguishing sex
plot(body, wing, col = rep(c("red", "blue", "green"), each = 3), xlim = c(6.5, 9.5),
    ylim = c(10, 14), cex = 1.5, lwd = 2, frame.plot = FALSE, las = 1, pch = 16,
    xlab = "Body length", ylab = "Wing span")

summary(lm <- lm(wing ~ pop-1 + body))     # Same as fm2

library(lme4)
summary(lmm1 <- lmer(wing ~ (1|pop) + body))  # Fit the model
ranef(lmm1)                                   # Print random effects

alpha_j <- fixef(lmm1)[1]+ranef(lmm1)$pop[,1]
cbind(fixed = coef(lm)[1:3], random = alpha_j)

op <- par(lwd = 3)
abline(lm$coef[1], lm$coef[4], col = "red", lty = 2)
abline(lm$coef[2], lm$coef[4], col = "blue", lty = 2)
abline(lm$coef[3], lm$coef[4], col = "green", lty = 2)
abline(alpha_j[1], fixef(lmm1)[2], col = "red")
abline(alpha_j[2], fixef(lmm1)[2], col = "blue")
abline(alpha_j[3], fixef(lmm1)[2], col = "green")
abline(fixef(lmm1), col = "black")
legend(6.5, 14, c("Catalonia", "Aragon", "Navarra"), col=c("blue", "green", "red"),
    lty = 1, pch = 16, bty = "n", cex = 1.5)
par(op)

summary(lmm2 <- lmer(wing ~ body + (1|pop) + (0+body|pop)))

# 3.4.2 Random effects for a Poisson data distribution: normal-Poisson generalised linear mixed model (GLMM)
# ------------------------------------------------------------------------
summary(glmm <- glmer(mites ~ body + (1|pop), family = poisson))

# 3.5 Summary and outlook (no code)

# 3.6 Exercises
# =============

# Define and plot data (10 times larger data set than the toy data set)
clone.size <- 10               # clone size
pop <- factor(rep(c(rep("Navarra", 3), rep("Aragon", 3), rep("Catalonia", 3)),
    levels = c("Navarra", "Aragon", "Catalonia"), clone.size))
wing <- rep(c(10.5, 10.6, 11.0, 12.1, 11.7, 13.5, 11.4, 13.0, 12.9), clone.size)
body <- rep(c(6.8, 8.3, 9.2, 6.9, 7.7, 8.9, 6.9, 8.2, 9.2), clone.size)
sex <- rep(factor(c("M","F","M","F","M","F","M","F","M"), levels = c("M", "F")), clone.size)
mites <- rep(c(0, 3, 2, 1, 0, 7, 0, 9, 6), clone.size)
color <- rep(c(0.45, 0.47, 0.54, 0.42, 0.54, 0.46, 0.49, 0.42, 0.57), clone.size)
damage <- rep(c(0,2,0,0,4,2,1,0,1), clone.size)

