#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
#
# Chapter 2. What are hierarchical models and how do we analyze them?
# =========================================================================


# 2.2 Random variables, probability density functions, statistical models,
#    probability, and statistical inference
# ==========================================================================

dbinom(0:5, size = 5, prob = 0.2)

pnorm(200, mean = 190, sd = 10)

pnorm(200, mean = 190, sd = 10) - pnorm(180, mean = 190, sd = 10)

f <- function(x, mu, sigma){
 (1 / sqrt(2*pi*sigma^2)) * exp( -((x-mu)^2)/(2*sigma^2))
}

integrate(f, lower = 180, upper = 200, mu = 190, sigma = 10)


# 2.2.1 Statistical models
# ------------------------------------------------------------------------

# 2.2.2 Joint, marginal, and conditional distributions
# ------------------------------------------------------------------------
Y <- 0:5   # Possible values of Y (# surveys with peregrine sightings)
X <- 0:5   # Possible values of X (# fledged young)
p <- plogis(-1.2 + 2*X) # p as function of X
round(p, 2)


# Joint distribution [Y, X]
lambda <- 0.4
joint <- matrix(NA, length(Y), length(X))
rownames(joint) <- paste("y=", Y, sep="")
colnames(joint) <- paste("x=", X, sep="")
for(i in 1:length(Y)) {
  joint[,i] <- dbinom(Y, 5, p[i]) * dpois(X[i], lambda)
}

round(joint, 3)

margX <- colSums(joint)
round(margX, 4)

margY <- rowSums(joint)
round(margY, 4)

YgivenX <- joint / matrix(margX, nrow(joint), ncol(joint), byrow=TRUE)
round(YgivenX, 2)


# 2.2.3 Statistical inference (no code)

