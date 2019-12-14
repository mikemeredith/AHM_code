#   Applied hierarchical modeling in ecology
#   Modeling distribution, abundance and species richness using R and BUGS
#   Volume 1: Prelude and Static models
#   Marc Kéry & J. Andy Royle
# Chapter 6. Modeling abundance with counts of unmarked individuals
#    in closed populations: binomial N-mixture models
# =========================================================================

# 6.16 Exercises
# ------------------------------------------------------------------------

# ~~~~~~ This code snippet will not run as-is ~~~~~~~~~~~~~~~~~~~~

# Solution A:
# range(mhbdata[,12:14], na.rm = TRUE)
# day.mean <- mean(as.matrix(mhbdata[,12:14]), na.rm = TRUE)
# day.sd <- sd(c(as.matrix(mhbdata[,12:14])), na.rm = TRUE)
# original.pred.day <- 15:110
# pred.day <- (original.pred.day - day.mean) / day.sd
# new<- data.frame(day=pred.day)
# pred<-predict(fm31,type="det",newdata=new,appendData=TRUE)
# head(pred)

# plot(Predicted ~ original.pred.day, pred,type="l",xlab="Date (1 = 1 April)", ylab="Expected detection prob",ylim=c(0,1), lwd = 2)
# lines(lower ~ original.pred.day, pred,type="l",col="red", lwd = 2)
# lines(upper ~ original.pred.day, pred,type="l",col="red", lwd = 2)

