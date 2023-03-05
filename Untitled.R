library(RCurl)
data_url <- "https://raw.githubusercontent.com/kentranz/socialMobilityCOVID/master/data/raw/TorontoCovid.csv"
toronto.covid <- read.csv(text = getURL(data_url),
                          header=TRUE)

############################# FUNCTIONS #############################
getmuhat <- function(sampleXY, complexity = 1){
  formula <- paste0("y ~ ",
                    if (complexity==0) {
                      "1"
                    } else {
                      paste0("poly(x, ", complexity, ", raw = FALSE)") 
                    }
  )
  
  fit <- lm(as.formula(formula), data = sampleXY)
  tx = sampleXY$x
  ty = fit$fitted.values
  
  range.X = range(tx)
  val.rY  = c( mean(ty[tx == range.X[1]]), 
               mean(ty[tx == range.X[2]]) )
  
  ## From this we construct the predictor function
  muhat <- function(x){
    if ("x" %in% names(x)) {
      ## x is a dataframe containing the variate named
      ## by xvarname
      newdata <- x
    } else 
      ## x is a vector of values that needs to be a data.frame
    { newdata <- data.frame(x = x) }
    ## The prediction
    ## 
    val = predict(fit, newdata = newdata)
    val[newdata$x < range.X[1]] = val.rY[1]
    val[newdata$x > range.X[2]] = val.rY[2]
    val
  }
  ## muhat is the function that we need to calculate values 
  ## at any x, so we return this function from getmuhat
  muhat
}

getXYpop <- function(xvarname, yvarname, pop){
  popData <- pop[, c(xvarname, yvarname)]
  names(popData) <- c("x", "y")
  popData
}

getXYSample <- function(xvarname, yvarname, samp, pop){
  sampData <- pop[samp, c(xvarname, yvarname)]
  names(sampData) <- c("x", "y")
  sampData
}

getSampleComp <- function(pop, size, replace=FALSE) {
  N <- nrow(as.data.frame(pop))
  samp <- rep(FALSE, N)
  samp[sample(1:N, size, replace = replace)] <- TRUE
  samp
}

sampSize <- function(samp) {popSize(samp)}

popSize <- function(pop) {nrow(as.data.frame(pop))}

getmubar <- function(muhats) {
  function(x) {
    Ans <- sapply(muhats, FUN=function(muhat){muhat(x)})
    apply(Ans, MARGIN=1, FUN=mean)
  }
}

gettauFun <- function(pop, xvarname, yvarname){
  pop   = na.omit(pop[, c(xvarname, yvarname)])
  
  # rule = 2 means return the nearest y-value when extrapolating, same as above.
  # ties = mean means that repeated x-values have their y-values averaged, as above.
  tauFun = approxfun(pop[,xvarname], pop[,yvarname], rule = 2, ties = mean)
  return(tauFun)
}

apse_all <- function(Ssamples, Tsamples, complexity, tau){
  ## average over the samples S
  ##
  N_S <- length(Ssamples)
  muhats <- lapply(Ssamples, 
                   FUN=function(sample) getmuhat(sample, complexity)
  )
  ## get the average of these, mubar
  mubar <- getmubar(muhats)
  
  rowMeans(sapply(1:N_S, 
                  FUN=function(j){
                    T_j <- Tsamples[[j]]
                    S_j <- Ssamples[[j]]
                    muhat <- muhats[[j]]
                    ## Take care of any NAs
                    T_j <- na.omit(T_j)
                    y <- c(S_j$y, T_j$y)
                    x <- c(S_j$x, T_j$x)
                    
                    tau_x    <- tau(x)
                    muhat_x <- muhat(x)
                    mubar_x <- mubar(x)
                    
                    apse        <- (y - muhat_x)
                    bias2       <- (mubar_x - tau_x)
                    var_mutilde <-  (muhat_x - mubar_x)
                    var_y       <- (y - tau_x)
                    
                    squares <- rbind(apse, var_mutilde, bias2, var_y)^2
                    
                    ## return means
                    rowMeans(squares)
                  }
  ))
}
############################# FUNCTIONS #############################


### QUESTION 2 ###

# a)

date <- toronto.covid$ï..Episode.Date[2:307]
date <- rev(date)
days <- seq(1, 306)
cases <- toronto.covid$Case.Count[2:307]
cases <- rev(cases)
tor <- data.frame(date,days,cases)
colnames(tor) <- c("Date", "Day_Num", "New_Cases")
head(tor, 15)

# b) 


xydata <- data.frame(tor["Day_Num"], tor["New_Cases"])
colnames(xydata) <- c("x", "y")

par(mfrow=c(3,2))

# deg = 1
plot(xydata, main = "New COVID Cases in Toronto\n March 1 - December 31, 2020",
     xlab = "Day_Num", ylab = "New_Cases", cex = 0.8, cex.main = 0.8, 
     col = adjustcolor("darkgrey", 0.5), pch = 16)
muhat1 <- getmuhat(xydata,1)
curve(muhat1, lwd = 3, add = TRUE)
legend("topleft", legend = c("deg = 1"), lwd = 3, cex = 0.8, bty = "n")

# deg = 2
plot(xydata, main = "New COVID Cases in Toronto\n March 1 - December 31, 2020",
     xlab = "Day_Num", ylab = "New_Cases", cex = 0.8, cex.main = 0.8, 
     col = adjustcolor("darkgrey", 0.5), pch = 16)
muhat2 <- getmuhat(xydata,2)
curve(muhat2, add = TRUE, lwd = 3, col = 2)
legend("topleft", legend = c("deg = 2"), lwd = 3, cex = 0.8, bty = "n", col = 2)

# deg = 5
plot(xydata, main = "New COVID Cases in Toronto\n March 1 - December 31, 2020",
     xlab = "Day_Num", ylab = "New_Cases", cex = 0.8, cex.main = 0.8, 
     col = adjustcolor("darkgrey", 0.5), pch = 16)
muhat5 <- getmuhat(xydata,5)
curve(muhat5, add = TRUE, lwd = 3, col = 3)
legend("topleft", legend = c("deg = 5"), lwd = 3, cex = 0.8, bty = "n", col = 3)

# deg = 10
plot(xydata, main = "New COVID Cases in Toronto\n March 1 - December 31, 2020",
     xlab = "Day_Num", ylab = "New_Cases", cex = 0.8, cex.main = 0.8, 
     col = adjustcolor("darkgrey", 0.5), pch = 16)
muhat10 <- getmuhat(xydata,10)
curve(muhat10, add = TRUE, lwd = 3, col = 4)
legend("topleft", legend = c("deg = 10"), lwd = 3, cex = 0.8, bty = "n", col = 4)

# deg = 15
plot(xydata, main = "New COVID Cases in Toronto\n March 1 - December 31, 2020",
     xlab = "Day_Num", ylab = "New_Cases", cex = 0.8, cex.main = 0.8, 
     col = adjustcolor("darkgrey", 0.5), pch = 16)
muhat15 <- getmuhat(xydata,15)
curve(muhat15, add = TRUE, lwd = 3, col = 5)
legend("topleft", legend = c("deg = 15"), lwd = 3, cex = 0.8, bty = "n", col = 5)

# deg = 20
plot(xydata, main = "New COVID Cases in Toronto\n March 1 - December 31, 2020",
     xlab = "Day_Num", ylab = "New_Cases", cex = 0.8, cex.main = 0.8, 
     col = adjustcolor("darkgrey", 0.5), pch = 16)
muhat20 <- getmuhat(xydata,15)
curve(muhat20, add = TRUE, lwd = 3, col = 6)
legend("topleft", legend = c("deg = 20"), lwd = 3, cex = 0.8, bty = "n", col = 6)


# c)

M <- 50
n <- 100
set.seed(341)
samps <- lapply(1:M, FUN= function(i){getSampleComp(tor, n)})
Ssamples <- lapply(samps, FUN= function(Si){getXYSample("Day_Num", "New_Cases", Si, tor)})
Tsamples <- lapply(samps, FUN= function(Si){getXYSample("Day_Num", "New_Cases", !Si, tor)})

# Fit a polynomial of the given complexity/degree to every sample
# and save the results in a list
muhats1 <- lapply(Ssamples, getmuhat, complexity = 1)
muhats2 <- lapply(Ssamples, getmuhat, complexity = 2)
muhats5 <- lapply(Ssamples, getmuhat, complexity = 5)
muhats10 <- lapply(Ssamples, getmuhat, complexity = 10)
muhats15 <- lapply(Ssamples, getmuhat, complexity = 15)
muhats20 <- lapply(Ssamples, getmuhat, complexity = 20)


# d)

par(mfrow=c(3,2))

# deg = 1
plot(xydata, main = "New COVID Cases in Toronto\n March 1 - December 31, 2020",
     xlab = "Day_Num", ylab = "New_Cases", cex = 0.8, cex.main = 0.8, 
     col = adjustcolor("darkgrey", 0.5), pch = 16)
for (i in 1:M) {
  curveFn <- muhats1[[i]]
  curve(curveFn, add = TRUE, lwd = 3, cex = 0.8, bty = "n", col = adjustcolor(1, 0.5))
}

# deg = 2
plot(xydata, main = "New COVID Cases in Toronto\n March 1 - December 31, 2020",
     xlab = "Day_Num", ylab = "New_Cases", cex = 0.8, cex.main = 0.8, 
     col = adjustcolor("darkgrey", 0.5), pch = 16)
for (i in 1:M) {
  curveFn <- muhats2[[i]]
  curve(curveFn, add = TRUE, lwd = 3, cex = 0.8, bty = "n", col = adjustcolor(2, 0.5))
}

# deg = 5
plot(xydata, main = "New COVID Cases in Toronto\n March 1 - December 31, 2020",
     xlab = "Day_Num", ylab = "New_Cases", cex = 0.8, cex.main = 0.8, 
     col = adjustcolor("darkgrey", 0.5), pch = 16)
for (i in 1:M) {
  curveFn <- muhats5[[i]]
  curve(curveFn, add = TRUE, lwd = 3, cex = 0.8, bty = "n", col = adjustcolor(3, 0.5))
}

# deg = 10
plot(xydata, main = "New COVID Cases in Toronto\n March 1 - December 31, 2020",
     xlab = "Day_Num", ylab = "New_Cases", cex = 0.8, cex.main = 0.8, 
     col = adjustcolor("darkgrey", 0.5), pch = 16)
for (i in 1:M) {
  curveFn <- muhats10[[i]]
  curve(curveFn, add = TRUE, lwd = 3, cex = 0.8, bty = "n", col = adjustcolor(4, 0.5))
}

# deg = 15
plot(xydata, main = "New COVID Cases in Toronto\n March 1 - December 31, 2020",
     xlab = "Day_Num", ylab = "New_Cases", cex = 0.8, cex.main = 0.8, 
     col = adjustcolor("darkgrey", 0.5), pch = 16)
for (i in 1:M) {
  curveFn <- muhats15[[i]]
  curve(curveFn, add = TRUE, lwd = 3, cex = 0.8, bty = "n", col = adjustcolor(5, 0.5))
}

# deg = 20
plot(xydata, main = "New COVID Cases in Toronto\n March 1 - December 31, 2020",
     xlab = "Day_Num", ylab = "New_Cases", cex = 0.8, cex.main = 0.8, 
     col = adjustcolor("darkgrey", 0.5), pch = 16)
for (i in 1:M) {
  curveFn <- muhats20[[i]]
  curve(curveFn, add = TRUE, lwd = 3, cex = 0.8, bty = "n", col = adjustcolor(6, 0.3))
}

# e)

tau.x <- gettauFun(xydata, "x", "y")
degrees <- 0:20
apse_vals <- sapply(degrees, FUN = function(complexity) {
  apse_all(Ssamples, Tsamples, complexity = complexity, tau = tau.x)
})

# print results
t(rbind(degrees, aspe = round(apse_vals, 5)))


# f)

par(mfrow=c(1,1))

plot(degrees, sqrt(apse_vals[1, ]), xlab = "Degree", 
     ylab = "", type = "l", col = "purple", lwd = 2, ylim = c(0, max(sqrt(apse_vals))))
lines(degrees, sqrt(apse_vals[2, ]), xlab = "Degree", 
      ylab = "", col = "blue", lwd = 2)
lines(degrees, sqrt(apse_vals[3, ]), xlab = "Degree", 
      ylab = "", col = "red", lwd = 2)
lines(degrees, sqrt(apse_vals[4, ]), xlab = "Degree", 
      ylab = "", col = "black", lwd = 2)
legend("topright", legend = c("APSE", "Var", "Bias2", "Var_y"),
       col = c("purple","blue","red","black"), lwd = 2, bty = "n", cex = 0.8)


# g)

plot(xydata, main = "New COVID Cases in Toronto\n March 1 - December 31, 2020",
     xlab = "Day_Num", ylab = "New_Cases", cex = 0.8, cex.main = 0.8, 
     col = adjustcolor("darkgrey", 0.5), pch = 16)
muhat.best <- getmuhat(xydata, 7)
curve(muhat.best, add = TRUE, lwd = 3, cex = 0.8, bty = "n", col = adjustcolor("purple", 0.5))



