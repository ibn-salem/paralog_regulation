########################################################################
#
# A script to find and test metrics for exclusive expression of gene pairs
#
########################################################################


x <- c(0, 1, 3, 2, 1, 0, 0, 15, 9, 3, 12, 8)
y <- c(9, 15, 8, 3, 11, 3, 0, 1, 2, 3, 0, 3)

# random
rand <- sapply(seq(100), function(i){
    x = runif(60) * 100
    y = runif(60) * 100
    ee(x,y)
})
summary(rand)
mean(rand)
sd(rand)

x = runif(60) * 10
y = runif(60) * 10
pl(x,y)

# exclusive expressed
x = c(runif(30)*2, runif(30)*10)
y = c(runif(30)*10, runif(30)*2)
pl(x,y)

# exclusive expressed all high
x = c(runif(30)*2, runif(30)*2+10)
y = c(runif(30)*2+10, runif(30)*2)
pl(x,y)


# exclusive expressed but lot of very low expressed
x = c(runif(10)*2, runif(10)*10, runif(40))
y = c(runif(10)*10, runif(10)*2, runif(40))
pl(x,y)


# exclusive expressed but lot of zeros
x = c(runif(10)*2, runif(10)*10, runif(10), rep(0, 30))
y = c(runif(10)*10, runif(10)*2, runif(10), rep(0, 30))
pl(x,y)


# all low
x <- runif(60)
y <- runif(60)
pl(x,y)



ee <- function(x, y){
    
    n = length(x)
    minVal <- mapply(min, x, y)
    maxVal <- mapply(max, x, y)
    
    ee <- (n/2) / sum(minVal / maxVal, na.rm=TRUE)
    return(ee)
}


ee.pc <- function(x, y, ps.count=0.01){
    
    n = length(x)
    minVal <- mapply(min, x, y)
    maxVal <- mapply(max, x, y) + ps.count
    
    ee <- (n/2) / sum(minVal / maxVal, na.rm=TRUE)
    return(ee)
}


pl <- function(x,y){
    plot(x,y, pch=20, cex=2, main=paste("EE =", signif(ee(x,y),3)))
}
