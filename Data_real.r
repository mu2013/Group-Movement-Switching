


library(mvtnorm)
library(Matrix)
set.seed(3809)

# Read data

raw <- read.csv("regular.scaled.njaarke.csv")  ## datafile
raw <- raw[-1,-1]

dates <- as.POSIXct(raw[,1], tz = "GMT")  
raw <- cbind(dates, raw[,-1])

# time difference between obs
diff_t <- diff(as.numeric(raw[,1]))/ 3600    # time in mins. Use /3600 to use diff in hours
dt=diff_t[1]

xy <- as.matrix(raw[,-1])
nai <- ncol(xy)/2
Nsamp = 
alen <- nrow(xy)

animaly<-xy[,2*(1:nai)] / 20          # All y coordinates for each individual
animalx<-xy[,2*(1:nai)-1] / 20        # All x coordinates for each individual 

# Coordinates for the black dot as an average of all other individuals. 
black_x <- mean(na.omit(animalx[1,]))
black_y <- mean(na.omit(animaly[1,])) 

# Initial parameter values
theta_x <- mean(animalx[1,][!is.na(animalx[1,])])
theta_y <- mean(animaly[1,][!is.na(animaly[1,])])  
theta <- c(theta_x,theta_y)

swlamda=c(0.1,0.5)

## rename the data for inferencing code 
ax = animalx
ay = animaly
all_t = seq(0,(alen-1)*2 ,dt)
all_s = matrix(c(1), nrow = alen, ncol = nai+1 ) ## nai+1 for animal and black dot
sw_ind = matrix(c("SP"),nrow = alen)

#black = rbind(black_x,black_y)
#mxsamp = max(all_t)

deltat = dt





