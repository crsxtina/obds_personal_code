vector <- 1:200
vector*123
vector[44]
x <- vector[1:15]
x
mean1 <- sum(x)/length(x)
mean2 <- mean(x)
mean2
median <- x[length(x)%%2]
median

vector2 <- c('actb', 100, 3.4)
vector2
z <- as.numeric(vector2[2])*4
z
vector3 <-  list('actb', 100, 3.4)
g <- vector3[[2]]*4
vector3
g
unlist(vector3)
