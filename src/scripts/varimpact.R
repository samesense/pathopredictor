library(varimpact)
Q_lib <- c("SL.mean", "SL.glmnet", "SL.ranger", "SL.rpartPrune", "SL.bayesglm")
g_lib <- c("SL.mean", "SL.glmnet")
args <- commandArgs(trailingOnly = TRUE)
xFile = args[1]
yFile = args[2]
outFile = args[3]

x = read.delim(xFile, sep='\t', header=TRUE)
y = as.list(read.delim(yFile, sep='\t', header=TRUE))
print( head(x) )
N = 200
num_normal = 5
X <- as.data.frame(matrix(rnorm(N * num_normal), N, num_normal))
print(head(X))
Y <- rbinom(N, 1, plogis(.2*X[, 1] + .1*X[, 2] - .2*X[, 3] + .1*X[, 3]*X[, 4] - .2*abs(X[, 4])))
vim <- varimpact(Y = y$y, data = x, minCell = 0L, verbose = FALSE) #, Q.library = Q_lib, g.library = g_lib)
#vim <- varimpact(Y = Y, data = X)
write.table(vim$results_all, outFile, row.names=TRUE, quote=FALSE, sep='\t')