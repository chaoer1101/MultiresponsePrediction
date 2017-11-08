
# test the term correspondence between R package mixOmics and Python sklearn in term of the partial least squares regression method

rm(list = ls())

library(mixOmics)

data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic

X = X[, 1:120]
Y = Y[, 1:5]

linn.pls <- pls(X, Y, ncomp = 2, mode = "regression")
linn.pls

sx=linn.pls$X
sy=linn.pls$Y

vx=linn.pls$variates$X
vy=linn.pls$variates$Y

# To test if the mixOmics loadings vectors are actually the weights (W and C in sklearn pls regression)

lx=linn.pls$loadings$X
ly=linn.pls$loadings$Y

w1 = lx[,1]
c1 = ly[,1]


# for the first component, check if the scores derived from two ways are equal
# sx and sy correspond to the X0 and Y0
t1=sx %*% w1
u1=sy %*% c1

# make sure the derived scores correspond to the mixOmics variates
all.equal(as.vector(vx[,1]), as.vector(t1))
all.equal(as.vector(vy[,1]), as.vector(u1))

# calculate the true loadings, reference the gender prediction paper 

p1 = t(sx) %*% t1
b1 = t(t1) %*% u1

# Partial Least Squares (PLS) methods for neuroimaging: A tutorial and review
# check if this equation is not correct
x1hat = t1 %*% t(p1)

y1hat = t1 %*% b1 %*% t(c1)
x1 = sx - x1hat
y1 = sy - y1hat


# check for the second component

w2 = lx[,2]
c2 = ly[,2]

t2=x1 %*% w2
u2=y1 %*% c2

# make sure the derived scores correspond to the mixOmics variates
all.equal(as.vector(vx[,2]), as.vector(t2))
all.equal(as.vector(vy[,2]), as.vector(u2))

pred = predict(linn.pls, X)





