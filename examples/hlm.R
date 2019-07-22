summary(colond) # colon cancer dataset

# fit model without censoring (i.e., ignoring it)
M <- hlm(log(time) ~ . | rx, data = colond)
M # brief display
summary(M) # full summary
vcov(M) # variance estimate

# fit same model but with censoring
hlm(cbind(log(time), status) ~ . | rx, data = colond)
