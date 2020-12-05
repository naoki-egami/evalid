set.seed(1234)

check <- c()
for(i in 1:1000){
  pv <- sample(x = c(0.001, 0.02, 0.04, 0.05, 0.1, 0.2, 0.003), size = 8, replace = TRUE)
  res <- pct(pv)
  check[i] <- all(pv[res$h_num] == sort(pv))
}
all(check)
