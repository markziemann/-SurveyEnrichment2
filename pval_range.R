
# First make a bunch of chisq tests
# n=number of tests to do: default=250
# acceptable test stat range: default=0.2
make_p_range <- function(n=250,bound=0.2){

l <- lapply(0:n, function(x) {
  tulip <- c(((2*n)+x), ((2*n)-x))
  res <- chisq.test(tulip, p = c(1/2, 1/2))
  c("stat"=res$stat,"p"=res$p.value)
})

res <- as.data.frame(do.call(rbind,l))
colnames(res) <- c("stat","p")
res$logp <- -log10(res$p)
res$lower <- res$stat - ( res$stat * bound )
res$upper <- res$stat + ( res$stat * bound )
res$logp_lower <- approx(x = res$stat, y =res$logp, xout = res$lower)[["y"]]
res$logp_upper <- approx(x = res$stat, y =res$logp, xout = res$upper)[["y"]]
res$p_lower <- 10^-res$logp_lower
res$p_upper <- 10^-res$logp_upper
return(res)
}

# show how it works
prange <- make_p_range(n=250,bound=0.2)
head(prange)


# this function will compare published and reproduction p-values,
# estimate the test statistics for both and return the variation
p_compare <- function(p_pub,p_repro,prange=prange) {
  logp_pub = -log10(p_pub)
  logp_repro = -log10(p_repro)
  y2=prange$logp[which(prange$logp > logp_pub)[1]]
  x2=prange$stat[which(prange$logp > logp_pub)[1]]
  y1=prange$logp[rev(which(prange$logp < logp_pub))[1]]
  x1=prange$stat[rev(which(prange$logp < logp_pub))[1]]
  ts_pub <- x1 + ( (logp_pub-y1) / (y2-y1) * (x2-x1) )
  y2=prange$logp[which(prange$logp > logp_repro)[1]]
  x2=prange$stat[which(prange$logp > logp_repro)[1]]
  y1=prange$logp[rev(which(prange$logp < logp_repro))[1]]
  x1=prange$stat[rev(which(prange$logp < logp_repro))[1]]
  ts_repro <- x1 + ( (logp_repro -y1) / (y2-y1) * (x2-x1) )
  v = abs(ts_repro - ts_pub) / ts_pub
  out <- c("test stat pub"=ts_pub,"test stat repro"=ts_repro,"variation"=v)
  return(out)
}

# test it out
p_compare(p_pub=4e-6,p_repro=6e-9,prange=prange)


# make up some data to test
pathways <- paste(rep("GO",20),1:20)
set.seed(100) ; p_pub <-sample(prange$p,20)
set.seed(200) ; p_repro <-sample(prange$p,20)
df <- data.frame(pathways,p_pub,p_repro)
df

# get the values
res <- lapply(1:nrow(df), function(i) {
  p_pub=df[i,"p_pub"]
  p_repro=df[i,"p_repro"]
  p_compare(p_pub=p_pub,p_repro=p_repro,prange=prange)
} )

# tidy up the data and check the 20% threshold
res <- do.call(rbind,res)
res2 <- cbind(df,res)
res2$pass <- res2$variation<0.2
table(res2$pass)
table(res2$pass)/nrow(res2)


