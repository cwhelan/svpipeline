sim.nodel <- function(coverage, insert.mean=200, insert.sd=30, max.insert=25000, num.noise=5) {
    num.true.reads <- rpois(1, coverage)
    true.inserts <- rnorm(num.true.reads, insert.mean, insert.sd)
    noise.inserts <- runif(num.noise, insert.mean, max.insert)
    c(true.inserts, noise.inserts)
}

sim.homdel <- function(coverage, del.size, insert.mean=200, insert.sd=30, max.insert=25000, num.noise=5) {
  num.true.reads <- rpois(1, coverage)
  true.inserts <- rnorm(num.true.reads, del.size + insert.mean, insert.sd)
  noise.inserts <- runif(num.noise, insert.mean, max.insert)
  c(true.inserts, noise.inserts)
}

sim.hetdel <- function(coverage, del.size, insert.mean=200, insert.sd=30, max.insert=25000, num.noise=5) {
  num.true.reads <- rpois(1, coverage)
  nodel.inserts <- rnorm(num.true.reads / 2, insert.mean, insert.sd)
  del.inserts <- rnorm(num.true.reads / 2, del.size + insert.mean, insert.sd)
  noise.inserts <- runif(num.noise, insert.mean, max.insert)
  c(nodel.inserts, del.inserts, noise.inserts)
}

likelihood <- function(y, w, mu, sigma) {
  m1.likelihoods <- dnorm(y, mu[1], sigma, log=TRUE)
  m2.likelihoods <- dnorm(y, mu[2], sigma, log=TRUE)
  
  weighted.likelihoods <- exp(w[1]) *  m1.likelihoods + exp(w[2]) * m2.likelihoods
  mean(weighted.likelihoods)
}

logsumexp <- function(x) {
  m <- max(x)
  s <- m + log(sum(exp(x - m)))
  s
}

gamma <- function(y, w, mu, sigma) {
  m1.likelihood <- w[1] + dnorm(y, mu[1], sigma, log=TRUE)
  print("m1")
  print(m1.likelihood)
  print(exp(m1.likelihood))
  m2.likelihood <- w[2] + dnorm(y, mu[2], sigma, log=TRUE)
  print("m2")
  print(m2.likelihood)
  print(exp(m2.likelihood))
  total.likelihood <- apply(cbind(m1.likelihood, m2.likelihood), 1, logsumexp)
  print("total likelihood")
  print(total.likelihood)
  gamma1 <- m1.likelihood - total.likelihood
  gamma2 <- m2.likelihood - total.likelihood
  cbind(gamma1,gamma2)
}

n.calc <- function(gamma.m) {
#  log(apply(exp(gamma.m), 2, sum))
  c(logsumexp(gamma.m[,1]), logsumexp(gamma.m[,2]))
}

w.update <- function(n, y) {
  n - log(length(y))
}

mu2.update <- function(gamma, y) {
   sum(exp(gamma[,2] + log(y))) / sum(exp(gamma[,2]))
}

em.step <- function(y, w, mu, sigma) {
  gamma.m <- gamma(y, w, mu, sigma)
  #print(gamma.m)
  n.m <- n.calc(gamma.m)
  #print(n.m)
  w.1 <- w.update(n.m, y)  
  mu2.1 <- mu2.update(gamma.m, y)
  list(w=w.1, mu2=mu2.1)
}

em.w <- function(y, w, mu, sigma, max.iter=10) {
  #y <- y[NNclean(y, ceiling(length(y)/5))$z == 1]
  #print(y)
  y <- y[mynnclean(y, sigma)]
  #print(mynnclean(y, sigma))
  if(length(y) == 0) {
    return(log(c(1,0)))
  }
  l <- likelihood(y, w, mu, sigma)
  #print(paste("l:", l))
  i <- 1
  repeat {  
    updates <- em.step(y,w,mu,sigma)
    w <- updates$w
    mu <- c(200, updates$mu2)
    l.prime <- likelihood(y, w, mu, sigma)  
    #print(paste("l:", l.prime))
    #print(paste("w:", w))
    #print(paste("mu:", mu))  
    i <- i+1
    if (abs(l.prime - l) < 0.00001 | i > max.iter) {
      break
    }
  l <- l.prime
  }
  print(mu)
  w
}

mynnclean <- function(x, sd) {
  m <- 3
  #print(x)
  d <- as.matrix(dist(x))
  ns = vector(mode="numeric", length=length(x))
  #print(paste("dim d",dim(d)))
  for (i in seq(1,length(x))) {
    ns[i] <- sort(d[i,])[m]
  }
  #print(ns)
  return(ns < 5 * sd)
}

testHomegrown <- function() {
  w <- log(c(.5,.5))
  mu <- c(200, 1000)
  sigma <- 30
  coverage <- 10
  tests <- 100
  noise <- 10
  delsize <- 300
  
  ws = vector(mode="numeric", length=tests)
  for (i in 1:tests) {
    #print(i)
    y <- sim.hetdel(coverage,delsize,num.noise=noise)
    #print(y)
    ws[i] <- em.w(y, w, mu, sigma)[1]
    #print(ws[i])
  }
  print(exp(ws))
  print(sum(exp(ws) > .25 & exp(ws) < .75))
  
  
  ws = vector(mode="numeric", length=tests)
  for (i in 1:tests) {
    #print(i)
    y <- sim.nodel(coverage, num.noise=noise)
    #print(y)
    ws[i] <- em.w(y, w, mu, sigma)[1]
    #print(ws[i])
  }
  print(exp(ws))
  print(sum(exp(ws) > .75))
  
  ws = vector(mode="numeric", length=tests)
  for (i in 1:tests) {
    #print(i)
    y <- sim.homdel(coverage,delsize, num.noise=noise)
    #print(y)
    ws[i] <- em.w(y, w, mu, sigma)[1]
    #print(ws[i])
  }
  print(exp(ws))
  print(sum(exp(ws) < .25))
}
  
# return 0 (no del), 1 (het), or 2 (hom)
call_with_mclust <- function(y, isize=200, isizeSD=30) {
  tryCatch ( {
    mclustout <- Mclust(y, 
                       modelNames=c('E'), 
                       G=1:2, # prior=priorControl(dof=199, scale=(isizeSD^2)*200),                       
                       initialization=list(noise=(NNclean(y, ceiling(length(y)/5))$z == 0)))
    components <- mclustout$G
    if (components == 1) {
      clustmean <- mclustout$parameters[['mean']]
      if (abs(clustmean - isize) < 1.5 * isizeSD) {
        return(0)
      } else {
        return(2)
      }
    } else {
      return(1)
    }
  }, error=function(err) return(NA))
}

# return 0 (no del), 1 (het), or 2 (hom)
call_with_mclust2 <- function(y, isize=200, isizeSD=30) {
  tryCatch ( {
    mclustout <- Mclust(y, 
                        modelNames=c('E'), 
                        G=2,  prior=priorControl(dof=199, scale=(isizeSD^2)*200),                       
                        initialization=list(noise=(NNclean(y, ceiling(length(y)/5))$z == 0)))    
      means <- mclustout$parameters[['mean']]
      if (abs(means[1] - means[2]) < 1.5 * isizeSD) {
        if (abs(mean(means) - isize) < 1.5 * isizeSD) {
          return(0)
        } else {
          return(2)
        }
      } else {
        # todo: check for mising proportion?
        pros <- mclustout$parameters[['pro']][1:2] / sum(mclustout$parameters[['pro']][1:2])
        if (pros[1] < .3 | pros[1] > .6) {
          if (abs(means[order(pros)[1]] - isize) < 1.5 * isizeSD) {
            return(0)
          } else {
            return(1)
          }
        } else {
          return(1)
        }
      }
  }, error=function(err) return(NA))
}

testmclust <- function(samples=300) {
  for (coverage in c(5,8,10,20)) {
    for (noise in c(1,10)) {
      for (dsize in c(75, 100, 200, 500, 1000)) {
        print(paste("Coverage =", coverage, "Number of noise points =", noise, "Deletion size =", dsize))
        truth = c(rep(0,(samples/3)), rep(1,(samples/3)), rep(2,(samples/3)))      
        nodelpred = vector(mode="numeric", length=(samples/3))
        for (i in 1:(samples/3)) {     
          y <- sim.nodel(coverage, num.noise=noise)
          nodelpred[i] <- call_with_mclust2(y)
        }
        hetpred = vector(mode="numeric", length=(samples/3))
        for (i in 1:(samples/3)) {        
          y <- sim.hetdel(coverage, dsize, num.noise=noise)
          hetpred[i] <- call_with_mclust(y)
        }
        hompred = vector(mode="numeric", length=(samples/3))
        for (i in 1:(samples/3)) {        
          y <- sim.homdel(coverage, dsize, num.noise=noise)
          hompred[i] <- call_with_mclust(y)
        }
        preds <- c(nodelpred,hetpred,hompred)
        print(table(data.frame(truth=truth, pred=preds)))
        print(paste("accuracy:", sum(ifelse(is.na(preds), -1, preds) == truth) / samples))
      }
    }
  }
}