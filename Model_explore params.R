rm(list = ls())
gc()
graphics.off()

quartz()

##http://rstudio-pubs-static.s3.amazonaws.com/159425_80725873417e42fdb13821c10a198281.html
##https://github.com/sokole/MCSim


#Sys.unsetenv("GITHUB_PAT")
#devtools::install_github('sokole/MCSim')
library(MCSim)
library(tidyverse)
library(fields)
library(vegan)
library(betapart)


################
##Explore parameters model##
################

set.seed(1234)

print(Sys.time())
time.ini <- print(Sys.time())

nsites <- 50

# Let's generate a data frame with coordinates for nsites sites, randomly distributed in a 10 x 10 environment
xy.coordinates <- data.frame(
  x = runif(nsites)*10,
  y = runif(nsites)*10
)
xy.coordinates <- xy.coordinates[order(xy.coordinates$x,xy.coordinates$y),]
print(xy.coordinates)
plot(y ~ x, xy.coordinates)

my.landscape <- MCSim::fn.make.landscape(
  site.coords = xy.coordinates, ## site coordinates
  m = rep(0.05,nsites),#rep(0.01,10), ##immigration rate for each site from the species pool
  #Ef = runif(nsites), ## Environmental variable value for each site
  Ef = sort(runif(nsites)), ## Environmental variable value for each site - since the sites are ordered on the x-y axise, using sort introduces an environmental gradient along the x axis
  JM = 1000000 ## metacommunity size (sensu Hubbell 2001) - species pool size
)
print(my.landscape$site.info)

plot(y ~ x, my.landscape$site.coords, pch=16, col=tim.colors(n=64)[ceiling(my.landscape$site.info$Ef*64)])

# disp.ab <- c(0,0.1,0.2,0.5,1,2)
# niche.width <- c(0.01,0.05,0.1,0.2,0.5,1,2,5,10,1000)
disp.ab <- c(0,2,5,10,1000,5e6)
niche.width <- c(0.01,0.1,1,5,1000)

Gamma <- Alpha <- Shannon <- beta.Sim <- beta.Sor <- beta.ab <- matrix(NA,length(niche.width),length(disp.ab))
MC.fin.mat.store <- list()

i.n <- 0
for(n in niche.width){
  
  #n=1000
  
  i.n <- i.n+1
  MC.fin.mat.store[[i.n]] <- list()
  
  
  # niche positions, niche breadths, and relative abundances for 10 species
  niche.positions <-  runif(10)
  niche.breadths <- rep(n,10)
  regional.rel.abund <- rep(1/10,10)
  #regional.rel.abund <- (1:100)/sum(1:100)
  
  
  
  ####Plot niches
  
  # -- function for plotting bell curves
  fn_norm_curve <- function(sigma=1, mu = 0,...) {
    curve(
      (1/sigma * sqrt(2 * pi)) * exp((-1  *(x - mu)^2) / (2 * sigma^2)), ...) #formula for bell curve
  }
  
  # -- Initialize plot of coenoclines
  plot(1,1,
       xlim = c(-0.2,1.2),
       ylim = c(0, (1/niche.breadths[1] * sqrt(2 * pi))),
       type = 'n',
       xlab = 'Environmental gradient',
       ylab = 'Prob. dens.',
       main = 'Niche positions and niche widths')
  
  #mypal <- rainbow(length(niche.positions))
  mypal <- tim.colors(length(niche.positions))
  
  # -- loop to plot each species' habitat preference
  for (i.spp in 1:length(niche.positions)){
    fn_norm_curve(
      mu = niche.positions[i.spp],
      sigma = niche.breadths[i.spp],
      add = TRUE,
      col = mypal[i.spp],
      n=1000)
  }
  
  # -- plot sites along the x-axis
  rug(my.landscape$site.info$Ef)
  
  legend(x = 'topleft',
         legend = 1:10,
         lty = 1,
         col = mypal,
         cex = .75,
         ncol = 4)
  
  
  ####Run simulations
  i.d <- 0
  for(d in disp.ab){
    
    #d=5e6
  
    i.d <- i.d+1
    
    print(paste("niche width =",n, ", niche index =", i.n, ";;; dispersal =", d, ",dispersal index =",i.d))
      
    # run a simulation with n.timestep generations
    n.timestep <- 100
    sim.result <- MCSim::fn.metaSIM(
      landscape = my.landscape,
      trait.Ef = niche.positions,
      trait.Ef.sd = niche.breadths, 
      gamma.abund = regional.rel.abund,
      W.r = d, ##slope of exponential dispersal kernel
      nu = 0,#0.001, ##probability that a novel species will appear during a recruitment event
      n.timestep = n.timestep, 
      sim.ID = "my_test_sim",
      output.dir.path = "my_sim_output_directory" 
    )
    
    print(sim.result$sim.result.name)
    
    head(sim.result$J.long)
    
    MC.fin <- sim.result$J.long[which(sim.result$J.long$timestep==n.timestep),]
    MC.fin.mat <- MC.fin %>% 
      pivot_wider(names_from=spp,values_from=c(count))
    MC.fin.mat <- MC.fin.mat[,-c(1,2)]
    MC.fin.mat.pa <- as.data.frame((MC.fin.mat>0)*1)
    
    #Gamma
    Gamma[i.n,i.d] <- length(which(colSums(MC.fin.mat)>0))
    
    #Alpha
    #rowSums(MC.fin.mat.pa)
    Alpha[i.n,i.d] <- mean(rowSums(MC.fin.mat.pa))
    
    #Shannon
    #diversity(MC.fin.mat)
    Shannon[i.n,i.d] <- mean(diversity(MC.fin.mat))
    
    #Beta
    #beta.pair(MC.fin.mat.pa)
    #beta.pair.abund(MC.fin.mat)$beta.bray
    beta <- lapply(beta.pair(MC.fin.mat.pa),mean)
    beta.Sim[i.n,i.d] <- beta$beta.sim
    beta.Sor[i.n,i.d] <- beta$beta.sor
    beta.ab[i.n,i.d] <- mean(beta.pair.abund(MC.fin.mat)$beta.bray)
    
    MC.fin.mat.store[[i.n]][[i.d]] <- MC.fin.mat
  }
}

print(time.ini)
print(Sys.time())

Gamma
Alpha
Shannon
beta.Sim
beta.Sor
beta.ab

save(MC.fin.mat.store,
     Gamma,
     Alpha,
     Shannon,
     beta.Sim,
     beta.Sor,
     beta.ab,file="Simulation_results.RData")



##Null model

mat.rand <- permatfull(MC.fin.mat.store[[2]][[3]],fixedmar="both")

Gamma.rand <- Alpha.rand <- Shannon.rand <- beta.Sim.rand <- beta.Sor.rand <- beta.ab.rand <- numeric(99)
for(i in 1:99){
  mat.rand.pa <- as.data.frame((mat.rand$perm[[i]]>0)*1)
  
  #Gamma
  Gamma.rand[i] <- length(which(colSums(mat.rand.pa)>0))
  
  #Alpha
  #rowSums(MC.fin.mat.pa)
  Alpha.rand[i] <- mean(rowSums(mat.rand.pa))
  
  #Shannon
  #diversity(MC.fin.mat)
  Shannon.rand[i] <- mean(diversity(mat.rand$perm[[i]]))
  
  #Beta
  #beta.pair(MC.fin.mat.pa)
  #beta.pair.abund(MC.fin.mat)$beta.bray
  beta.rand <- lapply(beta.pair(mat.rand.pa),mean)
  beta.Sim.rand[i] <- beta.rand$beta.sim
  beta.Sor.rand[i] <- beta.rand$beta.sor
  beta.ab.rand[i] <- mean(beta.pair.abund(mat.rand$perm[[i]])$beta.bray)
}

par(mfrow=c(2,3))
hist(Gamma.rand)
abline(v=Gamma[2,3],col="red")
hist(Alpha.rand)
abline(v=Alpha[2,3],col="red")
hist(Shannon.rand)
abline(v=Shannon[2,3],col="red")
hist(beta.Sor.rand)
abline(v=beta.Sor[2,3],col="red")
hist(beta.Sim.rand)
abline(v=beta.Sim[2,3],col="red")
hist(beta.ab.rand)
abline(v=beta.ab[2,3],col="red")
 








