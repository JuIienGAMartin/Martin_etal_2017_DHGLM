model{
#priors
## fixed effects: b, bv
for(h in 1:6){ 			   #for the mean part
	b[h] ~ dnorm(0, 0.001) 
	}
for(s in 1:3){ 			   #for the dispersion part
	bv[s] ~ dnorm(0, 0.001) 
	}
##random effects
###year
####for the mean part
sdyearmu ~ dunif(0,10) 		   #prior on standard deviation
tauyearmu <- 1/sdyearmu/sdyearmu   #conversion to a precision
s2yearmu <- sdyearmu * sdyearmu    #variance estimate
####for the dispersion part
sdyearV ~ dunif(0,10) 		   #prior on standard deviation
tauyearV <- 1/sdyearV/sdyearV	   #conversion to a precision
s2yearV <- sdyearV * sdyearV 	   #variance estimate
for (y in 1:n.years){		   #values for each year
	yearmu[y] ~ dnorm(0, tauyearmu)
	yearV[y] ~ dnorm(0, tauyearV)
	}
###permanent environmental effect(indivual identity)
taupe[1:2,1:2] ~ dwish(matpe[,],3) #wishart prior on the precision matrix
s2pe[1:2,1:2] <- inverse(taupe[,]) #conversion to a variance matrix
for (m in 1:n.id){ 		   #values for each individual
	pe[m,1:2] ~ dmnorm(zero[], taupe[,])
	} 
###Additive genetic effects
tauag[1:2,1:2] ~ dwish(matag[,],3) #wishart prior on the precision matrix
s2ag[1:2,1:2] <- inverse(tauag[,]) #conversion to a variance matrix 
####values for each founder (no known parents)
for (k in 1:n.fnd){  			
	ag[k,1:2] ~ dmnorm(matag[], tauag[,])
	}
####individuals with only father unknown
for (u in 1:nfunk){
	#unknown father breeding value		
	agfunk[u,1:2] ~ dmnorm(matag[], tauag[,])
	#mid-parent value for both mean and dispersion
	parag[funk[u],1] <- (ag[mid[funk[u]],1] + agfunk[u,1])/2
	parag[funk[u],2] <- (ag[mid[funk[u]],2] + agfunk[u,2])/2
	#breeding values 
	agi[funk[u],1:2] ~ dmnorm(parag[funk[u],1:2], tauag[,])
	ag[funk[u],1] <- agi[funk[u],1] * wsq[funk[u]]
	ag[funk[u],2] <- agi[funk[u],2] * wsq[funk[u]]
	}
####individuals with only mother unknown
for (uu in 1:n.munk){  
	#unknown mother breeding value	
	agmunk[uu,1:2] ~ dmnorm(matag[], tau.ag[,]) 
	#mid-parent value for both mean and dispersion 
	parag[munk[uu],1] <- (ag[fid[munk[uu]],1] + agmunk[uu,1])/2
	parag[munk[uu],2] <- (ag[fid[munk[uu]],2] + agmunk[uu,2])/2
	#breeding values
	agi[munk[uu],1:2] ~ dmnorm(parag[munk[uu],1:2], tauag[,])
	ag[munk[uu],1] <- agi[munk[uu],1] * wsq[munk[uu]]
	ag[munk[uu],2] <- agi[munk[uu],2] * wsq[munk[uu]]
	}
####individuals with both parents known
for (d in 1:n.des){
	#mid-parent value for both mean and dispersion
	parag[des[d],1] <- (ag[fid[des[d]],1] + ag[mid[des[d]],1])/2
	parag[des[d],2] <- (ag[fid[des[d]],2] + ag[mid[des[d]],2])/2
	#breeding values
	agi[des[d],1:2] ~ dmnorm(parag[des[d],1:2],tau.ag[,])
	ag[des[d],1] <- agi[des[d],1] * wsq[des[d]]
	ag[des[d],2] <- agi[des[d],2] * wsq[des[d]]
	}
#Likelihood
for(i in 1:n){
	#each observation has its own mean and precision
	docility[i] ~ dnorm(muy[i], tauerr[i])
	#express the model in terms of the variance 
	tauerr[i] <- 1/s2y[i]
  ## Model for the mean (Mean part)
	muy[i] <- b[1] + b[2]*trial[i] + b[3]*day[i] + b[4]*time[i] + b[5]*age1[i] + b[6]*age2[i] + yearmu[year[i]] + pe[id[i],1] + ag[id[i],1]
  ## Model for the variance (Dispersion part)
	log(s2y[i]) <- bv[1] + bv[2]*age1[i] + bv[3]*age2[i] + yearV[year[i]] + pe[id[i],2] + ag[id[i],2]
	}
}
