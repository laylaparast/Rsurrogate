delta.estimate= function(yone,yzero, var = FALSE, conf.int = FALSE, weight = NULL, weight.perturb=NULL) {
	if(is.null(weight)){weight = rep(1,length(yone)+length(yzero))}
	delta = sum(weight[(1:length(yone))]*yone)/sum(weight[(1:length(yone))]) - sum(weight[(length(yone)+1):(length(yone)+length(yzero))]*yzero)/(sum(weight[(length(yone)+1):(length(yone)+length(yzero))]))
	if(var | conf.int) { 
		if(is.null(weight.perturb)) {weight.perturb = matrix(rexp(500*(length(yone)+length(yzero)), rate=1), ncol = 500)}
		delta.p.vec = apply(weight.perturb, 2, function(x) sum(x[(1:length(yone))]*yone)/sum(x[(1:length(yone))]) - sum(x[(length(yone)+1):(length(yone)+length(yzero))]*yzero)/(sum(x[(length(yone)+1):(length(yone)+length(yzero))])) )
		}
	if(conf.int){
		conf.l.normal = delta - 1.96*sd(delta.p.vec)
		conf.u.normal = delta + 1.96*sd(delta.p.vec)
		conf.l.quantile = quantile(delta.p.vec, 0.025)
		conf.u.quantile = quantile(delta.p.vec, 0.975)
	}	
	if(!var & !conf.int) {return(list("delta" = delta))}
	if(var & !conf.int) {return(list("delta" = delta, "delta.var" = var(delta.p.vec)))}
	if(conf.int) {return(list("delta" = delta, "delta.var" = var(delta.p.vec), "conf.int.normal" = c(conf.l.normal, conf.u.normal), "conf.int.quantile" = as.vector(c(conf.l.quantile, conf.u.quantile))))}
}


delta.s.estimate = function(sone,szero,yone,yzero, weight.perturb = NULL, number="single",  type="robust", warn.te = FALSE, warn.support = FALSE, extrapolate = FALSE, transform = FALSE) {
	if(is.null(number)) {number = "single"}
	if(is.null(type)) {type = "robust"}
	if(substr(number,1,1) == "s") {number = "single"}
	if(substr(number,1,1) == "m") {number = "multiple"}
	if(substr(type,1,1) == "r") {type = "robust"}
	if(substr(type,1,1) == "m") {type = "model"}
	if(!(substr(type, 1,1) %in% c("r", "m"))) {print("Invalid type, default `robust' being used"); type = "robust"}
	if(!(substr(number, 1,1) %in% c("s", "m"))) {print("Invalid number, default `single' being used"); number = "single"}
	if(number == "multiple") {
		if(dim(as.matrix(sone))[2] <2) {print("Invalid number selection, there does not appear to be multiple surrogate markers; default `single' being used"); number = "single" }
	}
	if(number == "single" & dim(as.matrix(sone))[2] > 1) {print("Single setting being used but dimension of surrogate is greater than one; using first surrogate"); sone = as.matrix(sone)[,1]; szero = as.matrix(szero)[,1] }
	p.test = wilcox.test(x=yone, y=yzero, alternative = "two.sided")$p.value
	if(p.test > 0.05 & is.null(weight.perturb) & !warn.te) {print("Warning: it looks like the treatment effect is not significant; may be difficult to interpret the residual treatment effect in this setting")}
	if(number == "single" & type == "robust") {
	range.1 = range(sone)
	range.0 = range(szero)
	range.ind = (range.1[1] > range.0[1]) | (range.1[2] < range.0[2])
	if( range.ind & is.null(weight.perturb) & !warn.support & !extrapolate & !transform) {
		print("Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation")
	}
	}
	if(is.null(weight.perturb)) {weight = rep(1,length(yone)+length(yzero))}
	if(!is.null(weight.perturb)) {weight = weight.perturb}	
	if(type == "model" & number == "single") {
		reg.y = lm(c(yone,yzero)~c(sone,szero) + c(rep(1, length(sone)), rep(0,length(szero))) + c(rep(1, length(sone)), rep(0,length(szero)))*c(sone,szero), weights = weight)
		beta0 = reg.y$coeff[1]
		beta1 = reg.y$coeff[2]
		beta2 = reg.y$coeff[3]
		beta3 = reg.y$coeff[4]
		reg.s = lm(c(sone,szero) ~ c(rep(1, length(sone)), rep(0,length(szero))), weights = weight)
		alpha0 = reg.s$coeff[1]
		alpha1 = reg.s$coeff[2]
		m = beta0 + (beta1+beta3)*szero + beta2 
		delta.s = sum(weight[(length(yone)+1):(length(yone)+length(yzero))]*(m-yzero))/sum(weight[(length(yone)+1):(length(yone)+length(yzero))])
	}
	if(type == "robust"  & number == "single") {
		h.select = bw.nrd(sone)*(length(sone)^(-0.25))
		if(transform){	 	
    	mean.o= mean(c(szero, sone))
  		sd.o = sd(c(szero, sone))
    	s0.new = pnorm((szero- mean.o)/sd.o)
    	s1.new = pnorm((sone - mean.o)/sd.o)
		} 
		if(!transform){
			s0.new = szero
			s1.new = sone
		}
		mu.1.s0 = sapply(s0.new,pred.smooth,zz=s1.new, bw=h.select, y1=yone, weight = weight[(1:length(yone))])
  		if(sum(is.na(mu.1.s0))>0 & extrapolate){
  			print(paste("Note: ", sum(is.na(mu.1.s0)), " values extrapolated."))
    		c.mat = cbind(szero, mu.1.s0)
    		for(o in 1:length(mu.1.s0)) {
    			if(is.na(mu.1.s0[o])){
    				distance = abs(s0.new - s0.new[o])
    				c.temp = cbind(c.mat, distance)
    				c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where mean is not na
    				new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
    				mu.1.s0[o] = new.est[1]   #in case there are multiple matches
    				}
    		}
		}
		delta.s = sum(weight[(length(yone)+1):(length(yone)+length(yzero))]*(mu.1.s0-yzero))/sum(weight[(length(yone)+1):(length(yone)+length(yzero))])	
	}
	if(type == "model"  & number == "multiple") {
		working = as.matrix(lm(yone~sone, weights = weight[(1:length(yone))] )$coeff)
		s.tilde.0 = cbind(rep(1,length(szero[,1])),szero)%*%working
		s.tilde.1 = cbind(rep(1,length(sone[,1])),sone)%*%working
		delta.s = sum(weight[(length(yone)+1):(length(yone)+length(yzero))]*(s.tilde.0-yzero))/sum(weight[(length(yone)+1):(length(yone)+length(yzero))])
	}
	if(type == "robust" & number == "multiple") {
		working = as.matrix(lm(yone~sone, weights = weight[c(1:length(yone))])$coeff)
		s.tilde.0 = cbind(rep(1,length(szero[,1])),szero)%*%working
		s.tilde.1 = cbind(rep(1,length(sone[,1])),sone)%*%working
		if(transform){	 	
    	mean.o= mean(c(s.tilde.0, s.tilde.1))
  		sd.o = sd(c(s.tilde.0, s.tilde.1))
    	s0.new = pnorm((szero- mean.o)/sd.o)
    	s1.new = pnorm((sone - mean.o)/sd.o)
		} 
		if(!transform){
			s0.new = s.tilde.0
			s1.new = s.tilde.1
		}
		mu.1.s0 = sapply(s0.new,pred.smooth,zz=s1.new, y1=yone, weight = weight[(1:length(yone))])
		if(sum(is.na(mu.1.s0))>0 & extrapolate){
			print(paste("Note: ", sum(is.na(mu.1.s0)), " values extrapolated."))
    		c.mat = cbind(s0.new, mu.1.s0)
    		for(o in 1:length(mu.1.s0)) {
    			if(is.na(mu.1.s0[o])){
    				distance = abs(s0.new - s0.new[o])
    				c.temp = cbind(c.mat, distance)
    				c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where mean is not na
    				new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
    				mu.1.s0[o] = new.est[1]   #in case there are multiple matches
    				}
    		}
		}
		delta.s = sum(weight[(length(yone)+1):(length(yone)+length(yzero))]*(mu.1.s0-yzero))/sum(weight[(length(yone)+1):(length(yone)+length(yzero))])
	}
	return(delta.s)
}

R.s.estimate = function(sone,szero,yone,yzero, var = FALSE, conf.int = FALSE, weight.perturb = NULL,number = "single", type = "robust", extrapolate = FALSE, transform = FALSE) {
	warn.te = FALSE
	warn.support = FALSE
	if(substr(number,1,1) == "s") {number = "single"}
	if(substr(number,1,1) == "m") {number = "multiple"}
	if(substr(type,1,1) == "r") {type = "robust"}
	if(substr(type,1,1) == "m") {type = "model"}
	if(substr(type,1,1) == "f") {type = "freedman"}
	if(!(substr(type, 1,1) %in% c("r", "m", "f"))) {print("Warning: Invalid type, default `robust' being used"); type = "robust"}
	if(!(substr(number, 1,1) %in% c("s", "m"))) {print("Warning: Invalid number, default `single' being used"); number = "single"}
	if(number == "multiple") {
		if(dim(as.matrix(sone))[2] <2) {print("Invalid number selection, there does not appear to be multiple surrogate markers; default `single' being used"); number = "single" }
	}
	if(number == "single" & dim(as.matrix(sone))[2] > 1) {print("Single setting being used but dimension of surrogate is greater than one; using first surrogate"); sone = as.matrix(sone)[,1]; szero = as.matrix(szero)[,1] }
	p.test = wilcox.test(x=yone, y=yzero, alternative = "two.sided")$p.value
	if(p.test > 0.05) {print("Warning: it looks like the treatment effect is not significant; may be difficult to interpret the proportion of treatment effect explained in this setting")
		warn.te = TRUE}
	if(mean(yone)-mean(yzero) < 0) {print("Warning: it looks like you need to switch the treatment groups")}
	if(number == "single" & type == "robust") {
	range.1 = range(sone)
	range.0 = range(szero)
	range.ind = (range.1[1] > range.0[1]) | (range.1[2] < range.0[2])
	if(range.ind & !extrapolate & !transform) {
		print("Warning: observed supports do not appear equal, may need to consider a transformation or extrapolation")
		warn.support = TRUE
	}
	}
	if(type != "freedman"){
		delta = delta.estimate(yone = yone, yzero=yzero)$delta	
		delta.s = delta.s.estimate(sone=sone, szero=szero, yone=yone, yzero=yzero, number = number, type = type, warn.te = warn.te, warn.support = warn.support, extrapolate = extrapolate, transform = transform)
		R.s = 1-delta.s/delta	
	}
	if(type == "freedman"  & number == "single") {
		reg.y = lm(c(yone,yzero)~c(sone, szero) + c(rep(1, length(sone)),rep(0,length(szero))))
		beta.treat.s = reg.y$coeff[3]
		reg.y = lm(c(yone,yzero)~c(rep(1, length(sone)),rep(0,length(szero))))
		beta.treat.nos = reg.y$coeff[2]
		R.s = 1-beta.treat.s/beta.treat.nos
	}
	if(type == "freedman"  & number == "multiple") {
		treat.ind = c(rep(1,length(yone)), rep(0, length(yzero)))
		beta.treat.s = summary(lm(c(yone, yzero)~cbind(treat.ind, rbind(sone, szero))))$coeff[2]
		beta.treat.nos = summary(lm(c(yone, yzero)~treat.ind))$coeff[2]
		R.s = 1-beta.treat.s/beta.treat.nos
	}
	if(var | conf.int){
		if(is.null(weight.perturb)) {
	weight.perturb = matrix(rexp(500*(length(yone)+length(yzero)), rate=1), ncol = 500)}
		if(type != "freedman"){
			delta.s.p.vec = apply(weight.perturb, 2, delta.s.estimate, sone=sone, szero=szero, yone = yone, yzero = yzero,number = number, type = type, extrapolate = extrapolate, transform = transform)
			delta.p.vec = as.vector(unlist(apply(weight.perturb, 2, delta.estimate, yone = yone, yzero=yzero, var= FALSE, conf.int = FALSE)))
			R.p = 1-delta.s.p.vec/delta.p.vec
		}
		if(type == "freedman"  & number == "single") {
			beta.treat.s.p<- apply( weight.perturb, 2 , function(x) (lm(c(yone,yzero)~c(sone, szero) + c(rep(1, length(sone)),rep(0,length(szero))), weights = x) )$coeff[3])
			beta.treat.nos.p = apply( weight.perturb, 2 , function(x) (lm(c(yone,yzero)~ c(rep(1, length(sone)),rep(0,length(szero))), weights = x) )$coeff[2])
			R.p = 1-beta.treat.s.p/beta.treat.nos.p
		}
		if(type == "freedman"  & number == "multiple") {
		beta.treat.s.p<- apply(weight.perturb, 2 , function(x) (lm(c(yone,yzero)~cbind(treat.ind, rbind(sone, szero)), weights = x))$coeff[2])
		beta.treat.nos.p = apply( weight.perturb, 2 , function(x) (lm(c(yone, yzero)~treat.ind,weights = x))$coeff[2])
			R.p = 1-beta.treat.s.p/beta.treat.nos.p
		}
	if(conf.int & type != "freedman")	{
		conf.l.normal.delta = delta - 1.96*sd(delta.p.vec)
		conf.u.normal.delta = delta + 1.96*sd(delta.p.vec)
		conf.l.quantile.delta = quantile(delta.p.vec, 0.025)
		conf.u.quantile.delta = quantile(delta.p.vec, 0.975)
		
		conf.l.normal.delta.s = delta.s - 1.96*sd(delta.s.p.vec)
		conf.u.normal.delta.s = delta.s + 1.96*sd(delta.s.p.vec)
		conf.l.quantile.delta.s = quantile(delta.s.p.vec, 0.025)
		conf.u.quantile.delta.s = quantile(delta.s.p.vec, 0.975)
	}
	if(conf.int) {	
		conf.l.normal.R.s = R.s - 1.96*sd(R.p)
		conf.u.normal.R.s = R.s + 1.96*sd(R.p)
		conf.l.quantile.R.s = quantile(R.p, 0.025)
		conf.u.quantile.R.s = quantile(R.p, 0.975)
		if(type!= "freedman") {fieller.ci.calc = as.vector(fieller.ci(delta.s.p.vec, delta.p.vec, delta.s , delta))}
		if(type== "freedman") {fieller.ci.calc=as.vector(fieller.ci(beta.treat.s.p, beta.treat.nos.p, beta.treat.s, beta.treat.nos))}
	}	
	}
	if(type != "freedman") {
		if(!var & !conf.int) {return(list("delta" = delta, "delta.s" =delta.s, "R.s" = R.s))}
		if(var & !conf.int) {return(list("delta" = delta, "delta.s" =delta.s, "R.s" = R.s, "delta.var" = var(delta.p.vec), "delta.s.var" = var(delta.s.p.vec), "R.s.var" = var(R.p)))}
		if(conf.int) {return(list("delta" = delta, "delta.s" =delta.s, "R.s" = R.s, "delta.var" = var(delta.p.vec), "delta.s.var" = var(delta.s.p.vec), "R.s.var" = var(R.p), "conf.int.normal.delta" = c(conf.l.normal.delta, conf.u.normal.delta), "conf.int.quantile.delta" = as.vector(c(conf.l.quantile.delta, conf.u.quantile.delta)), "conf.int.normal.delta.s" = c(conf.l.normal.delta.s, conf.u.normal.delta.s), "conf.int.quantile.delta.s" = as.vector(c(conf.l.quantile.delta.s, conf.u.quantile.delta.s)),
"conf.int.normal.R.s" = c(conf.l.normal.R.s, conf.u.normal.R.s), "conf.int.quantile.R.s" = as.vector(c(conf.l.quantile.R.s, conf.u.quantile.R.s)), "conf.int.fieller.R.s" = fieller.ci.calc))}	
	}
	if(type == "freedman") {
		if(!var & !conf.int) {return(list("R.s" = as.vector(R.s)))}
		if(var & !conf.int) {return(list("R.s" = as.vector(R.s), "R.s.var" = var(R.p)))}
		if(conf.int) {return(list("R.s" = as.vector(R.s),  "R.s.var" = as.vector(var(R.p)), 
"conf.int.normal.R.s" = as.vector(c(conf.l.normal.R.s, conf.u.normal.R.s)), "conf.int.quantile.R.s" = as.vector(as.vector(c(conf.l.quantile.R.s, conf.u.quantile.R.s))), "conf.int.fieller.R.s" = fieller.ci.calc))}
	}
}

fieller.ci = function(perturb.delta.s, perturb.delta, delta.s, delta) {
	num = (perturb.delta.s-(delta.s/delta)*perturb.delta)^2
	sigma11 = var(perturb.delta.s)
	sigma22 = var(perturb.delta)
	sigma12 = cov(perturb.delta.s, perturb.delta)
	den = sigma11 - 2*(delta.s/delta)*sigma12 + (delta.s/delta)^2*sigma22
	c.alpha = quantile(num/den,probs=0.95)
	#quadratic
	a = delta^2-sigma22*c.alpha
	b = -2*delta.s*delta + c.alpha*sigma12*2
	c = delta.s^2 - c.alpha*sigma11
	root.1 = (-b+sqrt(b^2 - 4*a*c))/(2*a)
	root.2 = (-b-sqrt(b^2 - 4*a*c))/(2*a)
	return(c(1-root.1,1-root.2))
	
}

VTM<-function(vc, dm){
     matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
    }
    
 
Kern.FUN <- function(zz,zi,bw) 
  { 
    out = (VTM(zz,length(zi))- zi)/bw
	dnorm(out)/bw
           
  }
  

pred.smooth <-function(zz,zi.one, bw=NULL,y1, weight = NULL) { 	
  	if(is.null(bw)) { bw = bw.nrd(zz)/((length(zz))^(0.25))}
  	if(is.null(weight)) {weight = rep(1,length(y1))}
    return(sum(weight*Kern.FUN(zz,zi.one,bw=bw)*y1)/sum(weight*Kern.FUN(zz,zi.one,bw=bw)))
 
  }

 delta.surv.estimate= function(xone,xzero, deltaone, deltazero, t, var = FALSE, conf.int = FALSE, weight = NULL, weight.perturb = NULL) {
	if(is.null(weight)){weight = rep(1,length(xone)+length(xzero))}
	censor1.t = censor.weight(xone, deltaone, t, weight = weight[1:length(xone)])
	censor0.t = censor.weight(xzero, deltazero, t, weight = weight[(1+length(xone)):(length(xone)+length(xzero))])
	delta = sum(1*(xone > t)*weight[1:length(xone)])/sum(weight[1:length(xone)]*censor1.t) - sum(1*(xzero > t)*weight[(1+length(xone)):(length(xone)+length(xzero))])/sum(weight[(1+length(xone)):(length(xone)+length(xzero))]* censor0.t)
	if(var | conf.int) { 
		if(is.null(weight.perturb)) {weight.perturb = matrix(rexp(500*(length(xone)+length(xzero)), rate=1), ncol = 500)}
		delta.p.vec = apply(weight.perturb, 2, function(x) {  	censor1.t.p = censor.weight(xone, deltaone, t, weight = x[1:length(xone)]); censor0.t.p = censor.weight(xzero, deltazero, t, weight = x[(1+length(xone)):(length(xone)+length(xzero))]);
sum(1*(xone > t)*x[1:length(xone)])/sum(x[1:length(xone)]*censor1.t) - sum(1*(xzero > t)*x[(1+length(xone)):(length(xone)+length(xzero))])/sum(x[(1+length(xone)):(length(xone)+length(xzero))]* censor0.t)})
		}
	if(conf.int){
		conf.l.normal = delta - 1.96*sd(delta.p.vec)
		conf.u.normal = delta + 1.96*sd(delta.p.vec)
		conf.l.quantile = quantile(delta.p.vec, 0.025)
		conf.u.quantile = quantile(delta.p.vec, 0.975)
	}	
	if(!var & !conf.int) {return(list("delta" = delta))}
	if(var & !conf.int) {return(list("delta" = delta, "delta.var" = var(delta.p.vec)))}
	if(conf.int) {return(list("delta" = delta, "delta.var" = var(delta.p.vec), "conf.int.normal" = c(conf.l.normal, conf.u.normal), "conf.int.quantile" = as.vector(c(conf.l.quantile, conf.u.quantile))))}
}
 
censor.weight = function(data.x, data.delta, t, weight=NULL) {
	if(is.null(weight)) {weight = rep(1,length(data.x))}
	S.KM = survfit(Surv(data.x,1-data.delta)~1, weights = weight)
	S.t.KM = approx(S.KM$time,S.KM$surv,t)$y
	return(S.t.KM)
}



delta.s.surv.estimate = function(xone,xzero, deltaone, deltazero, sone, szero, t, weight.perturb = NULL, landmark, extrapolate = FALSE, transform = FALSE) {
	
	delta.check = delta.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, conf.int = TRUE)
	if(delta.check$conf.int.quantile[1] < 0 & delta.check$conf.int.quantile[2] > 0 & is.null(weight.perturb)) {print("Warning: it looks like the treatment effect is not significant; may be difficult to interpret the residual treatment effect in this setting")}
	if(delta.check$delta < 0 & is.null(weight.perturb)) {print("Warning: it looks like you need to switch the treatment groups")}
	range.1 = range(sone, na.rm = T)
	range.0 = range(szero, na.rm = T)
	range.ind = (range.1[1] > range.0[1]) & (range.1[2] < range.0[2])
	if( range.ind & is.null(weight.perturb)) {
		print("Warning: observed supports to not appear equal, may need to consider a transformation or extrapolation")
	}
	if(is.null(weight.perturb)) {weight = rep(1,length(xone)+length(xzero))}
	if(!is.null(weight.perturb)) {weight = weight.perturb}	
	weight.group1 = weight[1:length(xone)]
	weight.group0 = weight[(1+length(xone)):(length(xone)+length(xzero))]
	mu.1.s0 = pred.smooth.surv(xone.f=xone[xone>landmark], deltaone.f = deltaone[xone>landmark], sone.f=sone[xone>landmark], szero.one = szero[xzero>landmark], myt=t, weight.pred = weight.group1[xone>landmark], extrapolate = extrapolate, transform = transform)
	censor0.t = censor.weight(xzero, deltazero, t, weight = weight.group0)
	censor0.landmark = censor.weight(xzero, deltazero, landmark, weight = weight.group0)
	delta.s = sum(weight.group0[xzero>landmark]*mu.1.s0)/sum(weight.group0*censor0.landmark) - sum(weight.group0*1*(xzero>t))/sum(weight.group0*censor0.t)
	return(delta.s)
}


  

R.s.surv.estimate = function(xone,xzero, deltaone, deltazero, sone, szero, t, weight.perturb = NULL, landmark, extrapolate = FALSE, transform = FALSE, conf.int =FALSE, var = FALSE, incremental.value = FALSE) {
	delta = delta.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t)$delta	
	delta.s = delta.s.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, sone=sone, szero=szero, t=t, landmark=landmark, extrapolate = extrapolate, transform = transform)
	R.s = 1-delta.s/delta	
	
	if(var | conf.int){
		if(is.null(weight.perturb)) {
			
	weight.perturb = matrix(rexp(500*(length(xone)+length(xzero)), rate=1), ncol = 500)}
		delta.s.p.vec = apply(weight.perturb, 2, delta.s.surv.estimate, xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, sone=sone, szero=szero, t=t, landmark=landmark, extrapolate = extrapolate, transform = transform)
		print(weight.perturb[1:10, 1:10])
		delta.p.vec = as.vector(unlist(apply(weight.perturb, 2, delta.surv.estimate, xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, var= FALSE, conf.int = FALSE)))
		R.p = 1-delta.s.p.vec/delta.p.vec
		if(conf.int)	{
		conf.l.normal.delta = delta - 1.96*sd(delta.p.vec)
		conf.u.normal.delta = delta + 1.96*sd(delta.p.vec)
		conf.l.quantile.delta = quantile(delta.p.vec, 0.025)
		conf.u.quantile.delta = quantile(delta.p.vec, 0.975)
		
		conf.l.normal.delta.s = delta.s - 1.96*sd(delta.s.p.vec)
		conf.u.normal.delta.s = delta.s + 1.96*sd(delta.s.p.vec)
		conf.l.quantile.delta.s = quantile(delta.s.p.vec, 0.025)
		conf.u.quantile.delta.s = quantile(delta.s.p.vec, 0.975)

		conf.l.normal.R.s = R.s - 1.96*sd(R.p)
		conf.u.normal.R.s = R.s + 1.96*sd(R.p)
		conf.l.quantile.R.s = quantile(R.p, 0.025)
		conf.u.quantile.R.s = quantile(R.p, 0.975)
		
		fieller.ci.calc = as.vector(fieller.ci(delta.s.p.vec, delta.p.vec, delta.s , delta))
	}
	}
	if(incremental.value) {
		R.t.est = R.t.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, landmark = landmark,var = var, weight.perturb = weight.perturb)
		delta.t = R.t.est$delta.t
		R.t = R.t.est$R.t
		IV = R.s - R.t
		if(var | conf.int){
		delta.t.p.vec = apply(weight.perturb, 2, delta.t.surv.estimate, xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, landmark=landmark)
		R.p.t = 1-delta.t.p.vec/delta.p.vec
		IV.p = R.p - R.p.t
		if(conf.int)	{		
		conf.l.normal.delta.t = delta.t - 1.96*sd(delta.t.p.vec)
		conf.u.normal.delta.t = delta.t + 1.96*sd(delta.t.p.vec)
		conf.l.quantile.delta.t = quantile(delta.t.p.vec, 0.025)
		conf.u.quantile.delta.t = quantile(delta.t.p.vec, 0.975)

		conf.l.normal.R.t = R.t - 1.96*sd(R.p.t)
		conf.u.normal.R.t = R.t + 1.96*sd(R.p.t)
		conf.l.quantile.R.t = quantile(R.p.t, 0.025)
		conf.u.quantile.R.t = quantile(R.p.t, 0.975)
		
		conf.l.normal.iv = IV - 1.96*sd(IV.p)
		conf.u.normal.iv = IV + 1.96*sd(IV.p)
		conf.l.quantile.iv = quantile(IV.p, 0.025)
		conf.u.quantile.iv = quantile(IV.p, 0.975)
		
		fieller.ci.calc.t = as.vector(fieller.ci(delta.t.p.vec, delta.p.vec, delta.t , delta))
	}
	}

	}
	r.list = list("delta" = delta, "delta.s" =delta.s, "R.s" = R.s)
	if(var & !conf.int) {r.list = c(r.list, list("delta.var" = var(delta.p.vec), "delta.s.var" = var(delta.s.p.vec), "R.s.var" = var(R.p)))}
	if(conf.int) {r.list = c(r.list, list("delta.var" = var(delta.p.vec), "delta.s.var" = var(delta.s.p.vec), "R.s.var" = var(R.p), "conf.int.normal.delta" = c(conf.l.normal.delta, conf.u.normal.delta), "conf.int.quantile.delta" = as.vector(c(conf.l.quantile.delta, conf.u.quantile.delta)), "conf.int.normal.delta.s" = c(conf.l.normal.delta.s, conf.u.normal.delta.s), "conf.int.quantile.delta.s" = as.vector(c(conf.l.quantile.delta.s, conf.u.quantile.delta.s)),
"conf.int.normal.R.s" = c(conf.l.normal.R.s, conf.u.normal.R.s), "conf.int.quantile.R.s" = as.vector(c(conf.l.quantile.R.s, conf.u.quantile.R.s)), "conf.int.fieller.R.s" = fieller.ci.calc))}
	if(incremental.value) {
		r.list = c(r.list, list("delta.t" =delta.t, "R.t" = R.t, "incremental.value" = IV ))
		if(var & !conf.int) {r.list = c(r.list, list("delta.t.var" = var(delta.t.p.vec), "R.t.var" = var(R.p), "incremental.value.var" = var(IV.p)))}
		if(conf.int) { r.list = c(r.list, list("delta.t.var" = var(delta.t.p.vec), "R.t.var" = var(R.p), "incremental.value.var" = var(IV.p),  "conf.int.normal.delta.t" = c(conf.l.normal.delta.t, conf.u.normal.delta.t), "conf.int.quantile.delta.t" = as.vector(c(conf.l.quantile.delta.t, conf.u.quantile.delta.t)),
"conf.int.normal.R.t" = c(conf.l.normal.R.t, conf.u.normal.R.t), "conf.int.quantile.R.t" = as.vector(c(conf.l.quantile.R.t, conf.u.quantile.R.t)), "conf.int.fieller.R.t" = fieller.ci.calc.t, "conf.int.normal.iv" = c(conf.l.normal.iv, conf.u.normal.iv), "conf.int.quantile.iv" = as.vector(c(conf.l.quantile.iv, conf.u.quantile.iv))))			
		}
		}
	return(r.list)	

}


 delta.t.surv.estimate = function(xone,xzero, deltaone, deltazero, t, weight.perturb = NULL, landmark) {
 		delta.check = delta.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, conf.int = TRUE)
	if(delta.check$conf.int.quantile[1] < 0 & delta.check$conf.int.quantile[2] > 0) {print("Warning: it looks like the treatment effect is not significant; may be difficult to interpret the residual treatment effect in this setting")}
	if(delta.check$delta < 0) {print("Warning: it looks like you need to switch the treatment groups")}
	if(is.null(weight.perturb)) {weight = rep(1,length(xone)+length(xzero))}
	if(!is.null(weight.perturb)) {weight = weight.perturb}	
	weight.group1 = weight[1:length(xone)]
	weight.group0 = weight[(1+length(xone)):(length(xone)+length(xzero))]
	censor0.t = censor.weight(xzero, deltazero, t, weight = weight.group0)
	censor0.landmark = censor.weight(xzero, deltazero, landmark, weight = weight.group0)
	censor1.t = censor.weight(xone, deltaone, t, weight = weight.group1)
	censor1.landmark = censor.weight(xone, deltaone, landmark, weight = weight.group1)
	cond.prop = (sum(weight.group1*1*(xone>t))/sum(weight.group1*censor1.t))/(sum(weight.group1*1*(xone>landmark))/sum(weight.group1*censor1.landmark))
	delta.t = sum(weight.group0*1*(xzero>landmark))/sum(weight.group0*censor0.landmark)*cond.prop - sum(weight.group0*1*(xzero>t))/sum(weight.group0*censor0.t)
	return(delta.t)
}


R.t.surv.estimate = function(xone,xzero, deltaone, deltazero, t, weight.perturb = NULL, landmark, var = FALSE, conf.int = FALSE) {	
	delta = delta.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t)$delta	
	delta.t = delta.t.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero,  t=t, landmark=landmark)
	R.t = 1-delta.t/delta	
	if(var | conf.int){
		if(is.null(weight.perturb)) {
	weight.perturb = matrix(rexp(500*(length(xone)+length(xzero)), rate=1), ncol = 500)}
		delta.t.p.vec = apply(weight.perturb, 2, delta.t.surv.estimate, xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, landmark=landmark)
		delta.p.vec = as.vector(unlist(apply(weight.perturb, 2, delta.surv.estimate, xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, var= FALSE, conf.int = FALSE)))
		R.p = 1-delta.t.p.vec/delta.p.vec
		if(conf.int)	{
		conf.l.normal.delta = delta - 1.96*sd(delta.p.vec)
		conf.u.normal.delta = delta + 1.96*sd(delta.p.vec)
		conf.l.quantile.delta = quantile(delta.p.vec, 0.025)
		conf.u.quantile.delta = quantile(delta.p.vec, 0.975)
		
		conf.l.normal.delta.t = delta.t - 1.96*sd(delta.t.p.vec)
		conf.u.normal.delta.t = delta.t + 1.96*sd(delta.t.p.vec)
		conf.l.quantile.delta.t = quantile(delta.t.p.vec, 0.025)
		conf.u.quantile.delta.t = quantile(delta.t.p.vec, 0.975)

		conf.l.normal.R.t = R.t - 1.96*sd(R.p)
		conf.u.normal.R.t = R.t + 1.96*sd(R.p)
		conf.l.quantile.R.t = quantile(R.p, 0.025)
		conf.u.quantile.R.t = quantile(R.p, 0.975)
		
		fieller.ci.calc = as.vector(fieller.ci(delta.t.p.vec, delta.p.vec, delta.t , delta))
	}
	}
	if(!var & !conf.int) {return(list("delta" = delta, "delta.t" =delta.t, "R.t" = R.t))}
	if(var & !conf.int) {return(list("delta" = delta, "delta.t" =delta.t, "R.t" = R.t, "delta.var" = var(delta.p.vec), "delta.t.var" = var(delta.t.p.vec), "R.t.var" = var(R.p)))}
	if(conf.int) {return(list("delta" = delta, "delta.t" =delta.t, "R.t" = R.t, "delta.var" = var(delta.p.vec), "delta.t.var" = var(delta.t.p.vec), "R.t.var" = var(R.p), "conf.int.normal.delta" = c(conf.l.normal.delta, conf.u.normal.delta), "conf.int.quantile.delta" = as.vector(c(conf.l.quantile.delta, conf.u.quantile.delta)), "conf.int.normal.delta.t" = c(conf.l.normal.delta.t, conf.u.normal.delta.t), "conf.int.quantile.delta.t" = as.vector(c(conf.l.quantile.delta.t, conf.u.quantile.delta.t)),
"conf.int.normal.R.t" = c(conf.l.normal.R.t, conf.u.normal.R.t), "conf.int.quantile.R.t" = as.vector(c(conf.l.quantile.R.t, conf.u.quantile.R.t)), "conf.int.fieller.R.t" = fieller.ci.calc))}	
}




pred.smooth.surv <- function(xone.f, deltaone.f, sone.f, szero.one, myt, weight.pred, extrapolate, transform)
  { 
    if(transform){	 	
    	mean.o= mean(c(sone.f, szero.one))
  		sd.o = sd(c(szero.one, sone.f))
    	sone.f.new = pnorm((sone.f - mean.o)/sd.o)
    	szero.one.new = pnorm((szero.one - mean.o)/sd.o)
		sone.f = sone.f.new
		szero.one = szero.one.new
	}
    bwini = bw.nrd(sone.f)
    n.s = length(sone.f)
    bw <- bwini/(n.s^0.11)
    kerni.ss = Kern.FUN(zz=sone.f,zi=szero.one,bw)           
    tmpind = (xone.f<=myt)&(deltaone.f==1); tj = xone.f[tmpind]; 
    kerni.1 = t(weight.pred*t(kerni.ss))
    pihamyt0.tj.ss = helper.si(tj, "<=", xone.f, Vi=t(kerni.1)) ## n.tj x n.ss matrix ##   
    dLamhat.tj.ss = t((kerni.1[,tmpind]))/pihamyt0.tj.ss; 
    ret = apply(dLamhat.tj.ss,2,sum)
    Phat.ss  =exp(-ret)
    if(sum(is.na(Phat.ss))>0 & extrapolate){
    	print(paste("Note: ", sum(is.na(Phat.ss)), " values extrapolated."))
    	c.mat = cbind(szero.one, Phat.ss)
    	for(o in 1:length(Phat.ss)) {
    		if(is.na(Phat.ss[o])){
    			distance = abs(szero.one - szero.one[o])
    			c.temp = cbind(c.mat, distance)
    			c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where predication is not na
    			new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
    			Phat.ss[o] = new.est[1]   #in case there are multiple matches
    	}
    }
	}
    return(Phat.ss)
    }

cumsum2 <- function(mydat)  
  {
    if(is.null(dim(mydat))) return(cumsum(mydat))
    else{
      out <- matrix(cumsum(mydat), nrow=nrow(mydat))
      out <- out - VTM(c(0, out[nrow(mydat), -ncol(mydat)]), nrow(mydat))
      return(out)
    }
  }

helper.si <- function(yy,FUN,Yi,Vi=NULL)   ## sum I(yy FUN Yi)Vi
  {  
    if(FUN=="<"|FUN=="<=") { yy <- -yy; Yi <- -Yi}
    if(substring(FUN,2,2)=="=") yy <- yy + 1e-8 else yy <- yy - 1e-8
    pos <- rank(c(yy,Yi))[1:length(yy)] - rank(yy)
    if(is.null(Vi)){return(pos)}else{
      Vi <- cumsum2(as.matrix(Vi)[order(Yi),,drop=F])
      out <- matrix(0, nrow=length(yy), ncol=dim(as.matrix(Vi))[2])
      out[pos!=0,] <- Vi[pos,]
      if(is.null(dim(Vi))) out <- c(out)
      return(out) ## n.y x p
    }
  } 


augment.est.vector = function(point.delta, perturb.delta, treat.ind, basis, weights) {
	pi.a = sum(treat.ind==1)/length(treat.ind)
	nu =  t(basis)%*%(1*(treat.ind==1) - pi.a)
	nu.perturb = t(apply(weights, 2, perturb.nu.vector, mat = cbind(treat.ind, basis)))	
	cov.est = cov(perturb.delta, nu.perturb)
	var.est = var(nu.perturb)
	eps = chol2inv(var.est)%*%t(cov(perturb.delta, nu.perturb))
	aug = point.delta - t(eps)%*%nu  
	aug.var = var(perturb.delta)+t(eps)%*%var(nu.perturb)%*%eps - 2*cov.est%*%eps
	aug.perturb = perturb.delta - nu.perturb%*%eps
	return(list("aug.estimate" = aug, "eps" = eps, "nu" = nu, "nu.perturb" = nu.perturb, "aug.var" = aug.var, "aug.perturb" = aug.perturb))
}

perturb.nu.vector = function(mat, weights) {
	#mat should be treat.ind in first column and basis in others
	pi.a = sum(mat[,1]==1)/length(mat[,1])
    nu.p = t(as.matrix(mat[,-1]))%*%(weights*(1*(mat[,1]==1) - pi.a)) 
	return(nu.p)
}


Aug.R.s.surv.estimate = function(xone,xzero, deltaone, deltazero, sone, szero, t, weight.perturb = NULL, landmark, extrapolate = FALSE, transform = FALSE, basis.delta.one, basis.delta.zero,basis.delta.s.one = NULL, basis.delta.s.zero = NULL, incremental.value = FALSE) {
	if(is.null(basis.delta.s.one)) {basis.delta.s.one = basis.delta.one}
	if(is.null(basis.delta.s.zero)) {basis.delta.s.zero = basis.delta.zero}
	delta = delta.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t)$delta	
	delta.s = delta.s.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, sone=sone, szero=szero, t=t, landmark=landmark, extrapolate = extrapolate, transform = transform)
	R.s = 1-delta.s/delta	
 	if(is.null(weight.perturb)) {
	weight.perturb = matrix(rexp(500*(length(xone)+length(xzero)), rate=1), ncol = 500)}
	delta.s.p.vec = apply(weight.perturb, 2, delta.s.surv.estimate, xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, sone=sone, szero=szero, t=t, landmark=landmark, extrapolate = extrapolate, transform = transform)
	delta.p.vec = as.vector(unlist(apply(weight.perturb, 2, delta.surv.estimate, xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, var= FALSE, conf.int = FALSE)))
	R.p = 1-delta.s.p.vec/delta.p.vec
	
	#Augmented estimators
	if(dim(as.matrix(basis.delta.s.one))[2] == 1 ) {basis.delta.s.use =  rbind(cbind(basis.delta.s.one,basis.delta.s.one^2), cbind(basis.delta.s.zero, basis.delta.s.zero^2))}
	if(dim(as.matrix(basis.delta.s.one))[2] > 1 ) {basis.delta.s.use = rbind(basis.delta.s.one, basis.delta.s.zero) }
	if(dim(as.matrix(basis.delta.one))[2] == 1 ) {basis.delta.use =  rbind(cbind(basis.delta.one, basis.delta.one^2), cbind(basis.delta.zero, basis.delta.zero^2))}
	if(dim(as.matrix(basis.delta.one))[2] > 1 ) {basis.delta.use = rbind(basis.delta.one, basis.delta.zero) }
	treat.vector = c(rep(1, length(xone)), rep(0,length(xzero)))
	aug.delta.total = augment.est.vector(point.delta=delta, perturb.delta = delta.p.vec, treat.ind = treat.vector, basis = basis.delta.use, weights = weight.perturb) 
	aug.delta = aug.delta.total$aug.estimate
	sd.aug.delta= sd(aug.delta.total$aug.perturb) 
	conf.normal.aug.delta = c(aug.delta.total$aug.estimate - 1.96*sqrt(aug.delta.total$aug.var) , aug.delta.total$aug.estimate + 1.96*sqrt(aug.delta.total$aug.var) )
	conf.quantile.aug.delta = as.vector(c(quantile(aug.delta.total$aug.perturb, 0.025), quantile(aug.delta.total$aug.perturb, 0.975)))
	
	aug.delta.s.total = augment.est.vector(point.delta=delta.s, perturb.delta = delta.s.p.vec, treat.ind = treat.vector, basis  = basis.delta.s.use, weights = weight.perturb)
	aug.delta.s = aug.delta.s.total$aug.estimate
	sd.aug.delta.s= sd(aug.delta.s.total$aug.perturb)  
	conf.normal.aug.delta.s = c(aug.delta.s.total$aug.estimate - 1.96*sqrt(aug.delta.s.total$aug.var) , aug.delta.s.total$aug.estimate + 1.96*sqrt(aug.delta.s.total$aug.var) )
	conf.quantile.aug.delta.s = as.vector(c(quantile(aug.delta.s.total$aug.perturb, 0.025), quantile(aug.delta.s.total$aug.perturb, 0.975)))

	aug.R = 1-aug.delta.s.total$aug.estimate/aug.delta.total$aug.estimate
	sd.aug.R = sd(1- aug.delta.s.total$aug.perturb/aug.delta.total$aug.perturb)
	conf.normal.aug.R = c(aug.R - 1.96*sd.aug.R, aug.R + 1.96*sd.aug.R)
	conf.quantile.aug.R = as.vector(c(quantile(1- aug.delta.s.total$aug.perturb/aug.delta.total$aug.perturb, 0.025), quantile(1- aug.delta.s.total$aug.perturb/aug.delta.total$aug.perturb, 0.975)))
	conf.fieller.aug.R = fieller.ci(as.vector(aug.delta.s.total$aug.perturb), as.vector(aug.delta.total$aug.perturb), aug.delta.s.total$aug.estimate, aug.delta.total$aug.estimate)

	if(incremental.value) {
		R.t.est = R.t.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, landmark = landmark,var = TRUE, weight.perturb = weight.perturb)
		delta.t = R.t.est$delta.t
		R.t = R.t.est$R.t
		delta.t.p.vec = apply(weight.perturb, 2, delta.t.surv.estimate, xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, landmark=landmark)
		R.p.t = 1-delta.t.p.vec/delta.p.vec
		
		aug.delta.t.total = augment.est.vector(point.delta=delta.t, perturb.delta = delta.t.p.vec, treat.ind = treat.vector, basis = rbind(as.matrix(basis.delta.s.one), as.matrix(basis.delta.s.zero)), weights = weight.perturb)
		aug.delta.t = aug.delta.t.total$aug.estimate
		sd.aug.delta.t= sd(aug.delta.t.total$aug.perturb)  
		conf.normal.aug.delta.t = c(aug.delta.t.total$aug.estimate - 1.96*sqrt(aug.delta.t.total$aug.var) , aug.delta.t.total$aug.estimate + 1.96*sqrt(aug.delta.t.total$aug.var) )
		conf.quantile.aug.delta.t = as.vector(c(quantile(aug.delta.t.total$aug.perturb, 0.025), quantile(aug.delta.t.total$aug.perturb, 0.975)))

		aug.R.t = 1-aug.delta.t.total$aug.estimate/aug.delta.total$aug.estimate
		sd.aug.R.t = sd(1- aug.delta.t.total$aug.perturb/aug.delta.total$aug.perturb)
		conf.normal.aug.R.t = c(aug.R.t - 1.96*sd.aug.R.t, aug.R.t + 1.96*sd.aug.R.t)
		conf.quantile.aug.R.t = as.vector(c(quantile(1- aug.delta.t.total$aug.perturb/aug.delta.total$aug.perturb, 0.025), quantile(1- aug.delta.t.total$aug.perturb/aug.delta.total$aug.perturb, 0.975)))
		conf.fieller.aug.R.t = fieller.ci(as.vector(aug.delta.t.total$aug.perturb), as.vector(aug.delta.total$aug.perturb), aug.delta.t.total$aug.estimate, aug.delta.total$aug.estimate)
		
		aug.IV = aug.R - aug.R.t
		IV.p = (1- aug.delta.s.total$aug.perturb/aug.delta.total$aug.perturb) - (1- aug.delta.t.total$aug.perturb/aug.delta.total$aug.perturb)
				
		conf.l.normal.iv = aug.IV - 1.96*sd(IV.p)
		conf.u.normal.iv = aug.IV + 1.96*sd(IV.p)
		conf.l.quantile.iv = quantile(IV.p, 0.025)
		conf.u.quantile.iv = quantile(IV.p, 0.975)		
	}

	r.list = list("aug.delta" = aug.delta, "aug.delta.s" =aug.delta.s, "aug.R.s" = aug.R, "aug.delta.var" = var(aug.delta.total$aug.perturb), "aug.delta.s.var" = var(aug.delta.s.total$aug.perturb), "aug.R.s.var" = var(1- aug.delta.s.total$aug.perturb/aug.delta.total$aug.perturb), "conf.int.normal.aug.delta" = conf.normal.aug.delta, "conf.int.quantile.aug.delta" = conf.quantile.aug.delta, "conf.int.normal.aug.delta.s" = conf.normal.aug.delta.s, "conf.int.quantile.aug.delta.s" = conf.quantile.aug.delta.s,
"conf.int.normal.aug.R.s" = conf.normal.aug.R, "conf.int.quantile.aug.R.s" = conf.quantile.aug.R, "conf.int.fieller.aug.R.s" = conf.fieller.aug.R)
	if(incremental.value) {
		r.list = c(r.list, list("aug.delta.t" =aug.delta.t, "aug.R.t" = aug.R.t, "aug.incremental.value" = aug.IV,  "aug.delta.t.var" = var(aug.delta.t.total$aug.perturb), "aug.R.t.var" = var(1- aug.delta.t.total$aug.perturb/aug.delta.total$aug.perturb), "aug.incremental.value.var" = var(IV.p), "conf.int.normal.aug.delta.t" = conf.normal.aug.delta.t, "conf.int.quantile.aug.delta.t" = conf.quantile.aug.delta.t,
"conf.int.normal.aug.R.t" = conf.normal.aug.R.t, "conf.int.quantile.aug.R.t" = conf.quantile.aug.R.t, "conf.int.fieller.aug.R.t" = conf.fieller.aug.R.t, "conf.int.normal.aug.iv" = c(conf.l.normal.iv, conf.u.normal.iv), "conf.int.quantile.aug.iv" = as.vector(c(conf.l.quantile.iv, conf.u.quantile.iv))))			
}
	return(r.list)	
}

