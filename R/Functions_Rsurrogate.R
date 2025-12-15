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
    	h.select = bw.nrd(s1.new)*(length(s1.new)^(-0.25))
		} 
		if(!transform){
			s0.new = szero
			s1.new = sone
		}
		mu.1.s0 = sapply(s0.new,pred.smooth,zz=s1.new, bw=h.select, y1=yone, weight = weight[(1:length(yone))])
  		if(sum(is.na(mu.1.s0))>0 & extrapolate){
  			if(!warn.support){print(paste("Note: ", sum(is.na(mu.1.s0)), " values extrapolated."))}
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
			if(!warn.support){print(paste("Note: ", sum(is.na(mu.1.s0)), " values extrapolated."))}
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

R.s.estimate = function(sone,szero,yone,yzero, var = FALSE, conf.int = FALSE, weight.perturb = NULL,number = "single", type = "robust", extrapolate = FALSE, transform = FALSE, warn.te = FALSE, warn.support = FALSE) {
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
	if(p.test > 0.05 & !warn.te) {print("Warning: it looks like the treatment effect is not significant; may be difficult to interpret the proportion of treatment effect explained in this setting")
		warn.te = TRUE}
	if(number == "single" & type == "robust") {
	range.1 = range(sone)
	range.0 = range(szero)
	range.ind = (range.1[1] > range.0[1]) | (range.1[2] < range.0[2])
	if(range.ind & !extrapolate & !transform & !warn.support) {
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
			delta.s.p.vec = apply(weight.perturb, 2, delta.s.estimate, sone=sone, szero=szero, yone = yone, yzero = yzero,number = number, type = type, extrapolate = extrapolate, transform = transform, warn.te = warn.te, warn.support = warn.support)
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
	num = (as.vector(perturb.delta.s)-as.vector(rep((delta.s/delta), length(perturb.delta)))*as.vector(perturb.delta))^2
	sigma11 = var(perturb.delta.s)
	sigma22 = var(perturb.delta)
	sigma12 = cov(perturb.delta.s, perturb.delta)
	den = sigma11 - 2*(delta.s/delta)*sigma12 + (delta.s/delta)^2*sigma22
	c.alpha = quantile(as.vector(num)/as.vector(den),probs=0.95)
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

  delta.surv.estimate= function(xone,xzero, deltaone, deltazero, t, var = FALSE, conf.int = FALSE, weight = NULL, weight.perturb = NULL, approx=T) {
	if(is.null(weight)){weight = rep(1,length(xone)+length(xzero))}
	censor1.t = censor.weight(xone, deltaone, t, weight = weight[1:length(xone)], approx=approx)
	censor0.t = censor.weight(xzero, deltazero, t, weight = weight[(1+length(xone)):(length(xone)+length(xzero))], approx=approx)
	delta = sum(1*(xone > t)*weight[1:length(xone)])/sum(weight[1:length(xone)]*censor1.t) - sum(1*(xzero > t)*weight[(1+length(xone)):(length(xone)+length(xzero))])/sum(weight[(1+length(xone)):(length(xone)+length(xzero))]* censor0.t)
	if(var | conf.int) { 
		if(is.null(weight.perturb)) {weight.perturb = matrix(rexp(500*(length(xone)+length(xzero)), rate=1), ncol = 500)}
		delta.p.vec = apply(weight.perturb, 2, function(x) {  	censor1.t.p = censor.weight(xone, deltaone, t, weight = x[1:length(xone)],approx=approx); censor0.t.p = censor.weight(xzero, deltazero, t, weight = x[(1+length(xone)):(length(xone)+length(xzero))],approx=approx);
sum(1*(xone > t)*x[1:length(xone)])/sum(x[1:length(xone)]*censor1.t.p) - sum(1*(xzero > t)*x[(1+length(xone)):(length(xone)+length(xzero))])/sum(x[(1+length(xone)):(length(xone)+length(xzero))]* censor0.t.p)})
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
 
censor.weight = function(data.x, data.delta, t, weight=NULL, approx = T) {
	if(is.null(weight)) {weight = rep(1,length(data.x))}
	S.KM = survfit(Surv(data.x,1-data.delta)~1, weights = weight)
	if(approx) {S.t.KM = approx(S.KM$time,S.KM$surv,t, rule = 2)$y}
	if(!approx) {S.t.KM = summary(S.KM, times = t)$surv}
	return(S.t.KM)
}



delta.s.surv.estimate = function(xone,xzero, deltaone, deltazero, sone, szero, t, weight.perturb = NULL, landmark, extrapolate = FALSE, transform = FALSE, approx = T, warn.te = FALSE, warn.support = FALSE) {
	
	delta.check = delta.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, conf.int = TRUE,approx=approx)
	if(delta.check$conf.int.quantile[1] < 0 & delta.check$conf.int.quantile[2] > 0 & is.null(weight.perturb) & !warn.te) {print("Warning: it looks like the treatment effect is not significant; may be difficult to interpret the residual treatment effect in this setting")}
	range.1 = range(sone, na.rm = T)
	range.0 = range(szero, na.rm = T)
	range.ind = (range.1[1] > range.0[1]) & (range.1[2] < range.0[2])
	if( range.ind & is.null(weight.perturb) & !warn.support) {
		print("Warning: observed supports to not appear equal, may need to consider a transformation or extrapolation")
	}
	if(is.null(weight.perturb)) {weight = rep(1,length(xone)+length(xzero))}
	if(!is.null(weight.perturb)) {weight = weight.perturb}	
	weight.group1 = weight[1:length(xone)]
	weight.group0 = weight[(1+length(xone)):(length(xone)+length(xzero))]
	mu.1.s0 = pred.smooth.surv(xone.f=xone[xone>landmark], deltaone.f = deltaone[xone>landmark], sone.f=sone[xone>landmark], szero.one = szero[xzero>landmark], myt=t, weight.pred = weight.group1[xone>landmark], extrapolate = extrapolate, transform = transform, warn.support = warn.support)
	censor0.t = censor.weight(xzero, deltazero, t, weight = weight.group0, approx=approx)
	censor0.landmark = censor.weight(xzero, deltazero, landmark, weight = weight.group0, approx=approx)
	delta.s = sum(weight.group0[xzero>landmark]*mu.1.s0)/sum(weight.group0*censor0.landmark) - sum(weight.group0*1*(xzero>t))/sum(weight.group0*censor0.t)
	return(delta.s)
}


  

R.s.surv.estimate = function(xone,xzero, deltaone, deltazero, sone, szero, t, weight.perturb = NULL, landmark, extrapolate = FALSE, transform = FALSE, conf.int =FALSE, var = FALSE, incremental.value = FALSE, approx = T) {
	delta = delta.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t,approx=approx)$delta	
	delta.s = delta.s.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, sone=sone, szero=szero, t=t, landmark=landmark, extrapolate = extrapolate, transform = transform,approx=approx)
	R.s = 1-delta.s/delta	
	
	if(var | conf.int){
		if(is.null(weight.perturb)) {
			
	weight.perturb = matrix(rexp(500*(length(xone)+length(xzero)), rate=1), ncol = 500)}
		delta.s.p.vec = apply(weight.perturb, 2, delta.s.surv.estimate, xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, sone=sone, szero=szero, t=t, landmark=landmark, extrapolate = extrapolate, transform = transform,approx=approx, warn.te = TRUE, warn.support = TRUE)
		delta.p.vec = as.vector(unlist(apply(weight.perturb, 2, delta.surv.estimate, xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, var= FALSE, conf.int = FALSE,approx=approx)))
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
		R.t.est = R.t.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, landmark = landmark,var = var, weight.perturb = weight.perturb,approx=approx)
		delta.t = R.t.est$delta.t
		R.t = R.t.est$R.t
		IV = R.s - R.t
		if(var | conf.int){
		delta.t.p.vec = apply(weight.perturb, 2, delta.t.surv.estimate, xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, landmark=landmark,approx=approx)
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


 delta.t.surv.estimate = function(xone,xzero, deltaone, deltazero, t, weight.perturb = NULL, landmark, approx = T) {
 		delta.check = delta.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, conf.int = TRUE,approx=approx)
	if(delta.check$conf.int.quantile[1] < 0 & delta.check$conf.int.quantile[2] > 0) {print("Warning: it looks like the treatment effect is not significant; may be difficult to interpret the residual treatment effect in this setting")}
	if(is.null(weight.perturb)) {weight = rep(1,length(xone)+length(xzero))}
	if(!is.null(weight.perturb)) {weight = weight.perturb}	
	weight.group1 = weight[1:length(xone)]
	weight.group0 = weight[(1+length(xone)):(length(xone)+length(xzero))]
	censor0.t = censor.weight(xzero, deltazero, t, weight = weight.group0,approx=approx)
	censor0.landmark = censor.weight(xzero, deltazero, landmark, weight = weight.group0,approx=approx)
	censor1.t = censor.weight(xone, deltaone, t, weight = weight.group1,approx=approx)
	censor1.landmark = censor.weight(xone, deltaone, landmark, weight = weight.group1,approx=approx)
	cond.prop = (sum(weight.group1*1*(xone>t))/sum(weight.group1*censor1.t))/(sum(weight.group1*1*(xone>landmark))/sum(weight.group1*censor1.landmark))
	delta.t = sum(weight.group0*1*(xzero>landmark))/sum(weight.group0*censor0.landmark)*cond.prop - sum(weight.group0*1*(xzero>t))/sum(weight.group0*censor0.t)
	return(delta.t)
}


R.t.surv.estimate = function(xone,xzero, deltaone, deltazero, t, weight.perturb = NULL, landmark, var = FALSE, conf.int = FALSE, approx = T) {	
	delta = delta.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t,approx=approx)$delta	
	delta.t = delta.t.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero,  t=t, landmark=landmark,approx=approx)
	R.t = 1-delta.t/delta	
	if(var | conf.int){
		if(is.null(weight.perturb)) {
	weight.perturb = matrix(rexp(500*(length(xone)+length(xzero)), rate=1), ncol = 500)}
		delta.t.p.vec = apply(weight.perturb, 2, delta.t.surv.estimate, xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, landmark=landmark,approx=approx)
		delta.p.vec = as.vector(unlist(apply(weight.perturb, 2, delta.surv.estimate, xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, var= FALSE, conf.int = FALSE,approx=approx)))
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




pred.smooth.surv <- function(xone.f, deltaone.f, sone.f, szero.one, myt, weight.pred, extrapolate, transform, ps.weight = NULL, warn.support = FALSE)
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
    if(is.null(ps.weight)) {kerni.1 = t(weight.pred*t(kerni.ss))}
    if(!is.null(ps.weight)) {kerni.1 = t(ps.weight*weight.pred*t(kerni.ss))}
    pihamyt0.tj.ss = helper.si(tj, "<=", xone.f, Vi=t(kerni.1)) ## n.tj x n.ss matrix ##   
    dLamhat.tj.ss = t((kerni.1[,tmpind]))/pihamyt0.tj.ss; 
    ret = apply(dLamhat.tj.ss,2,sum)
    Phat.ss  =exp(-ret)
    if(sum(is.na(Phat.ss))>0 & extrapolate){
    	if(!warn.support) {print(paste("Note: ", sum(is.na(Phat.ss)), " values extrapolated."))}
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


Aug.R.s.surv.estimate = function(xone,xzero, deltaone, deltazero, sone, szero, t, weight.perturb = NULL, landmark, extrapolate = FALSE, transform = FALSE, basis.delta.one, basis.delta.zero,basis.delta.s.one = NULL, basis.delta.s.zero = NULL, incremental.value = FALSE, approx = T) {
	delta = delta.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t,approx=approx)$delta	
	delta.s = delta.s.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, sone=sone, szero=szero, t=t, landmark=landmark, extrapolate = extrapolate, transform = transform,approx=approx)
	R.s = 1-delta.s/delta	
 	if(is.null(weight.perturb)) {
	weight.perturb = matrix(rexp(500*(length(xone)+length(xzero)), rate=1), ncol = 500)}
	delta.s.p.vec = apply(weight.perturb, 2, delta.s.surv.estimate, xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, sone=sone, szero=szero, t=t, landmark=landmark, extrapolate = extrapolate, transform = transform,approx=approx)
	delta.p.vec = as.vector(unlist(apply(weight.perturb, 2, delta.surv.estimate, xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, var= FALSE, conf.int = FALSE,approx=approx)))
	R.p = 1-delta.s.p.vec/delta.p.vec
	
	#Augmented estimators
	if(dim(as.matrix(basis.delta.one))[2] == 1 ) {basis.delta.use =  rbind(cbind(basis.delta.one, basis.delta.one^2), cbind(basis.delta.zero, basis.delta.zero^2))}
	if(dim(as.matrix(basis.delta.one))[2] > 1 ) {basis.delta.use = rbind(basis.delta.one, basis.delta.zero) }
	basis.delta.s.use = basis.delta.use
	if(!is.null(basis.delta.s.one) | !is.null(basis.delta.s.zero)) {print("Note: basis.delta.s has been set to be equal to basis.delta.")}
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
		R.t.est = R.t.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, landmark = landmark,var = TRUE, weight.perturb = weight.perturb,approx=approx)
		delta.t = R.t.est$delta.t
		R.t = R.t.est$R.t
		delta.t.p.vec = apply(weight.perturb, 2, delta.t.surv.estimate, xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, landmark=landmark,approx=approx)
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

fieller.calculate.me = function(a, B, var.mat){
	Z.squared = qnorm(0.975)^2
	A1 = - (2*B*a - Z.squared*var.mat[1,2]*2)
	A2 = (2*B*a - Z.squared*2*var.mat[1,2])^2 -4*(B^2 - Z.squared*var.mat[2,2])*(a^2 - Z.squared*var.mat[1,1])
	A3 = 2*(B^2 - Z.squared*var.mat[2,2])
	int.1 = (A1 - sqrt(A2)) / A3
	int.2 = (A1 + sqrt(A2)) / A3
	#interval for a/b is int.1, int.2
	#interval for R is 1-int.2, 1- int.1
	return(c(1+int.1, 1+int.2))
}

R.s.estimate.me = function(sone,szero,yone,yzero,parametric = FALSE, estimator = "n", me.variance, extrapolate = TRUE, transform = FALSE, naive = FALSE, Ronly = TRUE) {
	if(!(estimator %in% c("d", "q", "n"))) {print("Warning: Invalid estimator, default `nonlinear' being used"); estimator = "n"}
	if(!parametric & estimator == "d") {print("Warning: Invalid estimator, there is no disattentuated estimator for nonparametric approach; default `nonlinear' being used"); estimator = "n"}
	outcome.vec = c(yone, yzero)
	treat.vec = c(rep(1,length(yone)),rep(0,length(yzero)))
	s.vec = c(sone, szero)	
	
	#IF PARAMETRIC
	if(parametric) {
		#Naive estimate
		m1 = lm(outcome.vec~treat.vec)
		m2.error = lm(outcome.vec~treat.vec + s.vec)
		Z.model1 = m1$coef[2]
		Z.model2.error = m2.error$coef[2]
		R.vec.error = 1-Z.model2.error/Z.model1
		beta.var.error.vec= summary(m2.error)$coeff[2,2]^2
		residual.var.2 = summary(m2.error)$sigma^2
		design.mat = cbind(rep(1,length(treat.vec)), treat.vec)
		rr = t(design.mat)%*%design.mat
		sigma.2 = residual.var.2 * solve(rr)
		term.2 = sigma.2[2,1]  
		R.var.error.vec = (Z.model2.error^2/ Z.model1^2)*(beta.var.error.vec/Z.model2.error^2 -2*term.2/(Z.model2.error*Z.model1) + summary(m1)$coeff[2,2]^2/Z.model1^2)
		var.mat.fieller = matrix(c(beta.var.error.vec, term.2, term.2,summary(m1)$coeff[2,2]^2), nrow = 2, ncol = 2, byrow = 2)
		R.var.error.vec.fieller = fieller.calculate.me(Z.model2.error, Z.model1, var.mat.fieller)
		R.error.normal = c(R.vec.error - 1.96*sqrt(R.var.error.vec),R.vec.error + 1.96*sqrt(R.var.error.vec) )
		beta.error.normal = c(Z.model2.error - 1.96*sqrt(beta.var.error.vec),Z.model2.error + 1.96*sqrt(beta.var.error.vec) )
		fullnaive.list = list("R.naive" = as.numeric(R.vec.error), "R.naive.var" = as.numeric(R.var.error.vec), "R.naive.CI.normal" = as.numeric(R.error.normal), "R.naive.CI.fieller" = as.numeric(R.var.error.vec.fieller),"B1star.naive" = as.numeric(Z.model2.error), "B1star.naive.var" = as.numeric(beta.var.error.vec), "B1star.naive.CI.normal" = as.numeric(beta.error.normal))
		Ronlynaive.list = list("R.naive" = as.numeric(R.vec.error), "R.naive.var" = as.numeric(R.var.error.vec), "R.naive.CI.normal" = as.numeric(R.error.normal), "R.naive.CI.fieller" = as.numeric(R.var.error.vec.fieller))
		if(estimator == "d") {
			S.model2.mm = m2.error$coef[3]
			var.w = var(s.vec)
			var.error = me.variance
			var.g = var(treat.vec)
			cov.gw = cov(treat.vec, s.vec)
			num = var.w*cov.gw - cov.gw*(var.w - var.error)
			den = var.g*(var.w -var.error) - cov.gw*cov.gw
			rho = num/den
			beta1.fixed = Z.model2.error - rho*S.model2.mm
			R.fixed = 1-beta1.fixed/Z.model1
			####variation attenuation   
			beta.var.mm.vec = summary(m2.error)$coeff[2,2]^2 + rho^2*summary(m2.error)$coeff[3,2]^2 - 2*rho*vcov(m2.error)[2,3]
			R.var.mm.vec = (beta1.fixed^2/ Z.model1^2)*(beta.var.mm.vec/beta1.fixed^2 -2*term.2/(beta1.fixed*Z.model1) + summary(m1)$coeff[2,2]^2/Z.model1^2)
			var.mat.fieller = matrix(c(beta.var.mm.vec, term.2, term.2,summary(m1)$coeff[2,2]^2), nrow = 2, ncol = 2, byrow = 2)
			R.var.mm.vec.fieller = fieller.calculate.me(beta1.fixed, Z.model1, var.mat.fieller)
			R.mm.normal = c(R.fixed - 1.96*sqrt(R.var.mm.vec),R.fixed + 1.96*sqrt(R.var.mm.vec) )
			beta.mm.normal = c(beta1.fixed - 1.96*sqrt(beta.var.mm.vec),beta1.fixed + 1.96*sqrt(beta.var.mm.vec) )
			#summary
			full.list = list("R.corrected.dis" = as.numeric(R.fixed), "R.corrected.var.dis" = as.numeric(R.var.mm.vec), "R.corrected.CI.normal.dis" = as.numeric(R.mm.normal), "R.corrected.CI.fieller.dis" = as.numeric(R.var.mm.vec.fieller),"B1star.corrected.dis" = as.numeric(beta1.fixed), "B1star.corrected.var.dis" = as.numeric(beta.var.mm.vec), "B1star.corrected.CI.normal.dis" = as.numeric(beta.mm.normal))
			Ronly.list = list("R.corrected.dis" = as.numeric(R.fixed), "R.corrected.var.dis" = as.numeric(R.var.mm.vec), "R.corrected.CI.normal.dis" = as.numeric(R.mm.normal), "R.corrected.CI.fieller.dis" = as.numeric(R.var.mm.vec.fieller))	
	
		}
		if(estimator == "q" | estimator == "n") {
		n.pseudo = 50
		lambda.vec = c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5,1.75, 2.0)
		Z.model1.temp = lm(outcome.vec~treat.vec)$coef[2]
		intercept.model1.temp = lm(outcome.vec~treat.vec)$coef[1]
		#naive, lambda=0
		Z.model2.temp = lm(outcome.vec~treat.vec+s.vec)$coef[2]
		R.vec.temp = 1-Z.model2.temp/Z.model1.temp
		Z.model2.pseudo = matrix(nrow = length(lambda.vec), ncol = n.pseudo)
		S.model2.pseudo = matrix(nrow = length(lambda.vec), ncol = n.pseudo)
		intercept.model2.pseudo = matrix(nrow = length(lambda.vec), ncol = n.pseudo)
		#pseudo-errors
		pe = matrix(rnorm(length(outcome.vec)*n.pseudo,0,1), ncol = n.pseudo)
		mat.lambda.for.C = c()
		mat.A = c()
		
		for(k in 1:length(lambda.vec))	{
			mat.for.C.ee1 = matrix(nrow = length(outcome.vec), ncol = n.pseudo)
			mat.for.C.ee2 = matrix(nrow = length(outcome.vec), ncol = n.pseudo)
			mat.for.C.ee3 = matrix(nrow = length(outcome.vec), ncol = n.pseudo)
			mat.for.C.ee4 = matrix(nrow = length(outcome.vec), ncol = n.pseudo)
			mat.for.C.ee5 = matrix(nrow = length(outcome.vec), ncol = n.pseudo)
			mat.for.A.1 = matrix(nrow = length(outcome.vec), ncol = n.pseudo)
			mat.for.A.2 = matrix(nrow = length(outcome.vec), ncol = n.pseudo)
			mat.for.A.3 = matrix(nrow = length(outcome.vec), ncol = n.pseudo)
			mat.for.A.4 = matrix(nrow = length(outcome.vec), ncol = n.pseudo)
			mat.for.A.5 = matrix(nrow = length(outcome.vec), ncol = n.pseudo)
			mat.for.A.6 = matrix(nrow = length(outcome.vec), ncol = n.pseudo)
			mat.for.A.7 = matrix(nrow = length(outcome.vec), ncol = n.pseudo)
			mat.for.A.8 = matrix(nrow = length(outcome.vec), ncol = n.pseudo)
			mat.for.A.9 = matrix(nrow = length(outcome.vec), ncol = n.pseudo)
			mat.for.A.10 = matrix(nrow = length(outcome.vec), ncol = n.pseudo)
			mat.for.A.11 = matrix(nrow = length(outcome.vec), ncol = n.pseudo)
			mat.for.A.12 = matrix(nrow = length(outcome.vec), ncol = n.pseudo)
			mat.for.A.13 = matrix(nrow = length(outcome.vec), ncol = n.pseudo)
			
			for(j in 1:50)	{
				pseudo.vals = s.vec+sqrt(lambda.vec[k])*sqrt(me.variance)*pe[,j]
				model.simex = lm(outcome.vec~treat.vec+pseudo.vals)
				Z.model2.pseudo[k,j] = model.simex$coef[2]
				S.model2.pseudo[k,j] = model.simex$coef[3]
				intercept.model2.pseudo[k,j] = model.simex$coef[1]
				mat.for.C.ee1[,j] = 2*(outcome.vec - intercept.model2.pseudo[k,j] - Z.model2.pseudo[k,j]*treat.vec - S.model2.pseudo[k,j]*pseudo.vals)*(-1)
				mat.for.C.ee2[,j] = 2*(outcome.vec - intercept.model2.pseudo[k,j] - Z.model2.pseudo[k,j]*treat.vec - S.model2.pseudo[k,j]*pseudo.vals)*(-treat.vec)
				mat.for.C.ee3[,j] = 2*(outcome.vec - intercept.model2.pseudo[k,j] - Z.model2.pseudo[k,j]*treat.vec - S.model2.pseudo[k,j]*pseudo.vals)*(-pseudo.vals)
				mat.for.C.ee4[,j] = 2*(outcome.vec - intercept.model1.temp  - Z.model1.temp*treat.vec)*(-1)
				mat.for.C.ee5[,j] = 2*(outcome.vec - intercept.model1.temp  - Z.model1.temp*treat.vec)*(-treat.vec)
				mat.for.A.1[,j] = rep(2,length(outcome.vec))
				mat.for.A.2[,j] = 2*treat.vec
				mat.for.A.3[,j] = 2*pseudo.vals
				mat.for.A.4[,j] = 2*treat.vec
				mat.for.A.5[,j] = 2*treat.vec^2
				mat.for.A.6[,j] = 2*treat.vec*pseudo.vals
				mat.for.A.7[,j] = 2*pseudo.vals
				mat.for.A.8[,j] = 2*treat.vec*pseudo.vals
				mat.for.A.9[,j] = 2*pseudo.vals^2
				mat.for.A.10[,j] = rep(2,length(outcome.vec))
				mat.for.A.11[,j] = 2*treat.vec
				mat.for.A.12[,j] = 2*treat.vec
				mat.for.A.13[,j] = 2*treat.vec^2
			}
			chi.ee1 = apply(mat.for.C.ee1, 1, mean)
			chi.ee2 = apply(mat.for.C.ee2, 1, mean)
			chi.ee3 = apply(mat.for.C.ee3, 1, mean)
			chi.ee4 = apply(mat.for.C.ee4, 1, mean)
			chi.ee5 = apply(mat.for.C.ee5, 1, mean)
			mat.lambda.for.C = cbind(mat.lambda.for.C,chi.ee1, chi.ee2, chi.ee3, chi.ee4, chi.ee5)
			A.mat.hold = matrix(c(mean(mat.for.A.1), mean(mat.for.A.2), mean(mat.for.A.3), 0, 0, mean(mat.for.A.4), mean(mat.for.A.5), mean(mat.for.A.6), 0,0,mean(mat.for.A.7), mean(mat.for.A.8), mean(mat.for.A.9),0,0,0,0,0,mean(mat.for.A.10), mean(mat.for.A.11), 0,0,0,mean(mat.for.A.12), mean(mat.for.A.13)), nrow = 5, ncol = 5, byrow = T)
			if(k==1) {mat.A = A.mat.hold}
			if(k > 1) {mat.A = as.matrix(bdiag(mat.A, A.mat.hold))}
		}
		C11 =cov(mat.lambda.for.C)
		A11 = mat.A
		mean.ps = apply(Z.model2.pseudo,1,mean)
		mean.Sbeta = apply(S.model2.pseudo,1,mean)
		mean.intercept = apply(intercept.model2.pseudo,1,mean)
		
		#quadratic estimate
		lambda.vec.sq = lambda.vec^2
		model.q = lm(mean.ps ~ lambda.vec + lambda.vec.sq)
		est.simex.q = model.q$coef%*%c(1,-1,1)
		R.vec.simex.q = 1-est.simex.q/Z.model1
		
		#quadratic variance
		G.gamma.lambda.p.q =  function(lam)	{return(matrix(c(1,lam,lam^2,0,0,0,0,0,0,0,0,0,0,0,1,lam,lam^2,0,0,0,0,0,0,0,0,0,0,0,1,lam,lam^2,0,0,0,0,0,0,0,0,0,0,0,1,0,rep(0,10),1), nrow = 11, ncol = 5, byrow = F))}
		G.gamma.neg1 = G.gamma.lambda.p.q(-1)
		S.gamma = c()
		for(h in 1:length(lambda.vec)){
			S.gamma = cbind(S.gamma, G.gamma.lambda.p.q(lambda.vec[h]))
		}
		D.gamma = S.gamma %*% t(S.gamma)
		Sigma = solve(A11)%*%C11%*%t(solve(A11))
		#see carroll book, don't need C.mat
		Sigma.gamma = solve(D.gamma)%*%S.gamma%*% Sigma %*%t(S.gamma)%*% solve(D.gamma)
		Simex.variance = t(G.gamma.neg1)%*%Sigma.gamma %*% G.gamma.neg1 
		sample.size = length(treat.vec)
		beta.var.simex.quad.vec = Simex.variance[2,2]/sample.size
		cov.simex = Simex.variance[2,5]/sample.size
		R.var.simex.q.vec = (est.simex.q^2/ Z.model1.temp^2)*(beta.var.simex.quad.vec/est.simex.q^2 -2*cov.simex/(est.simex.q*Z.model1.temp) + summary(lm(outcome.vec~treat.vec))$coef[2,2]^2/Z.model1.temp^2)
		var.mat.fieller = matrix(c(Simex.variance[2,2], Simex.variance[2,5], Simex.variance[2,5],Simex.variance[5,5]), nrow = 2, ncol = 2, byrow = 2)/sample.size
		R.var.simex.q.vec.fieller = fieller.calculate.me(est.simex.q, Z.model1.temp, var.mat.fieller)
		R.simex.q.normal = c(R.vec.simex.q - 1.96*sqrt(R.var.simex.q.vec), R.vec.simex.q + 1.96*sqrt(R.var.simex.q.vec))
		beta.simex.q.normal = c(est.simex.q - 1.96*sqrt(beta.var.simex.quad.vec), est.simex.q + 1.96*sqrt(beta.var.simex.quad.vec))
		if(estimator == "q") {
			full.list = list("R.corrected.q" = as.numeric(R.vec.simex.q), "R.corrected.var.q" = as.numeric(R.var.simex.q.vec), "R.corrected.CI.normal.q" = as.numeric(R.simex.q.normal), "R.corrected.CI.fieller.q" = as.numeric(R.var.simex.q.vec.fieller),"B1star.corrected.q" = as.numeric(est.simex.q), "B1star.corrected.var.q" = as.numeric(beta.var.simex.quad.vec), "B1star.corrected.CI.normal.q" = as.numeric(beta.simex.q.normal))
			Ronly.list = list("R.corrected.q" = as.numeric(R.vec.simex.q), "R.corrected.var.q" = as.numeric(R.var.simex.q.vec), "R.corrected.CI.normal.q" = as.numeric(R.simex.q.normal), "R.corrected.CI.fieller.q" = as.numeric(R.var.simex.q.vec.fieller))
		}  #end quadratic return
		
		if(estimator == "n") {
		#nonlinear est and variance
		parameters = c(NA,NA,NA)
		Z.model2.simex.nl = NA
		R.vec.simex.nl = NA
		beta.var.simex.nl.vec = NA
		R.var.simex.nl.vec = NA
		R.var.simex.nl.vec.fieller = c(NA,NA)
		R.simex.nl.normal = c(NA,NA)
		beta.simex.nl.normal = c(NA,NA)
		tryCatch({
		model.nl = nls(mean.ps ~ a.param + b.param/(c.param+lambda.vec), start = list(a.param=1, b.param=1, c.param = 1), control = nls.control(maxiter = 1000,  warnOnly = TRUE, minFactor = 1/100000))
		parameters = summary(model.nl)$parameters[,1]
		alpha_0 = parameters[1]
		alpha_1 = parameters[2]
		alpha_2 = parameters[3]
		est.simex.nl = parameters[1] + parameters[2]/(parameters[3] -1)
		parameters = c(NA,NA,NA)
		model.nl.intercept = nls(mean.intercept ~ a.param + b.param/(c.param+lambda.vec), start = list(a.param=1, b.param=1, c.param = 1), control = nls.control(maxiter = 1000,  warnOnly = TRUE, minFactor = 1/100000))
		parameters = summary(model.nl.intercept)$parameters[,1]
		alpha_0_star = parameters[1]
		alpha_1_star = parameters[2]
		alpha_2_star = parameters[3]
		model.nl.Sbeta = nls(mean.Sbeta ~ a.param + b.param/(c.param+lambda.vec), start = list(a.param=1, b.param=1, c.param = 1), control = nls.control(maxiter = 1000,  warnOnly = TRUE, minFactor = 1/100000))
		parameters = summary(model.nl.Sbeta)$parameters[,1]
		alpha_0_starstar = parameters[1]
		alpha_1_starstar = parameters[2]
		alpha_2_starstar = parameters[3]
	
		#variance of nonlinear simex
		param.nl.vector = c(	alpha_0, alpha_1, alpha_2, alpha_0_star, alpha_1_star, alpha_2_star, alpha_0_starstar, alpha_1_starstar, alpha_2_starstar)
		G.gamma.lambda.p.nl =  function(lam){return(matrix(c(1,(alpha_2 + lam)^(-1),-1*alpha_1*(alpha_2 +lam)^(-2),0,0,0,0,0,0,0,0,0,0,0,1,(alpha_2_star + lam)^(-1),-1*alpha_1_star*(alpha_2_star +lam)^(-2),0,0,0,0,0,0,0,0,0,0,0,1, (alpha_2_starstar + lam)^(-1),-1*alpha_1_starstar*(alpha_2_starstar +lam)^(-2),0,0,rep(0,9),1,0, rep(0,10),1), nrow = 11, ncol = 5, byrow = F))}
		G.gamma.neg1 = G.gamma.lambda.p.nl(-1)
		S.gamma = c()
		for(h in 1:length(lambda.vec)){
			S.gamma = cbind(S.gamma, G.gamma.lambda.p.nl(lambda.vec[h]))
		}
		D.gamma = S.gamma %*% t(S.gamma)
		Sigma = solve(A11)%*%C11%*%t(solve(A11))
		#see carroll book, don't need C.mat
		Sigma.gamma = solve(D.gamma)%*%S.gamma%*% Sigma %*%t(S.gamma)%*% solve(D.gamma)
		Simex.variance = t(G.gamma.neg1)%*%Sigma.gamma %*% G.gamma.neg1 
		Z.model2.simex.nl = est.simex.nl
		R.vec.simex.nl = 1-Z.model2.simex.nl/Z.model1
		beta.var.simex.nl.vec = Simex.variance[2,2]/sample.size 
		cov.simex = Simex.variance[2,5]/sample.size
		R.var.simex.nl.vec = (Z.model2.simex.nl^2/ Z.model1.temp^2)*(beta.var.simex.nl.vec/Z.model2.simex.nl^2 -2*cov.simex/(Z.model2.simex.nl*Z.model1.temp) + summary(lm(outcome.vec~treat.vec))$coef[2,2]^2/Z.model1.temp^2)
		var.mat.fieller = matrix(c(Simex.variance[2,2], Simex.variance[2,5], Simex.variance[2,5],Simex.variance[5,5]), nrow = 2, ncol = 2, byrow = 2)/sample.size
		R.var.simex.nl.vec.fieller = fieller.calculate.me(Z.model2.simex.nl, Z.model1.temp, var.mat.fieller) 
		R.simex.nl.normal = c(R.vec.simex.nl - 1.96*sqrt(R.var.simex.nl.vec), R.vec.simex.nl + 1.96*sqrt(R.var.simex.nl.vec))
		beta.simex.nl.normal = c(Z.model2.simex.nl - 1.96*sqrt(beta.var.simex.nl.vec), Z.model2.simex.nl + 1.96*sqrt(beta.var.simex.nl.vec))
		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

		full.list = list("R.corrected.nl" = as.numeric(R.vec.simex.nl), "R.corrected.var.nl" = as.numeric(R.var.simex.nl.vec), "R.corrected.CI.normal.nl" = as.numeric(R.simex.nl.normal), "R.corrected.CI.fieller.nl" = as.numeric(R.var.simex.nl.vec.fieller),"B1star.corrected.nl" = as.numeric(Z.model2.simex.nl), "B1star.corrected.var.nl" = as.numeric(beta.var.simex.nl.vec), "B1star.corrected.CI.normal.nl" = as.numeric(beta.simex.nl.normal))
		Ronly.list = list("R.corrected.nl" = as.numeric(R.vec.simex.nl), "R.corrected.var.nl" = as.numeric(R.var.simex.nl.vec), "R.corrected.CI.normal.nl" = as.numeric(R.simex.nl.normal), "R.corrected.CI.fieller.nl" = as.numeric(R.var.simex.nl.vec.fieller))		
		} #end nonlinear return	
		}  #end simex i.e. q or nl
	
	}  #close parametric if statement
	
	if(!parametric){
		
		#naive
		R.witherror = R.s.estimate(sone = sone, szero = szero, yone = yone, yzero = yzero, type = "robust", extrapolate = extrapolate, transform = transform)
		R.vec.error.NP = R.witherror$R.s
		deltas.vec.error.NP = R.witherror$delta.s
		delta.est = R.witherror$delta

	
		#estimating variance, closed form
		vec.temp = calculate.var.np(s1 = sone,s0=szero,y1 = yone,y0 = yzero)
		sample.size = length(treat.vec)
		multi.vec = cbind(rep(c(1,0),sample.size), rep(c(0,1),sample.size))
		vv.cov.mat = (vec.temp$total%*%multi.vec)/(sample.size^2)
		deltas.var.vec.error.NP = vv.cov.mat[1,1]
		delta.var.vec.error.NP = vv.cov.mat[2,2]
		theta1 = R.witherror$delta.s
		theta2 = R.witherror$delta 
		R.var.vec.error.NP = theta1^2/theta2^2*(vv.cov.mat[1,1]/theta1^2 - 2*vv.cov.mat[2,1]/(theta1*theta2) + vv.cov.mat[2,2]/theta2^2)
		var.mat.fieller = vv.cov.mat
		R.var.vec.error.NP.fieller = fieller.calculate.me(R.witherror$delta.s, R.witherror$delta, var.mat.fieller)
		R.naive.normal = c(R.vec.error.NP - 1.96*sqrt(R.var.vec.error.NP), R.vec.error.NP + 1.96*sqrt(R.var.vec.error.NP))
		deltas.naive.normal = c(deltas.vec.error.NP - 1.96*sqrt(deltas.var.vec.error.NP), deltas.vec.error.NP + 1.96*sqrt(deltas.var.vec.error.NP)) 
		fullnaive.list = list("R.naive" = as.numeric(R.vec.error.NP), "R.naive.var" = as.numeric(R.var.vec.error.NP), "R.naive.CI.normal" = as.numeric(R.naive.normal), "R.naive.CI.fieller" = as.numeric(R.var.vec.error.NP.fieller),"deltas.naive" = as.numeric(deltas.vec.error.NP), "deltas.naive.var" = as.numeric(deltas.var.vec.error.NP), "deltas.naive.CI.normal" = as.numeric(deltas.naive.normal))
		Ronlynaive.list = list("R.naive" = as.numeric(R.vec.error.NP), "R.naive.var" = as.numeric(R.var.vec.error.NP), "R.naive.CI.normal" = as.numeric(R.naive.normal), "R.naive.CI.fieller" = as.numeric(R.var.vec.error.NP.fieller))
		
		#SIMEX	
		n.pseudo = 50
		lambda.vec = c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5,1.75, 2.0)
		R.pseudo = matrix(nrow = length(lambda.vec), ncol = n.pseudo)
		#pseudo-errors
	
		pe = matrix(rnorm(length(outcome.vec)*n.pseudo,0,1), ncol = n.pseudo)
		mat.lambda.for.C = c()
		for(k in 1:length(lambda.vec))	{
			mat.for.C.ee1 = matrix(nrow = length(outcome.vec), ncol = n.pseudo)
			mat.for.C.ee2 = matrix(nrow = length(outcome.vec), ncol = n.pseudo)
		for(j in 1:50)	{
			pseudo.vals = s.vec+sqrt(lambda.vec[k])*sqrt(me.variance)*pe[,j]
			R.pseudo[k,j] = R.s.estimate(sone = pseudo.vals[treat.vec==1], szero = pseudo.vals[treat.vec==0], yone = outcome.vec[treat.vec==1], yzero = outcome.vec[treat.vec==0], type = "robust", extrapolate = extrapolate, transform=transform, warn.te = TRUE, warn.support = TRUE)$delta.s
			psi.hold = calculate.var.np(s1 = pseudo.vals[treat.vec==1],s0=pseudo.vals[treat.vec==0],y1 = outcome.vec[treat.vec==1],y0 = outcome.vec[treat.vec==0])$psionly 
			mat.for.C.ee1[,j] = t(psi.hold[1,])
			mat.for.C.ee2[,j] = t(psi.hold[2,])
		}
			chi.ee1 = apply(mat.for.C.ee1, 1, mean)
			chi.ee2 = apply(mat.for.C.ee2, 1, mean)
			mat.lambda.for.C = cbind(mat.lambda.for.C,chi.ee1, chi.ee2)
		}
	
		#note - even though this says R, see above that it is delta.s
		mean.ps = apply(R.pseudo,1,mean)
		#quadratic
		lambda.vec.sq = lambda.vec^2
		model.q = lm(mean.ps ~ lambda.vec + lambda.vec.sq)
		est.simex.q = model.q$coef%*%c(1,-1,1)
		R.vec.simex.q.NP = 1-est.simex.q/delta.est
		deltas.vec.simex.q.NP = est.simex.q
	
		#variance of quadratic simex
		G.gamma.lambda.np.q =  function(lam){return(matrix(c(1,lam,lam^2,0,0,0,0,1), nrow = 4, ncol = 2, byrow = F))}
		G.gamma.neg1 = G.gamma.lambda.np.q(-1)
		S.gamma = c()
		for(h in 1:length(lambda.vec)){
			S.gamma = cbind(S.gamma, G.gamma.lambda.np.q(lambda.vec[h]))
		}
		D.gamma = S.gamma %*% t(S.gamma)
		Sigma = cov(mat.lambda.for.C)
		Sigma.gamma = solve(D.gamma)%*%S.gamma%*% Sigma %*%t(S.gamma)%*% solve(D.gamma)
		Simex.variance = t(G.gamma.neg1)%*%Sigma.gamma %*% G.gamma.neg1 
		deltas.var.vec.simex.q.NP  = Simex.variance[1,1]/sample.size
		delta.var.vec.simex.q.NP  = Simex.variance[2,2]/sample.size
		cov.simex = Simex.variance[1,2]/sample.size
		delta.var.simex = Simex.variance[2,2]/sample.size
		theta1 = deltas.vec.simex.q.NP
		theta2 = R.witherror$delta 
		R.var.vec.simex.q.NP = theta1^2/theta2^2*(deltas.var.vec.simex.q.NP/theta1^2 - 2*cov.simex/(theta1*theta2) + delta.var.simex/theta2^2)
		var.mat.fieller = Simex.variance/sample.size
		R.var.vec.simex.q.NP.fieller = fieller.calculate.me(theta1, theta2, var.mat.fieller)
		R.simex.normal.q = c(R.vec.simex.q.NP - 1.96*sqrt(R.var.vec.simex.q.NP), R.vec.simex.q.NP + 1.96*sqrt(R.var.vec.simex.q.NP))
		deltas.simex.normal.q = c(deltas.vec.simex.q.NP - 1.96*sqrt(deltas.var.vec.simex.q.NP), deltas.vec.simex.q.NP + 1.96*sqrt(deltas.var.vec.simex.q.NP))
		if(estimator == "q") {
			full.list = list("R.corrected.q" = as.numeric(R.vec.simex.q.NP), "R.corrected.var.q" = as.numeric(R.var.vec.simex.q.NP), "R.corrected.CI.normal.q" = as.numeric(R.simex.normal.q ), "R.corrected.CI.fieller.q" = as.numeric(R.var.vec.simex.q.NP.fieller),"deltas.corrected.q" = as.numeric(deltas.vec.simex.q.NP), "deltas.corrected.var.q" = as.numeric(deltas.var.vec.simex.q.NP), "deltas.corrected.CI.normal.q" = as.numeric(deltas.simex.normal.q))
			Ronly.list = list("R.corrected.q" = as.numeric(R.vec.simex.q.NP), "R.corrected.var.q" = as.numeric(R.var.vec.simex.q.NP), "R.corrected.CI.normal.q" = as.numeric(R.simex.normal.q ), "R.corrected.CI.fieller.q" = as.numeric(R.var.vec.simex.q.NP.fieller))
		}  #end quadratic return
		
		if(estimator == "n") {
		#nonlinear
		parameters = c(NA,NA,NA)
		R.vec.simex.nl.NP = NA
		deltas.vec.simex.nl.NP = NA
		deltas.var.vec.simex.nl.NP = NA
		R.var.vec.simex.nl.NP = NA
		R.var.vec.simex.nl.NP.fieller = c(NA,NA)
		R.simex.normal.nl = c(NA,NA)
		deltas.simex.normal.nl = c(NA,NA)
		
		#put all nonlinear in trycatch
		tryCatch({
		model.nl = nls(mean.ps ~ a.param + b.param/(c.param+lambda.vec), start = list(a.param=1, b.param=1, c.param = 1), control = nls.control(warnOnly = TRUE, minFactor = 1/100000, maxiter = 1000))
		parameters = summary(model.nl)$parameters[,1]
		alpha_0 = parameters[1]
		alpha_1 = parameters[2]
		alpha_2 = parameters[3]
		est.simex.nl = parameters[1] + parameters[2]/(parameters[3] -1)
	
		#variance of nonlinear simex
		param.nl.vector = c(	alpha_0, alpha_1, alpha_2)
		G.gamma.lambda.np.nl =  function(lam){return(matrix(c(1,(alpha_2 + lam)^(-1),-1*alpha_1*(alpha_2 +lam)^(-2),0,0,0,0,1), nrow = 4, ncol = 2, byrow = F))}
		G.gamma.neg1 = G.gamma.lambda.np.nl(-1)
		S.gamma = c()
		for(h in 1:length(lambda.vec)){
			S.gamma = cbind(S.gamma, G.gamma.lambda.np.nl(lambda.vec[h]))
		}
		D.gamma = S.gamma %*% t(S.gamma)
		Sigma = cov(mat.lambda.for.C)
		Sigma.gamma = solve(D.gamma)%*%S.gamma%*% Sigma %*%t(S.gamma)%*% solve(D.gamma)
		Simex.variance = t(G.gamma.neg1)%*%Sigma.gamma %*% G.gamma.neg1 
		R.vec.simex.nl.NP = 1-est.simex.nl/delta.est
		deltas.vec.simex.nl.NP = est.simex.nl
		deltas.var.vec.simex.nl.NP  = Simex.variance[1,1]/sample.size
		delta.var.vec.simex.nl.NP  = Simex.variance[2,2]/sample.size
		cov.simex = Simex.variance[1,2]/sample.size
		delta.var.simex = Simex.variance[2,2]/sample.size
		theta1 = deltas.vec.simex.nl.NP
		theta2 = R.witherror$delta 
		R.var.vec.simex.nl.NP = theta1^2/theta2^2*(deltas.var.vec.simex.nl.NP/theta1^2 - 2*cov.simex/(theta1*theta2) + delta.var.simex/theta2^2)
		var.mat.fieller = Simex.variance/sample.size
		R.var.vec.simex.nl.NP.fieller = fieller.calculate.me(theta1, theta2, var.mat.fieller)
		R.simex.normal.nl = c(R.vec.simex.nl.NP - 1.96*sqrt(R.var.vec.simex.nl.NP), R.vec.simex.nl.NP + 1.96*sqrt(R.var.vec.simex.nl.NP))
		deltas.simex.normal.nl = c(deltas.vec.simex.nl.NP - 1.96*sqrt(deltas.var.vec.simex.nl.NP), deltas.vec.simex.nl.NP + 1.96*sqrt(deltas.var.vec.simex.nl.NP))}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
		full.list = list("R.corrected.nl" = as.numeric(R.vec.simex.nl.NP), "R.corrected.var.nl" = as.numeric(R.var.vec.simex.nl.NP), "R.corrected.CI.normal.nl" = as.numeric(R.simex.normal.nl), "R.corrected.CI.fieller.nl" = as.numeric(R.var.vec.simex.nl.NP.fieller ),"deltas.corrected.nl" = as.numeric(deltas.vec.simex.nl.NP), "deltas.corrected.var.nl" = as.numeric(deltas.var.vec.simex.nl.NP), "deltas.corrected.CI.normal.nl" = as.numeric(deltas.simex.normal.nl))
		Ronly.list = list("R.corrected.nl" = as.numeric(R.vec.simex.nl.NP), "R.corrected.var.nl" = as.numeric(R.var.vec.simex.nl.NP), "R.corrected.CI.normal.nl" = as.numeric(R.simex.normal.nl), "R.corrected.CI.fieller.nl" = as.numeric(R.var.vec.simex.nl.NP.fieller ))
	}#close nonlinear
	}#close nonparametric
	if(!naive & Ronly ) {return(Ronly.list)}
	if(!naive & !Ronly) {return(full.list)}
	if(naive & Ronly ) {return(c(Ronlynaive.list, Ronly.list))}
	if(naive & !Ronly ) {return(c(fullnaive.list, full.list))}
	
	}  #close function


calculate.var.np = function(s1,s0,y1,y0, extrapolate = TRUE) {
	big.mat.A = c()
	big.mat.B = c()
	big.mat.A.psionly = c()
	big.mat.B.psionly = c()
	pi1 = length(s1)/((length(s1) + length(s0)))
	pi0 = length(s0)/((length(s1) + length(s0)))
	#print("test")
	mu.1.s0 = 	sapply(s1,pred.smooth,zz=s1, y1=y1)
  	if(sum(is.na(mu.1.s0))>0 & extrapolate){
  		#print(paste("j = ",j))
  		#print("step 1")
    	c.mat = cbind(s1, mu.1.s0)
    	for(o in 1:length(mu.1.s0)) {
    		if(is.na(mu.1.s0[o])){
    			distance = abs(s1 - s1[o])
    			#print("step 2")
    			c.temp = cbind(c.mat, distance)
    			#print("step 3")
    			c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where mean is not na
    			new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
    			mu.1.s0[o] = new.est[1]   #in case there are multiple matches
    		}
    	}
		}
	density.s1 = density(s1)
	density.s0 = density(s0)
	r.est = approx(x = density.s0$x, y = density.s0$y, xout = s1, rule = 2)$y/approx(x = density.s1$x, y = density.s1$y, xout = s1)$y
	#STACK A's side by side
	mean1 = mean(y1)
	for(ww in 1:length(y1)) {
		dat.vec = y1
		A2 = (1/pi1)*(dat.vec[ww] - mean1)
		A1 = (1/pi1)*(dat.vec[ww] - mu.1.s0[ww])*r.est[ww]
		A.A.T = c(A1, A2)%*%t(c(A1,A2))
		big.mat.A = cbind(big.mat.A, A.A.T)
		big.mat.A.psionly = cbind(big.mat.A.psionly,c(A1, A2))
	}

	mu.1.s0 = 	sapply(s0,pred.smooth,zz=s1, y1=y1)
  	if(sum(is.na(mu.1.s0))>0 & extrapolate){
  		#print(paste("j = ",j))
  		#print("step 1")
    	c.mat = cbind(s0, mu.1.s0)
    	for(o in 1:length(mu.1.s0)) {
    		if(is.na(mu.1.s0[o])){
    			distance = abs(s0 - s0[o])
    			#print("step 2")
    			c.temp = cbind(c.mat, distance)
    			#print("step 3")
    			c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where mean is not na
    			new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
    			mu.1.s0[o] = new.est[1]   #in case there are multiple matches
    			}
    	}
	}

mean0 = mean(y0)
rr.temp = R.s.estimate(sone = s1, szero = s0, yone = y1, yzero = y0, type = "robust", extrapolate = TRUE, warn.te = TRUE, warn.support = TRUE)
first.term.deltas = rr.temp$delta.s + mean0
for(ww in 1:length(y0)) {
dat.vec = y0
B2 = -(1/pi0)*(dat.vec[ww] - mean0)
B1 = -(1/pi0)*(dat.vec[ww] - mean0 - (mu.1.s0[ww] - first.term.deltas))
B.B.T = c(B1, B2)%*%t(c(B1,B2))
big.mat.B = cbind(big.mat.B, B.B.T)
big.mat.B.psionly = cbind(big.mat.B.psionly, c(B1, B2))
}
total = cbind(big.mat.A,big.mat.B)

return(list("total"=total, "psionly" = cbind(big.mat.A.psionly, big.mat.B.psionly)))
}

me.variance.estimate  = function(replicates){	
	mean.i = apply(replicates,1,mean, na.rm = T)
	num.i = apply(replicates,1,function(x) sum(!is.na(x)))
	var.u = sum((replicates-mean.i)^2, na.rm = T)/sum(num.i)
	return(var.u)
}


R.multiple.surv = function(xone,xzero, deltaone, deltazero, sone, szero, type = 1, t, weight.perturb = NULL, landmark, extrapolate = FALSE, transform = FALSE, conf.int =FALSE, var = FALSE, incremental.value = FALSE, approx = T) {
	delta = delta.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t,approx=approx)$delta	
	delta.s = delta.multiple.surv(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, sone=sone, szero=szero, type = type, t=t, landmark=landmark, extrapolate = extrapolate, transform = transform,approx=approx)
	R.s = 1-delta.s/delta	
	
	if(var | conf.int){
		if(is.null(weight.perturb)) {
			
	weight.perturb = matrix(rexp(500*(length(xone)+length(xzero)), rate=1), ncol = 500)}
		delta.s.p.vec = apply(weight.perturb, 2, delta.multiple.surv, xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, sone=sone, szero=szero, type = type, t=t, landmark=landmark, extrapolate = extrapolate, transform = transform,approx=approx)
		print(weight.perturb[1:10, 1:10])
		delta.p.vec = as.vector(unlist(apply(weight.perturb, 2, delta.surv.estimate, xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, var= FALSE, conf.int = FALSE,approx=approx)))
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
		R.t.est = R.t.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, landmark = landmark,var = var, weight.perturb = weight.perturb,approx=approx)
		delta.t = R.t.est$delta.t
		R.t = R.t.est$R.t
		IV = R.s - R.t
		if(var | conf.int){
		delta.t.p.vec = apply(weight.perturb, 2, delta.t.surv.estimate, xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, landmark=landmark,approx=approx)
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

delta.multiple.surv = function(xone,xzero, deltaone, deltazero, sone, szero, type = 1, t, weight.perturb = NULL, landmark, extrapolate = FALSE, transform = FALSE, approx = T) {
	
	delta.check = delta.surv.estimate(xone=xone,xzero=xzero, deltaone=deltaone, deltazero=deltazero, t=t, conf.int = TRUE,approx=approx)
	if(delta.check$conf.int.quantile[1] < 0 & delta.check$conf.int.quantile[2] > 0 & is.null(weight.perturb)) {print("Warning: it looks like the treatment effect is not significant; may be difficult to interpret the residual treatment effect in this setting")}
	if(is.null(weight.perturb)) {weight = rep(1,length(xone)+length(xzero))}
	if(!is.null(weight.perturb)) {weight = weight.perturb}	
	weight.group1 = weight[1:length(xone)]
	weight.group0 = weight[(1+length(xone)):(length(xone)+length(xzero))]
	#type 1 = 2-stage smooth, 2 = 2-stage weighted, 3 = dr smoothed, 4 = 2-stage model, 5 = weighted only, 6 = dr model based
	if(type == 1) {
		cox.model = coxph(Surv(xone[xone>landmark], deltaone[xone>landmark])~sone[xone>landmark,], weights = 					weight.group1[xone>landmark])
		sone.sub = as.matrix(sone[xone>landmark,])
		szero.sub = as.matrix(szero[xzero>landmark,])
		score.one = sone.sub%*%cox.model$coef
		score.zero = szero.sub%*%cox.model$coef
		mu.1.s0 = pred.smooth.surv(xone.f=xone[xone>landmark], deltaone.f = deltaone[xone>landmark], 							sone.f=as.vector(score.one), szero.one = as.vector(score.zero), myt=t, weight.pred = weight.group1[xone>landmark], 		extrapolate = extrapolate, transform=transform)
		censor0.t = censor.weight(xzero, deltazero, t, weight = weight.group0)
		censor0.landmark = censor.weight(xzero, deltazero, landmark, weight = weight.group0)
		delta.s = sum(weight.group0[xzero>landmark]*mu.1.s0)/sum(weight.group0*censor0.landmark) - 								sum(weight.group0*1*(xzero>t))/sum(weight.group0*censor0.t)
	}
	if(type == 2) {
		treat.ind = c(rep(1,length(xone)), rep(0, length(xzero)))
		x.vector = c(xone, xzero)
		d.vector = c(deltaone, deltazero)
		s.matbig = rbind(sone, szero)
		pi.0 = sum(weight.group0)/(sum(weight.group1) + sum(weight.group0))
		ps.model = glm(1-treat.ind[x.vector > landmark] ~ s.matbig[x.vector>landmark,], family = binomial ,weights = 			weight[x.vector>landmark])
		ps.pred = predict(ps.model, type = "response")
		ps.pred.smooth = ps.pred
		treat.ind.reduce = treat.ind[x.vector>landmark]
		ps.pred.smooth.one = ps.pred[treat.ind.reduce == 1]	
		cox.model = coxph(Surv(xone[xone>landmark], deltaone[xone>landmark])~sone[xone>landmark,], weights = 					weight.group1[xone>landmark])
		sone.sub = as.matrix(sone[xone>landmark,])
		szero.sub = as.matrix(szero[xzero>landmark,])
		score.one = sone.sub%*%cox.model$coef
		score.zero = szero.sub%*%cox.model$coef
		mu.1.s0 = pred.smooth.surv(xone.f=xone[xone>landmark], deltaone.f = deltaone[xone>landmark], 					sone.f=as.vector(score.one), szero.one = as.vector(score.zero), myt=t, weight.pred = weight.group1[xone>landmark], 		extrapolate = extrapolate, transform=transform, ps.weight = ps.pred.smooth.one/(1-ps.pred.smooth.one))
		censor0.t = censor.weight(xzero, deltazero, t, weight = weight.group0)
		censor0.landmark = censor.weight(xzero, deltazero, landmark, weight = weight.group0)
		delta.s = sum(weight.group0[xzero>landmark]*mu.1.s0)/sum(weight.group0*censor0.landmark) - 								sum(weight.group0*1*(xzero>t))/sum(weight.group0*censor0.t)
	}
	if(type ==3 ){
		treat.ind = c(rep(1,length(xone)), rep(0, length(xzero)))
		x.vector = c(xone, xzero)
		d.vector = c(deltaone, deltazero)
		s.matbig = rbind(sone, szero)
		pi.0 = sum(weight.group0)/(sum(weight.group1) + sum(weight.group0))
		ps.model = glm(1-treat.ind[x.vector > landmark] ~ s.matbig[x.vector>landmark,], family = binomial ,weights = 			weight[x.vector>landmark])
		ps.pred = predict(ps.model, type = "response")
		ps.pred.smooth = ps.pred
		treat.ind.reduce = treat.ind[x.vector>landmark]
		censor0.t = censor.weight(xzero, deltazero, t, weight = weight.group0)
		censor0.landmark = censor.weight(xzero, deltazero, landmark, weight = weight.group0)
		censor1.t = censor.weight(xone, deltaone, t, weight = weight.group1)
		censor1.landmark = censor.weight(xone, deltaone, landmark, weight = weight.group1)
		censor.t.vector = censor1.t*treat.ind + censor0.t*(1-treat.ind)
		censor.landmark.vector = censor1.landmark*treat.ind + censor0.landmark*(1-treat.ind)
		#term in brackets
		bracket.term = (1*(treat.ind.reduce == 1)*ps.pred.smooth*censor1.landmark)/((1-ps.pred.smooth)*censor0.landmark) - 		1*(treat.ind.reduce == 0)
		first.term = 1*weight[x.vector>landmark]*(x.vector[x.vector>landmark] > t)/(censor.t.vector[x.vector>landmark]*pi.0)
		s.predictor = as.matrix(sone[xone>landmark,])
		x.adjust = xone[xone>landmark] - landmark
		cox.model = coxph(Surv(x.adjust, deltaone[xone>landmark])~s.predictor, weights = weight.group1[xone>landmark], 			model = TRUE)
		sone.sub = as.matrix(sone[xone>landmark,])
		szero.sub = as.matrix(szero[xzero>landmark,])
		sall.sub = as.matrix(s.matbig[x.vector>landmark,])
		score.one = sone.sub%*%cox.model$coef
		score.all = sall.sub%*%cox.model$coef
		mu.1.s0 = pred.smooth.surv(xone.f=xone[xone>landmark], deltaone.f = deltaone[xone>landmark], 							sone.f=as.vector(score.one), szero.one = as.vector(score.all), myt=t, weight.pred = weight.group1[xone>landmark], 		extrapolate = extrapolate, transform = transform)
		second.term = weight[x.vector>landmark]*mu.1.s0/(censor.landmark.vector[x.vector>landmark]*								pi.0)
		delta.s = (1/sum(weight))*sum(first.term*bracket.term - second.term*bracket.term)	
	}
	if(type == 4) {
		s.predictor = as.matrix(sone[xone>landmark,])
		x.adjust = xone[xone>landmark] - landmark
		cox.model = coxph(Surv(x.adjust, deltaone[xone>landmark])~s.predictor, weights = weight.group1[xone>landmark], 			model = TRUE)
		szero.sub = as.matrix(szero[xzero>landmark,])
		score.zero = exp(szero.sub%*%cox.model$coef)
		baseline.hazard = basehaz(cox.model, centered = FALSE)
		gap.time = t-landmark
		baseline.t <- approx(baseline.hazard$time,baseline.hazard$hazard,gap.time, rule = 2)$y
		mu.1.s0 = exp(-baseline.t*score.zero)
		censor0.t = censor.weight(xzero, deltazero, t, weight = weight.group0)
		censor0.landmark = censor.weight(xzero, deltazero, landmark, weight = weight.group0)
		delta.s = sum(weight.group0[xzero>landmark]*mu.1.s0)/sum(weight.group0*censor0.landmark) - 								sum(weight.group0*1*(xzero>t))/sum(weight.group0*censor0.t)
	}
	if(type == 5) {
		treat.ind = c(rep(1,length(xone)), rep(0, length(xzero)))
		x.vector = c(xone, xzero)
		s.matbig = rbind(sone, szero)
		pi.0 = sum(weight.group0)/(sum(weight.group0) + sum(weight.group1))
		ps.model = glm(1-treat.ind[x.vector > landmark] ~ s.matbig[x.vector>landmark,], family = binomial, weights = 			weight[x.vector>landmark])
		ps.pred = predict(ps.model, type = "response")
		ps.pred.smooth = ps.pred
		treat.ind.reduce = treat.ind[x.vector>landmark]
		censor0.t = censor.weight(xzero, deltazero, t, weight = weight.group0)
		censor0.landmark = censor.weight(xzero, deltazero, landmark, weight = weight.group0)
		censor1.t = censor.weight(xone, deltaone, t, weight = weight.group1)
		censor1.landmark = censor.weight(xone, deltaone, landmark, weight = weight.group1)
    	weight.censor.corrected = censor1.landmark/censor0.landmark
    	treat.term = sum(1*(xone[xone>landmark]>t)*ps.pred.smooth[treat.ind.reduce==1]* 										censor1.landmark*weight.group1[xone>landmark]/ (1-ps.pred.smooth[treat.ind.reduce==1]))/sum(censor1.t* 					pi.0*censor0.landmark*weight)
    	control.term = sum(1*(xzero>t)*weight.group0)/sum(censor0.t*pi.0*weight)
    	delta.s = treat.term-control.term
	}
	if(type == 6) {
		treat.ind = c(rep(1,length(xone)), rep(0, length(xzero)))
		x.vector = c(xone, xzero)
		d.vector = c(deltaone, deltazero)
		s.matbig = rbind(sone, szero)
		pi.0 = sum(weight.group0)/(sum(weight.group1) + sum(weight.group0))
		ps.model = glm(1-treat.ind[x.vector > landmark] ~ s.matbig[x.vector>landmark,], family = binomial ,weights = 			weight[x.vector>landmark])
		ps.pred = predict(ps.model, type = "response")
		ps.pred.smooth = ps.pred
		treat.ind.reduce = treat.ind[x.vector>landmark]
		censor0.t = censor.weight(xzero, deltazero, t, weight = weight.group0)
		censor0.landmark = censor.weight(xzero, deltazero, landmark, weight = weight.group0)
		censor1.t = censor.weight(xone, deltaone, t, weight = weight.group1)
		censor1.landmark = censor.weight(xone, deltaone, landmark, weight = weight.group1)
		censor.t.vector = censor1.t*treat.ind + censor0.t*(1-treat.ind)
		censor.landmark.vector = censor1.landmark*treat.ind + censor0.landmark*(1-treat.ind)
		#term in brackets
		bracket.term = (1*(treat.ind.reduce == 1)*ps.pred.smooth*censor1.landmark)/((1-ps.pred.smooth)*censor0.landmark) - 		1*(treat.ind.reduce == 0)
		first.term = 1*weight[x.vector>landmark]*(x.vector[x.vector>landmark] > t)/(censor.t.vector[x.vector>landmark]*			pi.0)
	
		s.predictor = as.matrix(sone[xone>landmark,])
		x.adjust = xone[xone>landmark] - landmark
		cox.model = coxph(Surv(x.adjust, deltaone[xone>landmark])~s.predictor, weights = weight.group1[xone>landmark], 			model = TRUE)
		sall.sub = as.matrix(s.matbig[x.vector>landmark,])
		score.all = sall.sub%*%cox.model$coef
		baseline.hazard = basehaz(cox.model, centered = FALSE)
		gap.time = t-landmark
		baseline.t <- approx(baseline.hazard$time,baseline.hazard$hazard,gap.time, rule = 2)$y
		cox.prediction = exp(-baseline.t*exp(score.all))
		second.term = weight[x.vector>landmark]*cox.prediction/(censor.landmark.vector[x.vector>landmark]*pi.0)
		delta.s = (1/sum(weight))*sum(first.term*bracket.term - second.term*bracket.term)
	}
	return(delta.s)
}

