#########################################################
## [1] Function cal_hier_score26 -- calculate the hierarchical scores for L=2,3,4,5,6

cal_hier_score26<-function(data, kmax=10000, ptim=100, anneal.coeff=1e-6,  myoutf = "result.txt")
{
	mytf = unique(c(data[,1], data[,2]))

 	conOut = file(myoutf, "w")
 	close(conOut)

	for(lev in 2:6)
	{

	 	cat("\n")
 		res = matrix(0, length(mytf), lev)
		row.names(res) = mytf
 		colnames(res) = 1:lev
 		hs = rep(0, ptim)


	 	for(p in 1:ptim)
 		{
	## Initiation
 			cat("\r\r\r", lev, "->", p)
 			tfs = sample(1:lev, length(mytf), replace = T)
 			names(tfs) = mytf
 			tmp1 = tfs[data[,1]]
 			tmp2 = tfs[data[,2]]
 			xx = sum(tmp1==tmp2)
 			emax = (sum(tmp1>tmp2)+xx)/(sum(tmp1<tmp2)+xx)
 			if(emax==Inf)
 			{
 				emax = 999999
 			}


 			for(k in 1:kmax)
 			{
 				while(1)
 				{
 					tmp.tf = sample(mytf, 1)
 					if(sum(tfs==tfs[tmp.tf])>1) break
 				}
 				xx 	= 1:lev
 				old.s = tfs[tmp.tf]
 				xx 	= xx[xx!=old.s]
 				tfs[tmp.tf] = sample(xx, 1)
 				tmp1 = tfs[data[,1]]
 				tmp2 = tfs[data[,2]]
 				xx = sum(tmp1==tmp2)
 				etmp = (sum(tmp1>tmp2)+xx)/(sum(tmp1<tmp2)+xx)
 				if(etmp==Inf)
 				{
 					etmp = 999999
 				}
 				if(etmp>emax)
 				{
 					emax = etmp
 				}else
				{
					tt = anneal.coeff*k
					prob = exp((etmp-emax)/tt)
					if(prob<0.9)
					{
						tfs[tmp.tf] = old.s
					}else
					{
						emax = etmp
					}
				}
 			}
 			for(i in 1:length(mytf))
 			{
 				res[i, tfs[i]]= res[i, tfs[i]] + 1
 			}
 			hs[p] = emax
 		}
 		res = res/ptim

 		max.lev = apply(res, 1, which.max)
 		tmp1 = max.lev[data[,1]]
 		tmp2 = max.lev[data[,2]]
 		xx = sum(tmp1==tmp2)
 		hs1 = (sum(tmp1>tmp2)+xx)/(sum(tmp1<tmp2)+xx)
 		xx = 0
 		for(i in 1:lev)
 		{
 			tmp = sum(max.lev==i)
 			xx = xx+ tmp*(tmp-1)
 		}
 		tmp = length(mytf)
 		ww = xx/tmp/(tmp-1)
 		xx = sum(tmp1==tmp2)*(1-ww)/ww/2
 		hs2 = (sum(tmp1>tmp2)+xx)/(sum(tmp1<tmp2)+xx)

 		sco = matrix(0, nrow(data), 3)
 		xx1 = res[data[,1], ]
 		xx2 = res[data[,2], ]
 		for(ss in 2:lev)
 			for(tt in 1:(ss-1))
 			{
 				sco[, 1] = sco[, 1] + xx1[,ss]*xx2[,tt]
 			}
 		for(ss in 1:(lev-1))
 			for(tt in (ss+1):lev)
 			{
 				sco[, 3] = sco[, 3] + xx1[,ss]*xx2[,tt]
 			}
 		for(ss in 1:lev)
 		{
 				sco[, 2] = sco[, 2] + xx1[,ss]*xx2[,ss]
 		}
 		xx = apply(sco, 2, sum)
 		hs3 = (xx[1]+xx[2])/(xx[3]+xx[2])
 		conOut = file(myoutf, "a")
 		curLine = paste("Lev=", lev, sep="")
 		curLine = paste(curLine, "\tHS=", hs1, "\tCHS=", hs2, "\tPHS=", hs3, sep="")
 		writeLines(curLine, conOut)
 		close(conOut)
 		write.table(res, myoutf, sep="\t", quote=F, append=T)
	}
}


#########################################################
## [2] Function Opt_anneal.coeff3 -- Plot HS vs. anneal.coeff

Opt_anneal.coeff3<-function(data, kmax=10000, anneal.coeff=(1:10)*1e-7, lev=3, myoutf = "result.txt")
{
	mytf = unique(c(data[,1], data[,2]))
	hs = rep(0, length(anneal.coeff))
	for(p in 1:length(anneal.coeff))
	{
		## Initiation
		cat("\r\r", anneal.coeff[p])
		tfs = sample(1:lev, length(mytf), replace = T)
		names(tfs) = mytf
		tmp1 = tfs[data[,1]]
		tmp2 = tfs[data[,2]]
		xx = sum(tmp1==tmp2)
		emax = (sum(tmp1>tmp2)+xx)/(sum(tmp1<tmp2)+xx)
		if(emax==Inf)
		{
			emax = 999999
		}

		for(k in 1:kmax)
		{
			while(1)
			{
				tmp.tf = sample(mytf, 1)
				if(sum(tfs==tfs[tmp.tf])>1) break
			}
			xx 	= 1:lev
			old.s = tfs[tmp.tf]
			xx 	= xx[xx!=old.s]
			tfs[tmp.tf] = sample(xx, 1)
			tmp1 = tfs[data[,1]]
			tmp2 = tfs[data[,2]]
			xx = sum(tmp1==tmp2)
			etmp = (sum(tmp1>tmp2)+xx)/(sum(tmp1<tmp2)+xx)
			if(etmp==Inf)
			{
				etmp = 999999
			}
			if(etmp>emax)
			{
				emax = etmp
			}else
			{
				tt = anneal.coeff[p]*k
				prob = exp((etmp-emax)/tt)
				if(prob<0.9)
				{
					tfs[tmp.tf] = old.s
				}else
				{
					emax = etmp
				}
			}
		}
		hs[p] = emax
	}
	res = cbind(anneal.coeff, hs)
	write.table(res, myoutf, sep="\t")
	plot(-log10(res[,1]), res[,2], ylim=c(0, max(hs)), xlab="-log10(anneal.coeff)", ylab="Hierarchical score (HS)", main="", type="b")
}

#########################################################
## [3] Function Opt_Kmax3 -- Plot HS vs. kmax

Opt_kmax3<-function(data, kmax=10000, anneal.coeff=1e-6, lev=3, myoutf = "result.txt")
{
	mytf = unique(c(data[,1], data[,2]))
	hs = rep(0, kmax)

	## Initiation
	tfs = sample(1:lev, length(mytf), replace = T)
	names(tfs) = mytf
	tmp1 = tfs[data[,1]]
	tmp2 = tfs[data[,2]]
	xx = sum(tmp1==tmp2)
	emax = (sum(tmp1>tmp2)+xx)/(sum(tmp1<tmp2)+xx)
	if(emax==Inf)
	{
		emax = 999999
	}

	for(k in 1:kmax)
	{
		if(k%%100==0)
		{
			cat("\r\r", k)
		}
		while(1)
		{
			tmp.tf = sample(mytf, 1)
			if(sum(tfs==tfs[tmp.tf])>1) break
		}
		xx 	= 1:lev
		old.s = tfs[tmp.tf]
		xx 	= xx[xx!=old.s]
		tfs[tmp.tf] = sample(xx, 1)
		tmp1 = tfs[data[,1]]
		tmp2 = tfs[data[,2]]
		xx = sum(tmp1==tmp2)
		etmp = (sum(tmp1>tmp2)+xx)/(sum(tmp1<tmp2)+xx)
		if(etmp==Inf)
		{
			etmp = 999999
		}
		if(etmp>emax)
		{
			emax = etmp
		}else
		{
			tt = anneal.coeff*k
			prob = exp((etmp-emax)/tt)
			if(prob<0.9)
			{
				tfs[tmp.tf] = old.s
			}else
			{
				emax = etmp
			}
		}
		hs[k] = emax
	}

	res = cbind(1:kmax, hs)
	write.table(res, myoutf, sep="\t")
	plot(log10(res[,1]), res[,2], ylim=c(0, max(hs)), xlab="log10(kmax)", ylab="Hierarchical score (HS)", main="", type="b")
}
