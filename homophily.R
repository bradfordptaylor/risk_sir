library(tidyverse)
library(scales)

#plot pars
ggsaveheight <- 10
ggsavewidth <- 1.5*ggsaveheight
ggsaveunits <- "cm"

#pars for
betaLL <- 1
betaHH <- 2
betaHL_wellmixed <- sqrt(betaHH*betaLL)
gamma <- 1
N <- 1
SL <- N
SH <- N
#r1 <- betaLL*SL+betaHL*SH
#r2 <- betaHL*SL+betaHH*SH

#FUNCTIONS
fHstar <- function(betaLL,betaHL,betaHH,gamma,SL,SH){
  a <- betaLL*SL-betaHH*SH + betaHL*(SH-SL)
  b <- betaHH*SH-betaLL*SL-betaHL*(SH+SL)
  c <- betaHL*SH
  fH <- (-b-sqrt(b^2-4*a*c))/(2*a)
  return(fH)
}

betaHLfromfH <- function(fH,betaLL,betaHH,SL,SH){
  betaHLstar <- -fH*(1-fH)*(betaHH*SH-betaLL*SL)/(SH*(1-fH)^2-SL*fH^2)
  return(betaHLstar)
}
betaHHfromfH <- function(fH,betaLL,betaHL,SL,SH){
  betaHH <- (((SH*(1-fH)^2-SL*fH^2)*betaHLstar)/(-fH*(1-fH))+betaLL*SL)/SH
  return(betaHLstar)
}

betaLLrescale <- function(betaLL,betaHH,betaHL,SL,SH){
  currdelta <- betaHL-sqrt(betaLL*betaHH)
  betaLLout <- betaLL-currdelta*SH/SL
  return(betaLLout)
}

betaHHrescale <- function(betaLL,betaHH,betaHL,SL,SH){
  currdelta <- betaHL-sqrt(betaLL*betaHH)
  betaHHout <- betaHH-currdelta*SL/SH
  return(betaHHout)
}

#solve for next population
PyI_full <- function(betaLL,betaHL,betaHH,gamma,SL,SH,IL,IH){
  fH <- IH/(IL+IH)
  #L recover, H recover, L infected, H infected
  Py <- c(gamma*(1-fH),gamma*fH,SL*(fH*betaHL+(1-fH)*betaLL),SH*(fH*betaHH+(1-fH)*betaHL))
  return(Py)
}

#remove extinct epidemics
PyI <- function(betaLL,betaHL,betaHH,gamma,SL,SH,IL,IH){
  Py <- PyI_full(betaLL,betaHL,betaHH,gamma,SL,SH,IL,IH)
  Py <- Py[c((IL==1 & IH>0) | IL>1,(IH==1 & IL>0) | IH>1,TRUE,TRUE)]
  return(Py)
}

PxIL_full <- function(IL){
  currIL <- c(IL-1,IL,IL+1,IL)
  return(currIL)
}
PxIL <- function(IL,IH){
  currIL <- PxIL_full(IL)
  currIL <- currIL[c((IL==1 & IH>0) | IL>1,(IH==1 & IL>0) | IH>1,TRUE,TRUE)]
  return(currIL)
}
PxIH_full <- function(IH){
  currIH <- c(IH,IH-1,IH,IH+1)
  return(currIH)
}
PxIH <- function(IL,IH){
  currIH <- PxIH_full(IH)
  currIH <- currIH[c((IL==1 & IH>0) | IL>1,(IH==1 & IL>0) | IH>1,TRUE,TRUE)]
  return(currIH)
}

Ptbl2tplus1 <- function(lasttbl,betaLL,betaHL,betaHH,gamma,SL,SH,IL,IH){
  currPyfun <- function(IL,IH,py){
    currprepy <- PyI(betaLL,betaHL,betaHH,gamma,SL,SH,IL,IH)
    currPy <- py*currprepy/sum(currprepy)
    return(currPy)
  }
  currPylist <- mapply(currPyfun,lasttbl %>% pull(pIL),lasttbl %>% pull(pIH),lasttbl %>% pull(py))
  currILlist <- mapply(PxIL,lasttbl %>% pull(pIL),lasttbl %>% pull(pIH))
  currIHlist <- mapply(PxIH,lasttbl %>% pull(pIL),lasttbl %>% pull(pIH))
  currtbl <- tibble(pIL=as.vector(unlist(currILlist)),pIH=as.vector(unlist(currIHlist)),py=as.vector(unlist(currPylist))) %>% 
    group_by(pIL,pIH) %>% mutate(py=sum(py)) %>% distinct() %>% ungroup()
  return(currtbl)
}

fHprobdist <- function(t,betaLL,betaHL,betaHH,gamma,SL,SH,IL,IH){
  currPyfun <- function(IL,IH,py){
    currprepy <- PyI(betaLL,betaHL,betaHH,gamma,SL,SH,IL,IH)
    currPy <- py*currprepy/sum(currprepy)
    return(currPy)
  }
  
  lasttbl <- tibble(t=0,pIL=IL,pIH=IH,py=1)
  alltbl <- list(lasttbl)
  for (currn in seq(1,t)){
    #print(currn)
    currPylist <- mapply(currPyfun,lasttbl %>% pull(pIL),lasttbl %>% pull(pIH),lasttbl %>% pull(py))
    currILlist <- mapply(PxIL,lasttbl %>% pull(pIL),lasttbl %>% pull(pIH))
    currIHlist <- mapply(PxIH,lasttbl %>% pull(pIL),lasttbl %>% pull(pIH))
    currtbl <- tibble(t=currn,pIL=as.vector(unlist(currILlist)),pIH=as.vector(unlist(currIHlist)),py=as.vector(unlist(currPylist))) %>% 
      group_by(pIL,pIH) %>% mutate(py=sum(py)) %>% distinct() %>% ungroup()
    alltbl <- append(alltbl,list(currtbl))
    lasttbl <- currtbl
  }
  finaltbl <- bind_rows(alltbl)
  return(finaltbl)
}


getPatn <- function(N,IH0,n0){
  IL0 <- N-IH0
  currfH <- IH0/(IL0+IH0)
  currbetaHL <- betaHLfromfH(currfH,betaLL,betaHH,SL,SH)
  currbetaLL <- betaLLrescale(betaLL,betaHH,currbetaHL,SL,SH)
  currbetaHH <- betaHHrescale(betaLL,betaHH,currbetaHL,SL,SH)
  test <- fHprobdist(n0,currbetaLL,currbetaHL,currbetaHH,gamma,SL,SH,IL0,IH0)
  testfH <- test %>% mutate(fH=pIH/(pIL+pIH)) %>% group_by(t,fH) %>% summarize(py=sum(py))
  return(testfH)
}

getPatn_notequil <- function(N,IH0,n0,startIL,startIH){
  IL0 <- N-IH0
  currfH <- IH0/(IL0+IH0)
  currbetaHL <- betaHLfromfH(currfH,betaLL,betaHH,SL,SH)
  currbetaLL <- betaLLrescale(betaLL,betaHH,currbetaHL,SL,SH)
  currbetaHH <- betaHHrescale(betaLL,betaHH,currbetaHL,SL,SH)
  test <- fHprobdist(n0,currbetaLL,currbetaHL,currbetaHH,gamma,SL,SH,startIL,startIH)
  testfH <- test %>% mutate(fH=pIH/(pIL+pIH)) %>% group_by(t,fH) %>% summarize(py=sum(py))
  return(testfH)
}

checkfH <- function(N,IH0){
  IL0 <- N-IH0
  currfH <- IH0/(IL0+IH0)
  currbetaHL <- betaHLfromfH(currfH,betaLL,betaHH,SL,SH)
  currbetaLL <- betaLLrescale(betaLL,betaHH,currbetaHL,SL,SH)
  currbetaHH <- betaHHrescale(betaLL,betaHH,currbetaHL,SL,SH)
  newfH <- fHstar(currbetaLL,currbetaHL,currbetaHH,gamma,SL,SH)
  print(currfH==newfH)
}
checkfH(1,.7)
mysmoothcount <- function(x,y,sigma0,xdom){
  smoothfun<-function(xfoc,x,y,sigma0){
    newy <- sum(y*exp(-.5*((xfoc-x)/sigma0)^2))
  }
  smoothed <- unlist(lapply(xdom,smoothfun,x,y,sigma0))
  smoothed <- smoothed/sum(smoothed)
  return(smoothed)
}
mymean <- function(fHtbl){
  currmean <- fHtbl %>% group_by(t) %>% summarize(meanpy=sum(fH*py))#sum((fHtbl %>% pull(fH))*(fHtbl %>% pull(py)))
  return(currmean)
}
myvar <- function(fHtbl){
  currvar <- fHtbl %>% group_by(t) %>% summarize(secondpy=sum(fH^2*py)) %>% left_join(mymean(fHtbl)) %>% mutate(varpy=secondpy-meanpy^2) %>% select(-secondpy)#sum((fHtbl %>% pull(fH))^2*(fHtbl %>% pull(py)))-mymean(fHtbl)^2
  return(currvar)
}
#####
#Generate data starting at a certain fraction
Nscale <- 1#10000
maxt <- 200#500#500
allPatn <- lapply(Nscale*seq(11,20),function(x) getPatn(Nscale*20,x,maxt))
allPatn_startIL <- lapply(seq(11,20),function(x) getPatn_notequil(20,x,maxt,1,0))
allPatn_startIH <- lapply(seq(11,20),function(x) getPatn_notequil(20,x,maxt,0,1))

#Distributions at certain periods
currsmoothx <- seq(0,1,by=.001)
currsigma <- .0075
currt <- 200#maxt
currPatn <- allPatn
currsmooth1 <- currPatn[[2]] %>% filter(t==currt)
currsmoothy1 <- mysmoothcount(currsmooth1 %>% pull(fH),currsmooth1 %>% pull(py),currsigma,currsmoothx)
currsmooth2 <- currPatn[[5]] %>% filter(t==currt)
currsmoothy2 <- mysmoothcount(currsmooth2 %>% pull(fH),currsmooth2 %>% pull(py),currsigma,currsmoothx)
currsmooth3 <- currPatn[[7]] %>% filter(t==currt)
currsmoothy3 <- mysmoothcount(currsmooth3 %>% pull(fH),currsmooth3 %>% pull(py),currsigma,currsmoothx)
smoothplot <- bind_rows(list(tibble(fH=currsmoothx,Py=currsmoothy1,homophily=12/20)),list(tibble(fH=currsmoothx,Py=currsmoothy2,homophily=15/20)),list(tibble(fH=currsmoothx,Py=currsmoothy3,homophily=17/20)))

#BLOP
hex <- hue_pal()(3)
ggplot(smoothplot,aes(x=fH,y=Py,color=as.character(homophily))) + 
  geom_line(size=1) + geom_vline(xintercept = 12/20,color=hex[1],linetype = "dashed",size=1,alpha=.5) +
  #geom_vline(xintercept= fHstar_stoch_fracR0(2,12/20,betaLL,betaHH,SL,SH,gamma),color=hex[1],alpha=.5) +
  #geom_vline(xintercept= fHstar_stoch(betaLL,betaHLfromfH(12/20,betaLL,betaHH,SL,SH),betaHH,gamma,SL,SH),color=hex[1],alpha=.5) +
  geom_vline(xintercept = 15/20,color=hex[2],linetype = "dashed",size=1,alpha=.5) +
  #geom_vline(xintercept= fHstar_stoch_fracR0(2,15/20,betaLL,betaHH,SL,SH,gamma),color=hex[2],alpha=.5) +
  geom_vline(xintercept = 17/20,color=hex[3],linetype = "dashed",size=1,alpha=.5) +
  #geom_vline(xintercept= fHstar_stoch_fracR0(2,17/20,betaLL,betaHH,SL,SH,gamma),color=hex[3],alpha=.5) +
  labs(x=expression(paste("equilibrium risk distribution, ",frac(I[H], I[L]+I[H]))),y="kernel density probability",color="homophily") + 
  lims(x=c(.4,1))+ 
  theme_classic() +
  theme(legend.position = "none")
#+ legend()
ggsave("~/projects/covid_spread/figures/masteratt.pdf",width=ggsavewidth,height=ggsaveheight,units = ggsaveunits)
ggsave("~/projects/covid_spread/figures/masteratt.png",width=ggsavewidth,height=ggsaveheight,units = ggsaveunits)

#dyn
ggplot(probreltbl %>% filter(t%%2==0),aes(x=t,y=probincrease/(probincrease+probdecrease),color=as.character(fHstar))) + geom_line(size=1) +
  geom_hline(yintercept=.5,size=1,linetype = "dashed") +
  #labs(x="time (population changes)",y=expression(frac(P(f[2]>f[2]^"*"), P(f[2]>f[2]^"*")+P(f[2]<f[2]^"*"))) ,colour="homophily") +
  labs(x="time (population changes)",y=expression(paste("relative increased sociality, ",rho)) ,colour="homophily") +
  theme_classic() +
  lims(x=c(0,200))
ggsave("~/projects/covid_spread/figures/increasevsdecrease_initfH.pdf",width=ggsavewidth,height=ggsaveheight,units = ggsaveunits)
ggsave("~/projects/covid_spread/figures/increasevsdecrease_initfH.png",width=ggsavewidth,height=ggsaveheight,units = ggsaveunits)

probincreasetbl <- bind_rows(lapply(seq(1,9,by=1),function(x){allPatn[[x]] %>% mutate(fHstar=(10+x)/20) %>% filter(fH>fHstar) %>% group_by(t) %>% summarize(probincrease=sum(py)) %>% mutate(fHstar=(10+x)/20)}))
probdecreasetbl <- bind_rows(lapply(seq(1,9,by=1),function(x){allPatn[[x]] %>% mutate(fHstar=(10+x)/20) %>% filter(fH<fHstar) %>% group_by(t) %>% summarize(probdecrease=sum(py)) %>% mutate(fHstar=(10+x)/20)}))
probincreasetbl_startIL <- bind_rows(lapply(seq(1,9,by=1),function(x){allPatn_startIL[[x]] %>% mutate(fHstar=(10+x)/20) %>% filter(fH>fHstar) %>% group_by(t) %>% summarize(probincrease=sum(py)) %>% mutate(fHstar=(10+x)/20)}))
probdecreasetbl_startIL <- bind_rows(lapply(seq(1,9,by=1),function(x){allPatn_startIL[[x]] %>% mutate(fHstar=(10+x)/20) %>% filter(fH<fHstar) %>% group_by(t) %>% summarize(probdecrease=sum(py)) %>% mutate(fHstar=(10+x)/20)}))
probincreasetbl_startIH <- bind_rows(lapply(seq(1,9,by=1),function(x){allPatn_startIH[[x]] %>% mutate(fHstar=(10+x)/20) %>% filter(fH>fHstar) %>% group_by(t) %>% summarize(probincrease=sum(py)) %>% mutate(fHstar=(10+x)/20)}))
probdecreasetbl_startIH <- bind_rows(lapply(seq(1,9,by=1),function(x){allPatn_startIH[[x]] %>% mutate(fHstar=(10+x)/20) %>% filter(fH<fHstar) %>% group_by(t) %>% summarize(probdecrease=sum(py)) %>% mutate(fHstar=(10+x)/20)}))

probreltbl <- probincreasetbl %>% left_join(probdecreasetbl,by=c("t","fHstar")) %>% mutate(probdevdiff=probincrease-probdecrease)
probreltbl_startIL <- probincreasetbl_startIL %>% left_join(probdecreasetbl_startIL,by=c("t","fHstar")) %>% mutate(probdevdiff=probincrease-probdecrease)
probreltbl_startIH <- probincreasetbl_startIH %>% left_join(probdecreasetbl_startIH,by=c("t","fHstar")) %>% mutate(probdevdiff=probincrease-probdecrease)

#Init L
ggplot(probreltbl_startIL %>% filter(t%%2==0),aes(x=t,y=probincrease/(probincrease+probdecrease),color=as.character(fHstar))) + geom_line(size=1) +
  geom_hline(yintercept=.5,size=1,linetype = "dashed") +
  #labs(x="time (population changes)",y=expression(frac(P(f[2]>f[2]^"*"), P(f[2]>f[2]^"*")+P(f[2]<f[2]^"*"))) ,colour="homophily") +
  labs(x="time (population changes)",y=expression(paste("relative increased sociality, ",rho)) ,colour="homophily") +
  theme_classic() +
  lims(x=c(0,200))
ggsave("~/projects/covid_spread/figures/rho_initnotrisky.png",width=ggsavewidth,height=ggsaveheight,units = ggsaveunits)

#Init L
ggplot(probreltbl_startIL %>% filter(t%%2==0),aes(x=t,y=probincrease/(probincrease+probdecrease),color=as.character(fHstar))) + geom_line(size=1) +
  geom_hline(yintercept=.5,size=1,linetype = "dashed") +
  #labs(x="time (population changes)",y=expression(frac(P(f[2]>f[2]^"*"), P(f[2]>f[2]^"*")+P(f[2]<f[2]^"*"))) ,colour="homophily") +
  labs(x="time (population changes)",y=expression(paste("relative increased sociality, ",rho)) ,colour="homophily") +
  theme_classic() +
  lims(x=c(0,500),y=c(.45,.55))
ggsave("~/projects/covid_spread/figures/rho_initnotrisky_inset.png",width=ggsavewidth,height=ggsaveheight,units = ggsaveunits)

#Init H
ggplot(probreltbl_startIH %>% filter(t%%2==0),aes(x=t,y=probincrease/(probincrease+probdecrease),color=as.character(fHstar))) + geom_line(size=1) +
  geom_hline(yintercept=.5,size=1,linetype = "dashed") +
  #labs(x="time (population changes)",y=expression(frac(P(f[2]>f[2]^"*"), P(f[2]>f[2]^"*")+P(f[2]<f[2]^"*"))) ,colour="homophily") +
  labs(x="time (population changes)",y=expression(paste("relative increased sociality, ",rho)) ,colour="homophily") +
  theme_classic() +
  lims(x=c(0,500))
ggsave("~/projects/covid_spread/figures/rho_initrisky.png",width=ggsavewidth,height=ggsaveheight,units = ggsaveunits)

ggplot(probreltbl_startIL %>% filter(t%%2==0),aes(x=t,y=probincrease/(probincrease+probdecrease),color=as.character(fHstar))) + geom_line(size=1) +
  geom_line(data=probreltbl_startIH %>% filter(t%%2==0),aes(x=t,y=probincrease/(probincrease+probdecrease),color=as.character(fHstar)),linetype = "dashed",size=1) +
  geom_hline(yintercept=.5,size=1,linetype = "dashed") +
  #labs(x="time (population changes)",y=expression(frac(P(f[2]>f[2]^"*"), P(f[2]>f[2]^"*")+P(f[2]<f[2]^"*"))) ,colour="homophily") +
  labs(x="time (population changes)",y=expression(paste("relative increased sociality, ",rho)) ,colour="homophily") +
  theme_classic() +
  lims(x=c(0,200))
ggsave("~/projects/covid_spread/figures/rho_initboth.pdf",width=ggsavewidth,height=ggsaveheight,units = ggsaveunits)
ggsave("~/projects/covid_spread/figures/rho_initboth.png",width=ggsavewidth,height=ggsaveheight,units = ggsaveunits)

#supplementary
#fH equilibrium relative to homophily 
betaHLvec <- seq(0,10,by=.01)
ggplot(tibble(betaHL=betaHLvec,fH=fHstar(betaLL,betaHLvec,betaHH,gamma,SL,SH)),aes(x=betaHL,y=fH)) + geom_line(size=2) + 
  #geom_vline(xintercept=sqrt(betaLL*betaHH),size=1,linetype = "dashed") +
  labs(y=expression(paste("Equilbrium risk distribution, ",frac(I[H], I[L]+I[H]))),x=expression(paste("cross-transmission, ",beta[HL]))) + lims(x=c(0,5),y=c(.5,1))+
  theme_classic(base_size = 14)
ggsave("~/projects/covid_spread/figures/prosocialdistribution.pdf",width=ggsavewidth,height=ggsaveheight,units = ggsaveunits)
ggsave("~/projects/covid_spread/figures/prosocialdistribution.png",width=ggsavewidth,height=ggsaveheight,units = ggsaveunits)

#long dyn
longdyntbl <- getPatn(20,14,1250)
write_tsv(longdyntbl,"~/projects/covid_spread/data/longdyntbl.tsv")
longdyntbl <- read_tsv("~/projects/covid_spread/data/longdyntbl.tsv")
longdyntbl <- longdyntbl %>% rename(fH=f2)
longdyntbl_stats <- myvar(longdyntbl)
ggplot(longdyntbl_stats,aes(x=t,y=varpy)) + geom_line(size=1) +
  labs(y="risk distribution probability distribution variance",x="time (population changes)") +theme_classic()
ggsave("~/projects/covid_spread/figures/longdyn_variance.pdf",width=ggsavewidth,height=ggsaveheight,units = ggsaveunits)
ggsave("~/projects/covid_spread/figures/longdyn_variance.png",width=ggsavewidth,height=ggsaveheight,units = ggsaveunits)

Nbase <- 5
IH0base <- 3
trun <- 75
Nprodvec <- seq(1,30)
getPatn(Nbase,IH0base,trun)
stochextlist <- lapply(Nprodvec,function(x) getPatn(x*Nbase,x*IH0base,trun))
stochextfhstar <- IH0base/Nbase
stochext_grow_t <- lapply(stochextlist,function(x) x %>% group_by(t) %>% summarise(rho = sum(py[fH>stochextfhstar])/sum(py[fH!=stochextfhstar])))
stochext_grow_end <- unlist(lapply(stochextlist,function(x) x %>% filter(t==trun-1) %>% summarise(rho = sum(py[fH>stochextfhstar])/sum(py[fH!=stochextfhstar])) %>% pull(rho)))

ggplot(tibble(N0 = Nbase*Nprodvec,rho=stochext_grow_end),aes(x=N0,y=rho)) + 
  geom_vline(xintercept=trun,linewidth=1,linetype="dashed") +
  geom_point() +
  theme_classic()#scale_y_continuous(trans='log10') + 
ggsave("~/projects/covid_spread/figures/stochextinctionbias_full.pdf",width=ggsavewidth,height=ggsaveheight,units = ggsaveunits)
ggsave("~/projects/covid_spread/figures/stochextinctionbias_full.png",width=ggsavewidth,height=ggsaveheight,units = ggsaveunits)

ggplot(tibble(N0 = Nbase*Nprodvec,rho=stochext_grow_end) %>% filter(N0>75),aes(x=N0,y=rho)) + geom_point() +
  theme_classic()#scale_y_continuous(trans='log10') + 
ggsave("~/projects/covid_spread/figures/stochextinctionbias_noextinction.pdf",width=ggsavewidth,height=ggsaveheight,units = ggsaveunits)
ggsave("~/projects/covid_spread/figures/stochextinctionbias_noextinction.png",width=ggsavewidth,height=ggsaveheight,units = ggsaveunits)

ggplot(tibble(N0 = Nbase*Nprodvec,rho=stochext_grow_end) %>% filter(N0<75),aes(x=N0,y=rho)) + geom_point() +
  theme_classic()#scale_y_continuous(trans='log10') + 
ggsave("~/projects/covid_spread/figures/stochextinctionbias_withextinction.pdf",width=ggsavewidth,height=ggsaveheight,units = ggsaveunits)
ggsave("~/projects/covid_spread/figures/stochextinctionbias_withextinction.png",width=ggsavewidth,height=ggsaveheight,units = ggsaveunits)