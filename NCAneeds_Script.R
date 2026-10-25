
############################# BIENVENUE #############################

############# I SOLEMNLY SWEAR THAT I AM UP TO NO GOOD ##############

#################################
## Loading packages

library(MASS)
library(lavaan)

################################
## The NCArigour function
## Estimates expected necessity effect due to correlation between X and Y
## The following values are plugged into the function:
## n = sample size
## rxy = empirical correlation between X and Y
## dnecc = empirical necessity effect
## ci.req = requested confidence interval (default = 0.95)
## cr_fdh = T if estimating CR-FDH (default) and cr_fdh = F if estimating CE-FDH
## R = number of simulated datasets (default = 1000)
## silent = F if printing results and silent = T if not

NCArigour <- function(n,rxy,dnecc,ci.req=0.95,cr_fdh=T,R=1000,silent=F){
  
  rm <- matrix(c(1,rxy,rxy,1),nrow=2) ## correlation matrix
  
  simres <- unlist(sapply(1:R, function(i){ ## for each simulation
    set.seed(i) ## for reproducibility
    df <- data.frame(mvrnorm(n=n,mu=c(0,0),Sigma=rm,empirical=T)) ## the data
    
    names(df) <- c("x","y") ## naming variables in data
    
    x.max <- max(df$x)
    y.max <- max(df$y)
    x.min <- min(df$x)
    y.min <- min(df$y)
    
    scope <- (x.max-x.min)*(y.max-y.min)
    
    sx <- df$x[order(df$x)] ## sorted X
    sy <- df$y[order(df$x)] ## Y in the order of sorted X
    
    cmy <- cummax(sy) ## cumulative max on the Y-variable
    cp <- ifelse(sy>=cmy,1,0) ## if the Y-value is a ceiling point
    
    cx <- sx[cp==1] ## X-values for ceiling points
    cy <- sy[cp==1] ## Y-values for ceiling points
    
    minx <- which(df$x==min(df$x))
    maxy <- which(df$y==max(df$y))
    same <- sum(minx %in% maxy)
    
    if(same>0) nec <- 0
    
    if(same==0){
      
      if(cr_fdh==F){
        
        numc <- length(cx)
        x.base <- cx[2:numc] - cx[1:(numc-1)]
        y.height <- y.max - cy[1:(numc-1)]
        nec <- sum(x.base*y.height) / scope
      }
      
      if(cr_fdh==T){
        
        fitc <- lm(cy ~ cx)
        ylim <- predict(fitc,newdata=data.frame(cx=x.min))
        xlim <- (y.max-fitc$coefficients[1])/fitc$coefficients[2]
        nec <- (0.5 * (y.max-ylim) * (xlim-x.min)) / scope
      }
    }
    
    return(nec) ## returns the necessity effect
  }))
  
  d <- mean(simres, na.rm=T) ## mean/expected necessity effect across the R effects
  p1 <- (1-ci.req)/2 ## lower percentile
  p2 <- 1-p1 ## upper percentile
  ci <- quantile(simres, probs=c(p1,p2)) ## required percentile confidence interval
  pv <- round(sum(simres>dnecc, na.rm=T)/R,3) ## frequency of effects stronger than the empirical effect
  tx1 <- ifelse(pv<0.05,"null rejected","null retained") ## text if p < 0.05 or not
  tx2 <- ifelse(ci[2]<dnecc,"null rejected","null retained") ## text if empirical effect above range of spur. or not
  bm <- round(max(simres,na.rm=T),3) ## max spurious effect
  
  if(silent==F){ ## if printing the results
    
    cat("\n\n*****",
        "\nd(obs)=",round(dnecc,3),
        "\nmax(sim)=",bm,
        "\np(sim>d(obs))=",pv,
        "\n",tx1,
        "\n**",
        "\nd(sim)=",round(d,3),
        "\nci(percentile)=[",round(ci[1],3),",",round(ci[2],3),"]",
        "\n",tx2,
        "\n*****\n")
  }
  
  res <- c(d,ci,pv) ## what the function returns
  return(res)
}

#################################
## Data

orvar1 <- c("A.a","C.a","R.a","A","C","R","P","S") ## labels, study 1

orvar2 <- orvar3 <- c("A","C","R") ## labels, study 2 and 3

olist <- list(orvar1,orvar2,orvar3) ## list of labels

## js = job satisfaction; wm = work meaningfulness;
## pa = positive affect; we = work engagement

## effects, study 1

d1.js <- c(0.22,0.20,0.16,0.15,0.29,0.20,0,0)
d1.wm <- c(0.13,0.24,0.17,0.15,0.28,0.15,0,0)
d1.pa <- c(0.17,0.30,0.20,0.17,0.31,0.19,0,0)
d1.we <- c(0.25,0.33,0.26,0.22,0.36,0.23,0.07,0.09)

## effects, study 2

d2.js <- c(0.26,0.28,0.29)
d2.wm <- c(0.16,0.20,0.19)
d2.pa <- c(0.15,0.33,0.29)
d2.we <- c(0.21,0.32,0.29)

## effects, study 3

d3.js <- c(0.35,0.40,0.32)
d3.wm <- c(0.19,0.24,0.12)
d3.pa <- c(0.30,0.33,0.23)
d3.we <- c(0.29,0.41,0.21)

## list with all effects

dlist <- list(list(d1.js,d1.wm,d1.pa,d1.we),
              list(d2.js,d2.wm,d2.pa,d2.we),
              list(d3.js,d3.wm,d3.pa,d3.we))

## correlations, study 1

r1.js <- c(0.65,0.36,0.52,0.62,0.46,0.53,0.18,0.37)
r1.wm <- c(0.55,0.35,0.50,0.63,0.44,0.54,0.22,0.42)
r1.pa <- c(0.54,0.45,0.52,0.60,0.49,0.51,0.29,0.43)
r1.we <- c(0.59,0.39,0.52,0.63,0.46,0.52,0.29,0.43)

## correlations, study 2

r2.js <- c(0.71,0.58,0.60)
r2.wm <- c(0.58,0.55,0.56)
r2.pa <- c(0.70,0.55,0.57)
r2.we <- c(0.67,0.54,0.54)

## correlations, study 3

r3.js <- c(0.79,0.73,0.74)
r3.wm <- c(0.82,0.69,0.67)
r3.pa <- c(0.73,0.67,0.66)
r3.we <- c(0.76,0.70,0.67)

## list with all correlations

rlist <- list(list(r1.js,r1.wm,r1.pa,r1.we),
              list(r2.js,r2.wm,r2.pa,r2.we),
              list(r3.js,r3.wm,r3.pa,r3.we))

nvar <- c(373,345,148) ## sample sizes

#################################
## Figure 2

panlab <- c("Job satisfaction", "Work meaningfulness",
            "Positive affect", "Work engagement") ## panel labels

country <- c("UK","China","Australia") ## countries

sf <- 0.8 ## sizing factor

if(dev.cur()==2) dev.off() ## removing earlier plots

par(mar=c(1,1,0,0), oma=c(3,2,2,2)) ## margins

## layout

layout_matrix <- matrix(c(1:12), nrow = 3, byrow = TRUE)

layout(layout_matrix, 
       widths = rep(1, 4), ## same width in all columns
       heights = c(2.4, 1, 1)) ## higher in first row

pl <- 1 ## keeping track of panels

for (j in 1:3){ ## three rows of panels
  
  n <- nvar[j] ## picking sample size
  orvar <- olist[[j]] ## picking labels
  duse <- dlist[[j]] ## picking effects
  ruse <- rlist[[j]] ## picking correlations
  
  yl <- rev(orvar) ## y-labels, reversed
  numef <- length(yl) ## number of effects
  upplim <- numef+0.5 ## adding some space
  
  for(i in 1:4){ ## four columns of panels
    
    dval <- duse[[i]] ## picking effects
    rval <- ruse[[i]] ## picking correlations
    
    plot(c(0,0.5), c(0.5,upplim), type="n", xlab="", ylab="", 
         xaxt="n", yaxt="n",cex.axis=sf, las=1) ## empty plot
    
    axis(1,seq(0,0.5,0.1), labels=F,las=2, cex.axis=sf) ## X-axis
    if(j==3) axis(1,seq(0,0.5,0.1), labels=seq(0,0.5,0.1),las=1, cex.axis=sf) ## X-axis
    axis(2,1:numef, labels=F,las=2, cex.axis=sf) ## Y-axis
    if(i==1) axis(2,1:numef, labels=yl,las=2, cex.axis=sf) ## Y-axis
    
    if(pl==10) mtext("Necessity Effect", side = 1, line=1.5, outer=T, cex=sf) ## adding X-label
    if(j==1) mtext(panlab[i],3, line=0.2, cex=sf) ## panel label
    if(i==4) mtext(country[j],4, line=0.2, cex=sf) ## country
    
    resemp <- t(sapply(1:numef,function(k){ ## for each necessity effect
      
      res <- NCArigour(n=n,rxy=rval[k],dnecc=dval[k],cr_fdh=F,silent=F) ## applying the NCArigour function
      
      yval <- (numef+1)-k ## Y-value (reversed to start on top)
      
      abline(h=yval, col="gray") ## a thin gray line
      
      arrows(res[2],yval,res[3],yval,angle=90,length=0.05, code=3,lwd=2) ## range of spuriousness
      points(dval[k], yval, pch=21, col="black",bg="black") ## the empirical necessity effect
      
      return(res) ## returning the results by the NCArigour function
      
    }))
    pl <- pl+1 ## next panel, please
  }
}

#################################
## Chinese data

n <- 345 ## sample size

rm <- matrix(c( ## correlations
  
  1.00, 0.67, 0.66, 0.71, 0.58, 0.70, 0.67,
  0.67, 1.00, 0.88, 0.58, 0.55, 0.55, 0.54,
  0.66, 0.88, 1.00, 0.60, 0.56, 0.57, 0.54,
  0.71, 0.58, 0.60, 1.00, 0.76, 0.78, 0.77,
  0.58, 0.55, 0.56, 0.76, 1.00, 0.71, 0.74,
  0.70, 0.55, 0.57, 0.78, 0.71, 1.00, 0.90,
  0.67, 0.54, 0.54, 0.77, 0.74, 0.90, 1.00), nrow=7)

rownames(rm) <- colnames(rm) <-
  c("AU","CO","RE","JS","WM","PA","WE")

#################################
## Alternative model

spurmod <- "

## Loadings

NS =~ CO+RE+1*AU
WB =~ WM+1*JS+PA+WE
POS =~ a*NS+a*WB

## Correlated residuals

CO ~~ RE
AU ~~ WM
PA ~~ WE

## (Error) variances

AU ~~ AU
CO ~~ CO
RE ~~ RE
JS ~~ JS
WM ~~ WM
PA ~~ PA
WE ~~ WE

NS ~~ NS
WB ~~ WB
POS ~~ 1*POS

"

fit.spur <- lavaan(spurmod, sample.cov=rm, ## fitting model to Chinese data 
              sample.nobs=n, meanstructure=F)

summary(fit.spur, fit.measures=T, ci=T, standardized=T, rsq=T) ## the results


########################## MISCHIEF MANAGED #########################

############################# AU REVOIR #############################



