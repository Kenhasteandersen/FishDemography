#
# Download and install the FishSizeSpectrum package:
#
library(devtools)
install_github("Kenhasteandersen/FishSizeSpectrum")
library(fishsizespectrum)

source("plottools.R")

colAgebased = "blue"


plotDemography <- function(p=baseparameters(W=5000)) {
  N = spectrum(p, nGrid=10000) # Calculate spectrum
  
  #age = seq(1,20)
  w_at_age = weightatage(p, seq(0, 200, length.out = 5000))
  age_at_w = interp1(w_at_age$w, w_at_age$age, N$w)
  # Recruitment at age 1:
  ix_age_1 = which(age_at_w>=1)[1]
  wRecruitment = N$w[ix_age_1]
  cat("Weight at recruitment: ", wRecruitment, '\n')
  N_at_age2 = N$NprR * N$g
  
  defaultplot(mfcol=c(2,2), oma=c(2,0.4,0.8,3))
  
  #
  # Age based:
  #
  #defaultplotvertical(2)
  # Numbers:
  semilogypanel(xlim=c(0,20), ylim=c(1,1e-3),
                ylab="Numbers (#)", xaxis=FALSE)
  Numbers = N$NprR*N$g
  lines(age_at_w, Numbers/Numbers[ix_age_1] , lw=3, col=colAgebased)
  vline(ageMaturation(p$W,p))
  mtext("Age based",side=top)
  makepanellabel('  a')
  # Biomass:
  semilogypanel(xlim=c(0,20), ylim=c(1,50),
                ylab="Biomass (g)", xlab="Age")
  Biomass = Numbers*N$w
  lines(age_at_w, Biomass/Biomass[ix_age_1], lw=3, col=colAgebased)
  vline(ageMaturation(p$W,p))
  makepanellabel()
  #
  # Make additional size-axis:
  #
  # Axis with ages:
  wa = weightatage(p, ages=seq(0,100,length.out = 100))
  weights = c(10,1000,3000,4000,4700)
  ages = interp1(wa$w, wa$age, weights)
  axis(side=bottom, line=2.3,
       at=ages,
       labels=FALSE)
  mtext("Weight (g)", side=bottom, line=3, at=13, adj=1)
  mtext(side=bottom, 
        at = ages,
        line=2.3,
        text=weights)
  #
  # Size-based:
  #
  # Number spectrum:
  #loglogpanel(xlim=c(wRecruitment, p$W), ylim=c(1,1e-4),
  #            ylab="Number spectrum (#/g)", xaxis=FALSE)
  ylim=c(1e-4,1)
  xlim=c(wRecruitment, p$W)
  plot(1, type='n', log='xy',
       ylim=ylim, 
       xlim=xlim, axes=FALSE, xlab='',ylab='', par(new=FALSE), xaxs='r')
  logaxes(bottom, lim=xlim, bExponential = TRUE, labels=FALSE, pow=NA)
  logaxes(right, lim=ylim, bExponential = TRUE, labels=TRUE, pow=NA)
  box(lwd=axis.lwd)
  mtext(side=right, line=1.5, "Number spectrum (#/g)")
  
  lines(N$w, N$NprR/N$NprR[ix_age_1], lw=3)
  vline(p$etaM*p$W)
  mtext("Size based",side=top)
  makepanellabel('     c')
  
  # Sheldon biomass:
  #loglogpanel(xlim=c(wRecruitment, p$W), ylim=c(1,1000),
  #            xlab="Weight (g)", ylab="Sheldon spectrum (g/g)")
  ylim=c(1,1000)
  xlim=c(wRecruitment, p$W)
  plot(1, type='n', log='xy',
       ylim=ylim, 
       xlim=xlim, axes=FALSE, xlab='',ylab='', par(new=FALSE), xaxs='r')
  logaxes(bottom, lim=xlim, bExponential = TRUE, labels=TRUE, pow=NA)
  logaxes(right, lim=ylim, bExponential = TRUE, labels=TRUE, pow=NA)
  box(lwd=axis.lwd)
  mtext(side=right, line=1.5, "Sheldon spectrum (g)")
  mtext(side=bottom, line=1.1, "Weight (g)")
  
  Sheldon = N$NprR*N$w^2
  lines(N$w, Sheldon/Sheldon[ix_age_1], lw=3)
  vline(p$etaM*p$W)
  mtext("Weight (g)     ", side=bottom, line=0, adj=1, at=1e-3)
  makepanellabel()
  #
  # Make additional age-axis:
  #
  # Axis with ages:
  wa = weightatage(p, ages=seq(0,20))
  axis(side=bottom, line=2.3,
       at=wa$w[2:length(wa$w)],
       labels=FALSE, col=colAgebased)
  mtext("Age      ", side=bottom, line=3, at=6e2, adj=1, col="blue")
  # Remove superfluous age-labels:
  wa$w10 = log10(wa$w)
  mindist = (max(wa$w10)-min(wa$w10))/30
  ix = 2
  while (ix<dim(wa)[1]) {
    if ( (wa$w10[ix]-wa$w10[ix-1]) < mindist) {
      wa = rbind( wa[1:ix-1,], wa[(ix+1):dim(wa)[1],] )
    } else
      ix = ix + 1
  }
  wa = wa[1:dim(wa)[1]-1,]
  mtext(side=bottom, 
        at = wa$w[2:length(wa$w)],
        line=2.3,
        text=wa$age[2:length(wa$w)], col=colAgebased)
}

plotGrowth = function() {
  # Rip data from Ursin (1967) "Mathematical model of fish growth", figure App. 13:
  l = 0.02*c(44, 53, 59, 66, 71, 78, 83, 87, 92, -1, -1, -1, 108,
             108,116,120, -1, 128, -1, 139, 141, 143, 148, 149, -1, 156, 157, 160, 173, 172, 172, 176, 178, 
             183, 192, 188, 193, 189) # In cm
  age = c(seq(0,26), seq(28,48,by = 2))
  ix = which(l!=-0.02) # Sort out missing data
  Ursin = data.frame(l=l[ix],age=age[ix])
  #
  # Set up plot
  #
  defaultplot()
  defaultpanel(xlim=c(0,48), ylim=c(0,4),
               xlab="Age (weeks)", ylab="Length (cm)")
  #
  # Fit and plot standard von B curve:
  #
  fit_vonB = nls(l ~ L*(1-exp(-K*(age-t0))), data=Ursin, start=list(K=0.1,L=4,t0=-4), 
                 lower=c(0,3,-10),upper=c(3,5,10), algorithm = "port")
  lines(Ursin$age, predict(fit_vonB, Ursin$age), lwd=3)
  #
  # Fit and plot bi-phasic growth:
  #
  # Biphasic growth function:
  laa = function(W,A,w0,age) {
    p = baseparameters(W=W)
    p$w0 = w0
    p$A = A
    return( weight2length( weightatage(p, ages=age)$w ))
  }
  fit_biphasic = nls(l ~ laa(W,A,0.008,age),
                     data=Ursin, start=list(W=1, A=0.1), 
                     lower=c(0,0),upper=c(2,0.3), algorithm = "port")
  lines(Ursin$age, predict(fit_biphasic, Ursin$age), lwd=3, col=colAgebased)
  # Plot data:
  points(Ursin$age, Ursin$l, ylim=c(0.5,4))
} 

plotPredictions = function(p=baseparameters(), n=25) {
  defaultplot(mfcol=c(1,3))
  W = 10^seq(0,6,length.out = n)  # Asymptotic sizes to run over
  A = c(0.5, 1, 2)*p$A # Growth rates
  #
  # Make predictions:
  #
  Fmsy = matrix(nrow=length(W), ncol=length(A))
  rmax = matrix(nrow=length(W), ncol=length(A))
  dwmat_dt = matrix(nrow=length(W), ncol=length(A))
  
  for (j in 1:length(A)) {
    for (i in 1:length(W)) {
      p = baseparameters(W[i], A[j])
      Fmsy[i,j] = calcRefpoints(p)$Fmsy
      rmax[i,j] = calc_rana1(W[i], p)
    }
    dwmat_dt[,j] = calcSelectionResponse(baseparamQG(0, NULL, p), W, F=0.3)$dwmdt
  }
  #
  # Plot F_msy
  #
  semilogxpanel(xlim=W, ylim=c(0,1.55),
                xlab="", ylab="\\textit{F}_{msy} (yr^{-1})")
  lines(W, Fmsy[,1], lw=1, col=stdgrey)
  lines(W, Fmsy[,2], lw=3)
  lines(W, Fmsy[,3], lw=2, col=stdgrey)
  makepanellabel()
  #
  # Plot r_max
  #
  semilogxpanel(xlim=W, ylim=c(0,1.55),
                xlab="Asymptotic weight (g)", ylab="Pop. growth, \\textit{r}_{max} (yr^{-1})")
  lines(W, rmax[,1], lw=1, col=stdgrey)
  lines(W, rmax[,2], lw=3)
  lines(W, rmax[,3], lw=2, col=stdgrey)
  makepanellabel()
  #
  # Plot FIE:
  #
  semilogxpanel(xlim=W, ylim=c(-0.0006, 0.0),
                xlab="", ylab="Evolution of \\textit{w}_{\\mathrm{mat}} (yr^{-1})")
  lines(W, dwmat_dt[,1], lw=1, col=stdgrey)
  lines(W, dwmat_dt[,2], lw=3)
  lines(W, dwmat_dt[,3], lw=2, col=stdgrey)
  makepanellabel()
  #
  # Add annotation
  #
  #arrows(3, -4, 3,-2)
}

plotTraitsFishbase <- function() {
  #
  # Download fishbase
  #
  #require(rfishbase)
  #require(taxize)
  #if (!file.exists('../data/fishbasedata.Rdata')) {
  #  fish_list <- species_list()
  #  flatfish_list <- species_list(Order="Pleuronectiformes")
  #  elasmobranch_list <- species_list(Class="Elasmobranchii")
  #  acti_list <- species_list(Class="Actinopterygii")
  #  clupeid_list <- species_list(Order="Clupeiformes")
  #  fish <- species(fish_list)
  #  growth <- popgrowth(fish_list)
  #  
  #  save(
  #    fish_list, flatfish_list, clupeid_list, elasmobranch_list, acti_list, fish, growth, 
  #    file='../data/fishbasedata.Rdata')
  #}
  #
  # .. or load it:
  #
  load('../data/fishbasedata.Rdata')
  # Keep only good quality data
  growth <- growth[which(growth$to>-1 & growth$to<1),]
  
  # create indices for: teleosts and elasmobranchs
  ixActi <- ( growth$sciname %in% acti_list )
  ixElasmo <- which( growth$sciname %in% elasmobranch_list )
  
  Loo <- growth$Loo
  Woo <- growth$Winfinity
  A <- calcA(growth$K, Loo)
  cat("Geom mean A: ", exp(mean(log(A[ixActi]))),".\n")
  
  popgrowth <- function(W,A,epsR,w0) {
    n = 3/4
    a = baseparameters()$a
    epsEgg = baseparameters()$epsEgg
    A*(1-n) * (W^(1-n) - w0^(1-n))^(-1) * ( (1-a)*log(W/w0)+log(epsEgg*epsR))
  }
  
  # Growth for fish assuming Woffspring = 1 mg
  r <- popgrowth(weight(Loo), A, baseparameters()$epsR, 0.001)
  # Growth for Elasmobranchs assuming Woffspring = 0.0044*Winf and eR=0.2:
  r[ixElasmo] <- popgrowth(weight(Loo[ixElasmo]), A[ixElasmo], 0.2, 0.0044*weight(Loo[ixElasmo]))
  # Set those with negative growth to a small number, just to have the points in:
  r[r<0] <- 0.1
  #
  # Trait space
  #
  defaultplot()
  factor = 1
  #loglogpanel(xlim=c(4,500), ylim=c(0.5, 150), bExponential=FALSE,
  Woo = weight(Loo)
  loglogpanel(xlim=c(.7,1e6), ylim=c(0.5, 150), bExponential=TRUE,
              xlab='Asymptotic weight $\\textit{W}_{\\infty}$ (g)', 
              ylab='Growth rate, $\\textit{A}$ (g$^{0.25}yr^{-1})')
  points(Woo[ixActi], A[ixActi], pch=dots, cex=factor*sqrt(r[ixActi]), col=gray(0.3, alpha=0.4))
  #points(Woo[ixElasmo], A[ixElasmo], pch=triangles, cex=2*factor*sqrt(r[ixElasmo]), col="white" )
  #points(Woo[ixElasmo], A[ixElasmo], pch=triangles, cex=factor*sqrt(r[ixElasmo]), col="black" )
  
  
  # # Highlight Coryphaena hippurus:
  # ixCoryphaena <- growth$sciname == "Coryphaena hippurus"
  # points(Woo[ixCoryphaena], A[ixCoryphaena], pch=dots, cex=factor*sqrt(r[ixCoryphaena]), col="blue")
  # 
  # # Highlight Cod
  # ixCod <- growth$sciname == "Gadus morhua"
  # points(Woo[ixCod], A[ixCod], pch=dots, cex=sqrt(r[ixCod]), col="orange")
  # 
  # # Highlight tiger shark
  # ixTiger <- growth$sciname == "Galeocerdo cuvier"
  # points(Woo[ixTiger], A[ixTiger], pch=dots, cex=sqrt(r[ixTiger]), col="yellow")
  # 
  # # Highlight herring
  # #ixHerring <- growth$sciname == "Clupea harengus"
  # #points(Loo[ixHerring], A[ixHerring], pch=21, cex=sqrt(r[ixHerring]), col="green", bg=rgb(0,1,0) )
  # 
  # # Highlight Gobies
  # goby_list <- species_list(Family="Gobiidae")
  # ixGoby <- which( growth$sciname %in% goby_list )
  # points(Woo[ixGoby], A[ixGoby], pch=dots, cex=factor*sqrt(r[ixGoby]), col="green" )
  # 
  # # Highlight Gasterosteidae (sticklebacks)
  # stickle_list <- species_list(Family="Gasterosteidae")
  # ixStickle <- which( growth$sciname %in% stickle_list )
  # points(Woo[ixStickle], A[ixStickle], pch=dots, cex=factor*sqrt(r[ixStickle]), col="magenta" )
  # 
  # # highlight rockfish
  # rockfish_list <- species_list(Family="Sebastidae")
  # ixRockfish <- which( growth$sciname %in% rockfish_list )
  # points(Woo[ixRockfish], A[ixRockfish], pch=dots, cex=factor*sqrt(r[ixGoby]), col="red" )
  
  
  # Distribution of growth:
  out=density(log(A[ixActi]))
  defaultpanel(log(c(0.7,1e6)), ylim=log(c(0.5, 150)), new=TRUE, xaxis = FALSE, yaxis=FALSE)
  polygon(x=log(1.75e6)-1.5*(c(out$y, 0*out$y)), y=c(out$x, out$x[seq(length(out$x),1,by=-1)]),
          border=NA, col=stdgrey)
  
  # Distribution of Winf
  #out=density(log(Woo[ixActi]))
  #polygon(x=(c(out$x,out$x[1])), y=2*c(out$y, out$y[1])-0.92, border=NA, col=stdgrey)
  
  loglogpanel(xlim=c(.7, 1e6), ylim=c(0.5, 150), bExponential=TRUE, new=TRUE)
}

plotEmergentDensityDependence = function(W=5000, n=25) {
  
  p_std = paramConsumerResourcemodel(W, 1e-8) # std. model driven by SR relation
  r_std = runspectrum(p_std)
  
  p_comp = paramConsumerResourcemodel(W, 100) # Fully driven by competition
  p_comp$thetaFish = 0 # No cannibalism
  r_comp = runspectrum(p_comp)
  
  p_can = paramConsumerResourcemodel(W, 100) # mainly driven by cannibalism
  r_can = runspectrum(p_can)
  
  defaultplot(mfcol=c(4,1))
  #
  # Resource
  #
  wrange = c(0.1,W)
  loglogpanel(xlim = wrange, 
              ylim=c(1000, 1e5),
              ylab = "Resource (g)",
              xlab = "Weight (g)", label=TRUE)

  wR = r$resource$wR
  NR = r$resource$NR/r$N[r$nSave,1,1]*(wR/p$w0)^2
  KR = p$KR * wR^(p$lambdaR)/r$N[r$nSave,1,1]*(wR/p$w0)^2
  #ribbon(wR, NR, KR)
  lines(wR, NR, lwd=2, col=black)
  lines(wR, KR, lwd=2, col=black, lty=dashed)
  
  #
  # Spectra
  #
  wrange = c(0.1,W)
  loglogpanel(xlim = wrange, 
              ylim=c(.5, 1e6),
              ylab = "Spectrum (g)",
              xlab = "Weight (g)", label=TRUE)
 
  
  w = r_std$fish$w
  Nstd = r_std$N[r$nSave,1,]/r_std$N[r$nSave,1,1]*(w/p_std$w0)^2
  N = r_comp$N[r_comp$nSave,1,]/r_comp$N[r$nSave,1,1]*(w/p_std$w0)^2
  #N_can = r_can$N[r$nSave,1,]/r_can$N[r$nSave,1,1]*(w/p_std$w0)^2
  #ribbon(w, ymax=N, ymin=Nstd)
  lines(w, Nstd, lty=dashed, lwd=2)
  lines(w, N, lwd=2)
  #lines(w, N_can, lwd=2, col=stdgrey)
  
  legend("bottomright", 
         legend=c("Early-life density dep.", "Competition density dep."),
                  #"Cannibalism density dep."), 
         lwd=2, lty=c(dashed,1,1),
         bty="n", col=c(black,black,stdgrey))
  #
  # Weight at age.
  #
  defaultpanel(xlim = c(0,30),
               ylim = c(0,W),
               xlab = "Age",
               ylab = "Weight (g)",
               label = TRUE)
  
  wa_std = calcWeightAtAge(p_std, r_std) # calc without DD
  wa_comp = calcWeightAtAge(p_comp, r)
  wa_can = calcWeightAtAge(p_can, r_can)
  
  #ribbon(wa$ages, wa_std$w, wa$w)
  lines(wa_std$ages, wa_std$w, lwd=2, lty=dashed)
  lines(wa_comp$ages, wa_comp$w, lwd=2)
  #lines(wa_can$ages, wa_can$w, lwd=2, col=stdgrey)
  #
  # Yield curves:
  #
  F = seq(0,1.5,length.out = n)
  Yield_std = 0*F
  Yield_comp = 0*F
  Yield_can = 0*F
  for (i in 2:length(F)) {
    # Early-life DD:
    p_std$F = F[i]
    r_std = runspectrum(p_std, r_std)
    Yield_std[i] = trapz(r_std$w, r_std$N[100,1,]*r_std$w*r_std$muF)
    
    # Competition DD:
    p_comp$F = F[i]
    r_comp = runspectrum(p_comp, r_comp)
    Yield_comp[i] = trapz(r_comp$w, r_comp$N[100,1,]*r_comp$w*r_comp$muF)
    
    # Competition DD:
    p_can$F = F[i]
    r_can = runspectrum(p_can, r_can)
    Yield_can[i] = trapz(r_can$w, r_can$N[100,1,]*r_can$w*r_can$muF)
  }
  
  defaultpanel(xlim = F,
               ylim = c(0,1),
               xlab = "Fishing mortality (yr$^{-1})",
               ylab = "Relative yield",
               label = TRUE)
  lines(F, Yield_std/max(Yield_std), lwd=2, lty=dashed)
  lines(F, Yield_comp/max(Yield_comp), lwd=2)
  #lines(F, Yield_can/max(Yield_can), lwd=2, col=stdgrey)
}



plotAll = function() {
  pdfplot('../demography.pdf',plotDemography, width=0.8*doublewidth, height=0.8*2*height)
  pdfplot('../growth.pdf', plotGrowth, width=singlewidth, height=height)
  pdfplot('../predictions.pdf', plotPredictions, width=doublewidth, height=0.8*height)
  pdfplot('../traitspace.pdf', plotTraitsFishbase, width=1.5*singlewidth, height=height)
  pdfplot('../emergent.pdf', plotEmergentDensityDependence, width=singlewidth, height=3*height)
}