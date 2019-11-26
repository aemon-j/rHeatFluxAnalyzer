#' Calculate sensible and latent heat flux
#'
#' @param ts surface temperature, degrees C.
#' @param  Uz wind speed, m/s.
#' @param  ta air temperature, degrees C.
#' @param  rh relative humidity, \%.
#' @param  hu height of wind measurement, m.
#' @param  ht height of air temperature measurement, m.
#' @param  hq heigh of humidity measurement, m.
#' @param  alt altitude of lake, m.
#' @export
#' @source Zeng, X., Zhao, M., Dickson, R.E. 1998. Intercomparison of Bulk
#'  aerodynamic algorithms for the computation of sea surface fluxes using
#'  TOGA COARE and TAO data. Journal of Climate 11: 2628-2644.


sens_latent <- function(ts,Uz,ta,rh,hu,ht,hq,alt,lat){

  # hardcoding limits seems critical to me...
  Uz[Uz < 0.2] <- 0.2

  # pre-define arrays
  tstar <- complex(length(Uz))
  qstar <- complex(length(Uz))
  u10 <- complex(length(Uz))
  t10 <- complex(length(Uz))
  q10 <- complex(length(Uz))
  wc <- complex(length(Uz))
  ts <- complex(real = ts)
  Uz <- complex(real = Uz)
  ta <- complex(real = ta)
  rh <- complex(real = rh)
  hu <- complex(real = hu)
  ht <- complex(real = ht)
  hq <- complex(real = hq)
  alt <- complex(real = alt)
  lat <- complex(real = lat)
  init_wnd <- Uz # initial wind speed

  # define constants
  const_vonKarman <- 0.41 # von Karman constant
  const_gas <- 287.1 # gas constant for dry air J kg-1 K-1
  const_SpecificHeatAir <- 1006 # Specific heat capacity of air, J kg-1 K-1
  const_Charnock <- 0.013 # charnock constant
  # gravitational acceleration, m/s2
  #const_Gravity <- 9.78033*(1+(0.0053*(sin(abs(lat))^2) - 5.8e-6*(sin(abs(lat))^2)))
  # gravitational acceleration, m/s2
  const_Gravity <- 9.780310*(1+(0.00530239*(sin(abs(lat))^2) -
                                  0.00000587*(sin(abs(2*lat))^2) - (31.55e-8)*alt))

  # calculate air pressure from altitude
  press <- 101325*(1 - 2.25577e-5*alt)^5.25588 # Pa
  press <- press/100 # mb

  # calculate humidity values
  e_s <- 6.11*exp(17.27*ta/(237.3 + ta)) # saturated vapour pressure at ta, mb
  e_a <- rh*e_s/100 # vapour pressure, mb
  q_z <- 0.622*e_a/press # specific humidity, kg kg-1
  e_sat <- 6.11*exp(17.27*ts/(237.3 + ts)) # saturated vapour pressure at ts, mb
  q_s <- 0.622*e_sat/press # humidity at saturation, kg kg-1

  # calculate gas constant for moist air, J kg-1 K-1
  R_a <- 287*(1 + 0.608*q_z)

  # calculate latent heat of vaporization, J kg-1
  xlv <- 2.501e6 - 2370*ts

  # calculate air density, kg/m3
  rho_a <- 100*press/(R_a*(ta + 273.16))

  # kinematic viscosity, m2 s-1
  KinV <- (1/rho_a)*(4.94e-8*ta + 1.7184e-5)

  # calculate virtual air temperature, K
  t_virt <- (ta + 273.16)*(1 + 0.61*q_z)

  # estimate initial values of u* and zo
  ustar <- Uz*sqrt(0.00104+0.0015/(1+exp((-Uz+12.5)/1.56)))
  zo <- (const_Charnock*ustar^2/const_Gravity) + (0.11*KinV/ustar)
  zo_prev <- zo*1.1
  for (i in 1:length(Uz)){
    while (abs(zo[i] - zo_prev[i])/abs(zo_prev[i]) > 0.00001){
      ustar[i] <- const_vonKarman*Uz[i]/(log(hu/zo[i]))
      dummy <- zo[i]
      zo[i] <- (const_Charnock*ustar[i]^2/const_Gravity) + (0.11*KinV[i]/ustar[i])
      zo_prev[i] <- dummy
    }
  }

  # calculate neutral transfer coefficients
  C_DN <- (ustar^2)/(Uz^2)
  re <- ustar*zo/KinV
  zot <- zo*exp(-2.67*(re)^(0.25) + 2.57)
  zot <- Re(zot)
  zoq <- zot
  C_HN <- const_vonKarman*sqrt(C_DN)/(log(hu/zot))
  C_EN <- C_HN

  # calculate neutral transfer coefficients at 10 m
  C_D10N <- (const_vonKarman/log(10/zo))*(const_vonKarman/log(10/zo))
  C_E10N <- (const_vonKarman*const_vonKarman)/(log(10/zo)*log(10/zoq))
  C_H10N <- C_E10N

  # calculate neutral latent and sensible heat fluxes, W/m2
  alhN <- rho_a*xlv*C_EN*Uz*(q_s-q_z)
  ashN <- rho_a*const_SpecificHeatAir*C_HN*Uz*(ts-ta)

  # calculate initial monin obukhov length scale, m
  obu <- (-rho_a*t_virt*(ustar*ustar*ustar))/(const_vonKarman*const_Gravity*
                                                (ashN/const_SpecificHeatAir +
                                                   0.61*(ta + 273.16)*alhN/xlv))

  # iteration to compute corrections for atmospheric stability
  zeta_thres <- 15
  zetam <- -1.574
  zetat <- -0.465

  for( i in 1:20){

    # calulate roughness lengths
    zo <- (0.013*((ustar^2)/const_Gravity)) + (0.11*(KinV/ustar))
    re <- (ustar*zo)/KinV
    xq <- 2.67*re^0.25 - 2.57
    xq[Re(xq) < 0] <- 0
    zoq <- zo/exp(xq)
    zot <- zoq

    # calculate ustar
    zeta <- hu/obu
    zeta[Re(zeta) < -zeta_thres] <- -zeta_thres
    zeta[Re(zeta) > zeta_thres] <- zeta_thres

    zi <- zeta_indices(Re(zeta),zetam,zetat)
    idx_m_vu <- zi$idx_m_vu
    idx_m_u <- zi$idx_m_u
    idx_s <- zi$idx_s
    idx_vs <- zi$idx_vs
    rm(zi)

    ustar[idx_m_vu] <- (Uz[idx_m_vu]*const_vonKarman)/((log((zetam*obu[idx_m_vu])/zo[idx_m_vu]) -
                          psi_zeng(1,zetam)) + 1.14*(((-zeta[idx_m_vu])^0.333) - ((-zetam)^0.333)))

    ustar[idx_m_u] <- (Uz[idx_m_u]*const_vonKarman)/(log(hu/zo[idx_m_u]) -
                                                       psi_zeng(1,zeta[idx_m_u]))

    ustar[idx_s] <- (Uz[idx_s]*const_vonKarman)/(log(hu/zo[idx_s]) + 5*zeta[idx_s])

    ustar[idx_vs] <- (Uz[idx_vs]*const_vonKarman)/((log(obu[idx_vs]/zo[idx_vs]) + 5) +
                                                     (5*log(zeta[idx_vs]) + zeta[idx_vs] - 1))

    # calculate tstar
    zeta <- ht/obu
    zeta[Re(zeta) < -zeta_thres] <- -zeta_thres
    zeta[Re(zeta) > zeta_thres] <- zeta_thres

    zi <- zeta_indices(Re(zeta),zetam,zetat)
    idx_t_vu <- zi$idx_t_vu
    idx_t_u <- zi$idx_t_u
    idx_s <- zi$idx_s
    idx_vs <- zi$idx_vs
    rm(zi)

    tstar[idx_t_vu] <- (const_vonKarman*(ta[idx_t_vu] - ts[idx_t_vu]))/
      ((log((zetat*obu[idx_t_vu])/zot[idx_t_vu]) - psi_zeng(2,zetat)) +
         0.8*((-zetat)^-0.333 - ((-zeta[idx_t_vu]))^-0.333))

    tstar[idx_t_u] <- (const_vonKarman*(ta[idx_t_u] - ts[idx_t_u]))/
      (log(ht/zot[idx_t_u]) - psi_zeng(2,zeta[idx_t_u]))

    tstar[idx_s] <- (const_vonKarman*(ta[idx_s] - ts[idx_s]))/
      (log(ht/zot[idx_s]) + 5*zeta[idx_s])

    tstar[idx_vs] <- (const_vonKarman*(ta[idx_vs] - ts[idx_vs]))/
      ((log(obu[idx_vs]/zot[idx_vs]) + 5) + (5*log(zeta[idx_vs]) + zeta[idx_vs] - 1))

    # calculate qstar
    zeta <- hq/obu
    zeta[Re(zeta) < -zeta_thres] <- -zeta_thres
    zeta[Re(zeta) > zeta_thres] <- zeta_thres

    zi <- zeta_indices(Re(zeta),zetam,zetat)
    idx_t_vu <- zi$idx_t_vu
    idx_t_u <- zi$idx_t_u
    idx_s <- zi$idx_s
    idx_vs <- zi$idx_vs
    rm(zi)

    qstar[idx_t_vu] <- (const_vonKarman*(q_z[idx_t_vu] - q_s[idx_t_vu]))/
      ((log((zetat*obu[idx_t_vu])/zoq[idx_t_vu]) - psi_zeng(2,zetat)) +
         0.8*((-zetat)^-0.333 - ((-zeta[idx_t_vu]))^-0.333))

    qstar[idx_t_u] <- (const_vonKarman*(q_z[idx_t_u] - q_s[idx_t_u]))/
      (log(hq/zoq[idx_t_u]) - psi_zeng(2,zeta[idx_t_u]))

    qstar[idx_s] <- (const_vonKarman*(q_z[idx_s] - q_s[idx_s]))/
      (log(hq/zoq[idx_s]) + 5*zeta[idx_s])

    qstar[idx_vs] <- (const_vonKarman*(q_z[idx_vs] - q_s[idx_vs]))/
      ((log(obu[idx_vs]/zoq[idx_vs]) + 5) + (5*log(zeta[idx_vs]) + zeta[idx_vs] - 1))

    # calculate zeta at 10 m
    zeta <- 10/obu
    zeta[Re(zeta) < -zeta_thres] <- -zeta_thres
    zeta[Re(zeta) > zeta_thres] <- zeta_thres

    zi <- zeta_indices(Re(zeta),zetam,zetat)
    idx_m_vu <- zi$idx_m_vu
    idx_m_u <- zi$idx_m_u
    idx_t_vu <- zi$idx_t_vu
    idx_t_u <- zi$idx_t_u
    idx_s <- zi$idx_s
    idx_vs <- zi$idx_vs
    rm(zi)

    # calculate wind speed at 10 m
    u10[idx_m_vu] <- (ustar[idx_m_vu]/const_vonKarman)*((log((zetam*obu[idx_m_vu])/zo[idx_m_vu]) -
                          psi_zeng(1,zetam)) + 1.14*(((-zeta[idx_m_vu])^0.333) - ((-zetam)^0.333)))

    u10[idx_m_u] <- (ustar[idx_m_u]/const_vonKarman)* (log(10/zo[idx_m_u]) -
                                                         psi_zeng(1,zeta[idx_m_u]))

    u10[idx_s] <- (ustar[idx_s]/const_vonKarman)*(log(10/zo[idx_s]) + 5*zeta[idx_s])

    u10[idx_vs] <- (ustar[idx_vs]/const_vonKarman)*((log(obu[idx_vs]/zo[idx_vs]) + 5) +
                                                      (5*log(zeta[idx_vs]) + zeta[idx_vs] - 1))

    # calcuate air temperature at 10 m
    t10[idx_t_vu] <- ((tstar[idx_t_vu]/const_vonKarman)*
                        ((log((zetat*obu[idx_t_vu])/zot[idx_t_vu]) - psi_zeng(2,zetat)) +
                           0.8*((-zetat)^-0.333 - ((-zeta[idx_t_vu]))^-0.333))) + ts[idx_t_vu]

    t10[idx_t_u] <- ((tstar[idx_t_u]/const_vonKarman)*
                       log(10/zot[idx_t_u]) - psi_zeng(2,zeta[idx_t_u])) + ts[idx_t_u]

    t10[idx_s] <- ((tstar[idx_s]/const_vonKarman)*(log(10/zot[idx_s]) + 5*zeta[idx_s])) + ts[idx_s]

    t10[idx_vs] <- ((tstar[idx_vs]/const_vonKarman)*((log(obu[idx_vs]/zot[idx_vs]) + 5) +
                                       (5*log(zeta[idx_vs]) + zeta[idx_vs] - 1))) + ts[idx_vs]

    # calcuate specific humidity at 10 m
    q10[idx_t_vu] <- ((qstar[idx_t_vu]/const_vonKarman)*
                        ((log((zetat*obu[idx_t_vu])/zoq[idx_t_vu]) - psi_zeng(2,zetat)) +
                           0.8*((-zetat)^-0.333 - ((-zeta[idx_t_vu]))^-0.333))) + q_s[idx_t_vu]

    q10[idx_t_u] <- ((qstar[idx_t_u]/const_vonKarman)*(log(10/zoq[idx_t_u]) -
                                                         psi_zeng(2,zeta[idx_t_u]))) + q_s[idx_t_u]

    q10[idx_s] <- ((qstar[idx_s]/const_vonKarman)*(log(10/zoq[idx_s]) + 5*zeta[idx_s])) + q_s[idx_s]

    q10[idx_vs] <- ((qstar[idx_vs]/const_vonKarman)*((log(obu[idx_vs]/zoq[idx_vs]) + 5) +
                                                       (5*log(zeta[idx_vs]) + zeta[idx_vs] - 1))) + q_s[idx_vs]

    # calculate transfer coefficients corrected for atmospheric stability
    C_H <- (-rho_a*const_SpecificHeatAir*ustar*tstar)/(rho_a*const_SpecificHeatAir*Uz*(ts - ta))
    C_E <- C_H
    C_D <- (ustar*ustar)/(Uz*Uz)

    # calculate tau and sensible and latent heat fluxes
    tau <- rho_a*ustar*ustar
    ash <- -rho_a*const_SpecificHeatAir*ustar*tstar
    alh <- -rho_a*xlv*ustar*qstar

    # calculate transfer coefficients at 10m
    C_H10 <- ash/(rho_a*const_SpecificHeatAir*u10*(ts - t10))
    C_E10 <- alh/(rho_a*xlv*u10*(q_s - q10))
    C_D10 <- (ustar*ustar)/(u10*u10)

    # calculate new monin obukhov length
    obu <- (-rho_a*t_virt*(ustar*ustar*ustar))/
      (const_Gravity*const_vonKarman*((ash/const_SpecificHeatAir) +
                                        (0.61*(ta + 273.16)*alh/xlv)))

    # alter zeta in stable cases ?? is this the VERY stable?
    zeta <- hu/obu
    zi <- zeta_indices(Re(zeta),zetam,zetat)
    idx_m_vu <- zi$idx_m_vu
    idx_m_u <- zi$idx_m_u
    idx_vs <- zi$idx_vs
    rm(zi)

    #Uz[idx_vs] <- max(c(Re(Uz[idx_vs]),0.1),na.rm = TRUE)
    Uz[idx_vs][Re(Uz[idx_vs])<0.1] <- 0.1

    # avoid singularity at um <- 0 for unstable conditions
    idx <- idx_m_vu | idx_m_u # unstable
    th <- (ta + 273.16)*(1000/press)^(const_gas/1004.67) # potential temperature
    thvstar <- tstar*(1 + 0.61*q_z/1000) + 0.61*th*qstar # temperature scaling parameter
    thv <- th*(1 + 0.61*q_z/1000) # virtual potential temperature
    wc[idx] <- 1*(-const_Gravity*ustar[idx]*thvstar[idx]/thv[idx])^0.333
    Uz[idx] <- sqrt(Uz[idx]*Uz[idx] + wc[idx]*wc[idx])

  }

  # if the measurements are taken at a height of 10 m, keep original input data
  if (Re(hu) == 10){
    u10 <- init_wnd
  }
  if (Re(ht) == 10){
    t10 <- ta
  }

  # convert specific humidity to relative humidity
  if (Re(hq) == 10){
    rh10 <- rh
  } else {
    es <- 6.1121*exp(17.502*t10/(t10 + 240.97))*(1.0007+3.46e-6*press)
    em <- q10*press/(0.378*q10 + 0.622)
    rh10 <- 100*em/es
    rh10[Re(rh10) > 100] <- 100
    rh10[Re(rh10) < 0] <- 0
  }

  # calculate evaporation [mm/day]
  rho_w <- 1000*(1-1.9549*0.00001*abs(ts-3.84)^1.68)
  Evap <- 86400*1000*alh/(rho_w*xlv)

  # define output matrix
  mm <- data.frame(tau=tau, alh=alh, ash=ash, ustar=ustar, tstar=tstar, qstar=qstar, u10=u10,
                   t10=t10, q10=q10,rh10=rh10, zo=zo, zot=zot, zoq=zoq, C_D=C_D, C_E=C_E, C_H=C_H,
                   C_D10=C_D10, C_E10=C_E10, C_H10=C_H10, C_D10N=C_D10N, C_E10N=C_E10N,
                   C_H10N=C_H10N, C_DN=C_DN, C_EN=C_EN, C_HN=C_HN, zeta=zeta, Evap=Evap,
                   rho_a=rho_a, q_s=q_s, ts=ts, ta=ta, q_z=q_z, rho_w=rho_w, xlv=xlv, obu=obu)
  mm <- as.data.frame(apply(mm,2,Re))

  return(mm)

}

zeta_indices <- function(zeta,zetam,zetat){

  idx_m_vu <- zeta < zetam            # very unstable conditions zetaM
  idx_m_u <- zeta < 0 & zeta >= zetam # unstable conditions zetaM
  idx_t_vu <- zeta < zetat            # very unstable conditions zetaT
  idx_t_u <- zeta < 0 & zeta >= zetat # unstable conditions zetaT
  idx_s <- zeta >= 0 & zeta <= 1      # stable conditions
  idx_vs <- zeta > 1                  # very stable conditions
  # reps stay undefined
  return(list(idx_m_vu=idx_m_vu,
              idx_m_u=idx_m_u,
              idx_t_vu=idx_t_vu,
              idx_t_u=idx_t_u,
              idx_s=idx_s,
              idx_vs=idx_vs))
}

psi_zeng <- function(k,zeta){
  chik <- (1 - 16*zeta)^0.25
  if( k == 1){
    psi_zeng <- 2*log((1 + chik)*0.5) + log((1+chik*chik)*0.5) -
      2*atan(chik) + (pi/2)
  } else {
    psi_zeng <- 2*log((1 + chik*chik)*0.5)
  }
  return(psi_zeng)
}


