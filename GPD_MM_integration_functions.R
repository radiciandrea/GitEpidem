# integration function 1

df <- function(t, x, parms) {
  
  # initial conditions and parameters
  with(parms, { 
    
    E = x[(1+n_r*0):(1*n_r)]
    J = x[(1+n_r*1):(2*n_r)]
    I = x[(1+n_r*2):(3*n_r)]
    A = x[(1+n_r*3):(4*n_r)]
    E_d = x[(1+n_r*4):(5*n_r)]
    
    #t_n = t[1]-t_s+1 # time of numerical integration to index matrix
    t_n = t[1]
    t_h = 24*(t - t_n) #should put t and not t[1]
    TM = temp_M[max(1,t_n-1),]*(t_h<t_sr[t_n, ]) + temp_M[t_n,]*(t_h>t_sr[t_n, ])
    Tm = temp_m[t_n, ]*(t_h<14) + temp_m[min(t_n+1, length(temp_m))]*(t_h>14)
    
    temp_h = ((TM+Tm)/2 + (TM-Tm)/2*cos(pi*(t_h+10)/(10+t_sr[t_n, ])))*(t_h<t_sr[t_n, ])+
      ((TM+Tm)/2 - (TM-Tm)/2*cos(pi*(t_h-t_sr[t_n, ])/(14-t_sr[t_n, ])))*(t_h>t_sr[t_n, ])*(t_h<14)+
      ((TM+Tm)/2 + (TM-Tm)/2*cos(pi*(t_h-14)/(10+t_sr[t_n, ])))*(t_h>14)
    
    delta_J = 1/(83.85 - 4.89*temp_h + 0.08*temp_h^2) #juvenile development rate (in SI: 82.42 - 4.87*temp_h + 0.08*temp_h^ 2)
    delta_I = 1/(50.1 - 3.574*temp_h + 0.069*temp_h^2) #first pre blood mean rate
    mu_E = -log(0.955 * exp(-0.5*((temp_h-18.8)/21.53)^6)) # egg mortality rate
    mu_J = -log(0.977 * exp(-0.5*((temp_h-21.8)/16.6)^6)) # juvenile mortality rate
    beta = (33.2*exp(-0.5*((temp_h-70.3)/14.1)^2)*(38.8 - temp_h)^1.5)*(temp_h<= 38.8) #fertility rate
    
    # ODE definition 
    dE = beta*(1-omega[t_n, ])*A - (h[t_n, ]*delta_E + mu_E)*E
    dJ = h[t_n, ]*(delta_E*E + sigma[t_n, ]*gamma*E_d) - (delta_J + mu_J + J/K[t_n, ])*J  
    dI = 0.5*delta_J*J - (delta_I + mu_A[t_n, ])*I
    dA = delta_I*I - mu_A[t_n, ]*A
    dE_d = beta*omega[t_n, ]*A -  h[t_n, ]*sigma[t_n, ]*E_d #I believe there should be an additional mortality due to winter
    
    dx <- c(dE, dJ, dI, dA, dE_d)
    
    # #ODE positive verification
    # dx <- ifelse(x + dx * dt < 0, 0, dx*dt)
    
    return(list(dx))})
}

# this has the log of the adults
# https://stackoverflow.com/questions/47401678/solving-odes-only-positive-solutions
# https://stackoverflow.com/questions/41648878/replacing-negative-values-in-a-model-system-of-odes-with-zero

df_log <- function(t, x, parms) {
  
  # initial conditions and parameters
  with(parms, { 
    
    logE = x[(1+n_r*0):(1*n_r)]
    logJ = x[(1+n_r*1):(2*n_r)]
    logI = x[(1+n_r*2):(3*n_r)]
    logA = x[(1+n_r*3):(4*n_r)]
    E_d = x[(1+n_r*4):(5*n_r)]
    
    #t_n = t[1]-t_s+1 # time of numerical integration to index matrix
    t_n = t[1]
    t_h = 24*(t - t_n) #should put t and not t[1]
    TM = temp_M[max(1,t_n-1),]*(t_h<t_sr[t_n, ]) + temp_M[t_n,]*(t_h>t_sr[t_n, ])
    Tm = temp_m[t_n, ]*(t_h<14) + temp_m[min(t_n+1, length(temp_m))]*(t_h>14)
    
    temp_h = ((TM+Tm)/2 + (TM-Tm)/2*cos(pi*(t_h+10)/(10+t_sr[t_n, ])))*(t_h<t_sr[t_n, ])+
      ((TM+Tm)/2 - (TM-Tm)/2*cos(pi*(t_h-t_sr[t_n, ])/(14-t_sr[t_n, ])))*(t_h>t_sr[t_n, ])*(t_h<14)+
      ((TM+Tm)/2 + (TM-Tm)/2*cos(pi*(t_h-14)/(10+t_sr[t_n, ])))*(t_h>14)
    
    delta_J = 1/(83.85 - 4.89*temp_h + 0.08*temp_h^2) #juvenile development rate (in SI: 82.42 - 4.87*temp_h + 0.08*temp_h^ 2)
    delta_I = 1/(50.1 - 3.574*temp_h + 0.069*temp_h^2) #first pre blood mean rate
    mu_E = -log(0.955 * exp(-0.5*((temp_h-18.8)/21.53)^6)) # egg mortality rate
    mu_J = -log(0.977 * exp(-0.5*((temp_h-21.8)/16.6)^6)) # juvenile mortality rate
    beta = (33.2*exp(-0.5*((temp_h-70.3)/14.1)^2)*(38.8 - temp_h)^1.5)*(temp_h<= 38.8) #fertility rate
    
    # ODE definition 
    dlogE = beta*(1-omega[t_n, ])*exp(logA)/exp(logE) - (h[t_n, ]*delta_E + mu_E)
    dlogJ = h[t_n, ]*(delta_E*exp(logE) + sigma[t_n, ]*gamma*E_d)/exp(logJ) - (delta_J + mu_J + exp(logJ)/K[t_n, ])  
    dlogI = 0.5*delta_J*exp(logJ)/exp(logI) - (delta_I + mu_A[t_n, ])
    dlogA = delta_I*exp(logI)/exp(logA) - mu_A[t_n, ]
    dE_d = beta*omega[t_n, ]*exp(logA) -  h[t_n, ]*sigma[t_n, ]*E_d #I believe there should be an additional mortality due to winter
    
    dx <- c(dlogE, dlogJ, dlogI, dlogA, dE_d)
    
    # #ODE positive verification
    # dx <- ifelse(x + dx * dt < 0, 0, dx*dt)
    
    return(list(dx))})
}