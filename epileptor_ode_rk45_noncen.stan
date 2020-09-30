functions {
        
   real[] epileptor_ode(real t, real[] y, real[] params, real[] x_r, int[] x_i) {
 
   real eta;  
   real dydt[size(y)];
   eta = params[1];        
      
   dydt[1] = 1.0 - y[1]*y[1]*y[1] - 2.0*y[1]*y[1]- y[2] + x_r[1];
   dydt[2] = (1.0/x_r[2])*(4*(y[1] - eta) - y[2] );
                
   return dydt;
    }

}
    
    
data {
    int<lower=1> nt;  
    real t0;
    real Ts[nt];  
    real dt;
    real eta_true;    
    real x_init;
    real z_init;
    real xlim[2];
    real zlim[2];
    real I1;
    real tau0;
    vector[nt] xs;
}

transformed data {
    real x_r[2];
    int x_i[0];
    real obs[nt,1];

    real std=1.;  
    
    x_r[1]=I1 ;        
    x_r[2]=tau0; 
    
    for (t in 1:nt) {
        obs[t,1] =xs[t];
        }
}


parameters {
    real x_init_star;
    real z_init_star;
    real eta_star;  
    real amplitude;
    real offset; 
    real<lower=0.0> eps;   
}


transformed parameters {
    real eta;
    real params[1];
    real y0[2];
    
    eta = eta_true + eta_star;

    params[1] = eta;
    y0[1] = x_init + std*x_init_star;
    y0[2] = z_init + std*z_init_star;

}


model {
    real yhat[nt,2];
    
    real xhat[nt,1];
    real zhat[nt,1];
    

    x_init_star ~ normal(0., 1.); 
    z_init_star ~ normal(0., 1.); 
    

    eta_star ~ normal(0., 1.0);  
    amplitude ~ normal(1., 1.0);
    offset ~ normal(0., 1.0); 
    eps ~ normal(0., 1.); 
    
    yhat = integrate_ode_rk45(epileptor_ode, y0, t0, Ts, params, x_r, x_i);
    


    for (t in 1:nt) {
        //obs[t,1] ~ normal(yhat[t,1], eps);                   
        xhat[t,1] = amplitude*yhat[t,1]+offset;
        zhat[t,1] = amplitude*yhat[t,2]+offset;
        target+=normal_lpdf(obs[t,1]| yhat[t,1], eps);   
  }
  
  
}


generated quantities {
    real yhat[nt,2];

    real x[nt,1];
    real z[nt,1];
    
    real xhat_q [nt,1];
    real zhat_q [nt,1];
    real x_ppc [nt,1];
    real z_ppc [nt,1];

    real log_lik [nt,1];

    
    yhat = integrate_ode_rk45(epileptor_ode, y0, t0, Ts, params, x_r, x_i);

    x[1,1]=y0[1];
    z[1,1]=y0[2];

    
    for (t in 1:(nt-1)) {
            x[t+1,1] =yhat[t+1,1];
            z[t+1,1] =yhat[t+1,2]; 
    } 
    
    
    for (t in 1:(nt)) {            
            //xhat_q[t,1] = amplitude*x[t,1]+offset;
            //zhat_q[t,1] = amplitude*z[t,1]+offset;
            xhat_q[t,1] = x[t,1];
            zhat_q[t,1] = z[t,1];
    } 
    
    
    for (i in 1:nt) {
         x_ppc[i,1] = normal_rng(xhat_q[i,1], eps);
         z_ppc[i,1] = normal_rng(zhat_q[i,1], eps);

  }
  
 
    for (i in 1:nt){
        log_lik[i,1]= normal_lpdf(obs[i,1]| xhat_q[i,1], eps);
      }

}