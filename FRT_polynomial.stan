functions {
  
  //// The ODE system. See https://mc-stan.org/docs/2_27/functions-reference/functions-ode-solver.html
  vector x_t(real time, vector x,
             real b2, real c2, real h2, real K2, real q2, real r2, real b1, real c1, real h1, real K1, real q1, real r1, real b0, real c0, real h0, real K0, real q0, real r0) {

    vector[2] dx_dt;
    
    // Prey:
    // dx_dt[1] = -b * pow(x[1], 1+q) / (1 + b * h * pow(x[1], 1+q)) * x[2] + r * x[1] * (1.0 - x[1]/K);
    // Predators:
    // dx_dt[2] = x[2] * c * b * pow(x[1], 1+q) / (1 + b * h * pow(x[1], 1+q));

    // prey:
    dx_dt[1] = - exp(b2 * pow(-0.5, 2) + b1 * -0.5 + b0) * pow(x[0],1+(q2 * pow(-0.5, 2) + q1 * -0.5 + q0)) / (1 + exp(b2 * pow(-0.5, 2) + b1 * -0.5 + b0) * exp(h2 * pow(-0.5, 2) + h1 * -0.5 + h0) * pow(x[0],1+(q2 * pow(-0.5, 2) + q1 * -0.5 + q0))) * x[1] + exp(r2 * pow(-0.5, 2) + r1 * -0.5 + r0) * x[0] * (1.0 - x[0]/exp(K2 * pow(-0.5, 2) + K1 * -0.5 + K0));
    // predator:
    dx_dt[2] = x[1] * exp(c2 * pow(-0.5, 2) + c1 * -0.5 + c0) * exp(b2 * pow(-0.5, 2) + b1 * -0.5 + b0) * pow(x[0],1+(q2 * pow(-0.5, 2) + q1 * -0.5 + q0)) / (1 + exp(b2 * pow(-0.5, 2) + b1 * -0.5 + b0) * exp(h2 * pow(-0.5, 2) + h1 * -0.5 + h0) * pow(x[0],1+(q2 * pow(-0.5, 2) + q1 * -0.5 + q0)));

    return dx_dt;
  }

}

data {
  
  int<lower=0> N_series;
  array[N_series] int<lower=0> rep_temp;

  array[N_series, 1] real time_end;
  array[N_series, 2] int x_start;
  array[N_series] vector[2] x_start_vector;
  array[N_series, 2] int x_end;
  
}

parameters {
  
  //// ODE parameters
  real b2_log; 
  real b1_log;
  real b0_log;
  
  real c2_log; 
  real c1_log;
  real c0_log;
  
  real h2_log; 
  real h1_log;
  real h0_log;
  
  real q2_log; 
  real q1_log;
  real q0_log;
  
  real K2_log; 
  real K1_log;
  real K0_log;
  
  real r2_log; 
  real r1_log;
  real r0_log;

  // array[N_series] vector[2] x_hat_start_log;
  
}

transformed parameters {
  
  array[N_series, 1] vector<lower=0>[2] x_hat_end;
    // ode integrator expects an array[N_data,N_times] vector[N_states]. This is why some arrays have the superficial indices here: array[N_series, 1]
  
  for(s in 1:N_series) {
      
      //// Version with latent initial state
      // x_hat_end[s,] = ode_rk45(x_t, exp(x_hat_start_log[s]), 0.0, time_end[s,],
      //                 exp(b_log), exp(c_log), exp(h_log), exp(K_log), q, exp(r_log));
      
      x_hat_end[s,] = ode_rk45(x_t, x_start_vector[s,], 0.0, time_end[s,],
                      b2_log,
                      b1_log, 
                      b0_log, 
                      c2_log, 
                      c1_log, 
                      c0_log, 
                      h2_log,
                      h1_log,
                      h0_log,
                      K2_log,
                      K1_log,
                      K0_log,
                      q2_log,
                      q1_log,
                      q0_log,
                      r2_log,
                      r1_log,
                      r0_log);
  }
}


model {
  // Priors
  b2_log ~ normal(0, 3); // b_log ~ normal(-4.106669, 1);
  b1_log ~ normal(0, 3); // b_log ~ normal(-4.106669, 1);
  b0_log ~ normal(0, 3); // b_log ~ normal(-4.106669, 1);
  c2_log ~ normal(0, 1);
  c1_log ~ normal(0, 1);
  c0_log ~ normal(0, 1);
  h2_log ~ normal(0, 3);
  h1_log ~ normal(0, 3);
  h0_log ~ normal(0, 3);
  K2_log ~ normal(5, 5);
  K1_log ~ normal(5, 5);
  K0_log ~ normal(5, 5);
  q2_log ~ normal(0, 5);
  q1_log ~ normal(0, 5);
  q0_log ~ normal(0, 5);
  r2_log ~ normal(0, 1);
  r1_log ~ normal(0, 1);
  r0_log ~ normal(0, 1);


  //// Statistical model
  for(s in 1:N_series) {

    //// Version with latent initinal state: Fitiing the unknown true initial state to the data,
    // x_start[s,] ~ poisson_log(x_hat_start_log[s]);

    //// Fitiing the final state from the ODE to the data.
    //// Fitiing the final state from the ODE to the data.
    if(x_start[s, 2] == 0) {
        x_end[s,1] ~ poisson(x_hat_end[s, 1, 1]); // index 1 is there to accomodate for data structure array[N_data,N_times] vector[N_states]
    } else {
       x_end[s,] ~ poisson(x_hat_end[s, 1]); // index 1 is there to accomodate for data structure array[N_data,N_times] vector[N_states]
    }
  }
}
