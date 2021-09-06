functions {
  
  //// The ODE system. See https://mc-stan.org/docs/2_27/functions-reference/functions-ode-solver.html
  vector x_t(real time, vector x,
             real b, real c, real h, real K, real q, real r) {

    vector[2] dx_dt;
    
    // Prey:
    dx_dt[1] = -b * pow(x[1], 1+q) / (1 + b * h * pow(x[1], 1+q)) * x[2] + r * x[1] * (1.0 - x[1]/K);
    // Predators:
    dx_dt[2] = x[2] * c * b * pow(x[1], 1+q) / (1 + b * h * pow(x[1], 1+q));

    return dx_dt;
  }

}

data {
  
  int<lower=0> N_series;
  array[N_series, 1] real time_end;
  array[N_series, 2] int x_start;
  array[N_series] vector<lower=0>[2] x_start_vector;
  array[N_series, 2] int x_end;
  
}

parameters {
  
  //// ODE parameters
  real b_log;
  real c_log;
  real h_log;
  real K_log;
  real q;
  real r_log;
  
  // array[N_series] vector[2] x_hat_start_log;
  
}

transformed parameters {
  
  array[N_series, 1] vector<lower=0>[2] x_hat_end;
    // ode integrator expects an array[N_data,N_times] vector[N_states]. This is why some arrays have the superficial indices here: array[N_series, 1]
  
  for(s in 1:N_series) {
      
      //// Version with latent initinal state
      // x_hat_end[s,] = ode_rk45(x_t, exp(x_hat_start_log[s]), 0.0, time_end[s,],
      //                 exp(b_log), exp(c_log), exp(h_log), exp(K_log), q, exp(r_log));
      
      x_hat_end[s,] = ode_rk45(x_t, x_start_vector[s,], 0.0, time_end[s,],
                      exp(b_log), exp(c_log), exp(h_log), exp(K_log), q, exp(r_log));
                      
  }
}


model {
  // Priors
  b_log ~ normal(-4.106669, 1);
  c_log ~ normal(-5.600194, 2);
  h_log ~ normal(-3.531813, 2);
  K_log ~ normal(7.901859, 2);
  q ~ normal(0.867401, 2);
  r_log ~ normal(-0.994875, 2);
  
  // sigma ~ normal(0, 1);
  

  //// Statistical model
  for(s in 1:N_series) {
    
    //// Version with latent initinal state: Fitiing the unknown true initial state to the data,
    // x_start[s,] ~ poisson_log(x_hat_start_log[s]);
    
    //// Fitiing the final state from the ODE to the data.
    x_end[s,] ~ poisson(x_hat_end[s, 1]); // index 1 is there to accomodate for data structure array[N_data,N_times] vector[N_states]
  }
}
