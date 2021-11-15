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
  array[N_series] int<lower=0> rep_temp;

  array[N_series, 1] real time_end;
  array[N_series, 2] int x_start;
  array[N_series] vector[2] x_start_vector;
  array[N_series, 2] int x_end;
  
}

parameters {
  
  //// ODE parameters
  real b_log; // grand mean parameter
  vector[3] b_log_temp_raw; // contrasts for the temperature levels
  real<lower=0> sigma_b; // standard deviation of the parameters
  
  real c_log;
  vector[3] c_log_temp_raw;
  real<lower=0> sigma_c;
  
  real h_log;
  vector[3] h_log_temp_raw;
  real<lower=0> sigma_h;

  real K_log; // grand mean parameter
  vector[3] K_log_temp_raw; // contrasts for the temperature levels
  real<lower=0> sigma_K; // standard deviation of the parameters
  
  real q;
  vector[3] q_temp_raw;
  real<lower=0> sigma_q;

  real r_log;
  vector[3] r_log_temp_raw;
  real<lower=0> sigma_r;

  // array[N_series] vector[2] x_hat_start_log;
  
}

transformed parameters {
  
  array[N_series, 1] vector<lower=0>[2] x_hat_end;
    // ode integrator expects an array[N_data,N_times] vector[N_states]. This is why some arrays have the superficial indices here: array[N_series, 1]
  
  vector[3] b_log_temp = b_log + sigma_b * b_log_temp_raw; // implies: b_log_temp ~ normal(b_log, sigma_v)
  vector[3] c_log_temp = c_log + sigma_c * c_log_temp_raw;
  vector[3] h_log_temp = h_log + sigma_h * h_log_temp_raw;
  vector[3] K_log_temp = K_log + sigma_K * K_log_temp_raw;
  vector[3] q_temp = q + sigma_q * q_temp_raw;
  vector[3] r_log_temp = r_log + sigma_r * r_log_temp_raw; // implies: b_log_temp ~ normal(b_log, sigma_v)
  
  for(s in 1:N_series) {
      
      //// Version with latent initial state
      // x_hat_end[s,] = ode_rk45(x_t, exp(x_hat_start_log[s]), 0.0, time_end[s,],
      //                 exp(b_log), exp(c_log), exp(h_log), exp(K_log), q, exp(r_log));
      
      int t = rep_temp[s];
      x_hat_end[s,] = ode_rk45(x_t, x_start_vector[s,], 0.0, time_end[s,],
                      exp(b_log_temp[t]), exp(c_log_temp[t]), exp(h_log_temp[t]), exp(K_log_temp[t]), q_temp[t], exp(r_log_temp[t]));
  }
}


model {
  // Priors
  b_log ~ normal(1.117, 1); // b_log ~ normal(-4.106669, 1);
  c_log ~ normal(-4.60517, 1);
  h_log ~ normal(-6.907755, 1);
  K_log ~ normal(8.477, 1);
  q ~ normal(-.56, 1);
  r_log ~ normal(0.1897936, 1);
  
  sigma_b ~ std_normal();
  sigma_c ~ std_normal();
  sigma_h ~ std_normal();
  sigma_K ~ std_normal();
  sigma_q ~ std_normal();
  sigma_r ~ std_normal();
  
  b_log_temp_raw ~ std_normal();
  c_log_temp_raw ~ std_normal();
  h_log_temp_raw ~ std_normal();
  K_log_temp_raw ~ std_normal();
  q_temp_raw ~ std_normal();
  r_log_temp_raw ~ std_normal();

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
