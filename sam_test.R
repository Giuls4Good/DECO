#' Implements the PaRIS algorithm
#'
#' @param N is the number of particles
#' @return Gives vectors X and Y representing the state space model
#' @export PARIS
PARIS = function(N = 10, N_tilde = 5, Y,  max_time, a = 0.7, b = 1, sigma_X = 0.2, sigma_Y = 1, h_tilde = function(x,y){x*y}, type = 'AR'){
  BPF_out = BPF(N = N, Y = Y, max_time = max_time, a = a, b = b, sigma_X = sigma_X, sigma_Y = sigma_Y)
  particles = BPF_out$particles
  weights   = BPF_out$weights
  
  tau = matrix(numeric(max_time*N), max_time, N)
  for (t in seq_len((max_time-1))){
    if (type != 'AR'){
      norm_particle_input = matrix(particles[(t+1),], N, N, byrow = T)  #Make the new particles like this so they can be inputted into the dnorm function
      q_weights = t(dnorm(norm_particle_input, mean = a*particles[t,], sd = sigma_X))  #test that this is right!!!
      
      J = matrix(numeric(N*N_tilde), N, N_tilde)
      for (i in 1:N){
        J[i,]  =  (1:N) %*% rmultinom(N_tilde, size = 1, prob = weights[t,]*q_weights[,i])
      }
    }else{
      J = accept_reject(N = N, N_tilde = N_tilde, old_weights = weights[t,], new_particles = particles[(t+1),], old_particles = particles[t,], sigma_X = sigma_X)
    }
    for (i in 1:N){
      inner = tau[t,][J[i,]] + h_tilde(particles[t,][J[i,]], rep(particles[(t+1),i], N_tilde) )
      tau[(t+1), i] = (1/N_tilde)*sum(inner)
    }
  }

  out = list(particles = particles, weights = weights, tau = tau)
  return(out)
}

