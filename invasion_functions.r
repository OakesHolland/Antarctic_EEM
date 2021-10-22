invasion_functions <- function(dist_mean,rep){
  
  intdata <- read.csv("interactions.csv") #Predator/prey interactions
  intdata2 <- read.csv("interactions2.csv") #Interactions where things eat each other 
  rdata <- read.csv("growth_rate.csv", header = FALSE) #Species growth rates, where known
  
  intdata = data.matrix(intdata)
  intdata2 = data.matrix(intdata2)
  intdata = intdata[1:nrow(intdata),2:ncol(intdata)]
  intdata2 = intdata2[1:nrow(intdata2),2:ncol(intdata2)]
  colnames(intdata) <- 1:ncol(intdata)
  colnames(intdata2) <- 1:ncol(intdata2)
  
  r_guess = rdata[1:nrow(rdata),2]
  r_data_pos = which(r_guess != "NA")
  r_guess = matrix(c(r_data_pos,r_guess[r_data_pos]),ncol=2)
  
  n_s = nrow(intdata)
  A_array = array(0,c(n_s,n_s,rep))
  r_array = matrix(0,n_s,rep)
  
  t_vec = seq(from=0,to=100,by=0.1) #Time vector, to = number of years)
  y_array = array(0,c(length(t_vec),n_s,rep))
  sol_array = c(list(),rep)
  
  outcomes = matrix(0,n_s-1,rep) #Vector to store pos/neg interactions
  outcomes_abund = outcomes #Vector to store relative change in abundance
  outcomes_dir = outcomes
  
  pred_prey_pair = generate_pred_prey_pairs(intdata)
  A_sign_original = intdata
  A_sign2 = intdata2
  ctr = 0 #Initialize a counter for number of simulations, final count of simulations is (ctr - reps)
  
  for (i in 1:rep){
    print(i/rep*100)
    flag2 = 0
    while (flag2 == 0){ # Only keep systems with reasonable relative change
      flag = 0
      while (flag == 0){ #Check that invasive-free systems are feasible
        ctr = ctr + 1
        A_sign = A_sign_original
        A_sign_2 = A_sign2
        
        
        A = runif(prod(dim(A_sign)),0,dist_mean)*A_sign #Generate random strengths between 0 and dist_mean with sign that corresponds to sign on matrix 
        #Second matrix for species that eat each other: generate a positive random strength (eats) & negative random strength (eaten) and 
        #and add together to get overall interaction strength
        A2 = (runif(prod(dim(A_sign_2)),0,dist_mean)* A_sign_2)+(runif(prod(dim(A_sign_2)),0,dist_mean)*(-A_sign_2)) 
        
        A = A + A2 #Add the matrices together to form an overall matrix of interaction strengths
        
        r = matrix(0,n_s,1)
        
        growth_mean = mean(r_guess[,2])
        growth_sd = sd(r_guess[,2])
        
        rlmean = log(growth_mean^2/sqrt(growth_sd^2 + growth_mean^2))
        rlsd = sqrt(log(1 + (growth_sd^2/growth_mean^2)))
        
        r[r_guess[,1]] = rlnorm(nrow(r_guess),meanlog = rlmean,sdlog = rlsd)
        
        r[r==0] = runif(sum(r==0))*max(r)
        
        
        
        for (j in 1:ncol(pred_prey_pair)){
          if (A[pred_prey_pair[j,1]] >= 0.2*abs(A[pred_prey_pair[j,2]])){
            A[pred_prey_pair[j,1]] = runif(1)*0.2*abs(A[pred_prey_pair[j,2]])
          }
        }
        
        eq_pred = solve(A,-r) #Solve for equation with invasive species
        eq_nopred = solve(A[2:nrow(A),2:ncol(A)],-r[2:nrow(r)]) #Solve for equations without invasive species
        
        jacobian_pred = matrix(0,nrow(A),ncol(A))
        for (k in 1:n_s){
          jacobian_pred[k,] = A[k,] * eq_pred[k]
          jacobian_pred[k,k] = jacobian_pred[k,k] + sum(A[k,]*eq_pred)
        }
        jacobian_pred[diag(1,nrow=n_s) == 1] = jacobian_pred[diag(1,nrow=n_s) == 1] + r
        
        As = A[2:ncol(A),2:nrow(A)]
        rs = r[2:nrow(r),]
        jacobian = matrix(0,nrow(As),ncol(As))
        for (k in 1:n_s-1){
          jacobian[k,] = As[k,]*eq_nopred[k]
          jacobian[k,k] = jacobian[k,k] + sum(As[k,]*eq_nopred)
        }
        jacobian[diag(1,nrow = n_s-1) == 1] = jacobian[diag(1,nrow = n_s-1) == 1] + rs
        
        
        eigs = eigen(jacobian)$values
        eigs_pred = eigen(jacobian_pred)$values
        
        
        # If non-invaded matrix is plausible, continue, else return to start of while loop
        if (all(eq_nopred>0) && max(Re(eigs))<0){
          flag = 1}
        
      }
      
      outcomes[,i] = eq_pred[-1]/eq_nopred > 1 #Does a species increase with an invasive species
      outcomes_abund[,i] = eq_pred[-1]/eq_nopred #Ratio of abundance with and without invasive species
      
      n_intro = matrix(0,prod(dim(eq_pred)),1)
      n_intro[1] = eq_pred[1]/10
      n_intro[-1] = eq_nopred
      
      init_dirs = n_intro*r + (A%*%n_intro)*n_intro
      outcomes_dir[,i] = init_dirs[-1] > 0
      
      if (max(outcomes_abund[,i])<100){ #Only allows systems with reasonable increases
        sol = ode_solve(function(t,y) return(derivs(t,y,A,r)),matrix(c(0,100)),n_intro)
        t = sol[1]
        y = sol[2]
        flag2 = 1
      }
    }
    t = matrix(unlist(t))
    y = matrix(unlist(y),ncol=n_s)
    y_int = interp1(t,y,t_vec)
    
    A_array[,,i] = A
    r_array[,i] = r
    y_array[,,i] = y_int
    
    
  }
  
  fname = paste('results_dmean_',toString(dist_mean),'_rep_',toString(rep),'.rda',sep='')
   save(file=fname,list = c("A_array","r_array","y_array","t_vec","outcomes","outcomes_abund","outcomes_dir", "ctr"))
  
  
}

#Code to extract pairs of interacting species

generate_pred_prey_pairs <- function(A){
  pred_prey_pair = matrix(0,nrow(A)*ncol(A),2)
  counter = 0
  n_s = nrow(A)
  
  for (j in 2:n_s){
    for (i in 1:(j-1)){
      if ((A[i,j] != 0) && (A[i,j] == - A[j,i])){
        counter = counter + 1
        if (A[i,j] == 1){
          pred_prey_pair[counter,1:2] = c(i,j)
        } else {
          pred_prey_pair[counter,1:2] = c(j,i)
        }
      }
    }
  }
  
  if (counter == 0){
    pred_prey_pair = matrix(nrow = 0,ncol = 0)
  } else {
    pred_prey_pair = pred_prey_pair[1:counter,1:2]
    for (i in 1:counter){
      c_pair = pred_prey_pair[i,1:2]
      pred_prey_pair[i,1] = sub2ind(dim(A),c_pair[1],c_pair[2])
      pred_prey_pair[i,2] = sub2ind(dim(A),c_pair[2],c_pair[1])
    }
  }
  return(pred_prey_pair)
  
}

sub2ind <- function(Adim,r,c){
  rows = Adim[1]
  return((c-1)*rows + r)
}

all <- function(r){
  out = 1
  for (i in 1:prod(dim(r))){
    if (r[i]){
      out = out*1
    } else {
      out = out*0
    }
  }
  if (out == 0){
    out = FALSE
  } else {
    out = TRUE
  }
  return(out)
}


derivs <- function(t,n,A,r){
  return(n*r + (A%*%n)*n)
  
}


ode_solve <- function(fun,t_window,IC){
  t = seq(from = t_window[1],to = t_window[2],by = 0.01)
  y = matrix(0,length(t),length(IC))
  y[1,] = IC
  for (i in 2:length(t)){
    y[i,] = pmax(y[i-1,] + (t[2]-t[1])*fun(t[i-1],y[i-1,]), 0)
  }
  return(list(t,y))
}


interp1 <- function(t,y,tv){
  out = matrix(0,length(tv),ncol(y))
  for (i in 1:ncol(y)){
    out[,i] = approx(t,y[,i],tv)$y
  }
  return(out)
}

