// Estimate MSY and PGY using TMB
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_SCALAR(spawner_init);
  DATA_SCALAR(deviance_init);
  DATA_INTEGER(SR); // 0: HS, 1: BH, 2: RI
  DATA_SCALAR(rec_par_a);
  DATA_SCALAR(rec_par_b);
  // DATA_SCALAR(rec_par_sd);
  DATA_SCALAR(rec_par_rho);
  DATA_SCALAR(bias_corrected_mean);
  DATA_ARRAY(rec_resid_mat);
  DATA_ARRAY(weight_mat);
  DATA_ARRAY(M_mat);
  DATA_INTEGER(Pope);
  DATA_INTEGER(sim_year);
  DATA_INTEGER(nsim);
  DATA_INTEGER(obj_catch); // 0: mean, 1: geomean
  DATA_INTEGER(objective); // 0: MSY, 1: PGY, 2: percentB0 or Bempirical
  DATA_SCALAR(objective_value); // Used for objective 1-2
  DATA_SCALAR(Fcurrent);
  DATA_INTEGER(Fcurrent_year);
  DATA_SCALAR(num_to_mass_scale);

  PARAMETER(x); //x = log(multiplier)
  
  // Type bias_corrected_mean = -Type(0.5)*pow(rec_par_sd,Type(2.0))/(Type(1.0)-pow(rec_par_rho, Type(2.0)));
  array<Type> F_mat(sim_year,nsim);
  array<Type> predN_mat(sim_year,nsim); 
  array<Type> N_mat(sim_year,nsim);
  array<Type> rec_deviance_mat(sim_year,nsim);
  array<Type> spawner_mat(sim_year,nsim);
  array<Type> catch_mat(sim_year,nsim);
  
  // Matrix of Fishing mortality
  for(int i=0; i<nsim; i++) { //replication of simulation 
    for(int t=0; t<sim_year; t++) {
      if(t<Fcurrent_year) {
        F_mat(t,i) = Fcurrent;
      }else{
        F_mat(t,i) = exp(x)*Fcurrent;
      }
    }
  }
  
  // Population dynamics
  for(int i=0; i<nsim; i++) {
    for(int t=0; t<sim_year; t++) {
      if(t==0) {
        if(SR == 0) { //Hockey-stick
          if(spawner_init < rec_par_b) {
            predN_mat(t,i) = spawner_init*rec_par_a;
          }else{
            predN_mat(t,i) = rec_par_b*rec_par_a;
          }
        }
        if(SR == 1) { //Beverton-Holt
          predN_mat(t,i) = rec_par_a*spawner_init/(1+rec_par_b*spawner_init);
        }
        if(SR == 2) { //Ricker
          predN_mat(t,i) = rec_par_a*spawner_init*exp(-rec_par_b*spawner_init);
        }
        N_mat(t,i) = predN_mat(t,i)*exp(rec_par_rho*deviance_init)*exp(rec_resid_mat(t,i));
        rec_deviance_mat(t,i) = log(N_mat(t,i)/predN_mat(t,i))-bias_corrected_mean;
        spawner_mat(t,i) = N_mat(t,i)*exp(-M_mat(t,i)-F_mat(t,i));
      }else{
        if(SR == 0) { //Hockey-stick
          vector<Type> rec_pred(2);
          rec_pred(0) = spawner_mat(t-1,i)*rec_par_a;
          rec_pred(1) = rec_par_b*rec_par_a;
          predN_mat(t,i) = min(rec_pred);
          // if(spawner_mat(t-1,i) < rec_par_b) {
          //   predN_mat(t,i) = spawner_mat(t-1,i)*rec_par_a;
          // }else{
          //   predN_mat(t,i) = rec_par_b*rec_par_a;
          // }
        }
        if(SR == 1) { //Beverton-Holt
          predN_mat(t,i) = rec_par_a*spawner_mat(t-1,i)/(1+rec_par_b*spawner_mat(t-1,i));
        }
        if(SR == 2) { //Ricker
          predN_mat(t,i) = rec_par_a*spawner_mat(t-1,i)*exp(-rec_par_b*spawner_mat(t-1,i));
        }
        N_mat(t,i) = predN_mat(t,i)*exp(rec_par_rho*rec_deviance_mat(t-1,i))*exp(rec_resid_mat(t,i));
        rec_deviance_mat(t,i) = log(N_mat(t,i)/predN_mat(t,i))-bias_corrected_mean;
        spawner_mat(t,i) = N_mat(t,i)*exp(-M_mat(t,i)-F_mat(t,i)); 
      }
    }
  }
  
  // Catch equation
  for(int i=0; i<nsim; i++) {
    for(int t=0; t<sim_year; t++) {
      if(Pope) {
        catch_mat(t,i) = num_to_mass_scale*weight_mat(t,i)*N_mat(t,i)*exp(-Type(0.5)*M_mat(t,i))*(1-exp(-F_mat(t,i)));
      }else{
        catch_mat(t,i) = num_to_mass_scale*weight_mat(t,i)*N_mat(t,i)*(1-exp(-M_mat(t,i)-F_mat(t,i)))*F_mat(t,i)/(M_mat(t,i)+F_mat(t,i));
      }
    }
  }
  
  // Get Catch or SSB in the final year
  Type obj = 0;
  for(int i=0; i<nsim; i++) {
    if(objective < 2) {
      if(obj_catch == 0) {
        obj += catch_mat(sim_year-1,i);
      }else{
        obj += log(catch_mat(sim_year-1,i));
      }
    }else{
      if(obj_catch == 0) {
        obj += num_to_mass_scale*spawner_mat(sim_year-1,i)*weight_mat(sim_year-1,i);
      }else{
        obj += num_to_mass_scale*log(spawner_mat(sim_year-1,i)*weight_mat(sim_year-1,i));
      }
    }
  }
  obj /= nsim;
  if(obj_catch == 1) { 
    obj = exp(obj); // geomean
  }
  
  if(objective == 0) { //  MSY
    obj = log(obj);
    obj *= -1; // negative log catch
  }else{ // PGY, percentB0, and Bempirical
    obj = pow(log(obj/objective_value), Type(2.0));
  }

  return obj;
  // 
  // ADREPORT(F_mat);
  // ADREPORT(N_mat);
  // ADREPORT(spawner_mat);
  // ADREPORT(catch_mat);
  // 
}
