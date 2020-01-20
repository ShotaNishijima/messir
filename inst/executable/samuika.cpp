// SAMUIKA: State-space Assessment Model Used for IKA

#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA //
  DATA_INTEGER(NYear);
  DATA_INTEGER(NStock);

  // Fishing process data
  DATA_ARRAY(M); //Natural morality coefficient (2 x NYear matrix)
  DATA_ARRAY(Weight); //2 x NYear matrix
  DATA_IVECTOR(SDlogF_key); //c(0,1) or c(0,0)
  DATA_INTEGER(logF_diff); // 0: 0 difference, 1:difference

  // Recruitment process data
  DATA_INTEGER(SR); //0:HS, 1:BH, 2:RI
  DATA_IARRAY(reca_key); //2 x NYear matrix
  DATA_IARRAY(recb_key); //2 x NYear matrix
  DATA_IARRAY(recSD_key); //2 x NYear matrix

  // Catch observation data
  DATA_INTEGER(NCatch); //# of catch data
  DATA_VECTOR(Catch); //vector of catch weights
  DATA_INTEGER(Pope); // Pope (1) or Baranov (0)?
  DATA_IARRAY(Catch_key); //# of rows = # of Catch (NCatch)
  // Column 0: Stock_ID, 1: iy (year from the first year), 2: SDlogC_key
  DATA_SCALAR(scale_num_to_mass);

  // Index observation data
  DATA_INTEGER(NIndex); //# of index values
  DATA_VECTOR(Index); //vector of index values
  DATA_IARRAY(Index_key); // # of rows = # of index (NIndex)
  // Column 0: Stock_ID, 1: iy, 2: q_key, 3:SDlogCPUE_key, 4:beta_key

  // Prior of log(Z)
  DATA_ARRAY(logZ_w); //2(Nstock) x NYear matrix (0 or 1)
  DATA_VECTOR(logZ_mean);
  DATA_VECTOR(logZ_sd);
  DATA_INTEGER(restrict_mean);

  // PARAMETERS //
  // Fishing process parameters
  PARAMETER_VECTOR(logSDlogF); // log(SD) for random walk of logF: process error
  PARAMETER(rho_SDlogF); // correlation coefficient for random walk of logF

  // Recruitment process parameters
  PARAMETER_VECTOR(rec_loga); // log(a) for SR relation
  PARAMETER_VECTOR(rec_logb); // log(b) for SR relation
  PARAMETER_VECTOR(rec_logSD); // log(SD) for SR relation: process error
  PARAMETER(trans_rho); // correlation coefficient for recruitment deviance

  // Catch observation parameter
  PARAMETER_VECTOR(logSDlogC); // log(SD) for log(Catch): observation error

  // Index observation parameters
  PARAMETER_VECTOR(logQ); //Catchability
  PARAMETER_VECTOR(logSDcpue); // log(SD) for CPUE: observation error
  PARAMETER_VECTOR(beta); //nonlinear coefficient

  // Random effects
  PARAMETER_ARRAY(logN); //matrix of 2 rows x NYear colums
  PARAMETER_ARRAY(logF); //matrix of 2 rows x NYear colums

  array<Type> F(NStock,NYear);
  array<Type> N(NStock,NYear);
  array<Type> SSN(NStock,NYear); //spawning stock number
  array<Type> predC(NStock,NYear); //predicted catch weight

  Type rec_rho = (exp(trans_rho)-Type(1.0))/(exp(trans_rho)+Type(1.0));

  // fishing process
  for(int i=0;i<NStock;i++){
    for(int j=0;j<NYear;j++){
      N(i,j)=exp(logN(i,j));
      F(i,j)=exp(logF(i,j));
      SSN(i,j)=N(i,j)*exp(-M(i,j)-F(i,j));
    }
  }

  // predicted catch
  if(Pope==0){
    for(int i=0;i<NStock;i++){
      for(int j=0;j<NYear;j++){
        predC(i,j)=(F(i,j)/(F(i,j)+M(i,j)))*N(i,j)*(1-exp(-M(i,j)-F(i,j)))*Weight(i,j);
        predC(i,j)=scale_num_to_mass*predC(i,j);
      }
    }
  }else{
    for(int i=0;i<NStock;i++){
      for(int j=0;j<NYear;j++){
        predC(i,j)=exp(-Type(0.5)*M(i,j))*N(i,j)*(1-exp(-F(i,j)))*Weight(i,j);
        predC(i,j)=scale_num_to_mass*predC(i,j);
        }
      }
  }

  matrix<Type> fvar(NStock,NStock);
  for(int i=0;i<NStock;i++){
    fvar(i,i)=pow(exp(logSDlogF(SDlogF_key(i))),Type(2.0));
  } // diagonal elements

  for(int i=0;i<NStock;i++){
    for(int j=0;j<NStock;j++){
      if(j!=i){
        fvar(i,j)=rho_SDlogF*exp(logSDlogF(SDlogF_key(i)))*exp(logSDlogF(SDlogF_key(j)));
      }
    }
  } // offdiagonal elements

  using namespace density;  // using multivariate normal distribution
  MVNORM_t<Type> neg_log_densityF(fvar);  // with the var-cov matirix "fvar" for MVN

  Type nll=0;

  if(logF_diff==0){
    for(int i=1;i<NYear;i++){
      nll+=neg_log_densityF(logF.col(i)-logF.col(i-1));
      SIMULATE {
        logF.col(i)=logF.col(i-1)+neg_log_densityF.simulate();
        F.col(i)=exp(logF.col(i));
      }
    } // F process likelihood
  }else{
    for(int i=2;i<NYear;i++){
      nll+=neg_log_densityF(logF.col(i)-Type(2.0)*logF.col(i-1)+logF.col(i-2));
      SIMULATE {
        logF.col(i)=Type(2.0)*logF.col(i-1)-logF.col(i-2)+neg_log_densityF.simulate();
        F.col(i)=exp(logF.col(i));
      }
    } // F process likelihood
  }

  vector<Type> sum_logZ(NStock);
  vector<Type> sum_logZ_w(NStock);
  sum_logZ.fill(0.0);
  sum_logZ_w.fill(0.0);
  for(int i=0;i<NStock;i++){
    for(int j=0;j<NYear;j++){
      // sum_logZ(i) += logZ_w(i,j)*log(F(i,j)+M(i,j));
      // sum_logZ(i) += logZ_w(i,j)*(F(i,j)+M(i,j));
      sum_logZ(i) += logZ_w(i,j)*(F(i,j));
      sum_logZ_w(i) += logZ_w(i,j);
      if(restrict_mean==0){
        if(logZ_w(i,j)>Type(0.0)){
          // nll-=dnorm(log(F(i,j)+M(i,j)),logZ_mean(i),logZ_sd(i),true);
          nll-=dnorm(log(F(i,j)),logZ_mean(i),logZ_sd(i),true);
        }
      }
    }
    if(restrict_mean==1){
      if(sum_logZ_w(i)>Type(0.0)){
        // nll-=dnorm(sum_logZ(i)/sum_logZ_w(i),logZ_mean(i),logZ_sd(i),true);
        nll-=dnorm(log(sum_logZ(i)/sum_logZ_w(i)),logZ_mean(i),logZ_sd(i),true);
      }
    }
  }

  // stock-recruitment process
  vector<Type> rec_a=exp(rec_loga);
  vector<Type> rec_b=exp(rec_logb);
  vector<Type> rec_SD=exp(rec_logSD);

  for(int j=1;j<NYear;j++){
    vector<Type> pred_logN(NStock); // logN predicted from stock-recruitment relationship
    if(SR==0){ //Hockey-Stick
      for(int i=0;i<NStock;i++){ // Stock_ID
          // vector<Type> rec_pred_HS(2);
          // rec_pred_HS(0)=rec_a(reca_key(i,j))*rec_b(recb_key(i,j));
          // rec_pred_HS(1)=rec_a(reca_key(i,j))*SSN(i,j-1);
        pred_logN(i)=CppAD::CondExpLt(rec_b(recb_key(i,j)),SSN(i,j-1),rec_b(recb_key(i,j)),SSN(i,j-1));
        pred_logN(i)*=rec_a(reca_key(i,j));
        pred_logN(i)=log(pred_logN(i));
      }
    }
    if(SR==1){ //Beverton-Holt
      for(int i=0;i<NStock;i++){ // Stock_ID
          pred_logN(i)=rec_a(reca_key(i,j))*SSN(i,j-1)/(1+rec_b(recb_key(i,j))*SSN(i,j-1));
          pred_logN(i)=log(pred_logN(i));
      }
    }
    if(SR==2){ //Ricker
      for(int i=0;i<NStock;i++){ // Stock_ID
          pred_logN(i)=rec_a(reca_key(i,j))*SSN(i,j-1)*exp(-rec_b(recb_key(i,j))*SSN(i,j-1));
          pred_logN(i)=log(pred_logN(i));
      }
    }

    matrix<Type> rec_var(NStock,NStock);
    for(int i=0;i<NStock;i++){
      rec_var(i,i)=pow(rec_SD(recSD_key(i,j)),Type(2.0));
    } // diagonal elements

    for(int k=0;k<NStock;k++){
      for(int i=0;i<NStock;i++){
        if(i!=k){
          rec_var(k,i)=rec_rho*rec_SD(recSD_key(i,j))*rec_SD(recSD_key(k,j));
        } //offdiagonal elements
      }
    }

    using namespace density;  // using multivariate normal distribution
    MVNORM_t<Type> neg_log_densityR(rec_var);  // with the var-cov matirix "rec_var" for MVN
    nll+=neg_log_densityR(logN.col(j)-pred_logN); // Recruit process likelihood
    SIMULATE{
      logN.col(j)=pred_logN+neg_log_densityR.simulate();
      N.col(j)=exp(logN.col(j));
    }
  }

  // Catch observation
  for(int i=0;i<NCatch;i++){
    nll-=dnorm(log(Catch(i)),log(predC(Catch_key(i,0),Catch_key(i,1))),exp(logSDlogC(Catch_key(i,2))),true);
    SIMULATE{
      Catch(i)=exp(rnorm(log(predC(Catch_key(i,0),Catch_key(i,1))),exp(logSDlogC(Catch_key(i,2)))));
    }
  }

  // Index observation
  vector<Type> predIndex(NIndex);
  for(int i=0;i<NIndex;i++){
    predIndex(i)=exp(logQ(Index_key(i,2)))*pow(N(Index_key(i,0),Index_key(i,1)),beta(Index_key(i,4)));
    nll-=dnorm(log(Index(i)),log(predIndex(i)),exp(logSDcpue(Index_key(i,3))),true);
    SIMULATE{
      Index(i)=exp(rnorm(log(predIndex(i)),exp(logSDcpue(Index_key(i,3)))));
    }
  }

  SIMULATE{
    REPORT(logF);
    REPORT(logN);
    REPORT(Catch);
    REPORT(Index);
    REPORT(F);
    REPORT(N);
    REPORT(predC);
    REPORT(SSN);
  }

  ADREPORT(logF);
  ADREPORT(logN);
  ADREPORT(F);
  ADREPORT(predC);
  ADREPORT(N);
  ADREPORT(SSN);
  // ADREPORT(rec_a);
  // ADREPORT(rec_b);
  // ADREPORT(rec_SD);

  return nll;
}
