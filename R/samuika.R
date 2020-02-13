
#' State-Space Assessment model Used for surume-IKA (samuika) ----
#'
#' @import TMB
#' @importFrom TMB MakeADFun
#' @param catch_data time-series of catch
#' @param weight_data time-series of weight
#' @param index_data time-series of indices
#' @param SR stock-recruit function ("BH", "RI", or "HS")
#' @param regime_year years when a regime shift occurred
#' @param regime_par parameters that changed by the regime shift (\code{c("a","b","sd")})
#' @param regime_key KEY representing regime (e.g., \code{c(0,1,0)} indicates the first and the third regimes were identical)
#' @param stock_shared_par parameter(s) that were shared between stocks (\code{c("a","b","sd","SDlogF","SDlogC")})
#' @param fixed_par parameter(s) fixed at the initial value(s) (\code{c("a","b","sd","rec_rho","SDlogF","rho_SDlogF","SDlogC","q","SDlogCPUE")})
#' @param add_cpue VECTOR of additional cpue time-series
#' @param add_cpue_info column 1: CPUE_ID, 2: Stock_ID, 3: Year_ID
#' @param add_cpue_tday VECTOR of days since the beginning of fishing season
#' @param add_cpue_covariate array used as covariate(s) for additional CPUE
#' @param add_cpue_error_type error distribution for fitting additional CPUE; 0: lonormal, 1:log-Laplace, 2: gamma, 3: normal, 4: Laplace
#' @param add_cpue_all fitting additional CPUE not separately per-year (0: FALSE, 1 (default): using mean(F+Z), 2: using mean(log(F+Z))
#' @param logZ_mean prior mean of logF (NOTICE: Not logZ!!)
#' @param logZ_sd prior SD of logF (NOTICE: Not logZ!!)
#' @param logZ_weight matrix for which years are included when using prior of \code{logZ_mean} and \code{logZ_sd} or \code{add_cpue}
#' @param fish_days duration of fishing season when using \code{add_cpue} (default: 180)
#' @param M natural mortality coefficient (default: 0.6)
#' @param Pope whether the Pope approximation is used (TRUE) or not (FALSE: default)
#' @param scale_num_to_mass scaling multiplier in conversion of number to mass
#' @param bias_correct bias correct option in \code{sdreport}
#' @param fixed_par which parameters among \code{c("a","b","sd","rec_rho","SDlogF","rho_SDlogF","SDlogC","q","SDlogCPUE")} are fixed at their initial values
#' @param map_add added factor to map other than fixed_par
#'
#' @encoding UTF-8
#' @export
samuika = function(
  catch_data,
  weight_data,
  index_data,
  SR = "BH",
  regime_year = NULL,
  regime_par = c("a","b","sd")[c(1:2)],
  regime_key = 0:length(regime_year),
  stock_shared_par = c("a","b","sd","SDlogF","SDlogC")[c(1:5)],
  SDlogCPUE_key = rep(0,length(unique(index_data$Index_ID))),
  beta_key = rep(0,length(unique(index_data$Index_ID))),
  beta_fix = 1,
  logZ_mean = NULL,
  logZ_sd = NULL,
  logZ_weight = NULL,
  restrict_mean = FALSE,
  add_cpue = NULL,
  add_cpue_info = NULL,
  add_cpue_tday = NULL,
  fish_days = 180,
  add_cpue_covariate = NULL,
  add_cpue_SDkey = NULL,
  add_cpue_error_type = 0,
  add_cpue_all = 1, # 0: FALSE, 1: modeling E(F+M), 2: modeling log(F+M)
  Fprocess_remove_year = NULL,
  # add_fixed_F = 0,
  # lambda = 0,
  p0_list = NULL,
  M = 0.6,
  Pope = FALSE,
  scale_num_to_mass = 0.1,
  bias_correct = TRUE,
  bias_correct_control = list(sd=FALSE),
  fixed_par = c("a","b","sd","rec_rho","SDlogF","rho_SDlogF","SDlogC","q","SDlogCPUE")[6],
  map_add = NULL,
  logF_diff = 0, #
  SDlogF_init = 0.2,
  rho_SDlogF_init = 0.0,
  reca_init = c(3.5),
  recb_init = ifelse(SR=="HS",14,0.04),
  recSD_init = 0.4,
  rec_rho_init = 0,
  SDlogC_init = 0.2,
  q_init = c(0.22,0.5),
  SDcpue_init = 0.2,
  N_init = 23.7,
  F_init = 0.38,
  silent = FALSE,
  maxrep = 100, # maximum replicate for SR="HS",
  HS_restrict = TRUE,
  nlminb_control = list(eval.max=2000,iter.max=2000)
) {

  argname <- ls()
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname

  # cat("One can backcast and/or forecast population dynamics if extending 'weight_data' \n")

  if (length(SDlogCPUE_key) != length(unique(index_data$Index_ID))) {
    stop("length(SDlogCPUE_key does not match the number of indices")
  }
  if (length(unique(beta_key)) != length(beta_fix)) {
    stop("length(beta_fix) does not match the number of beta_key")
  }

  index_data = index_data %>%
    mutate(Index_ID = Index_ID-min(Index_ID),Stock_ID = Stock_ID-min(Stock_ID))
  catch_data = catch_data %>%
    mutate(Stock_ID = Stock_ID-min(Stock_ID))

  weight_mat = weight_data %>% select(-Year) %>% as.matrix() %>% t()
  NYear=ncol(weight_mat)
  NStock=nrow(weight_mat)

  M_mat = matrix(M, nrow=nrow(weight_mat), ncol=ncol(weight_mat))
  if ("SDlogF" %in% stock_shared_par) {
    SDlogF_key = rep(0,NStock)
  } else {
    SDlogF_key = 1:NStock-1
  }

  SR_tmb=if_else(SR=="HS",0,if_else(SR=="BH",1,2))

  start_year = min(weight_data$Year)
  reca_key <- recb_key <- recSD_key <- matrix(0, nrow=nrow(weight_mat), ncol=ncol(weight_mat))
  if (!is.null(regime_year)) {
    if (length(regime_key) < length(regime_year)+1) {
      stop("'length(regime_key) < length(regime_year)+1' should be satisfied")
    }
    if (regime_key[1]!=0) {
      stop("'regime_key' must start at zero")
    }
    regime_iy = regime_year-start_year+1
    if ("a" %in% regime_par) {
      for (i in 1:length(regime_year)) {
        reca_key[,regime_iy[i]:ncol(reca_key)] <- regime_key[i+1]
      }
    }
    if ("b" %in% regime_par) {
      for (i in 1:length(regime_year)) {
        recb_key[,regime_iy[i]:ncol(recb_key)] <- regime_key[i+1]
      }
    }
    if ("sd" %in% regime_par) {
      for (i in 1:length(regime_year)) {
        recSD_key[,regime_iy[i]:ncol(recSD_key)] <- regime_key[i+1]
      }
    }
  }

  if (NStock > 1) {
    if (!("a" %in% stock_shared_par)) {
      for(i in 2:NStock) reca_key[i,] <- max(reca_key[i-1,])+1+reca_key[i,]
    }

    if (!("b" %in% stock_shared_par)) {
      for(i in 2:NStock) recb_key[i,] <- max(recb_key[i-1,])+1+recb_key[i,]
    }

    if (!("sd" %in% stock_shared_par)) {
      for(i in 2:NStock) recSD_key[i,] <- max(recSD_key[i-1,])+1+recSD_key[i,]
    }
  }

  NCatch = nrow(catch_data)
  Catch = catch_data$Catch_biomass %>% as.numeric()
  Catch_key = catch_data %>%
    mutate(iy = Year-start_year, key= 0) %>%
    select(Stock_ID, iy, key)
  if (!("SDlogC" %in% stock_shared_par)) {
    Catch_key = Catch_key %>% mutate(key = Stock_ID-min(Stock_ID))
  }
  Catch_key = as.matrix(Catch_key)

  NIndex = nrow(index_data)
  Index = index_data$Index %>% as.numeric()
  Index_key = index_data %>%
    mutate(iy = Year-start_year, q_key = Index_ID-min(Index_ID)) %>%
    mutate(sd_key = SDlogCPUE_key[Index_ID+1],beta_key = beta_key[Index_ID+1]) %>%
    select(Stock_ID, iy, q_key, sd_key, beta_key)
  Index_key = as.matrix(Index_key)

  if (is.null(logZ_mean)) {
    if (!is.null(logZ_sd)) {
      message("'logZ_sd' is not usable when 'logZ_mean' is null")
    }
    if (!is.null(logZ_weight)) {
      if (!is.null(add_cpue)) {
        message("'logZ_weight' is used for add_cpue")
      } else {
        message("'logZ_weight' is not usable when 'logZ_mean' is null")
      }
    }
    logZ_MEAN = rep(0,NStock)
    logZ_SD = rep(1,NStock)
    logZ_W = matrix(0,ncol=ncol(weight_mat),nrow=nrow(weight_mat))
  } else {
    if (is.null(logZ_sd)) {
      stop("Set 'logZ_sd' when using 'logZ_mean'!!!")
    }
    if (length(logZ_mean) != NStock) {
      stop("The length of 'logZ_mean' does not match the number of stocks")
    }
    logZ_MEAN = logZ_mean
    logZ_SD = logZ_sd
    logZ_MEAN[is.na(logZ_MEAN)] <- 0
    logZ_SD[is.na(logZ_SD)] <- 1

    if (is.null(logZ_weight)) {
      message("Restriction to logZ is applied to all years and all stocks")
      logZ_W = matrix(1,ncol=ncol(weight_mat),nrow=nrow(weight_mat))
    } else {
      logZ_W = logZ_weight %>% select(-Year) %>% as.matrix() %>% t()
      if (ncol(logZ_W) != NYear || ncol(logZ_W) != NYear) {
        stop("The size of logZ_weight is NOT appropriate")
      }
    }
  }
  trans_rho_init = log((1-rec_rho_init)/(1+rec_rho_init))

  if (is.null(add_cpue)) {
    use_add_cpue = 0
    add_cpue2 = c(1)
    add_cpue_info2 = matrix(0,ncol=3,nrow=1)
    add_cpue_tday = c(1)
    add_cpue_covariate = matrix(0,ncol=1,nrow=1)
    add_cpue_SD_key = c(1)
  } else {
    if (!is.null(logZ_mean)) {
      stop("Either logZ_mean or add_cpue should be selected!")
    }
    use_add_cpue = 1
    add_cpue2 = add_cpue
    if (is.null(add_cpue_info)) stop("add_cpue_info is neccessary when add_cpue is used")
    if (is.null(add_cpue_tday)) stop("add_cpue_tday is neccessary when add_cpue is used")
    if (nrow(add_cpue_info) != length(add_cpue)) {
      stop("nrow of add_cpue_info does not match length of add_cpue")
    }
    if (ncol(add_cpue_info) != 3) {
      stop("ncol of add_cpue_info should be three (CPUE_ID, Stock_ID, Year_ID)!!!")
    }
    if (length(add_cpue_tday) != length(add_cpue)) {
      stop("length of add_cpue_tday does not match length of add_cpue")
    }
    if (sum(!(add_cpue_info[,2] %in% Catch_key[,1])) != 0) {
      stop("Stock_ID (column 2) is not correct")
    }
    add_cpue_info2 = data.frame(CPUE_ID = add_cpue_info[,1]-min(add_cpue_info[,1]),
                                Stock_ID = add_cpue_info[,2],
                                iy = add_cpue_info[,3]-start_year) %>% as.matrix()
    if (is.null(add_cpue_covariate)) {
      add_cpue_covariate = matrix(0,ncol=1,nrow=length(add_cpue2))
    }
    if (is.null(add_cpue_SDkey)) add_cpue_SD_key = 1:(max(add_cpue_info2[,1]+1))-1
    if (add_cpue_all>0) {
      use_add_cpue = add_cpue_all + 1
      if (!is.null(logZ_weight)) {
        logZ_W = logZ_weight %>% select(-Year) %>% as.matrix() %>% t()
      } else {
        for (i in 1:nrow(logZ_W)) {
          for (j in 1:ncol(logZ_W)) {
            if ((i-1) %in% add_cpue_info2[,2]) {
              if ((j-1) %in% add_cpue_info2[,3]) logZ_W[i,j] <- 1
            }
          }
        }
      }
    }
  }
  F_incl_w = rep(1,nrow(weight_data))
  if (!is.null(Fprocess_remove_year)) {
    F_incl_w[weight_data$Year %in% Fprocess_remove_year] <- 0
  }

  data_list = list(NYear=NYear,NStock=NStock,M=M_mat,Weight=weight_mat,
                   SDlogF_key=SDlogF_key,logF_diff=logF_diff,
                   F_incl_w = F_incl_w,
                   # lambda = lambda,
                   # add_fixed_F=as.numeric(add_fixed_F),
                   SR=SR_tmb,reca_key=reca_key,recb_key=recb_key,recSD_key=recSD_key,
                   NCatch=NCatch,Catch=Catch,Pope=as.numeric(Pope),Catch_key=Catch_key,scale_num_to_mass=scale_num_to_mass,
                   NIndex=NIndex,Index=Index,Index_key=Index_key,
                   logZ_mean=logZ_MEAN,logZ_sd=logZ_SD,logZ_w=logZ_W,restrict_mean=as.numeric(restrict_mean),
                   use_add_cpue=use_add_cpue, fish_days=fish_days,add_cpue_info=add_cpue_info2,
                   add_cpue=add_cpue2,add_cpue_tday=add_cpue_tday,add_cpue_covariate=add_cpue_covariate,
                   add_cpue_SD_key = add_cpue_SD_key,cpue_add_error_type = add_cpue_error_type)

  if (is.null(p0_list)) {
    param_init = list(
      logSDlogF = rep_len(log(SDlogF_init),max(SDlogF_key)+1),
      rho_SDlogF = rho_SDlogF_init,
      # add_logF = rep(log(0),sum(F_incl_w==0)),
      rec_loga = rep_len(log(reca_init),max(reca_key)+1),
      rec_logb = rep_len(log(recb_init),max(recb_key)+1),
      rec_logSD = rep_len(log(recSD_init),max(recSD_key)+1),
      trans_rho = trans_rho_init,
      logSDlogC = rep_len(log(SDlogC_init),max(Catch_key[,3])+1),
      logQ = rep_len(log(q_init),max(Index_key[,3])+1),
      logSDcpue = rep_len(log(SDcpue_init), max(Index_key[,4])+1),
      beta = rep(1,length(beta_fix)),
      logN = matrix(log(23.7),nrow=dim(weight_mat)[1],ncol=dim(weight_mat)[2]),
      logF = matrix(log(0.38),nrow=dim(weight_mat)[1],ncol=dim(weight_mat)[2]),
      logQ_add = rep(log(1),length(unique(add_cpue_info2[,1]))),
      logSDcpue_add = rep(log(1),length(unique(add_cpue_SD_key))),
      # logSDcpue_add = rep(log(0.001),length(unique(add_cpue_SD_key))),
      alpha = matrix(0,ncol=ncol(add_cpue_covariate),nrow=length(unique(add_cpue_info2[,1])))
    )
  } else {
    param_init = p0_list
  }

  map = list()
  # map$logSDcpue_add = rep(factor(NA), length(param_init$logSDcpue_add))
  if (use_add_cpue==0) {
    map$logQ_add = rep(factor(NA),length(param_init$logQ_add))
    map$logSDcpue_add = rep(factor(NA),length(param_init$logSDcpue_add))
    map$alpha = rep(factor(NA),ncol(param_init$alpha)*nrow(param_init$alpha))
  }
  if (!is.null(fixed_par)) {
    message("Parameter in 'fixed_par' is fixed at the initial value")

    if ("a" %in% fixed_par) {
      map$rec_loga = rep(factor(NA),max(reca_key)+1)
      cat("'rec_loga' is fixed at", rep_len(log(reca_init),max(reca_key)+1),"\n")
    }
    if ("b" %in% fixed_par) {
      map$rec_logb = rep(factor(NA),max(recb_key)+1)
      cat("'rec_logb' is fixed at", rep_len(log(recb_init),max(recb_key)+1),"\n")
    }
    if ("sd" %in% fixed_par) {
      map$rec_logSD = rep(factor(NA),max(recSD_key)+1)
      cat("'rec_logSD' is fixed at", rep_len(log(recSD_init),max(recSD_key)+1),"\n")
    }
    if ("rec_rho" %in% fixed_par) {
      map$trans_rho = factor(NA)
      cat("'rec_rho' is fixed at", rec_rho_init,"\n")
    }
    if ("SDlogF" %in% fixed_par) {
      map$logSDlogF = rep(factor(NA),max(SDlogF_key)+1)
      cat("'logSDlogF' is fixed at", rep_len(log(SDlogF_init),max(SDlogF_key)+1),"\n")
    }
    if ("rho_SDlogF" %in% fixed_par) {
      map$rho_SDlogF = factor(NA)
      cat("'rho_SDlogF' is fixed at", rho_SDlogF_init,"\n")
    }
    if ("SDlogC" %in% fixed_par) {
      map$logSDlogC = rep(factor(NA),max(Catch_key[,3])+1)
      cat("'logSDlogC' is fixed at", rep_len(log(SDlogC_init),max(Catch_key[,3])+1),"\n")
    }
    if ("q" %in% fixed_par) {
      map$logQ = rep(factor(NA),max(Index_key[,3])+1)
      cat("'logQ' is fixed at", rep_len(log(q_init),max(Index_key[,3])+1),"\n")
    }
    if ("SDlogCPUE" %in% fixed_par) {
      map$logSDcpue = rep(factor(NA),max(Index_key[,4])+1)
      cat("'logSDcpue' is fixed at", rep_len(log(SDcpue_init),max(Index_key[,4])+1),"\n")
    }
  }
  if (NStock == 1) {
    map$trans_rho = factor(NA)
    cat("'rec_rho' is automatically fixed at", rec_rho_init,"because the number of stocks is one","\n")
  }

  if (any(is.na(beta_fix))) {
    map$beta = 0:(length(beta_fix)-1)
    map$beta[!is.na(beta_fix)] <- NA
    map$beta = factor(map$beta)
    } else {
    map$beta = rep(factor(NA),length(beta_fix))
    }

  if (!is.null(map_add)) {
    tmp = make_named_list(map_add)
    for(i in 1:length(tmp$map_add)) {
      map[[names(tmp$map_add)[i]]] <- map_add[[i]]
    }
  }
  # if (add_fixed_F==0) {
  #   map$add_logF = rep(factor(NA),length(param_init$add_logF))
  # }

  f = TMB::MakeADFun(data=data_list, parameters = param_init, random = c("logN","logF"), DLL="samuika",
                     map = map, silent=silent)
  fit = try(stats::nlminb(f$par, f$fn, f$gr,control = nlminb_control))
  if (class(fit)=="try-error") stop("ERROR in 'nlminb'")
  rep = try(TMB::sdreport(f, bias.correct = bias_correct, bias.correct.control = bias_correct_control))

  RES = list()
  RES$input = arglist
  RES$data_list = data_list
  RES$param_init = param_init
  RES$map = map
  RES$obj = f
  RES$opt = fit
  RES$par_list = f$env$parList(fit$par)
  RES$loglik <- loglik <- -(fit$objective)
  RES$N <- N <- NCatch+NIndex
  RES$npar <- npar <- length(fit$par)
  RES$AIC = -2*loglik + 2*npar
  RES$AICc = RES$AIC+2*npar*(npar+1)/(N-npar-1)
  RES$BIC = -2*loglik + npar*log(N)

  if (class(rep)=="try-error") {
    warning("ERROR in 'sdreport'")
  } else {
    RES$rep = rep
    rep_summary = summary(rep,"report")
    RES$rep_summary = rep_summary

    rec_pars = list()
    rec_pars$a = exp(f$env$parList(fit$par)$rec_loga)
    rec_pars$b = exp(f$env$parList(fit$par)$rec_logb)
    rec_pars$sd = exp(f$env$parList(fit$par)$rec_logSD)
    RES$rec_pars = rec_pars
    trans_rho = f$env$parList(fit$par)$trans_rho
    RES$rec_rho = (exp(trans_rho)-1)/(exp(trans_rho)+1)
    RES$SDlogF = exp(f$env$parList(fit$par)$logSDlogF)
    RES$rho_SDlogF = f$env$parList(fit$par)$rho_SDlogF
    RES$SDlogC = exp(f$env$parList(fit$par)$logSDlogC)
    RES$q = exp(f$env$parList(fit$par)$logQ)
    RES$SDcpue = exp(f$env$parList(fit$par)$logSDcpue)
    RES$beta = f$env$parList(fit$par)$beta
    RES$cor = stats::cov2cor(rep$cov.fixed)

    if (bias_correct) {
      VALUE = rep$unbiased$value
    } else {
      VALUE = rep$value
    }

    if (!is.null(add_cpue)) {
      RES$alpha = f$env$parList(fit$par)$alpha
      RES$q_add = exp(f$env$parList(fit$par)$logQ_add)
      RES$SDcpue_add = exp(f$env$parList(fit$par)$logSDcpue_add)
      RES$pred_cpue_add = as.numeric(exp(VALUE[names(VALUE)=="pred_log_cpue_add"]))
      RES$resid_cpue_add = as.numeric(VALUE[names(VALUE)=="resid_add"])
      add_cpue_table = add_cpue_info
      if (!is.null(add_cpue_covariate)) {
        add_cpue_table = bind_cols(add_cpue_table, tibble(cpue=add_cpue,tday = add_cpue_tday,covariate=add_cpue_covariate))
      }
      add_cpue_table = add_cpue_table %>%
        mutate(pred_cpue = RES$pred_cpue_add,resid = RES$resid_cpue_add)
      RES$add_cpue_table = add_cpue_table
      if (add_cpue_all > 0) {
        RES$meanZ = as.numeric(VALUE[names(VALUE)=="sum_logZ"]/VALUE[names(VALUE)=="sum_logZ_w"])
        if (add_cpue_all == 2) RES$meanZ = exp(RES$meanZ)
      }
    }

    N_est = matrix(VALUE[names(VALUE)=="N"],nrow=NStock)
    SSN_est = matrix(VALUE[names(VALUE)=="SSN"],nrow=NStock)
    F_est = matrix(VALUE[names(VALUE)=="F"],nrow=NStock)
    predC_est = matrix(VALUE[names(VALUE)=="predC"],nrow=NStock)
    predN_est = matrix(VALUE[names(VALUE)=="pred_N"],nrow=NStock)
    predN_est[predN_est==0] <- NA
    # rec_resid = matrix(VALUE[names(VALUE)=="rec_resid"],nrow=NStock)
    # rec_resid[rec_resid==0] <- NA

    Ncol = ifelse(bias_correct_control$sd,4,2)
    N_se = matrix(rep_summary[rownames(rep_summary)=="N",Ncol],nrow=NStock)
    SSN_se = matrix(rep_summary[rownames(rep_summary)=="SSN",Ncol],nrow=NStock)
    F_se = matrix(rep_summary[rownames(rep_summary)=="F",Ncol],nrow=NStock)
    predC_se = matrix(rep_summary[rownames(rep_summary)=="predC",Ncol],nrow=NStock)
    predN_se = matrix(rep_summary[rownames(rep_summary)=="pred_N",Ncol],nrow=NStock)

    B_est = N_est*weight_mat*scale_num_to_mass
    SSB_est = SSN_est*weight_mat*scale_num_to_mass
    B_se = N_se*weight_mat*scale_num_to_mass
    SSB_se = SSN_se*weight_mat*scale_num_to_mass

    colnames(N_est) <- colnames(SSN_est) <- colnames(F_est) <- colnames(predC_est) <-
      colnames(N_se) <- colnames(SSN_se) <- colnames(F_se) <- colnames(predC_se) <-
      colnames(B_est) <- colnames(SSB_est) <- colnames(B_se) <- colnames(SSB_se) <-
      colnames(predN_est) <- colnames(predN_se) <- 1:NYear+start_year-1
    rownames(N_est) <- rownames(SSN_est) <- rownames(F_est) <- rownames(predC_est) <-
      rownames(N_se) <- rownames(SSN_se) <- rownames(F_se) <- rownames(predC_se) <-
      rownames(B_est) <- rownames(SSB_est) <- rownames(B_se) <- rownames(SSB_se) <-
      rownames(predN_est) <- rownames(predN_se) <- 1:NStock-1

    RES$N_est = N_est
    RES$SSN_est = SSN_est
    RES$F_est = F_est
    RES$C_est = predC_est
    RES$B_est = B_est
    RES$SSB_est = SSB_est
    RES$predN_est = predN_est

    RES$N_se = N_se
    RES$SSN_se = SSN_se
    RES$F_se = F_se
    RES$C_se = predC_se
    RES$B_se = B_se
    RES$SSB_se = SSB_se
    RES$predN_se = predN_se

    reca_est <- recb_est <- recSD_est <-
      SDlogF_est <- SDlogC_est <- matrix(0, nrow=NStock, ncol=NYear)
    for (i in 1:NStock) {
      for (j in 1:NYear) {
        reca_est[i,j] <- rec_pars$a[reca_key[i,j]+1]
        recb_est[i,j] <- rec_pars$b[recb_key[i,j]+1]
        recSD_est[i,j] <- rec_pars$sd[recSD_key[i,j]+1]
        SDlogF_est[i,j] <- RES$SDlogF[SDlogF_key[i]+1]
        SDlogC_est[i,j] <- ifelse(nrow(filter(Catch_key %>% as_tibble(), Stock_ID==i-1, iy==j-1))==0,NA,
                                  RES$SDlogC[select(filter(Catch_key %>% as_tibble(), Stock_ID==i-1, iy==j-1),key) %>% unlist() +1])
      }
    }

    colnames(reca_est) <- colnames(recb_est) <- colnames(recSD_est) <-
      colnames(SDlogF_est) <- colnames(SDlogC_est) <- 1:NYear+start_year-1
    rownames(reca_est) <- rownames(recb_est) <- rownames(recSD_est) <-
      rownames(SDlogF_est) <- rownames(SDlogC_est) <- 1:NStock-1

    RES$reca_est = reca_est
    RES$recb_est = recb_est
    RES$recSD_est = recSD_est
    RES$SDlogF_est = SDlogF_est
    RES$SDlogC_est = SDlogC_est

    join_gather = function(X,Y) {
      full_join(X,
                mutate(t(Y) %>% as_tibble(), Year=as.numeric(colnames(Y))) %>%
                  gather(key=Stock_ID, value = value,-Year) %>%
                  mutate(Stock_ID = as.numeric(Stock_ID)),
                by = c("Stock_ID","Year"))
    }

    Summary_PopDyn = full_join(catch_data, Catch_key %>% as_tibble() %>% mutate(Year=iy+start_year), by = c("Stock_ID", "Year")) %>%
      as_tibble() %>% select(-iy, -key)
    Summary_PopDyn = join_gather(Summary_PopDyn, N_est) %>%
      rename(Stock_number = value)
    Summary_PopDyn = join_gather(Summary_PopDyn, SSN_est) %>%
      rename(Spawning_number = value)
    Summary_PopDyn = join_gather(Summary_PopDyn, F_est) %>%
      rename(F = value)
    Summary_PopDyn = join_gather(Summary_PopDyn, predC_est) %>%
      rename(Catch_est = value)
    Summary_PopDyn = join_gather(Summary_PopDyn, B_est) %>%
      rename(Stock_biomass = value)
    Summary_PopDyn = join_gather(Summary_PopDyn, SSB_est) %>%
      rename(Spawning_biomass = value)
    Summary_PopDyn = join_gather(Summary_PopDyn, reca_est) %>%
      rename(rec_a = value)
    Summary_PopDyn = join_gather(Summary_PopDyn, recb_est) %>%
      rename(rec_b = value)
    Summary_PopDyn = join_gather(Summary_PopDyn, recSD_est) %>%
      rename(rec_sd = value)
    Summary_PopDyn = join_gather(Summary_PopDyn, SDlogF_est) %>%
      rename(SDlogF = value)
    Summary_PopDyn = join_gather(Summary_PopDyn, predN_est) %>%
      rename(Stock_number_pred = value)
    Summary_PopDyn = Summary_PopDyn %>%
      mutate(rec_resid = log(Stock_number/Stock_number_pred)) %>%
      arrange(Stock_ID, Year)
    Summary_PopDyn = join_gather(Summary_PopDyn, SDlogC_est) %>%
      rename(SDlogC = value)
    Summary_PopDyn = Summary_PopDyn %>%
      mutate(Catch_resid = log(Catch_biomass/Catch_est)) %>%
      arrange(Stock_ID, Year)

    RES$Summary_PopDyn = Summary_PopDyn
    Summary_index = full_join(index_data %>% as_tibble(),
                              select(mutate(Index_key %>% as_tibble(), Year = iy+start_year),-iy),
                              by = c("Stock_ID","Year")) %>%
      mutate(q=RES$q[q_key+1],sd=RES$SDcpue[sd_key+1],beta=RES$beta[beta_key+1]) %>%
      select(-q_key,-sd_key,-beta_key) %>%
      left_join(select(Summary_PopDyn,Stock_ID,Year,Stock_number),
                by = c("Stock_ID","Year")) %>%
      mutate(pred_Index = q*Stock_number^beta) %>%
      mutate(residual = log(Index/pred_Index))

    RES$Summary_index = Summary_index
  }

  invisible(RES)
}


#' Retrospective analysis for samuika
#' @import tidyr
#' @importFrom tidyr %>%
#' @param samuika_res samuika object
#' @param n number of maximum removed years
#' @encoding UTF-8
#' @export
retro_samuika = function(samuika_res, n=5, first_remove_catch_year = NULL, first_remove_index_year = NULL,use_p0 = FALSE) {
  RES=list()
  RES$full <- res.c <- samuika_res

  Summary_PopDyn = samuika_res$Summary_PopDyn %>%
    mutate(retro_year = 0)

  for(i in 1:n) {
    input_tmp = res.c$input
    input_tmp$bias_correct_control = list(sd=FALSE)
    if (is.null(first_remove_catch_year)) {
      max_catch_year = max(input_tmp$catch_data$Year)
    } else {
      max_catch_year = first_remove_catch_year-(i-1)
    }
    if (is.null(first_remove_index_year)) {
      max_index_year = max(input_tmp$index_data$Year)
    } else {
      max_index_year = first_remove_index_year-(i-1)
    }

    catch_data2 = dplyr::filter(input_tmp$catch_data,Year<max_catch_year)
    index_data2 = dplyr::filter(input_tmp$index_data,Year<max_index_year)
    input_tmp$catch_data = catch_data2
    input_tmp$index_data = index_data2

    # if (!is.null(input_tmp$logZ_weight)) {
    #   logZ_weight2 = dplyr::filter(input_tmp$logZ_weight, Year < max(Year))
    #   input_tmp$logZ_weight = logZ_weight
    # }
    if (use_p0) {
      input_tmp$p0_list = res.c$par_list
    }

    res.c = try(do.call(samuika, input_tmp))
    if (class(res.c)=="try-error") {
      stop(paste0("Error in ",i,"th trial"))
    }
    RES$retro[[i]] = res.c
    Summary_PopDyn = Summary_PopDyn %>%
      bind_rows(res.c$Summary_PopDyn %>% mutate(retro_year = i))
  }
  RES$Summary_PopDyn = Summary_PopDyn

  invisible(RES)
}


#' Calculation of Mohn's rho from retrospective analysis of samuika
#' @import dplyr
#' @importFrom dplyr tibble
#' @param retro_res retro_samuika object
#' @param first_eval_year firstly-evaluated year
#' @param lag_year time-lag in evaluation of Mohn's rho
#' @encoding UTF-8
#' @export
get_mohn = function(retro_res, first_eval_year = NULL, lag_year = 1) {
  n = length(retro_res$retro)
  if (is.null(first_eval_year)) {
    first_eval_year = max(retro_res$full$Summary_PopDyn$Year)
  }
  N_res <- B_res <- SSN_res <- SSB_res <- F_res <- C_res <- c()
  for (i in lag_year:n) {
    eval_year = as.character(first_eval_year-(i-1))
    N_res = cbind(N_res,retro_res$retro[[i]]$N_est[,eval_year]/retro_res$full$N_est[,eval_year]-1)
    B_res = cbind(B_res,retro_res$retro[[i]]$B_est[,eval_year]/retro_res$full$B_est[,eval_year]-1)
    SSN_res = cbind(SSN_res,retro_res$retro[[i]]$SSN_est[,eval_year]/retro_res$full$SSN_est[,eval_year]-1)
    SSB_res = cbind(SSB_res,retro_res$retro[[i]]$SSB_est[,eval_year]/retro_res$full$SSB_est[,eval_year]-1)
    F_res = cbind(F_res,retro_res$retro[[i]]$F_est[,eval_year]/retro_res$full$F_est[,eval_year]-1)
    C_res = cbind(C_res,retro_res$retro[[i]]$C_est[,eval_year]/retro_res$full$C_est[,eval_year]-1)
  }
  Summary = tibble(N=rowMeans(N_res),B=rowMeans(N_res),
                   SSN=rowMeans(SSN_res),SSB=rowMeans(SSB_res),
                   F=rowMeans(F_res),Catch=rowMeans(C_res))
  Summary = Summary %>%
    mutate(Stock_ID = unique(retro_res$full$input$catch_data$Stock_ID)) %>%
    select(Stock_ID,everything())

  RES = list(Summary=Summary,N_res=N_res,B_res=B_res,SSN_res=SSN_res,
             SSB_res=SSB_res,F_res=F_res,C_res=C_res)

}
