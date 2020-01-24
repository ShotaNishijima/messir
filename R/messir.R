
##################################################################################
### Management and Evaluation for Sustainability of Surume-Ika with R (MESSIR) ###
##################################################################################

#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @importFrom dplyr tibble
#' @importFrom tidyr %>%
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr bind_rows
#' @importFrom dplyr bind_cols
#' @importFrom tibble as_tibble
#'
NULL

#' future simulation for Surume-Ika
#' @param assess_data assessment data
#' @param Fcurrent current fishing mortality coefficient
#' @encoding UTF-8
#' @export
future_sim = function(
  assess_data,
  Fcurrent,
  multi=1, #Fcurrentにかける乗数,
  sim_year=50,
  nsim=1000,
  Fcurrent_year=NULL,
  # pre_catch=NULL,
  SR="HS", # or "BH" or "RI"
  rec_arg=list(a=NULL,b=NULL,sd=0,rho=0,resid=NULL),
  bias_correct=TRUE,
  weight,
  num_to_mass_scale=1/1000,
  M=0.6,
  Pope=TRUE,
  seed=12345,
  HCR = FALSE,
  SBrefs = c("SBlimit","SBban"), #Biological Rererence points
  Ftarget = NULL, # Ftarget = beta*Fmsy
  array_to_tibble = FALSE,
  sim_rec_resid = NULL, # matrix of sim_year X nsim
  scenario_name = NULL
) {

  argname = ls()
  arglist = lapply(argname,function(x) eval(parse(text=x)))
  names(arglist) = argname

  # SR function
  if (SR == "HS") SRF = function(x,a=rec_arg$a,b=rec_arg$b) ifelse(x>b,b*a,x*a)
  if (SR == "BH") SRF = function(x,a=rec_arg$a,b=rec_arg$b) a*x/(1+b*x)
  if (SR == "RI") SRF = function(x,a=rec_arg$a,b=rec_arg$b) a*x*exp(-b*x)

  rec_pred_by_SR <- rec_deviance_to_SR <- NA
  for (i in 2:nrow(assess_data)) {
    rec_pred_by_SR <- c(rec_pred_by_SR,SRF(assess_data$Spawning_number[i-1]))
    rec_deviance_to_SR <- c(rec_deviance_to_SR,log(assess_data$Stock_number[i]/rec_pred_by_SR[i]))
  }

  rec_pred_incl_AR <- c(NA,NA)
  rec_resid_excl_AR <- rec_deviance_to_SR[1:2]
  for (i in 3:nrow(assess_data)) {
    rec_pred_incl_AR <- c(rec_pred_incl_AR,rec_pred_by_SR[i]*exp(rec_arg$rho*rec_deviance_to_SR[i-1]))
    rec_resid_excl_AR <- c(rec_resid_excl_AR,log(assess_data$Stock_number[i]/rec_pred_incl_AR[i]))
  }

  all_data = assess_data %>%
    mutate(Catch_number=Catch_biomass/(Weight*num_to_mass_scale),
           rec_pred_by_SR=rec_pred_by_SR,
           rec_deviance_to_SR=rec_deviance_to_SR,
           rec_pred_incl_AR=rec_pred_incl_AR,
           rec_resid_excl_AR=rec_resid_excl_AR) %>% as_tibble()

  if (HCR){
    if (is.null(SBrefs)) {
      stop("Set SBrefs = C(SBlimit, SBban)")
      } else {SNrefs=SBrefs/(weight*num_to_mass_scale)}
    if (is.null(Ftarget)) stop("Set Ftarget = beta*Fmsy")
  }

  F_func = function(x) {
    ifelse(!isTRUE(HCR),multi*Fcurrent,
           ifelse(x<SNrefs[2],0, #Bban
                  ifelse(x<SNrefs[1],Ftarget*(x-SNrefs[2])/(SNrefs[1]-SNrefs[2]), #Blimit
                         Ftarget)
           )
    )
  }

  # Catch equation
  if (Pope) {
    Catch_func = function(Stock,F,M) Stock*exp(-M/2)*(1-exp(-F))
  } else {
    Catch_func = function(Stock,F,M) (F/(F+M))*Stock*(1-exp(-M-F))
  }

  sim_array = array(0, dim=c(sim_year,ncol(all_data),nsim))
  dimnames(sim_array) = list(
    Future_Year = 1:sim_year,
    Category = colnames(all_data),
    Sim_ID=1:nsim
  )
  sim_array[,"Year",] = (max(assess_data$Year)+1):(max(assess_data$Year)+sim_year)
  sim_array[,"Weight",] = weight
  sim_array[,"M",] = M
  sigma = sqrt(rec_arg$sd^2/(1-rec_arg$rho^2)) #recruitment variance including autocorrelation
  bias_corrected_mean = ifelse(bias_correct, -0.5*sigma^2,0)

  if (is.null(sim_rec_resid)) {
    set.seed(seed)
    sim_array[,"rec_resid_excl_AR",] <- rnorm(prod(dim(sim_array[,"rec_resid_excl_AR",])),bias_corrected_mean,rec_arg$sd)
  } else {
    sim_array[,"rec_resid_excl_AR",] <- sim_rec_resid
  }

  for (y in 1:sim_year) {
    if (y == 1) {
      sim_array[y,"rec_pred_by_SR",] = SRF(unlist(dplyr::select(dplyr::filter(all_data,Year==max(Year)),Spawning_number)))
      # sim_array[y,"rec_deviance_to_SR",] = sim_array[y,"rec_resid_excl_AR",] + rec_arg$rho*unlist(dplyr::select(dplyr::filter(all_data,Year==max(Year)),rec_deviance_to_SR))
      sim_array[y,"rec_pred_incl_AR",] = sim_array[y,"rec_pred_by_SR",]*exp(rec_arg$rho*unlist(dplyr::select(dplyr::filter(all_data,Year==max(Year)),rec_deviance_to_SR)))
      sim_array[y,"Stock_number",] = sim_array[y,"rec_pred_incl_AR",]*exp(sim_array[y,"rec_resid_excl_AR",])
      sim_array[y,"rec_deviance_to_SR",] = log(sim_array[y,"Stock_number",]/sim_array[y,"rec_pred_by_SR",])-bias_corrected_mean
      sim_array[y,"F",] = ifelse(sim_array[y,"Year",] %in% Fcurrent_year,Fcurrent,
                                 F_func(unlist(dplyr::select(dplyr::filter(all_data,Year==max(Year)),Spawning_number))))
      sim_array[y,"Spawning_number",] = sim_array[y,"Stock_number",]*exp(-sim_array[y,"F",]-sim_array[y,"M",])
      sim_array[y,"Catch_number",] = Catch_func(sim_array[y,"Stock_number",],sim_array[y,"F",],sim_array[y,"M",])
    } else {
      sim_array[y,"rec_pred_by_SR",] = SRF(unlist(sim_array[y-1,"Spawning_number",]))
      # sim_array[y,"rec_deviance_to_SR",] = sim_array[y,"rec_resid_excl_AR",] + rec_arg$rho*(sim_array[y-1,"rec_deviance_to_SR",])
      sim_array[y,"rec_pred_incl_AR",] = sim_array[y,"rec_pred_by_SR",]*exp(rec_arg$rho*sim_array[y-1,"rec_deviance_to_SR",])
      sim_array[y,"Stock_number",] = sim_array[y,"rec_pred_incl_AR",]*exp(sim_array[y,"rec_resid_excl_AR",])
      sim_array[y,"rec_deviance_to_SR",] = log(sim_array[y,"Stock_number",]/sim_array[y,"rec_pred_by_SR",])-bias_corrected_mean
      sim_array[y,"F",] = ifelse(sim_array[y,"Year",] %in% Fcurrent_year,Fcurrent,
                                 F_func(sim_array[y-1,"Spawning_number",]))
      sim_array[y,"Spawning_number",] = sim_array[y,"Stock_number",]*exp(-sim_array[y,"F",]-sim_array[y,"M",])
      sim_array[y,"Catch_number",] = Catch_func(sim_array[y,"Stock_number",],sim_array[y,"F",],sim_array[y,"M",])
    }
  }

  sim_array[,"Stock_biomass",] = sim_array[,"Stock_number",]*sim_array[,"Weight",]*num_to_mass_scale
  sim_array[,"Spawning_biomass",] = sim_array[,"Spawning_number",]*sim_array[,"Weight",]*num_to_mass_scale
  sim_array[,"Catch_biomass",] = sim_array[,"Catch_number",]*sim_array[,"Weight",]*num_to_mass_scale
  sim_array[,"U",] = sim_array[,"Catch_biomass",]/sim_array[,"Stock_biomass",]

  Res = list()
  Res$input = arglist
  Res$sigma_for_bias_correction = sigma
  Res$sim_array = sim_array

  all_data = mutate(all_data, Status="Past")
  mean_table = as_tibble(apply(sim_array,c(1,2),mean,na.rm=TRUE)) %>%
    mutate(Status="Future")
  Res$mean_table = bind_rows(all_data,mean_table)

  median_table = as_tibble(apply(sim_array,c(1,2),median,na.rm=TRUE)) %>%
    mutate(Status="Future")
  Res$median_table = bind_rows(all_data,median_table)

  geomean = function(x) ifelse(min(x)<=0,median(x,na.rm=TRUE),exp(mean(log(x),na.rm=TRUE)))
  geomean_table = as_tibble(apply(sim_array,c(1,2),geomean)) %>%
    mutate(Status="Future")
  Res$geomean_table = bind_rows(all_data,geomean_table)

  H10_table = as_tibble(apply(sim_array,c(1,2),quantile,prob=0.9,na.rm=TRUE)) %>%
    mutate(Status="Future")
  Res$H10_table = bind_rows(all_data,H10_table)

  L10_table = as_tibble(apply(sim_array,c(1,2),quantile,prob=0.1,na.rm=TRUE)) %>%
    mutate(Status="Future")
  Res$L10_table = bind_rows(all_data,L10_table)

  Res$catch_rank = rank(sapply(1:nsim, function(i) sim_array[sim_year,"Catch_biomass",i]))

  if (array_to_tibble) {
    sim_tibble= mutate(all_data,Sim_ID=NA)
    for (i in 1:nsim) {
      sim_tibble = rbind(sim_tibble,mutate(as_tibble(sim_array[,,i]),Sim_ID=i))
    }
    Res$sim_tibble = sim_tibble
  }
  return(Res)
}


#' Estimating MSY-based management reference points for Surume-Ika
#' @import TMB
#' @importFrom TMB MakeADFun
#' @importFrom TMB sdreport
#' @importFrom dplyr filter
#' @inheritParams make_tmb_data
#' @param future_sim_res \code{future_sim} object
#' @encoding UTF-8
#' @export
est_MSY = function(
  future_sim_res, #result of future_sim()
  obj_catch = "mean", # other options: median or geomean
  nsim = NULL,
  MSY_year = NULL,
  PGY = NULL, # PGY = c(0.6,0.1)
  PGYlower_only = FALSE,
  percentB0 = NULL, #percentB0 = c(0.2,0.3,0.4)
  Bempirical = NULL, #特定のSBをtargetにする場合
  seed = NULL,
  # trace_multi = seq(0,10,by=0.1), # for Yield-curve
  TMB = FALSE,
  nlminb_control = list(eval.max=2000,iter.max=2000)
) {

  if (!is.null(nsim)) future_sim_res$input$nsim = nsim
  if (is.null(MSY_year)) MSY_year = max(future_sim_res$mean_table$Year)
  sim_year = max(MSY_year - max(future_sim_res$input$assess_data$Year),future_sim_res$input$sim_year)
  future_sim_res$input$sim_year = sim_year
  if (!is.null(seed)) future_sim_res$input$seed = seed
  if (future_sim_res$input$HCR) {
    warning("HCR option has been off automatically")
    future_sim_res$input$HCR = FALSE
  }

  future_sim_res = do.call(future_sim, future_sim_res$input)
  future_sim_res$input$sim_rec_resid = future_sim_res$sim_array[,"rec_resid_excl_AR",]

  x_grid = seq(-2.5,2.5,by=0.25)
  for_init = t(sapply(x_grid, function(x) {
      trace_input = future_sim_res$input
      trace_input$multi = exp(x)
      trace_input$nsim = min(trace_input$nsim,1000)
      trace_input$sim_rec_resid = NULL
      for_init_res = do.call(future_sim, trace_input)
      if (obj_catch == "mean") obj_list = for_init_res$mean_table
      if (obj_catch == "median") obj_list = for_init_res$median_table
      if (obj_catch == "geomean") obj_list = for_init_res$geomean_table

      obj_list = obj_list %>% dplyr::filter(Year==MSY_year) %>%
        dplyr::select(Catch_biomass,Stock_biomass,Spawning_biomass,Stock_number,Spawning_number,F) %>%
        mutate(F_to_Fcurrent=exp(x),x=x)
      return(obj_list)
  }))
  for_init = as.data.frame(for_init)
  for_init_catch = as.numeric(for_init$Catch_biomass)
  for_init_spawner = as.numeric(for_init$Spawning_biomass)
  for_init_x_grid = as.numeric(for_init$x)

  obj_msy = function(x) {
    future_sim_x = future_sim_res
    # future_sim_x$input$sim_year = sim_year
    future_sim_x$input$multi = exp(x)
    # future_sim_x$input$sim_rec_resid = future_sim_x$sim_array[,"rec_resid_excl_AR",]
    future_sim_x = do.call(future_sim, future_sim_x$input)
    if (obj_catch == "mean") catch_obj = future_sim_x$mean_table
    if (obj_catch == "median") catch_obj = future_sim_x$median_table
    if (obj_catch == "geomean") catch_obj = future_sim_x$geomean_table
    catch_obj = catch_obj %>% dplyr::filter(Year==MSY_year) %>% dplyr::select(Catch_biomass) %>% as.numeric()
    return(-log(catch_obj))
  }

  # x_grid = seq(-5,5,by=0.5)
  # obj_msy_grid = sapply(x_grid, function(x) obj_msy(x))
  # x_init = x_grid[which.min(obj_msy_grid)]

  x_init = for_init_x_grid[which.max(for_init_catch)]
  x_upper = ifelse(x_init == max(for_init_x_grid), 10, x_init+1)
  x_lower = ifelse(x_init == min(for_init_x_grid), -10, x_init-1)

  if (TMB) {
    if (obj_catch == "median") stop("TMB cannot be used when obj_catch=median")
    tmb_data_msy = make_tmb_data(future_sim_res,
                                 ifelse(obj_catch == "mean",0,1),
                                 0,0)
    objAD = TMB::MakeADFun(tmb_data_msy, list(x=x_init), DLL="est_MSY_tmb")
    msy_optim = stats::nlminb(objAD$par, objAD$fn, gr=objAD$gr,
                  lower=list(x=x_lower), upper=list(x=x_upper),contol=nlminb_control)
    multi_msy = as.numeric(exp(msy_optim$par))
    msy = exp(-msy_optim$objective)
  } else {
    msy_optim = optim(x_init, obj_msy, lower=x_lower, upper=x_upper, method="Brent")
    multi_msy = exp(msy_optim$par)
    msy = exp(-msy_optim$value)
  }
  cat("Fmsy to Fcurrent: ", multi_msy, "\n",sep="")

  Fmsy = multi_msy*future_sim_res$input$Fcurrent
  future_msy = future_sim_res
  future_msy$input$multi = multi_msy
  future_msy = do.call(future_sim, future_msy$input)

  if (obj_catch == "mean") msy_data = dplyr::filter(future_msy$mean_table,Year==MSY_year)
  if (obj_catch == "median") msy_data = dplyr::filter(future_msy$median_table,Year==MSY_year)
  if (obj_catch == "geomean") msy_data = dplyr::filter(future_msy$geomean_table,Year==MSY_year)

  ### B0 (F=0)
  future_F0 = future_sim_res
  future_F0$input$multi = 0.0
  future_F0 = do.call(future_sim, future_F0$input)

  if (obj_catch == "mean") B0_data = dplyr::filter(future_F0$mean_table,Year==MSY_year)
  if (obj_catch == "median") B0_data = dplyr::filter(future_F0$median_table,Year==MSY_year)
  if (obj_catch == "geomean") B0_data = dplyr::filter(future_F0$geomean_table,Year==MSY_year)

  RES = list()
  RES$msy = msy
  RES$multi_msy = multi_msy
  RES$future_msy = future_msy
  RES$future_F0 = future_F0
  RES$optim$msy = msy_optim

  SBcurrent = dplyr::filter(future_sim_res$input$assess_data, Year==max(Year)) %>%
    dplyr::select(Spawning_biomass) %>% as.numeric()
  Bcurrent = dplyr::filter(future_sim_res$input$assess_data, Year==max(Year)) %>%
    dplyr::select(Stock_biomass) %>% as.numeric()

  selected_values = c("Definition","Catch_biomass","Stock_biomass","Spawning_biomass","Stock_number","Spawning_number","F","U")

  msy_data = dplyr::mutate(msy_data,Definition="MSY") %>%
    dplyr::select(selected_values) %>%
    mutate(
      F_to_Fcurrent = multi_msy,
      B_to_B0 = Stock_biomass/B0_data$Stock_biomass,
      SB_to_SB0 = Spawning_biomass/B0_data$Spawning_biomass,
      B_to_Bcurrent = Stock_biomass/Bcurrent,
      SB_to_SBcurrent = Spawning_biomass/SBcurrent)

  summary_BRP = msy_data

  B0_data = dplyr::mutate(B0_data,Definition="B0") %>%
    dplyr::select(selected_values) %>%
    mutate(
      F_to_Fcurrent = 0,
      B_to_B0 = 1,
      SB_to_SB0 = 1,
      B_to_Bcurrent = Stock_biomass/Bcurrent,
      SB_to_SBcurrent = Spawning_biomass/SBcurrent)

  if (!is.null(PGY)) { # PGY lower
    # x_grid_PGYlower = seq(5,msy_optim$par,by=-0.5)
    # %1.2 = sapply(x_grid_PGYlower, function(x) obj_msy(x))

    pgy_lower_multi = NULL
    for (i in 1:length(PGY)) {
      # x_PGYlower_init =  x_grid_PGYlower[which((%1.2+msy*PGY[i])^2==min((%1.2+msy*PGY[i])^2))]
      for_init_x_grid_PGYlower = for_init_x_grid[for_init_x_grid>msy_optim$par]
      for_init_catch_PGYlower = for_init_catch[for_init_x_grid>msy_optim$par]
      x_PGYlower_init = for_init_x_grid_PGYlower[which(abs(for_init_catch_PGYlower-msy*PGY[i]) == min(abs(for_init_catch_PGYlower-msy*PGY[i])))]
      x_PGYlower_upper = ifelse(x_PGYlower_init == max(for_init_x_grid_PGYlower),10, x_PGYlower_init+1)
      x_PGYlower_lower = max(x_PGYlower_init-1,msy_optim$par)

      if (TMB) {
        tmb_data_pgy_lower = make_tmb_data(future_sim_res,
                                     ifelse(obj_catch == "mean",0,1),
                                     1,msy*PGY[i])
        objAD <- TMB::MakeADFun(tmb_data_pgy_lower,
                           list(x=x_PGYlower_init),
                           DLL="est_MSY_tmb")
        pgy_lower_optim = stats::nlminb(objAD$par, objAD$fn, gr=objAD$gr,
                            lower=list(x=x_PGYlower_lower), upper=list(x=x_PGYlower_upper),
                            contol=nlminb_control)
        if (pgy_lower_optim$objective > 1e-6) {
          init_vec = seq(x_PGYlower_upper-0.01,x_PGYlower_lower+0.01,length=10)
          for (ii in 1:length(init_vec)) {
            objAD <- TMB::MakeADFun(tmb_data_pgy_lower,
                               list(x=init_vec[ii]),
                               DLL="est_MSY_tmb")
            pgy_lower_optim = nlminb(objAD$par, objAD$fn, gr=objAD$gr,
                                     lower=list(x=x_PGYlower_lower), upper=list(x=x_PGYlower_upper),
                                     contol=nlminb_control)
            if (pgy_lower_optim$objective < 1e-6) break
          }
          if (pgy_lower_optim$objective > 1e-6) {
            stop("Estimation error in PGYlower: Try 'TMB=FALSE'")
          }
        }
      } else {
        pgy_lower_optim = optim(x_PGYlower_init, function(x) (obj_msy(x)+log(msy*PGY[i]))^2,
                                lower=x_PGYlower_lower,upper=x_PGYlower_upper,method="Brent")
      }
      pgy_lower_multi = c(pgy_lower_multi,as.numeric(exp(pgy_lower_optim$par)))
      RES$optim$pgy_lower[[i]] = pgy_lower_optim

      cat("F",PGY[i],"pgy_lower to Fcurrent: ", exp(pgy_lower_optim$par), "\n",sep="")
    }

    future_pgy_lower = lapply(1:length(PGY), function(i) {
      pgy_lower_input = future_sim_res$input
      pgy_lower_input$multi = pgy_lower_multi[i]
      do.call(future_sim, pgy_lower_input)
    })

    RES$future_pgy_lower = future_pgy_lower

    for (i in 1:length(PGY)) {
      if (obj_catch == "mean") summary_table = future_pgy_lower[[i]]$mean_table
      if (obj_catch == "median") summary_table = future_pgy_lower[[i]]$median_table
      if (obj_catch == "geomean") summary_table = future_pgy_lower[[i]]$geomean_table

      summary_table = dplyr::filter(summary_table,Year==MSY_year) %>%
        mutate(Definition=sprintf("%1.2fPGYlower",PGY[i])) %>%
        select(selected_values) %>%
        mutate(
          F_to_Fcurrent = pgy_lower_multi[i],
          B_to_B0 = Stock_biomass/B0_data$Stock_biomass,
          SB_to_SB0 = Spawning_biomass/B0_data$Spawning_biomass,
          B_to_Bcurrent = Stock_biomass/Bcurrent,
          SB_to_SBcurrent = Spawning_biomass/SBcurrent)
      summary_BRP = bind_rows(summary_BRP,summary_table)
    }

    if (!isTRUE(PGYlower_only)) { #PGY upper (if neccessary)
      # x_grid_PGYupper = seq(-5,msy_optim$par,by=0.5)
      # obj_PGYupper_grid = sapply(x_grid_PGYupper, function(x) obj_msy(x))

      pgy_upper_multi = NULL
      for (i in 1:length(PGY)) {
        # x_PGYupper_init =  x_grid_PGYupper[which((obj_PGYupper_grid+msy*PGY[i])^2==min((obj_PGYupper_grid+msy*PGY[i])^2))]
        for_init_x_grid_PGYupper = for_init_x_grid[for_init_x_grid<msy_optim$par]
        for_init_catch_PGYupper = for_init_catch[for_init_x_grid<msy_optim$par]
        x_PGYupper_init = for_init_x_grid_PGYupper[which(abs(for_init_catch_PGYupper-msy*PGY[i]) == min(abs(for_init_catch_PGYupper-msy*PGY[i])))]
        x_PGYupper_upper = min(x_PGYupper_init+1,as.numeric(msy_optim$par))
        x_PGYupper_lower = ifelse(x_PGYupper_init == min(for_init_x_grid_PGYupper),-10, x_PGYupper_init-1)

        if (TMB) {
          tmb_data_pgy_upper = make_tmb_data(future_sim_res,
                                             ifelse(obj_catch == "mean",0,1),
                                             1,msy*PGY[i])
          objAD = TMB::MakeADFun(tmb_data_pgy_upper, list(x=x_PGYupper_init), DLL="est_MSY_tmb")
          pgy_upper_optim = stats::nlminb(objAD$par, objAD$fn, gr=objAD$gr,
                                    lower=list(x=x_PGYupper_lower), upper=list(x=x_PGYupper_upper),
                                   contol=nlminb_control)
          if (pgy_upper_optim$objective > 1e-6) {
            init_vec = seq(x_PGYupper_lower+0.01,x_PGYupper_upper-0.01,length=10)
            for (ii in 1:length(init_vec)) {
              objAD <- TMB::MakeADFun(tmb_data_pgy_upper,list(x=init_vec[ii]),
                                 DLL="est_MSY_tmb")
              pgy_upper_optim = stats::nlminb(objAD$par, objAD$fn, gr=objAD$gr,
                                       lower=list(x=x_PGYupper_lower), upper=list(x=x_PGYupper_upper),
                                       contol=nlminb_control)
              if (pgy_upper_optim$objective < 1e-6) break
            }
            if (pgy_upper_optim$objective > 1e-6) {
              stop("Estimation error in PGYupper: Try 'TMB=FALSE'")
            }
          }
        } else {
          pgy_upper_optim = optim(x_PGYupper_init, function(x) (obj_msy(x)+log(msy*PGY[i]))^2,
                                  lower=x_PGYupper_lower,upper=x_PGYupper_upper,method="Brent")
        }

        pgy_upper_multi = c(pgy_upper_multi,as.numeric(exp(pgy_upper_optim$par)))
        RES$optim$pgy_upper[[i]] = pgy_upper_optim

        cat("F",PGY[i],"pgy_upper to Fcurrent: ", exp(pgy_upper_optim$par), "\n", sep="")
      }

      future_pgy_upper = lapply(1:length(PGY), function(i) {
        pgy_upper_input = future_sim_res$input
        pgy_upper_input$multi = pgy_upper_multi[i]
        do.call(future_sim, pgy_upper_input)
      })

      RES$future_pgy_upper = future_pgy_upper

      for (i in 1:length(PGY)) {
        if (obj_catch == "mean") summary_table = future_pgy_upper[[i]]$mean_table
        if (obj_catch == "median") summary_table = future_pgy_upper[[i]]$median_table
        if (obj_catch == "geomean") summary_table = future_pgy_upper[[i]]$geomean_table

        summary_table = dplyr::filter(summary_table,Year==MSY_year) %>%
          mutate(Definition=sprintf("%1.2fPGYupper",PGY[i])) %>%
          select(selected_values) %>%
          mutate(
            F_to_Fcurrent = pgy_upper_multi[i],
            B_to_B0 = Stock_biomass/B0_data$Stock_biomass,
            SB_to_SB0 = Spawning_biomass/B0_data$Spawning_biomass,
            B_to_Bcurrent = Stock_biomass/Bcurrent,
            SB_to_SBcurrent = Spawning_biomass/SBcurrent)
        summary_BRP = bind_rows(summary_BRP,summary_table)
      }
    }
  }

  if (!is.null(percentB0)) {
    obj_pB0 = function(x) {
      future_sim_x = future_sim_res
      # future_sim_x$input$sim_year = sim_year
      future_sim_x$input$multi = exp(x)
      # future_sim_x$input$sim_rec_resid = future_sim_x$sim_array[,"rec_resid_excl_AR",]
      future_sim_x = do.call(future_sim, future_sim_x$input)
      if (obj_catch == "mean") SB_obj = future_sim_x$mean_table
      if (obj_catch == "median") SB_obj = future_sim_x$median_table
      if (obj_catch == "geomean") SB_obj = future_sim_x$geomean_table
      SB_obj = SB_obj %>% dplyr::filter(Year==MSY_year) %>% dplyr::select(Spawning_biomass) %>% as.numeric()
      return(SB_obj)
    }

    # x_grid = seq(-5,5,by=0.5)
    # obj_pB0_grid = sapply(x_grid, function(x) obj_pB0(x))

    pB0_multi = NULL
    for (i in 1:length(percentB0)) {
      # x_B0_init = x_grid[which(abs(obj_pB0_grid-percentB0[i]*B0_data$Spawning_biomass) == min(abs(obj_pB0_grid-percentB0[i]*B0_data$Spawning_biomass)))]

      x_B0_init = for_init_x_grid[which(abs(for_init_spawner-percentB0[i]*B0_data$Spawning_biomass) == min(abs(for_init_spawner-percentB0[i]*B0_data$Spawning_biomass)))]
      x_B0_upper = ifelse(x_B0_init == max(for_init_x_grid),10, x_B0_init+1)
      x_B0_lower = ifelse(x_B0_init == min(for_init_x_grid),-10, x_B0_init-1)

      if (TMB) {
        tmb_data_pB0 = make_tmb_data(future_sim_res,
                                     ifelse(obj_catch == "mean",0,1),
                                     2,percentB0[i]*B0_data$Spawning_biomass)
        objAD <- TMB::MakeADFun(tmb_data_pB0, list(x=x_B0_init), DLL="est_MSY_tmb")
        B0_optim <- stats::nlminb(objAD$par, objAD$fn, gr=objAD$gr,
                                  lower=list(x=x_B0_lower), upper=list(x=x_B0_upper),contol=nlminb_control)
        if (B0_optim$objective > 1e-6) {
          stop("Estimation error in pB0: Try 'TMB=FALSE'")
        }
      } else {
        B0_optim = optim(x_B0_init, function(x) (obj_pB0(x)-percentB0[i]*B0_data$Spawning_biomass)^2,
                         lower=x_B0_lower,upper=x_B0_upper,method="Brent")
      }

      pB0_multi = c(pB0_multi,as.numeric(exp(B0_optim$par)))
      RES$optim$pB0[[i]] = B0_optim
      cat("F",percentB0[i],"B0 to Fcurrent: ", exp(B0_optim$par), "\n", sep="")
    }

    future_pB0 = lapply(1:length(percentB0), function(i) {
      pB0_input = future_sim_res$input
      pB0_input$multi = pB0_multi[i]
      do.call(future_sim, pB0_input)
    })

    RES$future_pB0 = future_pB0

    for (i in 1:length(percentB0)) {
      if (obj_catch == "mean") summary_table = future_pB0[[i]]$mean_table
      if (obj_catch == "median") summary_table = future_pB0[[i]]$median_table
      if (obj_catch == "geomean") summary_table = future_pB0[[i]]$geomean_table

      summary_table = dplyr::filter(summary_table,Year==MSY_year) %>%
        mutate(Definition=sprintf("%1.2fB0",percentB0[i])) %>%
        select(selected_values) %>%
        mutate(
          F_to_Fcurrent = pB0_multi[i],
          B_to_B0 = Stock_biomass/B0_data$Stock_biomass,
          SB_to_SB0 = Spawning_biomass/B0_data$Spawning_biomass,
          B_to_Bcurrent = Stock_biomass/Bcurrent,
          SB_to_SBcurrent = Spawning_biomass/SBcurrent)
      summary_BRP = bind_rows(summary_BRP,summary_table)
    }
  }

  if (!is.null(Bempirical)) {
    if (is.null(percentB0)) {
      obj_pB0 = function(x) {
        future_sim_x = future_sim_res
        # future_sim_x$input$sim_year = sim_year
        future_sim_x$input$multi = exp(x)
        # future_sim_x$input$sim_rec_resid = future_sim_x$sim_array[,"rec_resid_excl_AR",]
        future_sim_x = do.call(future_sim, future_sim_x$input)
        if (obj_catch == "mean") SB_obj = future_sim_x$mean_table
        if (obj_catch == "median") SB_obj = future_sim_x$median_table
        if (obj_catch == "geomean") SB_obj = future_sim_x$geomean_table
        SB_obj = SB_obj %>% dplyr::filter(Year==MSY_year) %>% dplyr::select(Spawning_biomass) %>% as.numeric()
        return(SB_obj)
      }

      # x_grid = seq(-5,5,by=0.5)
      # obj_pB0_grid = sapply(x_grid, function(x) obj_pB0(x))
    }

    Bempirical_multi = NULL
    for (i in 1:length(Bempirical)) {
      if (Bempirical[i]>as.numeric(B0_data$Spawning_biomass)) {
        warning("Bempirical is higher than B0: Should be reset")
      }
      # x_B0_init = x_grid[which(abs(obj_pB0_grid-Bempirical[i]) == min(abs(obj_pB0_grid-Bempirical[i])))]
      x_B0_init = for_init_x_grid[which(abs(for_init_spawner-Bempirical[i]) == min(abs(for_init_spawner-Bempirical[i])))]
      x_B0_upper = ifelse(x_B0_init == max(for_init_x_grid),10, x_B0_init+1)
      x_B0_lower = ifelse(x_B0_init == min(for_init_x_grid),-10, x_B0_init-1)
      if (TMB) {
        tmb_data_Bempirical = make_tmb_data(future_sim_res,
                                     ifelse(obj_catch == "mean",0,1),
                                     2,Bempirical[i])
        objAD = TMB::MakeADFun(tmb_data_Bempirical, list(x=x_B0_init), DLL="est_MSY_tmb")
        B0_optim = stats::nlminb(objAD$par, objAD$fn, gr=objAD$gr,
                           lower=list(x=x_B0_lower), upper=list(x=x_B0_upper),contol=nlminb_control)
        if (B0_optim$objective > 1e-6) {
          stop("Estimation error in Bempirical: Try 'TMB=FALSE'")
        }

      } else {
        B0_optim = optim(x_B0_init, function(x) (obj_pB0(x)-Bempirical[i])^2,
                         lower=x_B0_lower,upper=x_B0_upper,method="Brent")
      }
      Bempirical_multi = c(Bempirical_multi,as.numeric(exp(B0_optim$par)))
      RES$optim$Bempirical[[i]] = B0_optim
      cat("F_Bempirical_",Bempirical[i]," to Fcurrent: ", exp(B0_optim$par), "\n", sep="")
    }

    future_Bempirical = lapply(1:length(Bempirical), function(i) {
      Bempirical_input = future_sim_res$input
      Bempirical_input$multi = Bempirical_multi[i]
      do.call(future_sim, Bempirical_input)
    })

    RES$future_Bempirical = future_Bempirical

    for (i in 1:length(Bempirical)) {
      if (obj_catch == "mean") summary_table = future_Bempirical[[i]]$mean_table
      if (obj_catch == "median") summary_table = future_Bempirical[[i]]$median_table
      if (obj_catch == "geomean") summary_table = future_Bempirical[[i]]$geomean_table

      summary_table = dplyr::filter(summary_table,Year==MSY_year) %>%
        mutate(Definition=sprintf("Bempirical%4.0f",Bempirical[i])) %>%
        select(selected_values) %>%
        mutate(
          F_to_Fcurrent = Bempirical_multi[i],
          B_to_B0 = Stock_biomass/B0_data$Stock_biomass,
          SB_to_SB0 = Spawning_biomass/B0_data$Spawning_biomass,
          B_to_Bcurrent = Stock_biomass/Bcurrent,
          SB_to_SBcurrent = Spawning_biomass/SBcurrent)
      summary_BRP = bind_rows(summary_BRP,summary_table)
    }
  }

  summary_BRP = bind_rows(summary_BRP,B0_data)
  RES$summary_BRP = summary_BRP

  # if (!is.null(trace_multi)) {
  #   future_trace = lapply(1:length(trace_multi), function(i) {
  #     if (trace_multi[i]==0) { future_F0 } else {
  #       trace_input = future_sim_res$input
  #       trace_input$multi = trace_multi[i]
  #       do.call(future_sim, trace_input)
  #     }
  #   })
  #
  #   RES$future_trace = future_trace
  #
  #   trace_table = tibble()
  #   for (i in 1:length(trace_multi)) {
  #     if (obj_catch == "mean") summary_table = future_trace[[i]]$mean_table
  #     if (obj_catch == "median") summary_table = future_trace[[i]]$median_table
  #     if (obj_catch == "geomean") summary_table = future_trace[[i]]$geomean_table
  #
  #     summary_table = dplyr::filter(summary_table,Year==MSY_year) %>%
  #       mutate(Definition=sprintf("trace_multi%2.2f",trace_multi[i])) %>%
  #       select(selected_values) %>%
  #       mutate(
  #         F_to_Fcurrent = trace_multi[i],
  #         B_to_B0 = Stock_biomass/B0_data$Stock_biomass,
  #         SB_to_SB0 = Spawning_biomass/B0_data$Spawning_biomass,
  #         B_to_Bcurrent = Stock_biomass/Bcurrent,
  #         SB_to_SBcurrent = Spawning_biomass/SBcurrent)
  #     trace_table = bind_rows(trace_table,summary_table)
  #   }
  #   RES$trace_table = trace_table
  # }
  return(RES)
}

#' Function for making input data for est_MSY with TMB
#' @import dplyr
#' @importFrom dplyr filter
#' @param future_sim_res \code{future_sim} object
#' @param obj_catch objective of maximum catch: mean(0) or geomean(1)
#' @param objective # 0: MSY, 1: PGY, 2: percentB0
#' @param objective_value # for PGY and percentB0
#' @encoding UTF-8
#' @export
make_tmb_data = function(
  future_sim_res, #result of future_sim()
  obj_catch, # 0: mean, 1: geomean
  objective, # 0: MSY, 1: PGY, 2: percentB0,
  objective_value # for PGY and percentB0
) {
  spawner_init = dplyr::filter(future_sim_res$input$assess_data,Year==max(Year)) %>%
    dplyr::select(Spawning_number) %>% as.numeric()
  deviance_init = dplyr::filter(future_sim_res$mean_table,Year==max(future_sim_res$input$assess_data$Year)) %>%
    dplyr::select(rec_deviance_to_SR) %>% as.numeric()
  bias_corrected_mean = ifelse(future_sim_res$input$bias_correct, -0.5*future_sim_res$input$rec_arg$sd^2/(1-future_sim_res$input$rec_arg$rho^2),0)

  tmb_data = list(spawner_init=spawner_init,deviance_init=deviance_init,
                  SR = ifelse(future_sim_res$input$SR == "HS", 0, ifelse(future_sim_res$input$SR == "BH",1,2)),
                  rec_par_a = future_sim_res$input$rec_arg$a,
                  rec_par_b = future_sim_res$input$rec_arg$b,
                  # rec_par_sd = future_sim_res$input$rec_arg$sd,
                  rec_par_rho = future_sim_res$input$rec_arg$rho,
                  bias_corrected_mean = bias_corrected_mean,
                  rec_resid_mat = future_sim_res$sim_array[,"rec_resid_excl_AR",],
                  weight_mat = future_sim_res$sim_array[,"Weight",],
                  M_mat = future_sim_res$sim_array[,"M",],
                  Pope = as.numeric(future_sim_res$input$Pope),
                  sim_year = future_sim_res$input$sim_year,
                  nsim = future_sim_res$input$nsim,
                  obj_catch = obj_catch,
                  objective = objective,
                  objective_value = objective_value,
                  Fcurrent = future_sim_res$input$Fcurrent,
                  Fcurrent_year = ifelse(is.null(future_sim_res$input$Fcurrent_year), -1,
                                         max(future_sim_res$input$Fcurrent_year)-max(future_sim_res$input$assess_data$Year)),
                  num_to_mass_scale = future_sim_res$input$num_to_mass_scale
  )
  return(tmb_data)
}

#' Tracing for future_sim results with different F values
#' @import dplyr
#' @importFrom dplyr tibble
#' @importFrom dplyr filter
#' @importFrom dplyr bind_rows
#' @param future_sim_res \code{future_sim} object
#' @param trace_multi multiplier for Fcurrent to be traced
#' @param obj_catch objective value("mean", "median, or "geomean")
#' @param sim_year simulation years
#' @param nsim number of simulations
#' @param seed seed in \code{set.seed()}
#' @encoding UTF-8
#' @export
trace_future = function(
  future_sim_res,
  # trace_multi = c(0,exp(seq(-3,3,by=0.1))),
  trace_multi = seq(0,10,by=0.1),
  obj_catch = "mean",
  sim_year = NULL,
  nsim = NULL,
  seed = NULL
) {
  if (!is.null(nsim)) future_sim_res$input$nsim = nsim
  if (!is.null(sim_year)) future_sim_res$input$sim_year = sim_year
  if (!is.null(seed)) future_sim_res$input$seed = seed
  if (future_sim_res$input$HCR) {
    warning("HCR option has been off automatically")
    future_sim_res$input$HCR = FALSE
  }

  future_sim_res = do.call(future_sim, future_sim_res$input)

  trace_input = future_sim_res$input
  trace_input$sim_rec_resid = future_sim_res$sim_array[,"rec_resid_excl_AR",]

  future_trace = lapply(trace_multi, function(x) {
    trace_input$multi = x
    cat("simulating multiplier = ",x,"\n",sep="")
    do.call(future_sim, trace_input)
  })

  RES = list()
  RES$future_trace = future_trace

  selected_values = c("Definition","Catch_biomass","Stock_biomass","Spawning_biomass","Stock_number","Spawning_number","F","U")

  trace_table = tibble()
  for (i in 1:length(trace_multi)) {
    if (obj_catch == "mean") summary_table = future_trace[[i]]$mean_table
    if (obj_catch == "median") summary_table = future_trace[[i]]$median_table
    if (obj_catch == "geomean") summary_table = future_trace[[i]]$geomean_table

    summary_table = dplyr::filter(summary_table,Year==max(Year)) %>%
      mutate(Definition=sprintf("trace_multi%2.2f",trace_multi[i])) %>%
      select(selected_values)
    # %>%
    #   mutate(
    #     F_to_Fcurrent = trace_multi[i],
    #     B_to_B0 = Stock_biomass/B0_data$Stock_biomass,
    #     SB_to_SB0 = Spawning_biomass/B0_data$Spawning_biomass,
    #     B_to_Bcurrent = Stock_biomass/Bcurrent,
    #     SB_to_SBcurrent = Spawning_biomass/SBcurrent)
    trace_table = bind_rows(trace_table,summary_table)
  }

  RES$trace_table = trace_table
  # (RES$plot = g1)

  return (RES)
}

#' Function for visualizing yield curve from trace_future
#' @import ggplot2
#' @param trace_future_res \code{trace_future} object
#' @param title figure title
#' @encoding UTF-8
#' @export
visualize_yield_curve = function(trace_future_res,
                                 title = NULL) {
  alpha = 0.3
  trace_table = trace_future_res$trace_table
  g1 = ggplot(trace_table,aes(x=Spawning_biomass, y=Catch_biomass)) +
    geom_ribbon(aes(ymin=0,ymax=Catch_biomass),
                fill = "blue", alpha=0.4, colour = "blue", size = 1)+
      ylim(0,NA) +
      theme_bw(base_size=14)
  if (!is.null(title))g1 = g1 + ggtitle(title)
  (g1)
}

#' #' Theme for ggplot
#' #' @import ggplot2
#' #' @encoding UTF-8
#' #' @export
#' theme_SH <- function(){
#'   theme_bw(base_size=12) +
#'     theme(panel.grid = element_blank(),
#'           axis.text.x=element_text(size=11,color="black"),
#'           axis.text.y=element_text(size=11,color="black"),
#'           axis.line.x=element_line(size= 0.3528),
#'           axis.line.y=element_line(size= 0.3528),
#'           legend.position="none")
#' }

#' Function for plotting future_res
#' @import ggplot2
#' @param future_list list of \code{future_sim}
#' @param title figure title
#' @encoding UTF-8
#' @export
visualize_future = function(
  future_list,
  scenario_name = NULL,
  what_center = "mean", # or median or geomean
  CI_range = 0.8, # NULL if not neccessary
  title = NULL,
  # example_number = 0, # number of described simulations
  # seed = 12345, # seed for extracting examples
  BRP = NULL # under consideration...
) {

  if (is.null(scenario_name)) {
    for (i in 1:length(future_list)) {
      scenario_name = c(scenario_name, sprintf("scenario_%s",i))
    }
  }

  plot_data_all = tibble()
  for (i in 1:length(future_list)) {

    if (what_center == "mean") center_table = future_list[[i]]$mean_table
    if (what_center == "median") center_table = future_list[[i]]$median_table
    if (what_center == "geomean") center_table = future_list[[i]]$geomean_table

    center_table = center_table %>%
      mutate(Scenario = scenario_name[i],Range = "center")

    plot_data_all = bind_rows(plot_data_all,center_table)

    if (!is.null(CI_range)) {
      lower_data = as_tibble(apply(future_list[[i]]$sim_array,c(1,2),quantile,probs = (1-CI_range)/2)) %>%
        mutate(Status="Future",Scenario=scenario_name[i],Range="lower")
      plot_data_all = bind_rows(plot_data_all, lower_data)
      upper_data = as_tibble(apply(future_list[[i]]$sim_array,c(1,2),quantile,probs = 1-(1-CI_range)/2)) %>%
        mutate(Status="Future",Scenario=scenario_name[i],Range="upper")
      plot_data_all = bind_rows(plot_data_all, upper_data)
    }
  }

  plot_data_all = dplyr::select(plot_data_all,
                                Year,Stock_biomass,Spawning_biomass,Catch_biomass,F,
                                Status,Scenario,Range) %>%
    tidyr::gather(key=Y, value = value, Stock_biomass,Spawning_biomass,Catch_biomass,F)
  # dplyr::distinct()

  plot_data_center = mutate(plot_data_all,
                            Y_f = factor(plot_data_all$Y,
                                         levels=c("Stock_biomass","Spawning_biomass","Catch_biomass","F"))) %>%
    filter(Range == "center")

  plot_data_CI =  mutate(plot_data_all,
                         Y_f = factor(plot_data_all$Y,
                                      levels=c("Stock_biomass","Spawning_biomass","Catch_biomass","F"))) %>%
    filter(Range != "center") %>%
    spread(key = Range, value=value)

  alpha = 0.4
  g1 = ggplot(plot_data_center) +
    geom_line(data = dplyr::filter(plot_data_center),
              aes(x=Year,y=value,group=Scenario, colour=Scenario),size=1)+
    geom_line(data = dplyr::filter(plot_data_center,Status=="Past"&Scenario==scenario_name[1]),
              aes(x=Year,y=value),size=1, colour="black")+
    geom_ribbon(data = plot_data_CI,
                aes(x=Year,ymin=lower,ymax=upper,
                    group=Scenario, fill=Scenario),alpha=alpha) +
    facet_wrap(~Y_f, scales="free_y") +
    ylim(0,NA) +
    theme_bw(base_size=14)+
    theme(axis.title.y = element_blank())

  if (!is.null(title)) g1 = g1 + ggtitle(title)
  return(g1)
}
