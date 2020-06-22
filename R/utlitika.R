
### Ulility Functions in messir ### ----

#' output a table of samuika results
#' @param model \code{samuika} object
#' @encoding UTF-8
#' @export
out_summary_estimate <- function(model,CI=0.8) {
  samuika_model = model
  Summary_PopDyn = samuika_model$Summary_PopDyn %>%
    dplyr::select(Stock_ID, Year, Stock_biomass, Spawning_biomass, F, Catch_est) %>%
    gather(key=category, value=value, -Stock_ID, -Year) %>%
    mutate(category_f = factor(category,level=c("Stock_biomass","Spawning_biomass","F","Catch_est")))

  B_se = t(samuika_model$B_se) %>% as_tibble() %>%
    mutate(Year=as.numeric(colnames(samuika_model$B_se))) %>%
    gather(key=Stock_ID, value=SE, -Year) %>%
    mutate(Stock_ID=as.numeric(Stock_ID),
           category="Stock_biomass",category_f=factor("Stock_biomass"))

  SSB_se = t(samuika_model$SSB_se) %>% as_tibble() %>%
    mutate(Year=as.numeric(colnames(samuika_model$SSB_se))) %>%
    gather(key=Stock_ID, value=SE, -Year) %>%
    mutate(Stock_ID=as.numeric(Stock_ID),
           category="Spawning_biomass",category_f=factor("Spawning_biomass"))

  F_se = t(samuika_model$F_se) %>% as_tibble() %>%
    mutate(Year=as.numeric(colnames(samuika_model$F_se))) %>%
    gather(key=Stock_ID, value=SE, -Year) %>%
    mutate(Stock_ID=as.numeric(Stock_ID),
           category="F",category_f=factor("F"))

  Catch_se = t(samuika_model$C_se) %>% as_tibble() %>%
    mutate(Year=as.numeric(colnames(samuika_model$C_se))) %>%
    gather(key=Stock_ID, value=SE, -Year) %>%
    mutate(Stock_ID=as.numeric(Stock_ID),
           category="Catch_est",category_f=factor("Catch_est"))

  SE_data = bind_rows(B_se, SSB_se, F_se, Catch_se)

  Summary_PopDyn = left_join(Summary_PopDyn,SE_data, by=c("Stock_ID","Year","category","category_f"))

  Summary_PopDyn = Summary_PopDyn %>%
    mutate(CV = SE/value) %>%
    mutate(Cz = exp(qnorm(CI+(1-CI)/2)*sqrt(log(1+CV^2)))) %>%
    mutate(lower = value/Cz, upper = value*Cz) %>%
    mutate(type = "Estimate") %>%
    mutate(
      category2 = case_when(category == "Catch_est" ~ "Catch",
                            category == "Biomass" ~ "Stock_biomass",
                            category == "Spawning_Biomass" ~ "Spawning_biomass",
                            TRUE ~ category),
      category_f2 = case_when(category_f == "Catch_est" ~ "Catch",
                              category == "Biomass" ~ "Stock_biomass",
                              category == "Spawning_Biomass" ~ "Spawning_biomass",
                              TRUE ~ category_f)) %>%
    dplyr::select(-category,-category_f)

  Summary_PopDyn = Summary_PopDyn %>%
    rename(category=category2, category_f=category_f2)

  return(Summary_PopDyn)
}


#' output a table of samuika results (Stock number and Spawning_number)
#' @param model \code{samuika} object
#' @encoding UTF-8
#' @export
out_summary_estimate_n <- function(model,CI=0.8) {
  samuika_model = model
  Summary_PopDyn = samuika_model$Summary_PopDyn %>%
    dplyr::select(Stock_ID, Year, Stock_number, Spawning_number) %>%
    gather(key=category, value=value, -Stock_ID, -Year) %>%
    mutate(category_f = factor(category,level=c("Stock_number","Spawning_number")))

  N_se = t(samuika_model$N_se) %>% as_tibble() %>%
    mutate(Year=as.numeric(colnames(samuika_model$N_se))) %>%
    gather(key=Stock_ID, value=SE, -Year) %>%
    mutate(Stock_ID=as.numeric(Stock_ID),
           category="Stock_number",category_f=factor("Stock_number"))

  SSN_se = t(samuika_model$SSN_se) %>% as_tibble() %>%
    mutate(Year=as.numeric(colnames(samuika_model$SSN_se))) %>%
    gather(key=Stock_ID, value=SE, -Year) %>%
    mutate(Stock_ID=as.numeric(Stock_ID),
           category="Spawning_number",category_f=factor("Spawning_number"))

  SE_data = bind_rows(N_se, SSN_se)

  Summary_PopDyn = left_join(Summary_PopDyn,SE_data, by=c("Stock_ID","Year","category","category_f"))

  Summary_PopDyn = Summary_PopDyn %>%
    mutate(CV = SE/value) %>%
    mutate(Cz = exp(qnorm(CI+(1-CI)/2)*sqrt(log(1+CV^2)))) %>%
    mutate(lower = value/Cz, upper = value*Cz) %>%
    mutate(type = "Estimate")

  return(Summary_PopDyn)
}


#' Making a combined list with original list name(s)
#' @param ... list of \code{class("list")}
#' @encoding UTF-8
#' @export
make_named_list <- function(...){
  object_name <- as.character(eval(substitute(alist(...))))
  x <- list(...)
  vnames <- names(x)
  novnames <- !nzchar(vnames)
  if (length(novnames) > 0L){
    vnames[novnames] <- object_name[novnames]
  } else {
    vnames <- object_name
  }
  setNames(x, vnames)
}

#' Stock-recruitment function
#'
#' @encoding UTF-8
#' @export
SRF = function(a,b,x,SR="HS") {
  sapply(x, function(y) {
    if (SR=="HS") R = ifelse(y>b,b*a,y*a)
    if (SR=="BH") R = a*y/(1+b*y)
    if (SR=="RI") R = a*y*exp(-b*y)
    R
  })
}

#' Function to get stock-recruit data from samuika object
#' @import dplyr
#' @param model samuika object
#' @encoding UTF-8
#' @export
get_SRdata = function(model) {
  REC = t(model$N_est) %>% as_tibble() %>% mutate(Year = colnames(model$N_est)) %>%
    mutate(Year = as.numeric(as.character(Year))) %>%
    gather(key=Stock_ID, value = Stock_number, -Year) %>%
    filter(Year > min(Year))

  SSN = t(model$SSN_est) %>% as_tibble() %>% mutate(Year = colnames(model$SSN_est)) %>%
    mutate(Year = as.numeric(as.character(Year))) %>%
    gather(key=Stock_ID, value = Spawning_number, -Year) %>%
    filter(Year < max(Year)) %>%
    mutate(Year = Year +1)

  SRdata = full_join(REC, SSN) %>%
    dplyr::mutate(Stock_ID2 = as.numeric(Stock_ID)) %>%
    dplyr::select(-Stock_ID) %>%
    dplyr::rename(Stock_ID = Stock_ID2) %>%
    dplyr::select(Stock_ID, everything())

  SRdata2 = left_join(SRdata,
    select(model$Summary_PopDyn,Stock_ID,Year,Regime))

  SRdata2
}


#' plotting Kobe chart
#' @import ggplot2
#' @param assess_data assessment data
#' @param SBrefs values of SBtarget, SBlimit, and SBban
#' @param Fmsy value of Fmsy
#' @param labeled_year years to be labeled
#' @param title figure title
#' @encoding UTF-8
#' @export
visualize_kobe_chart = function(
  assess_data,
  SBrefs = c("SBtarget","SBlimit","SBban"),
  Fmsy,
  # HCR = TRUE,
  labeled_year = NULL,
  title = NULL
) {

  assess_data = dplyr::mutate(assess_data,
                              Fratio = F/Fmsy,
                              SBratio = Spawning_biomass/SBrefs[1])

  if (is.null(labeled_year)) {
    labeled_year = c(assess_data$Year[assess_data$Year%%5==0],max(assess_data$Year))
  }

  assess_data = dplyr::mutate(assess_data,
                              Fratio = F/Fmsy,
                              SBratio = Spawning_biomass/SBrefs[1]) %>%
    mutate(label=if_else(Year %in% labeled_year, as.character(Year), ""))


  xmax = max(assess_data$SBratio*1.1,2)
  ymax = max(assess_data$Fratio*1.2,2)

  g1 = ggplot(assess_data, aes(x=SBratio,y=Fratio)) +
    geom_polygon(data = tibble(SBratio=c(0,0,xmax,xmax), Fratio=c(0,ymax,ymax,0)),
                 fill="khaki1")+
    geom_polygon(data = tibble(SBratio=c(1,1,xmax,xmax), Fratio=c(0,1,1,0)),
                 fill="olivedrab2")+
    geom_polygon(data = tibble(SBratio=c(0,0,1,1), Fratio=c(1,ymax,ymax,1)),
                 fill="indianred1")+
    geom_vline(xintercept=c(1,SBrefs[2]/SBrefs[1],SBrefs[3]/SBrefs[1]),
               colour=c("#00533E","#edb918","#C73C2E"),size=1)+
    annotate("text",label="SBmsy",x=1+0.1,y=ymax*0.95,size=5)+
    annotate("text",label="SBlimit",x=SBrefs[2]/SBrefs[1]+0.1,y=ymax*0.95,size=5)+
    annotate("text",label="SBban",x=SBrefs[3]/SBrefs[1]+0.1,y=ymax*0.95,size=5)+
    annotate("text",label="Fmsy",x=xmax*0.95,y=1.1,size=5)+
    geom_path(size=1.2) +
    geom_point(size=3, colour="white") +
    ggrepel::geom_text_repel(aes(label=label),
                             size=6,box.padding=0.5,segment.color="gray") +
    theme(legend.position="none") +
    theme_bw(base_size=14)+
    coord_cartesian(xlim=c(0,xmax),ylim=c(0,ymax),expand=0) +
    ylab("F/Fmsy") + xlab("SB/SBmsy")

  if (!is.null(title)) g1 + ggtitle(title)

  (g1)
}


#' plotting (multiple) samuika results
#' @import ggplot2
#' @param list_samuika_res list of samuika objects
#' @encoding UTF-8
#' @export
plot_samuika_estimates = function(list_samuika_res,
                                  model_name = NULL,
                                  Stock_ID=0,
                                  CI = 0.8,
                                  CI_plot = TRUE,
                                  legend_position = "bottom",
                                  legend_text_size = 1.2,
                                  legend_key_size = 1
) {
  message(paste0("plotting Stock_ID=",Stock_ID," OK?"))
  NRES = length(list_samuika_res)
  name_list = NULL

  for (i in 1:NRES) {
    res = list_samuika_res[[i]]
    if (!is.null(model_name)) {
      name = model_name[i]
    } else {
      name = names(list_samuika_res)[i]
    }
    if (is.null(name)) {
      name = paste0("model",i)
    }
    name_list = c(name_list,name)
    Stock_ID_dbl = Stock_ID
    out_samuika = out_summary_estimate(res,CI=CI) %>%
      mutate(model0 = name) %>%
      dplyr::filter(Stock_ID == Stock_ID_dbl)
    if (i==1) {
      Summary_all = out_samuika
    } else {
      Summary_all = full_join(Summary_all,out_samuika)
    }
  }
  Summary_all = Summary_all %>%
    mutate(category_f2 = factor(category_f, level = c("Stock_biomass","Spawning_biomass","F","Catch"))) %>%
    mutate(model = factor(model0, level = name_list))

  g1 = ggplot(Summary_all, aes(x=Year))
  if (CI_plot) {
    g1 = g1+geom_ribbon(aes(ymin=lower,ymax=upper,group=model,fill=model),alpha=0.3)
  }
  g1 = g1 +
    geom_path(aes(y=value,group=model,colour=model,linetype=model),size=1.5)+
    facet_wrap(~category_f2,scales="free_y")+
    # scale_colour_manual(values = c("blue","red"))+
    # scale_fill_manual(values = c("blue","red"))+
    # scale_linetype_manual(values=c("solid","dashed"))+
    # # scale_size_manual(values=c(100,100,50))+
    theme_bw(base_size=18)+
    theme(legend.position=legend_position,legend.key.size=unit(legend_key_size, "cm"),
          legend.title=element_blank(),legend.text=element_text(size=rel(legend_text_size))
    )+ylim(0,NA)+
    # guides(linetype = guide_legend(override.aes = list(size = c(1.5,1.5))))+
    ylab("")
  g1
}

#' plotting retrospective analyses of samuika
#' @import ggplot2
#' @import dplyr
#' @importFrom  dplyr select
#' @importFrom  dplyr filter
#' @param samuika_retro samuika_retro object
#' @param is_number whether stock and spawning numbers are plotted instead of biomass
#' @param only_NF whether plot only N and F (remove S and C)
#' @encoding UTF-8
#' @export
plot_retro_samuika=function(samuika_retro,
                            mohn=NULL,Stock_ID=0,latest_year=NULL,
                            max_retro_year = NULL,forecast_year=0,
                            is_number = FALSE, only_NF = FALSE) {
  message(paste0("plotting Stock_ID=",Stock_ID," OK?"))
  if (is.null(latest_year)) {
    latest_year = max(samuika_retro$Summary_PopDyn$Year)
  }
  if (is.null(max_retro_year)) {
    max_retro_year = max(samuika_retro$Summary_PopDyn$retro_year)
  }
  Stock_ID_dbl = Stock_ID
  retro_table = samuika_retro$Summary_PopDyn %>%
    mutate(eval_year = latest_year+forecast_year-retro_year) %>%
    filter(Year < eval_year+forecast_year+1) %>%
    filter(Stock_ID == Stock_ID_dbl)
  if (is_number) {
    retro_table = retro_table %>%
      dplyr::select(Stock_ID, Year, Stock_number, Spawning_number, F, Catch_est, retro_year) %>%
      gather(key=variable, value=value, -Stock_ID, -Year, -retro_year)
    retro_table = retro_table %>%
      mutate(variable_f = factor(variable, levels=c("Stock_number","Spawning_number","F","Catch_est")),
             retro_year_f = factor(retro_year, levels=0:max_retro_year))
    if (only_NF) {
      retro_table = retro_table %>%
        filter(variable == "Stock_number"|variable == "F")
    }
  } else {
    retro_table = retro_table %>%
      dplyr::select(Stock_ID, Year, Stock_biomass, Spawning_biomass, F, Catch_est, retro_year) %>%
      gather(key=variable, value=value, -Stock_ID, -Year, -retro_year)
    retro_table = retro_table %>%
      mutate(variable_f = factor(variable, levels=c("Stock_biomass","Spawning_biomass","F","Catch_est")),
             retro_year_f = factor(retro_year, levels=0:max_retro_year))
    if (only_NF) {
      retro_table = retro_table %>%
        filter(variable == "Stock_biomass"|variable == "F")
    }
  }

  max_value = retro_table %>% group_by(variable_f) %>%
    summarise(value=max(value)*1.15,Year=min(Year))

  g4 = ggplot(data=filter(retro_table,retro_year>0), aes(x=Year, y=value)) +
    geom_path(data = filter(retro_table,retro_year==0),aes(group=retro_year_f,colour=retro_year_f),size=1.5,show.legend=FALSE,colour="black")+
    geom_path(aes(group=retro_year_f,colour=retro_year_f),size=1.2,show.legend=FALSE)+
    facet_wrap(~variable_f,ncol=2,scales="free_y")+
    geom_blank(data = max_value)+
    theme_bw(base_size=18)+ylab("")+
    ylim(0,NA)

  if (!is.null(mohn)) {
    if (is_number) {
      rho_data = mohn$Summary %>%
        dplyr::filter(Stock_ID==Stock_ID_dbl) %>%
        dplyr::select(-B,-SSB) %>%
        rename(Stock_number=N,Spawning_number=SSN,Catch_est=Catch) %>%
        gather(key=variable,value=rho,-Stock_ID) %>%
        mutate(variable_f = factor(variable, levels=c("Stock_number","Spawning_number","F","Catch_est")))
      if (only_NF) {
        rho_data = rho_data %>%
          filter(variable == "Stock_number"|variable == "F")
      }
    } else {
      rho_data = mohn$Summary %>%
        dplyr::filter(Stock_ID==Stock_ID_dbl) %>%
        dplyr::select(-N,-SSN) %>%
        rename(Stock_biomass=B,Spawning_biomass=SSB,Catch_est=Catch) %>%
        gather(key=variable,value=rho,-Stock_ID) %>%
        mutate(variable_f = factor(variable, levels=c("Stock_biomass","Spawning_biomass","F","Catch_est")))
      if (only_NF) {
        rho_data = rho_data %>%
          filter(variable == "Stock_biomass"|variable == "F")
      }
    }
    rho_data = rho_data %>%
      mutate(label=sprintf("rho == %.2f",rho)) %>%
      mutate(Year=min(retro_table$Year),value=max_value$value)
    g4 = g4 +
      geom_text(data=rho_data, parse=TRUE, aes(label=label,hjust=0,vjust=1),size=6)

  }
  g4
}

#' plotting stock-recruitment curve from samuika objects
#' @import ggplot2
#' @import dplyr
#' @importFrom  dplyr select
#' @importFrom  dplyr filter
#' @inheritParams get_SRdata
#' @param list_samuika_res list of samuika objects
#' @param model_name model names used in legend
#' @encoding UTF-8
#' @export
plot_SR_samuika = function(list_samuika_res,
                           model_name = NULL,
                           Stock_ID=0,
                           plot_point = TRUE,
                           draw_replace = TRUE,
                           Smsy = NULL,
                           xmax_ratio = 1.3,
                           ymax_ratio = 1.1,
                           path_size=1.2,
                           legend_position = "right",
                           legend_text_size = 1.2,
                           legend_key_size = 1,
                           palette_name = "Set1",
                           point_size = 2) {
  message(paste0("plotting Stock_ID=",Stock_ID," OK?"))
  NRES = length(list_samuika_res)
  name_list = NULL

  for (i in 1:NRES) {
    res = list_samuika_res[[i]]
    if (!is.null(model_name)) {
      name = model_name[i]
    } else {
      name = names(list_samuika_res)[i]
    }
    if (is.null(name)) {
      name = paste0("model",i)
    }
    name_list = c(name_list,name)

    if (i==1) {
      SRdata_est = get_SRdata(res) %>%
        mutate(model0 = name)
    } else {
      SRdata_est = full_join(SRdata_est,get_SRdata(res) %>%
                               mutate(model0 = name))
    }
  }
  Stock_ID_dbl = Stock_ID
  SRdata_est = SRdata_est %>%
    mutate(model = factor(model0, level = name_list)) %>%
    filter(Stock_ID == Stock_ID_dbl)

  xmax = max(SRdata_est$Spawning_number)*xmax_ratio
  ymax = max(SRdata_est$Stock_number)*ymax_ratio
  SSN = seq(0,xmax,length=1000)

  for (i in 1:NRES) {
    res = list_samuika_res[[i]]
    if (length(res$input$regime_key)!= 1 && !is.null(res$input$regime_year)) {
      warning("Use 'plot_SR_regime()' to plot different regimes; only the first regime SR curve is plotted")
    }
    Stock_ID_dbl = Stock_ID
    a = res$Summary_PopDyn %>% filter(Stock_ID==Stock_ID_dbl) %>% select(rec_a) %>% unique()
    a = a[[1]][1]
    b = res$Summary_PopDyn %>% filter(Stock_ID==Stock_ID_dbl) %>% select(rec_b) %>% unique()
    b = b[[1]][1]

    R_pred = SRF(a=a,b=b,x=SSN,SR=res$input$SR)
    if (i==1) {
      SRdata_pred = tibble(Stock_number=R_pred,Spawning_number=SSN) %>%
        mutate(model0 = name_list[i])
    } else {
      SRdata_pred = full_join(SRdata_pred,tibble(Stock_number=R_pred,Spawning_number=SSN) %>%
                                mutate(model0 = name_list[i]))
    }
  }
  Stock_ID_dbl = Stock_ID
  SRdata_pred = SRdata_pred %>%
    mutate(model = factor(model0, level = name_list)) %>%
    filter(Stock_ID == Stock_ID_dbl)

  g1 = ggplot(data=NULL,aes(x=Spawning_number,y=Stock_number))+
    geom_path(data = SRdata_pred, aes(group=model,colour=model),size=path_size)
  if (plot_point) {
    g1 = g1 + geom_point(data = SRdata_est, aes(group=model,colour=model),
                         size=point_size)
  }
  g1 = g1 + coord_cartesian(xlim=c(0,xmax),ylim=c(0,ymax),expand=0) +
    theme_bw(base_size=18)+
    scale_colour_brewer(palette=palette_name)+
    theme(legend.position=legend_position,legend.key.size=unit(legend_key_size, "cm"),
          legend.title=element_blank(),legend.text=element_text(size=rel(legend_text_size)))


  if (draw_replace) {
    g1 = g1 +
      geom_abline(intercept=0,slope = exp(list_samuika_res[[1]]$input$M[1]),colour="gray",size=1,linetype="dashed")
  }

  if (!is.null(Smsy)) {
    if (length(Smsy) != length(list_samuika_res)) {
      stop("The length of SBmsy does not macth with that of list_samuika_res")
    }
    SBmsy_tbl = tibble(Spawning_number = Smsy, model0 = name_list) %>%
      mutate(model = factor(model0, level = name_list))
    g1 = g1+
      geom_vline(data = SBmsy_tbl, aes(xintercept = Spawning_number, group=model, colour=model),
                 linetype="dotted",size=1.2)
  }
  g1
}


#' plotting stock-recruitment curve with regimes from a samuika object
#' @import ggplot2
#' @import dplyr
#' @import ggrepel
#' @importFrom  dplyr select
#' @importFrom  dplyr filter
#' @inheritParams get_SRdata
#' @param list_samuika_res list of samuika objects
#' @param model_name model names used in legend
#' @param regime_name regime names used in legend
#' @encoding UTF-8
#' @export
plot_SR_regime = function(list_samuika_res,
                           model_name = NULL,
                           Stock_ID=0,
                           regime_name=NULL,
                           plot_point = TRUE,
                           draw_replace = TRUE,
                           Smsy = NULL,
                           xmax_ratio = 1.3,
                           ymax_ratio = 1.1,
                           path_size=1.2,
                           legend_position = "right",
                           legend_text_size = 1.2,
                           legend_key_size = 1,
                          legend_off = TRUE,
                          palette_name = "Set1",
                          point_size=2,
                          connect_point = TRUE,
                          connect_point_size=1,
                          label_year=NULL) {
  message(paste0("plotting Stock_ID=",Stock_ID," OK?"))
  if (class(list_samuika_res) == "samuika") {
    list_samuika_res <- list(list_samuika_res)
  }
  NRES = length(list_samuika_res)
  name_list = NULL

  for (i in 1:NRES) {
    res = list_samuika_res[[i]]
    if (!is.null(model_name)) {
      name = model_name[i]
    } else {
      name = names(list_samuika_res)[i]
    }
    if (is.null(name)) {
      name = paste0("model",i)
    }
    name_list = c(name_list,name)

    if (i==1) {
      SRdata_est = get_SRdata(res) %>%
        mutate(model0 = name)
    } else {
      SRdata_est = full_join(SRdata_est,get_SRdata(res) %>%
                               mutate(model0 = name))
    }
  }
  Stock_ID_dbl = Stock_ID
  SRdata_est = SRdata_est %>%
    mutate(model = factor(model0, level = name_list)) %>%
    filter(Stock_ID == Stock_ID_dbl)

  xmax = max(SRdata_est$Spawning_number)*xmax_ratio
  ymax = max(SRdata_est$Stock_number)*ymax_ratio
  SSN = seq(0,xmax,length=1000)

  for (i in 1:NRES) {
    res = list_samuika_res[[i]]
    Stock_ID_dbl = Stock_ID
    if (is.null(regime_name)) {
      regime_name2 = res$input$regime_key %>% as.character() %>% unique()
    } else {
      if (class(regime_name) == "list") {
        if (length(regime_name) != NRES) {
          stop("'length(regime_name)' must satisfy the number of models")
          regime_name2 = regime_name[[i]]
        }
      } else {
        regime_name2 = regime_name
      }
    }

    ab = res$Summary_PopDyn %>% dplyr::filter(Stock_ID==Stock_ID_dbl) %>%
      dplyr::select(rec_a,rec_b,rec_sd) %>% distinct()
    ab = ab %>%
      dplyr::mutate(regime = regime_name2)

    for (j in 1:nrow(ab)) {
      R_pred = SRF(a=as.numeric(ab$rec_a[j]),b=as.numeric(ab$rec_b[j]),x=SSN,SR=res$input$SR)
      if (i==1 && j==1) {
        SRdata_pred = tibble(Stock_number=R_pred,Spawning_number=SSN) %>%
          mutate(model0 = name_list[i],regime=regime_name2[j])
      } else {
        SRdata_pred = full_join(SRdata_pred,tibble(Stock_number=R_pred,Spawning_number=SSN) %>%
                                  mutate(model0 = name_list[i],regime=regime_name2[j]))
      }
    }
  }
  unique_chr= SRdata_pred$regime %>% unique()
  # unique_dbl= SRdata_est$Regime %>% unique()
  SRdata_est = SRdata_est %>%
    mutate(regime = unique_chr[Regime+1])

  Stock_ID_dbl = Stock_ID
  SRdata_pred = SRdata_pred %>%
    mutate(model = factor(model0, level = name_list)) %>%
    filter(Stock_ID == Stock_ID_dbl)

  if (is.null(label_year)) {
    label_year = c(SRdata_est$Year[SRdata_est$Year %% 5 == 0],range(SRdata_est$Year)) %>%
      unique() %>% sort()
    }
  SRdata_est = SRdata_est %>%
    mutate(label = ifelse(Year %in% label_year, Year, NA))

  g1 = ggplot(data=NULL,aes(x=Spawning_number,y=Stock_number))+
    geom_path(data = SRdata_pred,
              aes(group=interaction(model,regime),colour=regime,linetype=model),
              size=path_size)
  if (plot_point) {
    g1 = g1 + geom_point(data = SRdata_est,
                         aes(group=interaction(model,regime),colour=regime,shape=model),
                         size=point_size)+
      scale_shape_discrete(solid=T)
    if (connect_point) {
      g1 = g1 + geom_path(data = SRdata_est,aes(group=model),colour="darkgray",size=connect_point_size)+
        ggrepel::geom_label_repel(data = SRdata_est,aes(label=label))
    }
  }
  g1 = g1 + coord_cartesian(xlim=c(0,xmax),ylim=c(0,ymax),expand=0) +
    theme_bw(base_size=18)+
    scale_colour_brewer(palette=palette_name)+
    theme(legend.position=legend_position,legend.key.size=unit(legend_key_size, "cm"),
          legend.title=element_blank(),legend.text=element_text(size=rel(legend_text_size)))
  if (legend_off) {
    g1 = g1 + theme(legend.position="none")
  }

  if (draw_replace) {
    g1 = g1 +
      geom_abline(intercept=0,slope = exp(list_samuika_res[[1]]$input$M[1]),colour="gray",size=1,linetype="dashed")
  }

  if (!is.null(Smsy)) {
    if (length(Smsy) != length(list_samuika_res)) {
      stop("The length of SBmsy does not macth with that of list_samuika_res")
    }
    SBmsy_tbl = tibble(Spawning_number = Smsy, model0 = name_list) %>%
      mutate(model = factor(model0, level = name_list))
    g1 = g1+
      geom_vline(data = SBmsy_tbl, aes(xintercept = Spawning_number, group=model, colour=model),
                 linetype="dotted",size=1.2)
  }
  g1

}




#' Function for plotting future_res
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom dplyr select
#' @importFrom tidyr gather
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
  BRP = NULL, # under consideration...
  line_size = 1,
  legend_position = "right",
  legend_text_size = 1.2,
  legend_key_size = 1,
  strip_text_size = 10,
  base_size = 16,
  BRP_line_size = 1,
  palette_name = "Set1",
  alpha = 0.4
) {

  if (is.null(scenario_name)) {
    for (i in 1:length(future_list)) {
      scenario_name = c(scenario_name, sprintf("scenario_%s",i))
    }
  } else {
    if (length(scenario_name) != length(future_list)) {
      stop("The length of scenario_name does not match with that of future_list")
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
  plot_data_all = plot_data_all %>%
    rename(Scenario0 = Scenario) %>%
    mutate(Scenario = factor(Scenario0, levels=scenario_name)) %>%
    select(-Scenario0)

  plot_data_center = mutate(plot_data_all,
                            Y_f = factor(plot_data_all$Y,
                                         levels=c("Stock_biomass","Spawning_biomass","Catch_biomass","F"))) %>%
    filter(Range == "center") %>%
    rename(Scenario0 = Scenario) %>%
    mutate(Scenario = factor(Scenario0, levels=scenario_name)) %>%
    select(-Scenario0)

  plot_data_CI =  mutate(plot_data_all,
                         Y_f = factor(plot_data_all$Y,
                                      levels=c("Stock_biomass","Spawning_biomass","Catch_biomass","F"))) %>%
    filter(Range != "center") %>%
    spread(key = Range, value=value) %>%
    rename(Scenario0 = Scenario) %>%
    mutate(Scenario = factor(Scenario0, levels=scenario_name)) %>%
    select(-Scenario0)

  g1 = ggplot(plot_data_center)

  if (!is.null(BRP)) {
    BRP_tbl = dplyr::select(BRP,Definition,Stock_biomass,Spawning_biomass,Catch_biomass,F) %>%
      tidyr::gather(key=Y,value=value,-Definition) %>%
      mutate(Y_f = factor(Y,levels=c("Stock_biomass","Spawning_biomass","Catch_biomass","F"))) %>%
      mutate(Def_f = factor(Definition,levels=BRP$Definition))
    g1 = g1 +
      geom_hline(data = BRP_tbl, aes(yintercept=value,linetype=Def_f),colour="darkgray",size=BRP_line_size)
    if (legend_position == "bottom") {
      g1 = g1 + guides(fill=guide_legend(nrow=2),colour = guide_legend(nrow = 2),linetype = guide_legend(nrow = 2))
    }
  }

  if (!is.null(CI_range)) {
    g1 = g1 + geom_ribbon(data = plot_data_CI,
                              aes(x=Year,ymin=lower,ymax=upper,
                                  group=Scenario, fill=Scenario),alpha=alpha)+
      scale_fill_brewer(palette=palette_name)
  }
  g1 = g1 + geom_line(data = dplyr::filter(plot_data_center),
              aes(x=Year,y=value,group=Scenario, colour=Scenario),size=line_size)+
    geom_line(data = dplyr::filter(plot_data_center,Status=="Past"&Scenario==scenario_name[1]),
              aes(x=Year,y=value),size=line_size, colour="black") +
    scale_colour_brewer(palette=palette_name)+
    facet_wrap(~Y_f, scales="free_y") +
    ylim(0,NA) +
    theme_bw(base_size=base_size)+
    theme(axis.title.y = element_blank())+
    theme(legend.position=legend_position,legend.key.size=unit(legend_key_size, "cm"),
          legend.title=element_blank(),legend.text=element_text(size=rel(legend_text_size)))+
    theme(strip.text.x = element_text(size = strip_text_size))

  if (!is.null(title)) g1 = g1 + ggtitle(title)

  return(g1)
}

#' Function for calculating probability of achiveing biological reference points
#' @import dplyr
#' @import tidyr
#' @param future_sim_res \code{future_sim} object
#' @param ref_value value of biological reference points
#' @param variable what variable ref_value indicates (default: "Spawning_biomass")
#' @encoding UTF-8
#' @export
get_prob_ref = function(future_sim_res, ref_value, variable="Spawning_biomass") {
  sim_table = future_sim_res$sim_array[,variable,]
  func = function(x) mean(x>ref_value)
  prob = apply(sim_table,1,func)
  names(prob) <- rev(rev(future_sim_res$mean_table$Year)[1:length(prob)])
  return(prob)
}


#' Function for converting samuika result to the data for est_MSY
#' @import dplyr
#' @param samuika_res \code{samuika} object
#' @param Stock_ID Stock_ID (default is 0)
#' @encoding UTF-8
#' @export
convert_samuika_estMSY = function(samuika_res, Stock_ID = 0, scale_num_to_mass = NULL) {
  message(paste0("Using Stock_ID=",Stock_ID," OK?"))
  model = samuika_res
  M = samuika_res$input$M
  if (is.null(scale_num_to_mass)) scale_num_to_mass = model$input$scale_num_to_mass
  Stock_ID_tbl = as.numeric(Stock_ID)
  assess_data= dplyr::filter(model$Summary_PopDyn, Stock_ID==Stock_ID_tbl) %>%
    mutate(Weight = Stock_biomass/(Stock_number*scale_num_to_mass),
           M = M,
           U = Catch_biomass/Stock_biomass)
  invisible(assess_data)
}

#' Function for plotting time series of catch biomass
#' @import ggplot2
#' @inheritParams out_summary_estimate
#' @param model \code{samuika} object
#' @encoding UTF-8
#' @export
plot_catch = function(model,CI=0.8,plot_CI=TRUE,Stock_name=NULL,plot_obs=TRUE,
                      base_size=18,path_size=1.5,colour="blue",point_size=2,
                      alpha=0.4,scales="fix",add_MSY=TRUE,
                      MSY_path_size=1,MSY_path_colour = "orange",MSY_path_linetype = "dashed"){
  data_est = out_summary_estimate(model,CI=CI) %>%
    filter(category == "Catch")
  if (is.null(Stock_name)) {
    data_est = data_est %>% mutate(Stock_name = as.character(Stock_ID))
  } else {
    data_est$Stock_name = Stock_name[data_est$Stock_ID+1]
  }
  gg = ggplot(data=data_est,aes(x=Year))
  if (plot_CI) {
    gg = gg + geom_ribbon(aes(ymin=lower,ymax=upper),alpha=alpha,fill=colour)
  }
  ymax = data_est$upper %>% max
  ymax = ymax*1.05
  xmin = data_est$Year %>% min -1
  xmax = data_est$Year %>% max + 1
  gg = gg +
    geom_path(aes(y=value),size=path_size,colour=colour)+
    # scale_colour_brewer(palette=palette_name)+
    theme_bw(base_size=base_size)+ylab("Catch")+
    ylim(0,NA)+xlim(xmin,xmax)+
    facet_wrap(vars(Stock_name),scales=scales)
  if (plot_obs) {
    data_obs = dplyr::select(model$Summary_PopDyn,Stock_ID,Year,Catch_biomass) %>% na.omit() %>%
      rename(value = Catch_biomass)
    if (is.null(Stock_name)) {
      data_obs = data_obs %>% mutate(Stock_name = as.character(Stock_ID))
    } else {
      data_obs$Stock_name = Stock_name[data_obs$Stock_ID+1]
    }
    gg = gg + geom_point(data = data_obs,aes(y=value),size=point_size)
  }
  if (add_MSY) {
    Summary_PopDyn = integrate_detBRF(model) %>%
      dplyr::rename(value = MSY) %>%
      dplyr::full_join(data_est %>% dplyr::select(Stock_ID,Year,Stock_name))
    gg = gg + geom_path(data=Summary_PopDyn,aes(x=Year,y=value),
                        size=MSY_path_size,colour=MSY_path_colour,linetype=MSY_path_linetype)
  }
  gg
}

#' Function for plotting time-series of stock biomass
#' @import ggplot2
#' @inheritParams out_summary_estimate
#' @param model \code{samuika} object
#' @encoding UTF-8
#' @export
plot_stock_biomass = function(model,CI=0.8,plot_CI=TRUE,Stock_name=NULL,
                      base_size=18,path_size=1.5,colour="blue",point_size=2,
                      alpha=0.4,scales="fix",add_MSY = FALSE,
                      MSY_path_size=1,MSY_path_colour = "orange",MSY_path_linetype = "dashed"){
  data_est = out_summary_estimate(model,CI=CI) %>%
    filter(category == "Stock_biomass")
  if (is.null(Stock_name)) {
    data_est = data_est %>% mutate(Stock_name = as.character(Stock_ID))
  } else {
    data_est$Stock_name = Stock_name[data_est$Stock_ID+1]
  }
  gg = ggplot(data=data_est,aes(x=Year))
  if (plot_CI) {
    gg = gg + geom_ribbon(aes(ymin=lower,ymax=upper),alpha=alpha,fill=colour)
  }
  ymax = data_est$upper %>% max
  ymax = ymax*1.05
  xmin = data_est$Year %>% min -1
  xmax = data_est$Year %>% max + 1
  gg = gg +
    geom_path(aes(y=value),size=path_size,colour=colour)+
    # scale_colour_brewer(palette=palette_name)+
    theme_bw(base_size=base_size)+ylab("Stock biomass")+
    ylim(0,NA)+xlim(xmin,xmax)+
    facet_wrap(vars(Stock_name),scales=scales)
  if (add_MSY) {
    Summary_PopDyn = integrate_detBRF(model) %>%
      dplyr::mutate(value = Nmsy*Weight*model$input$scale_num_to_mass) %>%
      dplyr::full_join(data_est %>% dplyr::select(Stock_ID,Year,Stock_name))
    gg = gg + geom_path(data=Summary_PopDyn,aes(x=Year,y=value),
                        size=MSY_path_size,colour=MSY_path_colour,linetype=MSY_path_linetype)
  }
  gg
}

#' Function for plotting time-series of stock number
#' @import ggplot2
#' @inheritParams out_summary_estimate_n
#' @param model \code{samuika} object
#' @encoding UTF-8
#' @export
plot_stock_number = function(model,CI=0.8,plot_CI=TRUE,Stock_name=NULL,plot_obs=TRUE,
                              base_size=18,path_size=1.5,colour="blue",point_size=2,
                              alpha=0.4,scales="fix",add_MSY=FALSE,
                             MSY_path_size=1,MSY_path_colour = "orange",MSY_path_linetype = "dashed"){
  data_est = out_summary_estimate_n(model,CI=CI) %>%
    filter(category == "Stock_number")
  if (is.null(Stock_name)) {
    data_est = data_est %>% mutate(Stock_name = as.character(Stock_ID))
  } else {
    data_est$Stock_name = Stock_name[data_est$Stock_ID+1]
  }
  gg = ggplot(data=data_est,aes(x=Year))
  if (plot_CI) {
    gg = gg + geom_ribbon(aes(ymin=lower,ymax=upper),alpha=alpha,fill=colour)
  }
  ymax = data_est$upper %>% max
  ymax = ymax*1.05
  xmin = data_est$Year %>% min -1
  xmax = data_est$Year %>% max + 1
  gg = gg +
    geom_path(aes(y=value),size=path_size,colour=colour)+
    # scale_colour_brewer(palette=palette_name)+
    theme_bw(base_size=base_size)+ylab("Stock number")+
    ylim(0,NA)+xlim(xmin,xmax)+
    facet_wrap(vars(Stock_name),scales=scales)

  if (plot_obs) {
    data_obs = dplyr::select(model$Summary_index,Index_ID,Stock_ID,Year,Index,q) %>% na.omit() %>%
      mutate(value = Index/q)
    if (is.null(Stock_name)) {
      data_obs = data_obs %>% mutate(Stock_name = as.character(Stock_ID))
    } else {
      data_obs$Stock_name = Stock_name[data_obs$Stock_ID+1]
    }
    gg = gg + geom_point(data = data_obs,aes(y=value),size=point_size)
  }
  if (add_MSY) {
    Summary_PopDyn = integrate_detBRF(model) %>%
      dplyr::rename(value = Nmsy) %>%
      dplyr::full_join(data_est %>% dplyr::select(Stock_ID,Year,Stock_name))
    gg = gg + geom_path(data=Summary_PopDyn,aes(x=Year,y=value),
                        size=MSY_path_size,colour=MSY_path_colour,linetype=MSY_path_linetype)
  }
  gg
}

#' Function for plotting time-series of spawning number
#' @import ggplot2
#' @inheritParams out_summary_estimate_n
#' @param model \code{samuika} object
#' @encoding UTF-8
#' @export
plot_spawning_number = function(model,CI=0.8,plot_CI=TRUE,Stock_name=NULL,
                             base_size=18,path_size=1.5,colour="blue",point_size=2,
                             alpha=0.4,scales="fix",add_MSY=FALSE,
                             MSY_path_size=1,MSY_path_colour = "orange",MSY_path_linetype = "dashed"){
  data_est = out_summary_estimate_n(model,CI=CI) %>%
    filter(category == "Spawning_number")
  if (is.null(Stock_name)) {
    data_est = data_est %>% mutate(Stock_name = as.character(Stock_ID))
  } else {
    data_est$Stock_name = Stock_name[data_est$Stock_ID+1]
  }
  gg = ggplot(data=data_est,aes(x=Year))
  if (plot_CI) {
    gg = gg + geom_ribbon(aes(ymin=lower,ymax=upper),alpha=alpha,fill=colour)
  }
  ymax = data_est$upper %>% max
  ymax = ymax*1.05
  xmin = data_est$Year %>% min -1
  xmax = data_est$Year %>% max + 1
  gg = gg +
    geom_path(aes(y=value),size=path_size,colour=colour)+
    # scale_colour_brewer(palette=palette_name)+
    theme_bw(base_size=base_size)+ylab("Spawning number")+
    ylim(0*ymax,NA)+xlim(xmin,xmax)+
    # coord_cartesian(xlim = c(xmin,xmax),ylim=c(0,ymax),expand=0)+
    facet_wrap(vars(Stock_name),scales=scales)
  if (add_MSY) {
    Summary_PopDyn = integrate_detBRF(model) %>%
      dplyr::rename(value = Smsy) %>%
      dplyr::full_join(data_est %>% dplyr::select(Stock_ID,Year,Stock_name))
    gg = gg + geom_path(data=Summary_PopDyn,aes(x=Year,y=value),
                        size=MSY_path_size,colour=MSY_path_colour,linetype=MSY_path_linetype)
  }
  gg
}

#' Function for plotting time-series of spawning biomass
#' @import ggplot2
#' @inheritParams out_summary_estimate
#' @param model \code{samuika} object
#' @encoding UTF-8
#' @export
plot_spawning_biomass = function(model,CI=0.8,plot_CI=TRUE,Stock_name=NULL,
                              base_size=18,path_size=1.5,colour="blue",point_size=2,
                              alpha=0.4,scales="fix",add_MSY,
                              MSY_path_size=1,MSY_path_colour = "orange",MSY_path_linetype = "dashed"){
  data_est = out_summary_estimate(model,CI=CI) %>%
    filter(category == "Spawning_biomass")
  if (is.null(Stock_name)) {
    data_est = data_est %>% mutate(Stock_name = as.character(Stock_ID))
  } else {
    data_est$Stock_name = Stock_name[data_est$Stock_ID+1]
  }
  gg = ggplot(data=data_est,aes(x=Year))
  if (plot_CI) {
    gg = gg + geom_ribbon(aes(ymin=lower,ymax=upper),alpha=alpha,fill=colour)
  }
  ymax = data_est$upper %>% max
  ymax = ymax*1.05
  xmin = data_est$Year %>% min -1
  xmax = data_est$Year %>% max + 1
  gg = gg +
    geom_path(aes(y=value),size=path_size,colour=colour)+
    # scale_colour_brewer(palette=palette_name)+
    theme_bw(base_size=base_size)+ylab("Stock biomass")+
    ylim(0,NA)+xlim(xmin,xmax)+
    facet_wrap(vars(Stock_name),scales=scales)
  if (add_MSY) {
    Summary_PopDyn = integrate_detBRF(model) %>%
      dplyr::mutate(value = Smsy*Weight*model$input$scale_num_to_mass) %>%
      dplyr::full_join(data_est %>% dplyr::select(Stock_ID,Year,Stock_name))
    gg = gg + geom_path(data=Summary_PopDyn,aes(x=Year,y=value),
                        size=MSY_path_size,colour=MSY_path_colour,linetype=MSY_path_linetype)
  }
  gg
}


#' Function for plotting time-series of fishing mortality coefficient
#' @import ggplot2
#' @inheritParams out_summary_estimate
#' @param model \code{samuika} object
#' @encoding UTF-8
#' @export
plot_F = function(model,CI=0.8,plot_CI=TRUE,Stock_name=NULL,
                                 base_size=18,path_size=1.5,colour="blue",point_size=2,
                                 alpha=0.4,scales="fix",add_MSY=FALSE,
                                 MSY_path_size=1,MSY_path_colour = "orange",MSY_path_linetype = "dashed"){
  data_est = out_summary_estimate(model,CI=CI) %>%
    filter(category == "F")
  if (is.null(Stock_name)) {
    data_est = data_est %>% mutate(Stock_name = as.character(Stock_ID))
  } else {
    data_est$Stock_name = Stock_name[data_est$Stock_ID+1]
  }
  gg = ggplot(data=data_est,aes(x=Year))
  if (plot_CI) {
    gg = gg + geom_ribbon(aes(ymin=lower,ymax=upper),alpha=alpha,fill=colour)
  }
  ymax = data_est$upper %>% max
  ymax = ymax*1.05
  xmin = data_est$Year %>% min -1
  xmax = data_est$Year %>% max + 1
  gg = gg +
    geom_path(aes(y=value),size=path_size,colour=colour)+
    # scale_colour_brewer(palette=palette_name)+
    theme_bw(base_size=base_size)+ylab("F")+
    ylim(0,NA)+xlim(xmin,xmax)+
    facet_wrap(vars(Stock_name),scales=scales)
  if (add_MSY) {
    Summary_PopDyn = integrate_detBRF(model) %>%
      dplyr::rename(value = Fmsy) %>%
      dplyr::full_join(data_est %>% dplyr::select(Stock_ID,Year,Stock_name))
    gg = gg + geom_path(data=Summary_PopDyn,aes(x=Year,y=value),
                        size=MSY_path_size,colour=MSY_path_colour,linetype=MSY_path_linetype)
  }
  gg
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


#' output summary table when changing beta
#' @import ggplot2
#' @import dplyr
#' @inheritParams get_prob_ref
#' @param future_HCR \code{future_sim} object with \code{HCR = TRUE}
#' @param SBmsy SBmsy as a target reference point
#' @param Fmsy Fmsy to be multiplied by beta
#' @param beta vector of beta (multipliers of Fmsy)
#' @param future_Fcurrent \code{future_sim} object with Fcurrent (optional)
#' @param filename file name of saving
#' @encoding UTF-8
#' @export
out_beta_table = function(future_HCR, SBmsy, Fmsy = NULL, beta = seq(0,1,by = 0.1), future_Fcurrent = NULL,
                          filename = "beta_table") {
  if (is.null(Fmsy)) stop("Set Fmsy to be multiplied by beta!!!")
  if (is.null(SBmsy)) stop("Set SBmsy for calculating achievement probability!!!")

  SBlimit = future_HCR$input$SBrefs[1]
  SBban = future_HCR$input$SBrefs[2]

  if (!is.null(future_Fcurrent)) {
    Ptarget_Fcurrent = get_prob_ref(future_Fcurrent, ref_value = SBmsy)
    Plimit_Fcurrent = get_prob_ref(future_Fcurrent, ref_value = SBlimit)
    Pban_Fcurrent = get_prob_ref(future_Fcurrent, ref_value = SBban)
    summary_table = future_Fcurrent$mean_table %>%
      dplyr::filter(Status == "Future") %>%
      dplyr::mutate(Ptarget = Ptarget_Fcurrent, Plimit = Plimit_Fcurrent, Pban = Pban_Fcurrent) %>%
      dplyr::mutate(scenario = "Fcurrent", beta = NA,scenario_no = 0)
    summary_table0 = summary_table
  }

  beta_seq = beta
  for (j in 1:length(beta_seq)) {
    cat("---",j,"---\n")
    beta2 = beta_seq[j]
    input_tmp = future_HCR$input
    input_tmp$Ftarget = beta2*Fmsy
    future_HCR2 = do.call(future_sim, input_tmp)

    Ptarget = get_prob_ref(future_HCR2, ref_value = SBmsy)
    Plimit = get_prob_ref(future_HCR2, ref_value = SBlimit)
    Pban = get_prob_ref(future_HCR2, ref_value = SBban)

    beta_table = future_HCR2$mean_table %>%
      dplyr::select(Year,Stock_biomass,Spawning_biomass,Catch_biomass,F,
                    Stock_number,Spawning_number,Status) %>%
      dplyr::filter(Status == "Future") %>%
      dplyr::mutate(Ptarget = Ptarget, Plimit = Plimit, Pban = Pban) %>%
      dplyr::mutate(scenario = paste0("beta=",beta2), beta = beta2,scenario_no = j)
    if (j==1 && is.null(future_Fcurrent)) {
      summary_table2 = beta_table
    } else {
      summary_table2 = full_join(summary_table, beta_table)
    }
    summary_table = summary_table2
  }
  summary_table2 = summary_table %>%
    dplyr::select(scenario_no,scenario, beta, Year,Stock_biomass,Spawning_biomass,Catch_biomass,Ptarget,Plimit,Pban,F)
  write.csv(summary_table2, row.names = FALSE,file = paste0(filename,"_all.csv"))

  catch_mean  = summary_table2 %>% dplyr::select(scenario_no,scenario, beta, Year,Catch_biomass) %>%
    tidyr::spread(key = Year, value = Catch_biomass) %>%
    dplyr::arrange(desc(scenario_no)) %>%
    dplyr::select(-scenario_no)
  write.csv(catch_mean, row.names = FALSE,file = paste0(filename,"_Catch_biomass_mean.csv"))

  stock_biomass_mean  = summary_table2 %>% dplyr::select(scenario_no,scenario, beta, Year,Stock_biomass) %>%
    tidyr::spread(key = Year, value = Stock_biomass) %>%
    dplyr::arrange(desc(scenario_no)) %>%
    dplyr::select(-scenario_no)
  write.csv(stock_biomass_mean, row.names = FALSE,file = paste0(filename,"_Stock_biomass_mean.csv"))

  spawning_biomass_mean  = summary_table2 %>% dplyr::select(scenario_no,scenario, beta, Year,Spawning_biomass) %>%
    tidyr::spread(key = Year, value = Spawning_biomass) %>%
    dplyr::arrange(desc(scenario_no)) %>%
    dplyr::select(-scenario_no)
  write.csv(spawning_biomass_mean, row.names = FALSE,file = paste0(filename,"_Spawning_biomass_mean.csv"))

  F_mean  = summary_table2 %>% dplyr::select(scenario_no,scenario, beta, Year,F) %>%
    tidyr::spread(key = Year, value = F) %>%
    dplyr::arrange(desc(scenario_no)) %>%
    dplyr::select(-scenario_no)
  write.csv(F_mean, row.names = FALSE,file = paste0(filename,"_F_mean.csv"))

  Ptarget_table  = summary_table2 %>% dplyr::select(scenario_no,scenario, beta, Year,Ptarget) %>%
    tidyr::spread(key = Year, value = Ptarget) %>%
    dplyr::arrange(desc(scenario_no)) %>%
    dplyr::select(-scenario_no)
  write.csv(Ptarget_table, row.names = FALSE,file = paste0(filename,"_Ptarget.csv"))

  Plimit_table  = summary_table2 %>% dplyr::select(scenario_no,scenario, beta, Year,Plimit) %>%
    tidyr::spread(key = Year, value = Plimit) %>%
    dplyr::arrange(desc(scenario_no)) %>%
    dplyr::select(-scenario_no)
  write.csv(Plimit_table, row.names = FALSE,file = paste0(filename,"_Plimit.csv"))

  Pban_table  = summary_table2 %>% dplyr::select(scenario_no,scenario, beta, Year,Pban) %>%
    tidyr::spread(key = Year, value = Pban) %>%
    dplyr::arrange(desc(scenario_no)) %>%
    dplyr::select(-scenario_no)
  write.csv(Pban_table, row.names = FALSE,file = paste0(filename,"_Pban.csv"))
  RES = list(summary_table=summary_table2,catch_mean=catch_mean, stock_biomass_mean=stock_biomass_mean, spawning_biomass_mean=spawning_biomass_mean,
             F_mean=F_mean,Ptarget=Ptarget_table,Plimit=Plimit_table,Pban=Pban_table)
  return(RES)
}


#' output summary table when changing Fmulti
#' @import ggplot2
#' @import dplyr
#' @import stringr
#' @import magrittr
#' @inheritParams get_prob_ref
#' @param future_Fcurrent \code{future_sim} object with Fcurrent
#' @param SBrefs Reference points for calculating achievement probability
#' @param filename file name of saving
#' @encoding UTF-8
#' @export
out_Fmulti_table = function(future_Fcurrent, SBrefs, Fmulti = seq(0,1,by = 0.1),scenario_name = NULL,filename = "Fmulti_table") {
  if (is.null(names(SBrefs))) {
    message("Use names(SBrefs) if one wants to give SBrefs names")
    names(SBrefs) <- stringr::str_c("SBref",1:length(SBrefs))
  }
  Pname = stringr::str_c("P",names(SBrefs))

  if (is.null(scenario_name)) {
    scenario_name = stringr::str_c("Fmulti=",round(Fmulti_seq,3))
  } else {
    if (length(scenario_name) != length(Fmulti)) stop("Lengths of 'scenario_name' and 'Fmulti' must be same")
  }

  Fmulti_seq = Fmulti
  for (j in 1:length(Fmulti_seq)) {
    cat("---",j,"---\n")
    input_tmp = future_Fcurrent$input
    input_tmp$multi = Fmulti_seq[j]
    future_Fmulti = do.call(future_sim, input_tmp)

    Fmulti_table = future_Fmulti$mean_table %>%
      dplyr::select(Year,Stock_biomass,Spawning_biomass,Catch_biomass,F,
                    Stock_number,Spawning_number,Status) %>%
      dplyr::filter(Status == "Future") %>%
      dplyr::mutate(scenario = scenario_name[j], Fmulti = Fmulti_seq[j],scenario_no = j)

    Pvalue = NULL
    for (k in 1:length(SBrefs)) {
      Pvalue = cbind(Pvalue,get_prob_ref(future_Fmulti, ref_value = SBrefs[k]))
    }
    Pyear = as.numeric(rownames(Pvalue))
    Pvalue = as.data.frame(Pvalue) %>%
      magrittr::set_colnames(value = Pname) %>%
      dplyr::mutate(Year = Pyear)
    Fmulti_table = full_join(Fmulti_table,Pvalue,by="Year")

    if (j==1) {
      summary_table = Fmulti_table
    } else {
      summary_table = full_join(summary_table, Fmulti_table)
    }
  }

  summary_table2 = summary_table %>%
    dplyr::select(-Stock_number, -Spawning_number, -Status) %>%
    dplyr::select(scenario_no,scenario, Fmulti, Year,Stock_biomass,Spawning_biomass,Catch_biomass,everything())
  write.csv(summary_table2, row.names = FALSE,file = paste0(filename,"_all.csv"))

  catch_mean  = summary_table2 %>% dplyr::select(scenario_no,scenario, Fmulti, Year,Catch_biomass) %>%
    tidyr::spread(key = Year, value = Catch_biomass) %>%
    dplyr::arrange(desc(scenario_no)) %>%
    dplyr::select(-scenario_no)
  write.csv(catch_mean, row.names = FALSE,file = paste0(filename,"_Catch_biomass_mean.csv"))

  stock_biomass_mean  = summary_table2 %>% dplyr::select(scenario_no,scenario, Fmulti, Year,Stock_biomass) %>%
    tidyr::spread(key = Year, value = Stock_biomass) %>%
    dplyr::arrange(desc(scenario_no)) %>%
    dplyr::select(-scenario_no)
  write.csv(stock_biomass_mean, row.names = FALSE,file = paste0(filename,"_Stock_biomass_mean.csv"))

  spawning_biomass_mean  = summary_table2 %>% dplyr::select(scenario_no,scenario, Fmulti, Year,Spawning_biomass) %>%
    tidyr::spread(key = Year, value = Spawning_biomass) %>%
    dplyr::arrange(desc(scenario_no)) %>%
    dplyr::select(-scenario_no)
  write.csv(spawning_biomass_mean, row.names = FALSE,file = paste0(filename,"_Spawning_biomass_mean.csv"))

  F_mean  = summary_table2 %>% dplyr::select(scenario_no,scenario, Fmulti, Year,F) %>%
    tidyr::spread(key = Year, value = F) %>%
    dplyr::arrange(desc(scenario_no)) %>%
    dplyr::select(-scenario_no)
  write.csv(F_mean, row.names = FALSE,file = paste0(filename,"_F_mean.csv"))
  RES = list(summary_table=summary_table2,catch_mean=catch_mean, stock_biomass_mean=stock_biomass_mean, spawning_biomass_mean=spawning_biomass_mean,
             F_mean=F_mean)

  for (k in 1:length(Pname)) {
    P_table  = summary_table2 %>% dplyr::select(scenario_no,scenario, Fmulti, Year,as.character(Pname[k])) %>%
      tidyr::spread(key = Year, value = as.character(Pname[k])) %>%
      dplyr::arrange(desc(scenario_no)) %>%
      dplyr::select(-scenario_no)
    write.csv(P_table, row.names = FALSE,file = paste0(filename,"_",Pname[k],".csv"))
    RES[[Pname[k]]] = P_table
  }
  return(RES)
}


#' Output summary table for bootstrapped samuika
#' @inheritParams samuika
#' @param boot_res \code{boot_samuika} object
#' @encoding UTF-8
#' @export
out_boot_table = function(boot_res, CI_range = 0.8) {
  RES = list()
  for(i in 1:length(boot_res)) {
    table_tmp = boot_res[[i]]$Summary_PopDyn %>% dplyr::mutate(boot_ID = i)
    if (i==1) {
      all_table = table_tmp
    } else {
      all_table = dplyr::bind_rows(all_table, table_tmp)
    }
  }
  tmp = dplyr::select(all_table,-Stock_number_pred,-rec_resid,-Catch_resid) %>%
    dplyr::select(boot_ID,everything()) %>%
    tidyr::gather(key = "Key", value = "Value", -boot_ID,-Stock_ID,-Year,-Regime)
  summary_table = tmp %>%
    group_by(Key,Stock_ID,Year,Regime) %>%
    summarise(Mean = mean(Value,na.rm=TRUE), Median = median(Value,na.rm=TRUE),
              Lower = quantile(Value, probs = (1-CI_range)/2,na.rm=TRUE),
              Upper = quantile(Value, probs = 1-(1-CI_range)/2,na.rm=TRUE))
  RES$all_table = all_table
  RES$summary_table = summary_table
  return(RES)
}

#' plot of bootstrapped stock number
#' @inheritParams out_boot_table
#' @inheritParams plot_stock_number
#' @param res \code{samuika} object
#' @param boot_res \code{boot_samuika} object
#' @encoding UTF-8
#' @export
plot_boot_stock_number = function(res, boot_res, CI=0.8, Stock_name = NULL,alpha=0.3,path_size=1.5,colour="red",...) {
  g1 = plot_stock_number(res,Stock_name=Stock_name,plot_CI=FALSE,plot_obs=FALSE,path_size=path_size,...)
  summary_boot = out_boot_table(boot_res,CI_range=CI)$summary_table %>%
    dplyr::filter(Key == "Stock_number")
  if (is.null(Stock_name)) {
    summary_boot = summary_boot %>% mutate(Stock_name = as.character(Stock_ID))
  } else {
    summary_boot$Stock_name = Stock_name[summary_boot$Stock_ID+1]
  }
  g2 = g1 + geom_ribbon(data=summary_boot, aes(x=Year,ymin=Lower,ymax=Upper),alpha=alpha,fill=colour)+
    geom_path(data=summary_boot, aes(x=Year,y=Median),size=path_size,colour=colour,linetype = "dashed")
  g2
}

#' plot of bootstrapped fishing mortality coefficient
#' @inheritParams out_boot_table
#' @inheritParams plot_F
#' @param res \code{samuika} object
#' @param boot_res \code{boot_samuika} object
#' @encoding UTF-8
#' @export
plot_boot_F = function(res, boot_res, CI=0.8, Stock_name = NULL,alpha=0.3,path_size=1.5,colour="red",...) {
  g1 = plot_F(res,Stock_name=Stock_name,plot_CI=FALSE,path_size=path_size,...)
  summary_boot = out_boot_table(boot_res,CI_range=CI)$summary_table %>%
    dplyr::filter(Key == "F")
  if (is.null(Stock_name)) {
    summary_boot = summary_boot %>% mutate(Stock_name = as.character(Stock_ID))
  } else {
    summary_boot$Stock_name = Stock_name[summary_boot$Stock_ID+1]
  }
  g2 = g1 + geom_ribbon(data=summary_boot, aes(x=Year,ymin=Lower,ymax=Upper),alpha=alpha,fill=colour)+
    geom_path(data=summary_boot, aes(x=Year,y=Median),size=path_size,colour=colour,linetype = "dashed")
  g2
}

#' plot of bootstrapped Catch
#' @inheritParams out_boot_table
#' @inheritParams plot_catch
#' @param res \code{samuika} object
#' @param boot_res \code{boot_samuika} object
#' @encoding UTF-8
#' @export
plot_boot_catch_obs = function(res, boot_res, CI=0.8, Stock_name = NULL,alpha=0.3,path_size=1.5,colour="red",...) {
  g1 = plot_catch(res,Stock_name=Stock_name,plot_CI=FALSE,path_size=path_size,plot_obs=TRUE,add_MSY=FALSE)
  summary_boot = out_boot_table(boot_res,CI_range=CI)$summary_table %>%
    dplyr::filter(Key == "Catch_biomass")
  if (is.null(Stock_name)) {
    summary_boot = summary_boot %>% mutate(Stock_name = as.character(Stock_ID))
  } else {
    summary_boot$Stock_name = Stock_name[summary_boot$Stock_ID+1]
  }
  g2 = g1 + geom_ribbon(data=summary_boot, aes(x=Year,ymin=Lower,ymax=Upper),alpha=alpha,fill=colour)+
    geom_path(data=summary_boot, aes(x=Year,y=Median),size=path_size,colour=colour,linetype = "dashed")
  g2
}
