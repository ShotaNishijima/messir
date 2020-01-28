
### Ulility Functions in messir ### ----

#' output a table of samuika results
#' @param model \code{samuika} object
#' @encoding UTF-8
#' @export
out_summary_estimate <- function(model,CI=0.8) {
  samuika_model = model
  Summary_PopDyn = samuika_model$Summary_PopDyn %>%
    select(Stock_ID, Year, Stock_biomass, Spawning_biomass, F, Catch_est) %>%
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
    select(-category,-category_f)

  Summary_PopDyn = Summary_PopDyn %>%
    rename(category=category2, category_f=category_f2)

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

#' Making a combined list with original list name(s)
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

  full_join(REC, SSN)
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
    out_samuika = out_summary_estimate(res,CI=CI) %>%
      mutate(model0 = name)
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
#' @param samuika_retro samuika_retro object
#' @encoding UTF-8
#' @export
plot_retro_samuika=function(samuika_retro,
                            mohn=NULL,Stock_ID=0,latest_year=NULL,
                            max_retro_year = NULL,forecast_year=1) {
  message(paste0("plotting Stock_ID=",Stock_ID," OK?"))
  if (is.null(latest_year)) {
    latest_year = max(samuika_retro$Summary_PopDyn$Year)
  }
  if (is.null(max_retro_year)) {
    max_retro_year = max(samuika_retro$Summary_PopDyn$retro_year)
  }
  retro_table = samuika_retro$Summary_PopDyn %>%
    mutate(eval_year = latest_year+forecast_year-retro_year) %>%
    filter(Year < eval_year+forecast_year) %>%
    filter(Stock_ID == Stock_ID) %>%
    select(Stock_ID, Year, Stock_biomass, Spawning_biomass, F, Catch_est, retro_year) %>%
    gather(key=variable, value=value, -Stock_ID, -Year, -retro_year) %>%
    mutate(variable_f = factor(variable, levels=c("Stock_biomass","Spawning_biomass","F","Catch_est")),
           retro_year_f = factor(retro_year, levels=0:max_retro_year))

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
    rho_data = mohn$Summary %>%
      select(-N,-SSN) %>%
      rename(Stock_biomass=B,Spawning_biomass=SSB,Catch_est=Catch) %>%
      gather(key=variable,value=rho) %>%
      mutate(variable_f = factor(variable, levels=c("Stock_biomass","Spawning_biomass","F","Catch_est"))) %>%
      mutate(label=sprintf("rho == %.2f",rho)) %>%
      mutate(Year=min(retro_table$Year),value=max_value$value)
    g4 = g4 +
      geom_text(data=rho_data, parse=TRUE, aes(label=label,hjust=0,vjust=1),size=6)

  }
  g4
}
