
### Ulility Functions in messir ### ----

#' output a table of samuika results
#' @param model \code{samuika} object
#' @encoding UTF-8
#' @export
out_summary_estimate <- function(model) {
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
    mutate(Cz = exp(qnorm(0.9)*sqrt(log(1+CV^2)))) %>%
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
