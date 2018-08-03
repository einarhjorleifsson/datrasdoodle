# https://chrisvoncsefalvay.com/defensive-programming-r

source("R/DATRAS_stripped-down-fitALK.R")

cpue_per_length_per_haul <- function(hh, hl, Latin) {

  by.haul.positive <-
    hl %>%
    filter(latin == Latin) %>%
    mutate(length = floor(length)) %>%
    # Note: we are collapsing sex and maturity
    group_by(id, latin, length) %>%
    summarise(n = sum(n)) %>%
    drop_na()

  all <-
    hh %>%
    filter(haulval == "V") %>%
    # Lets only carry forward variables "needed"
    select(id, year, quarter) %>%
    crossing(length = sort(unique(by.haul.positive$length)),
             latin = unique(by.haul.positive$latin))

  by.haul <-
    all %>%
    left_join(by.haul.positive, by = c("id", "length", "latin")) %>%
    mutate(n = replace_na(n, 0))

  return(by.haul)

}


fit_alk <- function(ca2, lengths = 1:200, model = "mlogit") {

  data.mlogit <-
    ca2 %>%
    dplyr::select(lngtclass = length,
                  age,
                  n) %>%
    # NOTE: Need to check interpetation of the n-variable in ca-data
    tidyr::uncount(n) %>%
    mlogit::mlogit.data(varying = NULL, choice = 'age', shape = 'wide')

  m <- mlogit::mlogit(age~1 | lngtclass, data = data.mlogit, reflevel = "1")

  x <- coefficients(m)

  p <-
    dplyr::data_frame(variable = names(x),
                      value = x) %>%
    tidyr::separate(variable, c("age", "parameter"), sep = ":", convert = TRUE) %>%
    dplyr::mutate(parameter = stringr::str_replace(parameter, "\\(", ""),
                  parameter = stringr::str_replace(parameter, "\\)", "")) %>%
    tidyr::spread(parameter, value)

  ages = c(min(ca2$age):max(ca2$age))

  alk <-
    tidyr::crossing(age = ages,
                    length = lengths) %>%
    dplyr::left_join(p, by = "age") %>%
    # This is for the specified reflev in the mlogit function call
    dplyr::mutate(intercept = tidyr::replace_na(intercept, 0),
                  lngtclass = tidyr::replace_na(lngtclass, 0),
                  p = exp(intercept + lngtclass * length)) %>%
    # Probabilities of all ages within a length class
    #   must sum to one
    dplyr::group_by(length) %>%
    dplyr::mutate(p = p / sum(p)) %>%
    dplyr::ungroup() %>%
    dplyr::select(age, length, p) %>%
    dplyr::as_tibble()

  return(alk)
}


