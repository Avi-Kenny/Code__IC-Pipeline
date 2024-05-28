# Processing CSVs
if (F) {
  
  library(magrittr)
  library(dplyr)
  
  tids <- c(1:3,13:15,25:27,37:39,45:47,49:51,59:64)
  folder <- paste0("../Figures + Tables/Janssen (partA) plots/Run 19 (added UL",
                   "OQ truncation)/CSVs")
  
  for (i in tids) {
    
    path_cve <- paste0(folder, "/cve_", i, ".csv")
    path_risk <- paste0(folder, "/risk_", i, ".csv")
    cve <- utils::read.csv(path_cve)
    risk <- utils::read.csv(path_risk)
    cve$which <- "cve"
    risk$which <- "risk"
    dat <- rbind(cve,risk)
    dat %<>% dplyr::filter(!is.na(y))
    dat %<>% dplyr::filter(x!=999)
    dat$i <- i
    dat %<>% dplyr::mutate(
      dataset = dplyr::case_when(
        i %in% c(1:3,45:47,59:61) ~ "Pooled",
        i %in% c(13:15) ~ "NA",
        i %in% c(25:27,49:51,62:64) ~ "LA",
        i %in% c(37:39) ~ "SA",
        TRUE ~ "Error"
      ),
      endpoint = dplyr::case_when(
        i %in% c(1:39) ~ "Moderate to Severe",
        i %in% c(45:51) ~ "Severe",
        i %in% c(59:64) ~ "Moderate",
        TRUE ~ "Error"
      ),
      marker = dplyr::case_when(
        i %in% c(1,13,25,37,45,49,59,62) ~ "Day29bindSpike",
        i %in% c(2,14,26,38,46,50,60,63) ~ "Day29bindRBD",
        i %in% c(3,15,27,39,47,51,61,64) ~ "Day29pseudoneutid50",
        TRUE ~ "Error"
      )
    )
    
    if (i==1) {
      dat_long <- dat
    } else {
      dat_long <- rbind(dat_long,dat)
    }
    
  }
  
  dat_long %<>% dplyr::relocate(dataset, .before=x)
  dat_long %<>% dplyr::relocate(marker, .before=x)
  dat_long %<>% dplyr::relocate(endpoint, .before=x)
  dat_long %<>% dplyr::relocate(which, .before=x)
  dat_long$overall <- NULL
  dat_long$i <- NULL
  write.table(dat_long, file="ENSEMBLE_partA.csv", sep=",", row.names=FALSE)
  
}