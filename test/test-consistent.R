source("common_functions.R")

gs_invalid <- load_gs("gs/LoveLab-InvalidNode")
consistency(gs_invalid)

gs_valid <- load_gs("gs/LoveLab-Debug")
consistency(gs_valid)
