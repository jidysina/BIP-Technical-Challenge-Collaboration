#read csv file
df <- read.csv("dataset/XXHqn.csv")
library(dplyr)
#selecting important variables
analysis_df <- select(df,QN12, QN13, QN16, QN24, QN25, QN86, QN22, QN19, QN20, QN21, QN88, QN89, QN90, QN91, QN100, QN101, QN102, QN23, QN85, QN87, QN99, QN103, QN104, QN27, QN28, QN29, QN30, QN43, QN48, qntb4, qnillict, QN26, QN84,qnowt, qnobese, weight, stratum, psu)
#writing into csv file
write.csv(
      analysis_df,
     "dataset/analysis_dataset.csv",
      row.names = FALSE)
