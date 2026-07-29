# DATA Preparation
# 2023 National YRBS — Co-occurring Adverse/Protective Factors and Adolescent Suicidality
#
# This script does three things, in order:
#   1. Load the recoded (QN#) dataset
#   2. Subset to only the variables we'd be using
#   3. Label the 1/2 codes as readable Yes/No factors


# Loading the required libraries using pacman
pacman::p_load(survey, tidyverse, readr, haven) # loading all packages at once


# Loading the data 
# XXHqn.csv is the CDC-recoded file: QN# = dichotomized (1=Yes/2=No) versions of each question,
# plus the derived composites and the survey design variables (weight, stratum, psu)
df <- read.csv("XXHqn.csv")


nrow(df) # checking num of rows for qc: should be 20,103 rows per the Data User's Guide.


# Selecting the important variables
analysis_df <- select(
  df, q1, q2, q3, raceeth, QN12, QN13, QN16, QN19, QN20, QN21, QN22, QN23, QN24, 
  QN25, QN26, QN27, QN28, QN29, QN30, QN43, QN48, QN84, QN85, QN86, QN87, QN88, 
  QN89, QN90, QN91, QN99, QN100, QN101, QN102, QN103, QN104, qntb4, qnillict, 
  qnowt, qnobese, weight, stratum, psu
  )


# Creating a binary variables vector of yes/no variables, 
# except the four demographics and three design variables, which aren't binary items.
binary_vars <- c(
  "QN12","QN13","QN16","QN19","QN20","QN21","QN22","QN23","QN24",
  "QN25","QN26","QN27","QN28","QN29","QN30","QN43","QN48","QN84","QN85",
  "QN86","QN87","QN88","QN89","QN90","QN91","QN99","QN100","QN101","QN102",
  "QN103","QN104","qntb4","qnillict","qnowt","qnobese"
)

# Relabelling the 1/2 binary variables as Yes/No factors for readability.
# New labelled columns ("_f") are created rather than overwriting the originals.
for (v in binary_vars) {                    # loops over every variable name in the list
  analysis_df[[paste0(v, "_f")]] <- factor( # creates a NEW column named e.g. "QN27_f"
    analysis_df[[v]],                       # takes the values from the original column, e.g. analysis_df$QN27
    levels = c(1, 2),                       # tells R the two non-missing raw codes to expect
    labels = c("Yes", "No")                 # maps them positionally: 1 as "Yes", 2 as "No"
  )
}


# Writing into csv file
write.csv(
      analysis_df,
     "analysis_dataset.csv",
      row.names = FALSE)


# Creating Survey Design Object
# Building the survey design object on the new analysis dataframe (with the factor columns)
analysisdes <- svydesign(
  id      = ~psu, # primary sampling unit, which defines the clustering
  strata  = ~stratum, # stratification: race/ethicity and metro status area (rural/urban areeas)
  weights = ~weight, # corrects for nonresponse + oversampling
  nest    = TRUE, # PSU codes repeat across strata; tells R to treat them as nested, not globally unique
  data    = analysis_df
)

# Summary of survey design
summary(analysisdes)