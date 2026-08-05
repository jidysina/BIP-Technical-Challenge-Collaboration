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




# STEP 3 - Building a comparison function
#
# For any outcome + predictor pair, this function returns: weighted % in the
# "Yes" group, weighted % in the "No" group, the percentage-point gap,
# the unweighted N in each group, and a survey design-adjusted significance test.

# binary_vars currently mixes outcomes and predictors together. It was
# built that way so as to label the 1/2 as yes/no
# They'd be split here, since the comparison function needs to
# treat them differently. setdiff() keeps everything in binary_vars
# that is NOT one of the four outcome names.
outcome_vars   <- c("QN27", "QN28", "QN29", "QN30")
predictor_vars <- setdiff(binary_vars, outcome_vars)


# THE FUNCTION 
# Takes one outcome, one predictor, and the survey design, and returns a
# single row describing how strongly they're associated.
test_cooccurrence <- function(outcome_var, predictor_var, design) {
  
  # design$variables is the actual dataframe living inside the design
  # object. [[predictor_var]] pulls out one column from it by name --
  # this is needed (instead of $) because predictor_var is a text string
  # that changes every time the function is called, not a fixed column
  des_yes <- subset(design, design$variables[[predictor_var]] == 1)  # everyone who answered "Yes" to the predictor
  des_no  <- subset(design, design$variables[[predictor_var]] == 2)  # everyone who answered "No"
  
  # svyciprop needs a formula, not just a column name e.g. it needs
  # ~I(QN27 == 1), not just the text "QN27". Since outcome_var is a changing text string,
  # paste0() builds that exact piece of text first (e.g. "~I(QN27 == 1)"),
  # and as.formula() then converts that text into a real formula R can use.
  form <- as.formula(paste0("~I(", outcome_var, " == 1)"))
  
  pct_yes <- svyciprop(form, des_yes, na.rm = TRUE, method = "xlogit")  # weighted % with the outcome, Yes group
  pct_no  <- svyciprop(form, des_no,  na.rm = TRUE, method = "xlogit")  # weighted % with the outcome, No group
  
  n_yes <- as.numeric(unwtd.count(form, des_yes))  # how many real (non-missing) responses back up pct_yes
  n_no  <- as.numeric(unwtd.count(form, des_no))   # same, for pct_no
  
  # Same paste0()/as.formula() trick as above, but building a formula for
  # the chi-square test instead, i.e. "~QN27 + QN24".
  chisq_form <- as.formula(paste0("~", outcome_var, " + ", predictor_var))
  
  # The design-adjusted association test (Rao-Scott correction) -- tests
  # whether outcome and predictor are related, accounting for the survey's
  # weighting and clustering.
  test <- svychisq(chisq_form, design, statistic = "F")
  
  # Package everything computed above into one tidy row to return.
  data.frame(
    outcome       = outcome_var,                                          # which outcome this row is about
    predictor     = predictor_var,                                        # which predictor this row is about
    pct_given_yes = round(as.numeric(pct_yes) * 100, 1),                  # weighted % with outcome, among "Yes" group
    pct_given_no  = round(as.numeric(pct_no) * 100, 1),                   # weighted % with outcome, among "No" group
    diff_pct      = round((as.numeric(pct_yes) - as.numeric(pct_no)) * 100, 1), # the gap between the two groups
    n_yes_group   = n_yes,                                                # unweighted sample size, "Yes" group
    n_no_group    = n_no,                                                 # unweighted sample size, "No" group
    p_value       = signif(unname(test$p.value), 3)                      # unname() removes a label R attaches automatically
  )
}


# Testing the function
test_cooccurrence("QN27", "QN24", analysisdes)



# THE LOOP
# Using a for loop to run test_cooccurrence() once for every combination of outcome and
# predictor, and collects every result into one big table.

# Table size - one row per combination, so number of outcomes * number of predictors.
n_rows <- length(outcome_vars) * length(predictor_vars)

# Building an empty results table first with the right size, and every
# column set to the right data type. Every value would be overridden when the loop runs.
cooccurrence_results <- data.frame(
  outcome       = character(n_rows),
  predictor     = character(n_rows),
  pct_given_yes = numeric(n_rows),
  pct_given_no  = numeric(n_rows),
  diff_pct      = numeric(n_rows),
  n_yes_group   = numeric(n_rows),
  n_no_group    = numeric(n_rows),
  p_value       = numeric(n_rows)
)

# row_i keeps track of the row we're about to fill in next. It starts
# at row 1, and 1 is added to it every time a combination is finished.
row_i <- 1

# Outer loop: go through each of the 4 outcomes, one at a time.
for (o in outcome_vars) {
  
  # Inner loop: for the current outcome, go through every predictor, one at a time.
  for (p in predictor_vars) {
    
    # Running the function on this one specific outcome-predictor pair.
    one_result <- test_cooccurrence(o, p, analysisdes)
    
    # Copying the one-row result into row number row_i of the results table.
    # The comma after row_i means "for this row, every column".
    cooccurrence_results[row_i, ] <- one_result
    
    # Moving on to the next row, ready for the next pair.
    row_i <- row_i + 1
  }
}

# Display the finished table
view(cooccurrence_results)