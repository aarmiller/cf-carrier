# cf-carrier
 scripts for cf-carrier study

The following files contain the pseudo code and R code that can be used to duplicate the findings of our 
paper "Carriers of cystic fibrosis". The Truven Marketscan Research databases can be purchased from IBM 
Watson Health. These databases are delivered as a series of SAS data files. At the University of Iowa,
these data are then converted to a relational database which is stored on a local drive. Because this 
data and computing structure is specific to the University or Iowa, we have provided the psuedo code 
(i.e., extraction procedures) that would need to be used to extract the same datasets from the raw Truven 
Marketscan files. We have noted what data the extracted files should contain.


Below is a brief summary of the file structure. The R scripts are located in two primary folders. The 
scripts folder contains the scripts for 

The scripts folder contains the following scripts:

  build_matched_cohorts.R - This script contains the procedures for identifying cohorts of carriers or cases of CF
                            along with the corresponding matched cohorts. This code specifically identifies the 
                            enrollee id's to extract. After running a this script you would then need to extract 
                            the data for these enrollees from the Truven database. Note: on our system we extract the
                            enrollee id's then store all the extracted data in a relational (SQLite) database.
                          
  main_analysis.R - This script contains the code for running the primary analysis in the study. Note: this script
                    uses the database of matched cohorts created from the enrollee id's identified in the 
                    build_mathced_cohorts.R script. This script also relies on the condition codesidentified in the 
                    condition_codes.R script and utilizes function from the analysis_functions.R and the 
                    small_db_functions.R scripts.
                    
  figures.R - This script contains all of the code that was used to generate the figures for the paper. This script
              utilizes the data containing odds ratio estimates that were generated in the main_analysis.R script. 
              This script also utilizes functions contained in the analysis_functions.R script.

  simulation.R - This script condatins the procedures and code that would need to be used to run the simulation analysis
                 described in the appendix of the paper.

The functions folder contains the following function scripts that were used for our analysis:

  small_db_functions.R - this file contains functions for working with the cohort (SQLite) databases that were generated
                         using the enrolids identified in the build_matched_cohorts.R script. These functions are used to
                         pull a longitudinal (timemap) of visits for all enrollee, and for returning the visits with
                         specific diagnosis codes.
                         
  analysis_functions.R - this file contains the primary functions that were utilized for the analysis in this study.
  
  condition_codes.R - this file contains the function for extracting all of the icd9 and icd10 codes that were used for
                      identifying the CF-related conditions.



