        ##################################################
        #         GBLUP-based multi-kernel models        #   
        #   Genomic prediction of single-cross hybrids   #
        ##################################################
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Removing variables:
        rm(list=ls())
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Directory, packages, and functions:
        setwd("~/Documents/tag-maize-ss-nss-gblup-multikernel/")
        source("scripts/0_functions.R")
        library("asreml")
        library("dplyr")
        library("tibble")
        library("magrittr")

        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Getting the blues:
        blues = readRDS("data/single_env_analyses/Final_blues_all_traits.rds")
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Getting the additive and dominance relationship matrices:
        add_sst = readRDS("data/genomic_matrices/all_groups/Add_sst.rds")
        add_nst = readRDS("data/genomic_matrices/all_groups/Add_nst.rds")
        dom = readRDS("data/genomic_matrices/all_groups/Dom_matrix.rds")
        sca = readRDS("data/genomic_matrices/all_groups/Sca_matrix.rds")
        
        add_inverse_sst = readRDS("data/genomic_matrices/all_groups/Add_inverse_sst.rds")
        add_inverse_nst = readRDS("data/genomic_matrices/all_groups/Add_inverse_nst.rds")
        dom_inverse = readRDS("data/genomic_matrices/all_groups/Dom_inverse.rds")
        sca_inverse = readRDS("data/genomic_matrices/all_groups/Sca_inverse.rds")
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Getting the BLUEs from the multi-environment trial analyses:
        pedigree = unique(blues$Pedigree)
        blues_pheno = readRDS("data/multi_env_analyses/Blues.rds")
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Defining the common maximum possible training set size across the four cv scenarios:
        hybrids = levels(blues$Pedigree)
        train_sizes = data.frame(Pedigree = character(), T2 = numeric(), T1F = numeric(), T1M = numeric(), T0 = numeric())
        
        for (i in hybrids) {
          
             print(i)
          
             genitors_i = unlist(strsplit(i, "/"))
             female_i = genitors_i[1]
             male_i = genitors_i[2]
          
             n_t2  = length(unique(blues$Pedigree[blues$Pedigree != i]))
             n_t1f = length(unique(blues$Pedigree[blues$Pedigree != i & blues$Pedigree1 != female_i]))
             n_t1m = length(unique(blues$Pedigree[blues$Pedigree != i & blues$Pedigree2 != male_i]))
             n_t0  = length(unique(blues$Pedigree[blues$Pedigree != i & blues$Pedigree1 != female_i & blues$Pedigree2 != male_i]))
             
             train_sizes = rbind(train_sizes, data.frame(Pedigree = i, T2 = n_t2, T1F = n_t1f, T1M = n_t1m, T0 = n_t0))
        }
        
        max_common_train_size = min(train_sizes[, c("T2", "T1F", "T1M", "T0")])
        cat("The common maximum possible training set size across the four scenarios:", max_common_train_size, "\n")
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Cross-validation and prediction accuracy estimation:
        results = list()       
        n_rep = 30
        trait = "G.Y" # define the trait to be analyzed (run the script separately for each trait)

        for (rep in 1:n_rep) {
             
             cat("Repetition:", rep, "\n")
        
             for (i in hybrids) {
                   
                   cat(sprintf("Repetition: %d | Hybrid: %s\n", rep, i))
                
                   genitors_i = unlist(strsplit(i, "/"))
                   female_i = genitors_i[1]
                   male_i = genitors_i[2]
                
                   # Training set:
                   tra_set = list(tra_set_t2 = blues[blues$Pedigree != i, ],
                                  tra_set_t1f = blues[blues$Pedigree != i & blues$Pedigree1 != female_i, ],
                                  tra_set_t1m = blues[blues$Pedigree != i & blues$Pedigree2 != male_i, ],
                                  tra_set_t0 = blues[blues$Pedigree != i & blues$Pedigree1 != female_i & blues$Pedigree2 != male_i, ])
                   
                   # Fitting the gs models to get the genomic estimated breeding values (gebv):
                   
                   for (j in names(tra_set)) {
                     
                   cat(sprintf("Repetition: %d | Hybrid: %s | Cross Validation Scenario: %s\n", rep, i, j))
                     
                   sampled_peds = sample(unique(tra_set[[j]]$Pedigree), 125) %>% droplevels()
                   tra_data = tra_set[[j]] %>% filter(Pedigree %in% sampled_peds) %>% droplevels()

                   model_1 = fit_design2_models(trait, tra_data, dom_inverse, type = 2) # model with D
                   model_2 = fit_design2_models(trait, tra_data, sca_inverse, type = 2) # model with S
                   
                   #-------------------------------------------------------------------------#
                   # The argument "type" defines the best-fitting model structure per trait: #
                   #   type = 1 → for P.H                                                    #
                   #   type = 2 → G.Y, E.H, and S.I                                          #
                   #   type = 3 → A.N                                                        #
                   # ------------------------------------------------------------------------#
                   
                   effects_1 = model_1$final_results
                   effects_2 = model_2$final_results
                   
                   # Getting the gebv of each untested hybrids:
                   gebv_u_1 = effects_1 %>% filter(Pedigree == i) %>% droplevels()
                   gebv_u_2 = effects_2 %>% filter(Pedigree == i) %>% droplevels()
                   
                   # Using different methods to predict single-cross performance:
                   
                   # Method 1:
                   y_hat_1_a = gebv_u_1[, "Intercept"] + gebv_u_1[, "GCA_1"] + gebv_u_1[, "GCA_2"]
                   y_hat_1_b = gebv_u_2[, "Intercept"] + gebv_u_2[, "GCA_1"] + gebv_u_2[, "GCA_2"]
                   
                   # Method 2:
                   y_hat_2_a = gebv_u_1[, "Intercept"] + gebv_u_1[, "GCA_1"] + gebv_u_1[, "GCA_2"] + gebv_u_1[, "SCA"]
                   y_hat_2_b = gebv_u_2[, "Intercept"] + gebv_u_2[, "GCA_1"] + gebv_u_2[, "GCA_2"] + gebv_u_2[, "SCA"]
                   
                   # Method 3 and 4:
                   y_all_1 = get_covariance_matrices(i, tra_data, model_1$model, add_sst, add_nst, dom, blues_pheno, trait)
                   y_hat_3_a = y_all_1$y_hat_3
                   y_hat_4_a = y_all_1$y_hat_4
                   
                   y_all_2 = get_covariance_matrices(i, tra_data, model_2$model, add_sst, add_nst, sca, blues_pheno, trait)
                   y_hat_3_b = y_all_2$y_hat_3
                   y_hat_4_b = y_all_2$y_hat_4

                   # Organizing the results:
                   results_dom = tibble(Method = c("1a", "2a", "3a", "4a"), Trait = c(y_hat_1_a, y_hat_2_a, y_hat_3_a, y_hat_4_a))
                   results_sca = tibble(Method = c("1b", "2b", "3b", "4b"), Trait = c(y_hat_1_b, y_hat_2_b, y_hat_3_b, y_hat_4_b))
                   results_all = dplyr::bind_rows(results_dom, results_sca)
                   colnames(results_all)[2] = trait 
                   
                   final_results = data.frame(Pedigree = i, C_V_Scenario = j, Rep = rep, results_all)
                   results = dplyr::bind_rows(results, final_results)
                   
                   rm(model_1, effects_1, gebv_u_1, y_all_1, y_hat_1_a, y_hat_2_a, y_hat_3_a, y_hat_4_a)
                   rm(model_2, effects_2, gebv_u_2, y_all_2, y_hat_1_b, y_hat_2_b, y_hat_3_b, y_hat_4_b)
                   gc()
                   
                   }
                   
              }
          
              # Saving results:
              saveRDS(results, file = paste0("", rep, "_", trait, ".rds")) # users should modify it according to their local directory

              saveRDS(results, file = paste0("results/Final_results_", rep, "_", trait, ".rds"))
              

              cat(sprintf("Results for repetition %d saved.\n", rep))
              
        }
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        
        
        