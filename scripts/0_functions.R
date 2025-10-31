        #########################################
        #    GBLUP-based multi-kernel models    #   
        #               Functions               #
        #########################################
      
        #--------------------------------------------------------------------------------------------------------------------------------------------------------------------#
        # Function to build the sca matrices:
        
        #-------------------------------------------------------------#
        # crosses: data.frame with columns: hybrid, p1 (SS), p2 (NSS) #
        # add_matrix_sst: additive relationship matrix for SS lines   #
        # add_matrix_nst: additive relationship matrix for NSS lines  #
        #-------------------------------------------------------------#
        
        build_s_matrix = function(crosses, add_matrix_sst, add_matrix_nst) {
          
              stopifnot(all(crosses$p1 %in% rownames(add_matrix_sst)))
              stopifnot(all(crosses$p2 %in% rownames(add_matrix_nst)))
          
              n_hybrids = nrow(crosses)
              s_mat = matrix(NA, nrow = n_hybrids, ncol = n_hybrids)
              rownames(s_mat) = crosses$hybrid
              colnames(s_mat) = crosses$hybrid
          
              for (i in 1:n_hybrids) {
              
                   for (j in 1:n_hybrids) {
              
                        sst_i = as.character(crosses$p1[i])  
                        sst_j = as.character(crosses$p1[j])  
                        nst_i = as.character(crosses$p2[i])  
                        nst_j = as.character(crosses$p2[j])  
                       
                        s_mat[i, j] = add_matrix_sst[sst_i, sst_j] * add_matrix_nst[nst_i, nst_j]
            
                    }
              }
          
              return(s_mat)
        
        }
        
        #--------------------------------------------------------------------------------------------------------------------------------------------------------------------#
        # Function to get the selected baseline models incorporating different variance structures:
        
        #---------------------------------------------------------------------------#
        # The argument "type" defines the best-fitting model structure per trait    #
        #                                                                           #
        #   type = 1 → for P.H (Plant Height)                                       #
        #              → uses the GEI = identity structure                          #
        #                                                                           #
        #   type = 2 → for G.Y, E.H, and S.I (Grain Yield, Ear Height, Silking)     #
        #              → uses the GEI = identity + diagonal structure               #
        #                                                                           #
        #   type = 3 → for A.N (Anthesis)                                           #
        #              → uses the GEI = diagonal structure                          #
        #                                                                           #
        # Two relationship matrices are tested for each trait (chosen_matrix):      #
        #   - "D" model → uses the dominance matrix (dom_inverse)                   #
        #   - "S" model → uses the SCA matrix (sca_inverse)                         #
        #---------------------------------------------------------------------------#
        
        fit_design2_models = function(trait, data, chosen_matrix, type) {
        
              fixed_effect = paste0("PV_", trait, " ~ Year_Location")
              w_var = data[[paste0("Weights_", trait)]]
              
              if (type == 1) {
                
              random_effect = "~ vm(Pedigree1, add_inverse_sst) + 
                                 vm(Pedigree2, add_inverse_nst) + 
                                 vm(Pedigree, chosen_matrix) + 
                                 Year_Location:Pedigree1 + 
                                 Year_Location:Pedigree2 + 
                                 Year_Location:Pedigree1:Pedigree2"
              
              } else if (type == 2) {
              
              random_effect = "~ vm(Pedigree1, add_inverse_sst) + 
                                 vm(Pedigree2, add_inverse_nst) + 
                                 vm(Pedigree, chosen_matrix) + 
                                 Year_Location:Pedigree1 + 
                                 Year_Location:Pedigree2 + 
                                 idh(Year_Location):Pedigree1:Pedigree2"
              
              } else if (type == 3) {
                
              random_effect = "~ vm(Pedigree1, add_inverse_sst) + 
                                 vm(Pedigree2, add_inverse_nst) + 
                                 vm(Pedigree, chosen_matrix) + 
                                 idh(Year_Location):Pedigree1 + 
                                 idh(Year_Location):Pedigree2 + 
                                 idh(Year_Location):Pedigree1:Pedigree2"
              
              }
              
              model = asreml(fixed = as.formula(fixed_effect),
                             random = as.formula(random_effect),
                             data = data,
                             weights = w_var,
                             family = asr_gaussian(dispersion = 1),
                             na.action = na.method(x = "include", y = "include"), 
                             workspace = 6e8,
                             maxit = 100)
              
              # Coefficients - Blues:
              blues_coef = as.data.frame(summary(model, coef = T)$coef.fixed)
              intercept = blues_coef["(Intercept)", "solution"]
              blues_rows = grepl("Year_Location_", rownames(blues_coef))
              mean_year_loc = mean(blues_coef[blues_rows, "solution"], na.rm = TRUE)
                  
              # Coefficients - Blups:
              blups_coef = as.data.frame(summary(model, coef = T)$coef.random)
                   
              # GCA - Parent 1 (Female):
              bp_coef_p1 = blups_coef[grepl("^(Pedigree1_|vm\\(Pedigree1, .*\\)_)[^:]+$", rownames(blups_coef)), ]
              bp_coef_p1$Pedigree1 = gsub(".*_", "", rownames(bp_coef_p1))
              rownames(bp_coef_p1) = NULL
                   
              # GCA - Parent 2 (Male): 
              bp_coef_p2 = blups_coef[grepl("^(Pedigree2_|vm\\(Pedigree2, .*\\)_)[^:]+$", rownames(blups_coef)), ]
              bp_coef_p2$Pedigree2 = gsub(".*_", "", rownames(bp_coef_p2))
              rownames(bp_coef_p2) = NULL
                   
              # SCA - Parent 1 x Parent 2:
              bp_coef_p3 = blups_coef[grepl("^vm\\(Pedigree, chosen_matrix\\)_[^/]+/[^/]+$", rownames(blups_coef)) | 
                                      grepl("^Pedigree1_[^:]+:Pedigree2_[^:]+$", rownames(blups_coef)), ]
              
              bp_coef_p3$Pedigree = ifelse(grepl("^Pedigree1_", rownames(bp_coef_p3)), 
                                           gsub("^Pedigree1_([^:]+):Pedigree2_(.+)$", "\\1/\\2", rownames(bp_coef_p3)),
                                           gsub(".*_", "", rownames(bp_coef_p3)))
              rownames(bp_coef_p3) = NULL
                   
              # Interactions:
              if (any(grepl(".*Pedigree1.*", rownames(blups_coef)) & grepl("Year_Location", rownames(blups_coef)))) {
                     
                       inter_1 = blups_coef[grepl("Year_Location.*Pedigree1", rownames(blups_coef)) & !grepl("Pedigree2", rownames(blups_coef)),]
                       inter_1 = inter_1 %>% tibble::rownames_to_column("Term") %>% 
                                 tidyr::separate(Term, into = c("Year_Location", "Pedigree1"), 
                                 sep = ":Pedigree1_") %>% 
                                 dplyr::group_by(Pedigree1) %>%
                                 dplyr::summarise(GxE_1 = mean(solution, na.rm = TRUE))

              } else {
                          
                       inter_1 = tibble::tibble(Pedigree1 = bp_coef_p1$Pedigree1, GxE_1 = 0)
                   
              }
                   
              if (any(grepl(".*Pedigree2.*", rownames(blups_coef)) & grepl("Year_Location", rownames(blups_coef)))) {
                     
                       inter_2 = blups_coef[grepl("Year_Location.*Pedigree2", rownames(blups_coef)) & !grepl("Pedigree1", rownames(blups_coef)),]
                       inter_2 = inter_2 %>% tibble::rownames_to_column("Term") %>% 
                                 tidyr::separate(Term, into = c("Year_Location", "Pedigree2"), 
                                 sep = ":Pedigree2_") %>% 
                                 dplyr::group_by(Pedigree2) %>%
                                 dplyr::summarise(GxE_2 = mean(solution, na.rm = TRUE))

              } else {
                     
                       inter_2 = tibble::tibble(Pedigree2 = bp_coef_p2$Pedigree2, GxE_2 = 0)
                     
              }
                   
              if (any(grepl("^Year_Location_.*:Pedigree1_.*:Pedigree2_.*$", rownames(blups_coef)))) {
                     
                       inter_3 = blups_coef[grepl("Year_Location.*[:_]Pedigree1.*[:_]Pedigree2.*", rownames(blups_coef)), ]
                       inter_3 = inter_3 %>%
                                 tibble::rownames_to_column("Term") %>%
                                 tidyr::separate(Term, into = c("Year_Location", "Pedigree1", "Pedigree2"), 
                                 sep = ":", remove = FALSE) %>%
                                 dplyr::mutate(
                                 Pedigree1 = gsub("Pedigree1_", "", Pedigree1),
                                 Pedigree2 = gsub("Pedigree2_", "", Pedigree2),
                                 Pedigree = paste0(Pedigree1, "/", Pedigree2)) %>%
                                 dplyr::group_by(Pedigree) %>%
                                 dplyr::summarise(GxE_3 = mean(solution, na.rm = TRUE))

              } else {
          
                       inter_3 = tibble::tibble(Pedigree = bp_coef_p3$Pedigree, GxE_3 = 0)
          
              }
        
              # Parent 1 (Female):
              bp_coef_p1 = bp_coef_p1 %>% mutate(Intercept = intercept, Mean_Env = mean_year_loc)
              bp_coef_p1$PV_Coef = bp_coef_p1[, "solution"] + bp_coef_p1[, "Intercept"] + bp_coef_p1[, "Mean_Env"]
              bp_coef_p1 = bp_coef_p1[, c(4, 5, 6, 1, 2, 3, 7)] %>% dplyr::arrange(desc(PV_Coef)) %>% dplyr::mutate(Rank_Coef = row_number())
              names(bp_coef_p1)[4:6] = c("GCA_1", "Std_Error_Coef", "Z_Ratio")
              gca_1 = bp_coef_p1
              
              # Parent 2 (Male): 
              bp_coef_p2 = bp_coef_p2 %>% mutate(Intercept = intercept, Mean_Env = mean_year_loc)
              bp_coef_p2$PV_Coef = bp_coef_p2[, "solution"] + bp_coef_p2[, "Intercept"] + bp_coef_p2[, "Mean_Env"]
              bp_coef_p2 = bp_coef_p2[, c(4, 5, 6, 1, 2, 3, 7)] %>% dplyr::arrange(desc(PV_Coef)) %>% dplyr::mutate(Rank_Coef = row_number())
              names(bp_coef_p2)[4:6] = c("GCA_2", "Std_Error_Coef", "Z_Ratio")
              gca_2 = bp_coef_p2
              
              # Hibryds:
              bp_coef_p3 = bp_coef_p3 %>% tidyr::separate(Pedigree, into = c("Pedigree1", "Pedigree2"), sep = "/", remove = FALSE)
              sca_h = bp_coef_p3[, c(4, 5, 6, 1, 2, 3)] %>% dplyr::arrange(desc(solution)) %>% dplyr::mutate(Rank_Coef = row_number())
              names(sca_h)[4:6] = c("SCA", "Std_Error_Coef", "Z_Ratio")
              bp_coef_p3 = left_join(bp_coef_p3, inter_1, by = "Pedigree1")
              bp_coef_p3 = left_join(bp_coef_p3, inter_2, by = "Pedigree2")
              bp_coef_p3 = left_join(bp_coef_p3, inter_3, by = "Pedigree")
              bp_coef_p3 = left_join(bp_coef_p3, bp_coef_p1[, c(1,4)], by = "Pedigree1")
              bp_coef_p3 = left_join(bp_coef_p3, bp_coef_p2[, c(1,4)], by = "Pedigree2")
              bp_coef_p3 = bp_coef_p3 %>% mutate(Intercept = intercept, Mean_Env = mean_year_loc)
              
              bp_coef_p3$PV_Coef_1 = bp_coef_p3[, "solution"] + bp_coef_p3[, "Intercept"] + bp_coef_p3[, "Mean_Env"] + 
                                     bp_coef_p3[, "GxE_1"] + bp_coef_p3[, "GxE_2"] + bp_coef_p3[, "GxE_3"] + 
                                     bp_coef_p3[, "GCA_1"] + bp_coef_p3[, "GCA_2"] 
              
              bp_coef_p3$PV_Coef_2 = bp_coef_p3[, "solution"] + bp_coef_p3[, "Intercept"] + bp_coef_p3[, "Mean_Env"] + 
                                     bp_coef_p3[, "GxE_1"] + bp_coef_p3[, "GxE_2"] + bp_coef_p3[, "GxE_3"] 
              
              bp_coef_p3 = (bp_coef_p3 %>% dplyr::mutate(
                            Rank_Coef_1 = rank(-PV_Coef_1, ties.method = "first"),
                            Rank_Coef_2 = rank(-PV_Coef_2, ties.method = "first")) %>%
                            dplyr::arrange(desc(PV_Coef_1))
                            )[, c(4, 5, 6, 12, 13, 7, 8, 9, 10, 11, 1, 2, 3, 14, 16, 15, 17)]
              
              names(bp_coef_p3)[11:13] = c("SCA", "Std_Error_Coef", "Z_Ratio")
              final_results = bp_coef_p3
                
              return(list(model = model, gca_m = gca_1, gca_p = gca_2, sca_h = sca_h, final_results = final_results))

        }
        
        #--------------------------------------------------------------------------------------------------------------------------------------------------------------------#
        # Function to computation of genetic covariance matrices for hybrid prediction:
        
        #----------------------------------------------------------------------------------------------#
        # Args:                                                                                        #
        #   i: name of the untested hybrid.                                                            #
        #   data: training dataset used to fit the model (contains tested hybrids and their genitors)  #
        #   model: ASReml model fitted with either the D or S structure                                #
        #   f_mat: additive relationship matrix for seed parent (SS) parental lines                    #
        #   m_mat: additive relationship matrix for pollen parent (NSS) parental lines                 #
        #   h_mat: relationship matrix for hybrids                                                     #
        #          (either dominance matrix D or SCA matrix S, depending on the model fitted)          #
        #   pheno: data frame with BLUEs or predicted phenotypic values across environments            #
        #   trait: trait name (character string, e.g., "G.Y", "E.H", "S.I", or "A.N")                  #
        #----------------------------------------------------------------------------------------------#
        
        get_covariance_matrices = function(i, data, model, f_mat, m_mat, h_mat, pheno, trait) {
          
              genitors_i = unlist(strsplit(i, "/"))
              female_i = genitors_i[1]
              male_i = genitors_i[2]
          
              hybrids_tt = levels(data$Pedigree)
              n_hybrids_tt = length(hybrids_tt)
              c_ut_m3 = numeric(n_hybrids_tt) 
              c_tt_m3 = matrix(0, n_hybrids_tt, n_hybrids_tt)
              c_ut_m4 = numeric(n_hybrids_tt) 
              c_tt_m4 = matrix(0, n_hybrids_tt, n_hybrids_tt)
             
              names(c_ut_m3) = hybrids_tt
              names(c_ut_m4) = hybrids_tt
             
              rownames(c_tt_m3) = hybrids_tt
              colnames(c_tt_m3) = hybrids_tt
              rownames(c_tt_m4) = hybrids_tt
              colnames(c_tt_m4) = hybrids_tt

              comp = summary(model)$varcomp
              sigma2_gca_f = comp["vm(Pedigree1, add_inverse_sst)", "component"]
              sigma2_gca_m = comp["vm(Pedigree2, add_inverse_nst)", "component"]
              sigma2_sca_h = comp["vm(Pedigree, chosen_matrix)", "component"]
          
              for (k in hybrids_tt) {
                
                   genitors_k = unlist(strsplit(k, "/"))
                   female_k = genitors_k[1]
                   male_k = genitors_k[2]
            
                   # Genetic covariance matrix of untested and tested hybrids:
                   gf_ik = ifelse(female_i %in% rownames(f_mat) & female_k %in% colnames(f_mat), f_mat[female_i, female_k], 0) 
                   gm_ik = ifelse(male_i %in% rownames(m_mat) & male_k %in% colnames(m_mat), m_mat[male_i, male_k], 0)
                   dh_ik = ifelse(i %in% rownames(h_mat) & k %in% colnames(h_mat), h_mat[i, k], 0)
                   c_ut_m3[k] = (gf_ik * sigma2_gca_f) + (gm_ik * sigma2_gca_m)
                   c_ut_m4[k] = (gf_ik * sigma2_gca_f) + (gm_ik * sigma2_gca_m) + (dh_ik * sigma2_sca_h)
                   
                   # Genetic covariance matrix of tested hybrids:
                   for (l in hybrids_tt) {
                    
                        genitors_l = unlist(strsplit(l, "/"))
                        female_l = genitors_l[1]
                        male_l = genitors_l[2]
                       
                        gf_kl = ifelse(female_k %in% rownames(f_mat) & female_l %in% colnames(f_mat), f_mat[female_k, female_l], 0)
                        gm_kl = ifelse(male_k %in% rownames(m_mat) & male_l %in% colnames(m_mat), m_mat[male_k, male_l], 0)
                        dh_kl = ifelse(k %in% rownames(h_mat) & l %in% colnames(h_mat), h_mat[k, l], 0)
                       
                        if (female_k == female_l && male_k == male_l) {
                         
                        st_err = (pheno %>% filter(Pedigree %in% l) %>% pull(paste0("Std_Error_", trait)))^2
                        c_tt_m3[k, l] = (gf_kl * sigma2_gca_f) + (gm_kl * sigma2_gca_m) + st_err
                        c_tt_m4[k, l] = (gf_kl * sigma2_gca_f) + (gm_kl * sigma2_gca_m) + (dh_kl * sigma2_sca_h) + st_err

                        } else {
                         
                        c_tt_m3[k, l] = (gf_kl * sigma2_gca_f) + (gm_kl * sigma2_gca_m)
                        c_tt_m4[k, l] = (gf_kl * sigma2_gca_f) + (gm_kl * sigma2_gca_m) + (dh_kl * sigma2_sca_h)
                       
                       }
                  
                   }
                   
              }
              
              y_t = pheno %>% 
                    filter(Pedigree %in% hybrids_tt) %>% 
                    mutate(Pedigree = factor(Pedigree, levels = hybrids_tt)) %>% 
                    arrange(Pedigree) %>% 
                    pull(paste0("Predicted_Values_", trait))

              y_hat_3 = as.numeric(c_ut_m3 %*% solve(c_tt_m3) %*% y_t)
              y_hat_4 = as.numeric(c_ut_m4 %*% solve(c_tt_m4) %*% y_t)
                   
              return(list(c_ut_m3 = c_ut_m3, c_ut_m4 = c_ut_m4, c_tt_m3 = c_tt_m3, c_tt_m4 = c_tt_m4, y_hat_3 = y_hat_3, y_hat_4 = y_hat_4))
        
        }
        
        #--------------------------------------------------------------------------------------------------------------------------------------------------------------------#
        
        

        