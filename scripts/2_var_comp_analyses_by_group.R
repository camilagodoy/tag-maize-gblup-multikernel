        #################################################
        #         BLUP-based multi-kernel models        #   
        #   Estimation of genetic variance components   #
        #################################################
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Removing variables:
        rm(list=ls()) 
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Directory, packages, and functions:
        setwd("~/Documents/tag-maize-ss-nss-gblup-multikernel/")
        library("tidyverse")
        library("asreml")

        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Getting the blues:
        blues_early = readRDS("data/single_env_analyses/Final_blues_early_all_traits.rds")
        blues_intermediate = readRDS("data/single_env_analyses/Final_blues_intermediate_all_traits.rds")
        blues_late = readRDS("data/single_env_analyses/Final_blues_late_all_traits.rds")
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Getting the additive and dominance relationship matrices:
        add_inverse_e_sst = readRDS("data/genomic_matrices/early/Add_inverse_er_sst.rds")
        add_inverse_e_nst = readRDS("data/genomic_matrices/early/Add_inverse_er_nst.rds")
        add_inverse_i_sst = readRDS("data/genomic_matrices/intermediate/Add_inverse_in_sst.rds")
        add_inverse_i_nst = readRDS("data/genomic_matrices/intermediate/Add_inverse_in_nst.rds")
        add_inverse_l_sst = readRDS("data/genomic_matrices/late/Add_inverse_la_sst.rds")
        add_inverse_l_nst = readRDS("data/genomic_matrices/late/Add_inverse_la_nst.rds")
        dom_inverse_e = readRDS("data/genomic_matrices/early/Dom_inverse_er.rds")
        dom_inverse_i = readRDS("data/genomic_matrices/intermediate/Dom_inverse_in.rds")
        dom_inverse_l = readRDS("data/genomic_matrices/late/Dom_inverse_la.rds")
        sca_inverse_e = readRDS("data/genomic_matrices/early/Sca_inverse_er.rds")
        sca_inverse_i = readRDS("data/genomic_matrices/intermediate/Sca_inverse_in.rds")
        sca_inverse_l = readRDS("data/genomic_matrices/late/Sca_inverse_la.rds")
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Removing invalid environments:        
        envs_drop_e = blues_early %>% 
                      dplyr::group_by(Year_Location) %>% 
                      dplyr::summarize(n_F = n_distinct(Pedigree1), 
                      n_M = n_distinct(Pedigree2), 
                      n_H = n_distinct(Pedigree)) %>% 
                      filter(n_F == 1 | n_M == 1)
        
        envs_drop_i = blues_intermediate %>% 
                      dplyr::group_by(Year_Location) %>% 
                      dplyr::summarize(n_F = n_distinct(Pedigree1), 
                      n_M = n_distinct(Pedigree2), 
                      n_H = n_distinct(Pedigree)) %>% 
                      filter(n_F == 1 | n_M == 1)
         
        envs_drop_l = blues_late %>% 
                      dplyr::group_by(Year_Location) %>% 
                      dplyr::summarize(n_F = n_distinct(Pedigree1), 
                      n_M = n_distinct(Pedigree2), 
                      n_H = n_distinct(Pedigree)) %>% 
                      filter(n_F == 1 | n_M == 1)
        
        blues_early = blues_early %>% filter(!Year_Location %in% envs_drop_e$Year_Location) %>% droplevels()
        blues_intermediate = blues_intermediate %>% filter(!Year_Location %in% envs_drop_i$Year_Location) %>% droplevels()
        blues_late = blues_late %>% filter(!Year_Location %in% envs_drop_l$Year_Location) %>% droplevels()
        
        nlevels(blues_early$Year_Location) # 18
        nlevels(blues_intermediate$Year_Location) # 36
        nlevels(blues_late$Year_Location) # 28
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Fitting the models by maturity group:
        
        # Early:
        
        # Grain yield:
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_gy_earl_d = asreml(fixed = PV_G.Y ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_e_sst) + 
                                         vm(Pedigree2, add_inverse_e_nst) + 
                                         vm(Pedigree, dom_inverse_e) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +                                
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_early,
                                         weights = Weights_G.Y,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_gy_earl_s = asreml(fixed = PV_G.Y ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_e_sst) + 
                                         vm(Pedigree2, add_inverse_e_nst) + 
                                         vm(Pedigree, sca_inverse_e) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +                                
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_early,
                                         weights = Weights_G.Y,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        model_gy_earl_s = update(model_gy_earl_s)
        model_gy_earl_s = update(model_gy_earl_s)
        
        # Plant height:
        
        # Considering GEI (identity) and assuming additive and dominance matrices for lines and hybrids:     
        model_ph_earl_d = asreml(fixed = PV_P.H ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_e_sst) + 
                                         vm(Pedigree2, add_inverse_e_nst) + 
                                         vm(Pedigree, dom_inverse_e) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         Year_Location:Pedigree1:Pedigree2,
                                         data = blues_early,
                                         weights = Weights_P.H,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        # Considering GEI (identity) and assuming additive and dominance matrices for lines and hybrids:     
        model_ph_earl_s = asreml(fixed = PV_P.H ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_e_sst) + 
                                         vm(Pedigree2, add_inverse_e_nst) + 
                                         vm(Pedigree, sca_inverse_e) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         Year_Location:Pedigree1:Pedigree2,
                                         data = blues_early,
                                         weights = Weights_P.H,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
     
        # Ear height:
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_eh_earl_d = asreml(fixed = PV_E.H ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_e_sst) + 
                                         vm(Pedigree2, add_inverse_e_nst) + 
                                         vm(Pedigree, dom_inverse_e) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_early,
                                         weights = Weights_E.H,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        model_eh_earl_d = update(model_eh_earl_d)
        model_eh_earl_d = update(model_eh_earl_d)
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_eh_earl_s = asreml(fixed = PV_E.H ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_e_sst) + 
                                         vm(Pedigree2, add_inverse_e_nst) + 
                                         vm(Pedigree, sca_inverse_e) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_early,
                                         weights = Weights_E.H,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        model_eh_earl_s = update(model_eh_earl_s)
        
        # Silking:
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_si_earl_d = asreml(fixed = PV_S.I ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_e_sst) + 
                                         vm(Pedigree2, add_inverse_e_nst) + 
                                         vm(Pedigree, dom_inverse_e) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_early,
                                         weights = Weights_S.I,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        model_si_earl_d = update(model_si_earl_d)
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_si_earl_s = asreml(fixed = PV_S.I ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_e_sst) + 
                                         vm(Pedigree2, add_inverse_e_nst) + 
                                         vm(Pedigree, sca_inverse_e) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,                                         
                                         data = blues_early,
                                         weights = Weights_S.I,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        # Anthesis:
        
        # Considering GEI (diagonal) and assuming additive and dominance matrices for lines and hybrids:     
        model_an_earl_d = asreml(fixed = PV_A.N ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_e_sst) + 
                                         vm(Pedigree2, add_inverse_e_nst) + 
                                         vm(Pedigree, dom_inverse_e) + 
                                         idh(Year_Location):Pedigree1 + idh(Year_Location):Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_early,
                                         weights = Weights_A.N,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        # Considering GEI (diagonal) and assuming additive and dominance matrices for lines and hybrids:     
        model_an_earl_s = asreml(fixed = PV_A.N ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_e_sst) + 
                                         vm(Pedigree2, add_inverse_e_nst) + 
                                         vm(Pedigree, sca_inverse_e) + 
                                         idh(Year_Location):Pedigree1 + idh(Year_Location):Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_early,
                                         weights = Weights_A.N,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Intermediate:
        
        # Grain yield:
        
        # Model 2 - Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_gy_inte_d = asreml(fixed = PV_G.Y ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_i_sst) + 
                                         vm(Pedigree2, add_inverse_i_nst) + 
                                         vm(Pedigree, dom_inverse_i) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +                                
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_intermediate,
                                         weights = Weights_G.Y,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_gy_inte_s = asreml(fixed = PV_G.Y ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_i_sst) + 
                                         vm(Pedigree2, add_inverse_i_nst) + 
                                         vm(Pedigree, sca_inverse_i) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +                                
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_intermediate,
                                         weights = Weights_G.Y,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        # Plant height:
        
        # Considering GEI (identity) and assuming additive and dominance matrices for lines and hybrids:     
        model_ph_inte_d = asreml(fixed = PV_P.H ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_i_sst) + 
                                         vm(Pedigree2, add_inverse_i_nst) + 
                                         vm(Pedigree, dom_inverse_i) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         Year_Location:Pedigree1:Pedigree2,
                                         data = blues_intermediate,
                                         weights = Weights_P.H,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        # Considering GEI (identity) and assuming additive and dominance matrices for lines and hybrids:     
        model_ph_inte_s = asreml(fixed = PV_P.H ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_i_sst) + 
                                         vm(Pedigree2, add_inverse_i_nst) + 
                                         vm(Pedigree, sca_inverse_i) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         Year_Location:Pedigree1:Pedigree2,
                                         data = blues_intermediate,
                                         weights = Weights_P.H,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        # Ear height:
        
        # Model 2 - Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_eh_inte_d = asreml(fixed = PV_E.H ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_i_sst) + 
                                         vm(Pedigree2, add_inverse_i_nst) + 
                                         vm(Pedigree, dom_inverse_i) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_intermediate,
                                         weights = Weights_E.H,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)

        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_eh_inte_s = asreml(fixed = PV_E.H ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_i_sst) + 
                                         vm(Pedigree2, add_inverse_i_nst) + 
                                         vm(Pedigree, sca_inverse_i) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_intermediate,
                                         weights = Weights_E.H,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        # Silking:
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_si_inte_d = asreml(fixed = PV_S.I ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_i_sst) + 
                                         vm(Pedigree2, add_inverse_i_nst) + 
                                         vm(Pedigree, dom_inverse_i) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_intermediate,
                                         weights = Weights_S.I,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_si_inte_s = asreml(fixed = PV_S.I ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_i_sst) + 
                                         vm(Pedigree2, add_inverse_i_nst) + 
                                         vm(Pedigree, sca_inverse_i) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_intermediate,
                                         weights = Weights_S.I,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        model_si_inte_s = update(model_si_inte_s)
        model_si_inte_s = update(model_si_inte_s)

        # Anthesis:
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_an_inte_d = asreml(fixed = PV_A.N ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_i_sst) + 
                                         vm(Pedigree2, add_inverse_i_nst) + 
                                         vm(Pedigree, dom_inverse_i) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_intermediate,
                                         weights = Weights_A.N,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        model_an_inte_d = update(model_an_inte_d)
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_an_inte_s = asreml(fixed = PV_A.N ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_i_sst) + 
                                         vm(Pedigree2, add_inverse_i_nst) + 
                                         vm(Pedigree, sca_inverse_i) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_intermediate,
                                         weights = Weights_A.N,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        model_an_inte_s = update(model_an_inte_s)
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Late:
        
        # Grain yield:
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_gy_late_d = asreml(fixed = PV_G.Y ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_l_sst) + 
                                         vm(Pedigree2, add_inverse_l_nst) + 
                                         vm(Pedigree, dom_inverse_l) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +                                
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_late,
                                         weights = Weights_G.Y,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_gy_late_s = asreml(fixed = PV_G.Y ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_l_sst) + 
                                         vm(Pedigree2, add_inverse_l_nst) + 
                                         vm(Pedigree, sca_inverse_l) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +                                
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_late,
                                         weights = Weights_G.Y,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        # Plant height:
        
        # Considering GEI (identity) and assuming additive and dominance matrices for lines and hybrids:     
        model_ph_late_d = asreml(fixed = PV_P.H ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_l_sst) + 
                                         vm(Pedigree2, add_inverse_l_nst) + 
                                         vm(Pedigree, dom_inverse_l) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         Year_Location:Pedigree1:Pedigree2,
                                         data = blues_late,
                                         weights = Weights_P.H,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        # Considering GEI (identity) and assuming additive and dominance matrices for lines and hybrids:     
        model_ph_late_s = asreml(fixed = PV_P.H ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_l_sst) + 
                                         vm(Pedigree2, add_inverse_l_nst) + 
                                         vm(Pedigree, sca_inverse_l) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         Year_Location:Pedigree1:Pedigree2,
                                         data = blues_late,
                                         weights = Weights_P.H,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)

        # Ear height:
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_eh_late_d = asreml(fixed = PV_E.H ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_l_sst) + 
                                         vm(Pedigree2, add_inverse_l_nst) + 
                                         vm(Pedigree, dom_inverse_l) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_late,
                                         weights = Weights_E.H,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        model_eh_late_d = update(model_eh_late_d)

        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_eh_late_s = asreml(fixed = PV_E.H ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_l_sst) + 
                                         vm(Pedigree2, add_inverse_l_nst) + 
                                         vm(Pedigree, sca_inverse_l) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_late,
                                         weights = Weights_E.H,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        # Silking:
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_si_late_d = asreml(fixed = PV_S.I ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_l_sst) + 
                                         vm(Pedigree2, add_inverse_l_nst) + 
                                         vm(Pedigree, dom_inverse_l) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_late,
                                         weights = Weights_S.I,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        # Model 5 - Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_si_late_s = asreml(fixed = PV_S.I ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_l_sst) + 
                                         vm(Pedigree2, add_inverse_l_nst) + 
                                         vm(Pedigree, sca_inverse_l) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_late,
                                         weights = Weights_S.I,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)

        # Anthesis:
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_an_late_d = asreml(fixed = PV_A.N ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_l_sst) + 
                                         vm(Pedigree2, add_inverse_l_nst) + 
                                         vm(Pedigree, dom_inverse_l) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_late,
                                         weights = Weights_A.N,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_an_late_s = asreml(fixed = PV_A.N ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_l_sst) + 
                                         vm(Pedigree2, add_inverse_l_nst) + 
                                         vm(Pedigree, sca_inverse_l) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues_late,
                                         weights = Weights_A.N,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        model_an_late_s = update(model_an_late_s)
        model_an_late_s = update(model_an_late_s)
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Extracting the variance components:
        
        # Early:
        summary(model_gy_earl_d)$varcomp
        summary(model_gy_earl_s)$varcomp
        summary(model_ph_earl_d)$varcomp
        summary(model_ph_earl_s)$varcomp
        summary(model_eh_earl_d)$varcomp
        summary(model_eh_earl_s)$varcomp
        summary(model_si_earl_d)$varcomp
        summary(model_si_earl_s)$varcomp
        summary(model_an_earl_d)$varcomp
        summary(model_an_earl_s)$varcomp
        
        # Intermediate:
        summary(model_gy_inte_d)$varcomp
        summary(model_gy_inte_s)$varcomp
        summary(model_ph_inte_d)$varcomp
        summary(model_ph_inte_s)$varcomp
        summary(model_eh_inte_d)$varcomp
        summary(model_eh_inte_s)$varcomp
        summary(model_si_inte_d)$varcomp
        summary(model_si_inte_s)$varcomp
        summary(model_an_inte_d)$varcomp
        summary(model_an_inte_s)$varcomp
        
        # Late:
        summary(model_gy_late_d)$varcomp
        summary(model_gy_late_s)$varcomp
        summary(model_ph_late_d)$varcomp
        summary(model_ph_late_s)$varcomp
        summary(model_eh_late_d)$varcomp
        summary(model_eh_late_s)$varcomp
        summary(model_si_late_d)$varcomp
        summary(model_si_late_s)$varcomp
        summary(model_an_late_d)$varcomp
        summary(model_an_late_s)$varcomp
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Checking the significance of each component:
        
        # Early:
        
        # Grain yield:
        model_gy_earl_redd_1 = update(model_gy_earl_d, random = ~ . - vm(Pedigree1, add_inverse_e_sst))
        model_gy_earl_redd_2 = update(model_gy_earl_d, random = ~ . - vm(Pedigree2, add_inverse_e_nst))
        model_gy_earl_redd_3 = update(model_gy_earl_d, random = ~ . - vm(Pedigree, dom_inverse_e))
        
        lrt(model_gy_earl_d, model_gy_earl_redd_1)
        lrt(model_gy_earl_d, model_gy_earl_redd_2)
        lrt(model_gy_earl_d, model_gy_earl_redd_3)
        
        model_gy_earl_reds_1 = update(model_gy_earl_s, random = ~ . - vm(Pedigree1, add_inverse_e_sst))
        model_gy_earl_reds_2 = update(model_gy_earl_s, random = ~ . - vm(Pedigree2, add_inverse_e_nst))
        model_gy_earl_reds_3 = update(model_gy_earl_s, random = ~ . - vm(Pedigree, sca_inverse_e))
        
        lrt(model_gy_earl_s, model_gy_earl_reds_1)
        lrt(model_gy_earl_s, model_gy_earl_reds_2)
        lrt(model_gy_earl_s, model_gy_earl_reds_3)
        
        # Plant height:
        model_ph_earl_redd_1 = update(model_ph_earl_d, random = ~ . - vm(Pedigree1, add_inverse_e_sst))
        model_ph_earl_redd_2 = update(model_ph_earl_d, random = ~ . - vm(Pedigree2, add_inverse_e_nst))
        model_ph_earl_redd_3 = update(model_ph_earl_d, random = ~ . - vm(Pedigree, dom_inverse_e))
        
        lrt(model_ph_earl_d, model_ph_earl_redd_1)
        lrt(model_ph_earl_d, model_ph_earl_redd_2)
        lrt(model_ph_earl_d, model_ph_earl_redd_3)
        
        model_ph_earl_reds_1 = update(model_ph_earl_s, random = ~ . - vm(Pedigree1, add_inverse_e_sst))
        model_ph_earl_reds_2 = update(model_ph_earl_s, random = ~ . - vm(Pedigree2, add_inverse_e_nst))
        model_ph_earl_reds_3 = update(model_ph_earl_s, random = ~ . - vm(Pedigree, sca_inverse_e))
        
        lrt(model_ph_earl_s, model_ph_earl_reds_1)
        lrt(model_ph_earl_s, model_ph_earl_reds_2)
        lrt(model_ph_earl_s, model_ph_earl_reds_3)
        
        # Ear height:
        model_eh_earl_redd_1 = update(model_eh_earl_d, random = ~ . - vm(Pedigree1, add_inverse_e_sst))
        model_eh_earl_redd_2 = update(model_eh_earl_d, random = ~ . - vm(Pedigree2, add_inverse_e_nst))
        model_eh_earl_redd_3 = update(model_eh_earl_d, random = ~ . - vm(Pedigree, dom_inverse_e))
        
        lrt(model_eh_earl_d, model_eh_earl_redd_1)
        lrt(model_eh_earl_d, model_eh_earl_redd_2)
        lrt(model_eh_earl_d, model_eh_earl_redd_3)
        
        model_eh_earl_reds_1 = update(model_eh_earl_s, random = ~ . - vm(Pedigree1, add_inverse_e_sst))
        model_eh_earl_reds_2 = update(model_eh_earl_s, random = ~ . - vm(Pedigree2, add_inverse_e_nst))
        model_eh_earl_reds_3 = update(model_eh_earl_s, random = ~ . - vm(Pedigree, sca_inverse_e))
        
        lrt(model_eh_earl_s, model_eh_earl_reds_1)
        lrt(model_eh_earl_s, model_eh_earl_reds_2)
        lrt(model_eh_earl_s, model_eh_earl_reds_3)
        
        # Silking:
        model_si_earl_redd_1 = update(model_si_earl_d, random = ~ . - vm(Pedigree1, add_inverse_e_sst))
        model_si_earl_redd_2 = update(model_si_earl_d, random = ~ . - vm(Pedigree2, add_inverse_e_nst))
        model_si_earl_redd_3 = update(model_si_earl_d, random = ~ . - vm(Pedigree, dom_inverse_e))
        
        lrt(model_si_earl_d, model_si_earl_redd_1)
        lrt(model_si_earl_d, model_si_earl_redd_2)
        lrt(model_si_earl_d, model_si_earl_redd_3)
        
        model_si_earl_reds_1 = update(model_si_earl_s, random = ~ . - vm(Pedigree1, add_inverse_e_sst))
        model_si_earl_reds_2 = update(model_si_earl_s, random = ~ . - vm(Pedigree2, add_inverse_e_nst))
        model_si_earl_reds_3 = update(model_si_earl_s, random = ~ . - vm(Pedigree, sca_inverse_e))
        
        lrt(model_si_earl_s, model_si_earl_reds_1)
        lrt(model_si_earl_s, model_si_earl_reds_2)
        lrt(model_si_earl_s, model_si_earl_reds_3)
        
        # Anthesis:
        model_an_earl_redd_1 = update(model_an_earl_d, random = ~ . - vm(Pedigree1, add_inverse_e_sst))
        model_an_earl_redd_2 = update(model_an_earl_d, random = ~ . - vm(Pedigree2, add_inverse_e_nst))
        model_an_earl_redd_3 = update(model_an_earl_d, random = ~ . - vm(Pedigree, dom_inverse_e))
        
        lrt(model_an_earl_d, model_an_earl_redd_1)
        lrt(model_an_earl_d, model_an_earl_redd_2)
        lrt(model_an_earl_d, model_an_earl_redd_3)
        
        model_an_earl_reds_1 = update(model_an_earl_s, random = ~ . - vm(Pedigree1, add_inverse_e_sst))
        model_an_earl_reds_2 = update(model_an_earl_s, random = ~ . - vm(Pedigree2, add_inverse_e_nst))
        model_an_earl_reds_3 = update(model_an_earl_s, random = ~ . - vm(Pedigree, sca_inverse_e))
        
        lrt(model_an_earl_s, model_an_earl_reds_1)
        lrt(model_an_earl_s, model_an_earl_reds_2)
        lrt(model_an_earl_s, model_an_earl_reds_3)
        
        # Intermediate:
        
        # Grain yield:
        model_gy_inte_redd_1 = update(model_gy_inte_d, random = ~ . - vm(Pedigree1, add_inverse_i_sst))
        model_gy_inte_redd_2 = update(model_gy_inte_d, random = ~ . - vm(Pedigree2, add_inverse_i_nst))
        model_gy_inte_redd_3 = update(model_gy_inte_d, random = ~ . - vm(Pedigree, dom_inverse_i))
        
        lrt(model_gy_inte_d, model_gy_inte_redd_1)
        lrt(model_gy_inte_d, model_gy_inte_redd_2)
        lrt(model_gy_inte_d, model_gy_inte_redd_3)
        
        model_gy_inte_reds_1 = update(model_gy_inte_s, random = ~ . - vm(Pedigree1, add_inverse_i_sst))
        model_gy_inte_reds_2 = update(model_gy_inte_s, random = ~ . - vm(Pedigree2, add_inverse_i_nst))
        model_gy_inte_reds_3 = update(model_gy_inte_s, random = ~ . - vm(Pedigree, sca_inverse_i))
        
        lrt(model_gy_inte_s, model_gy_inte_reds_1)
        lrt(model_gy_inte_s, model_gy_inte_reds_2)
        lrt(model_gy_inte_s, model_gy_inte_reds_3)
        
        # Plant height:
        model_ph_inte_redd_1 = update(model_ph_inte_d, random = ~ . - vm(Pedigree1, add_inverse_i_sst))
        model_ph_inte_redd_2 = update(model_ph_inte_d, random = ~ . - vm(Pedigree2, add_inverse_i_nst))
        model_ph_inte_redd_3 = update(model_ph_inte_d, random = ~ . - vm(Pedigree, dom_inverse_i))
        
        lrt(model_ph_inte_d, model_ph_inte_redd_1)
        lrt(model_ph_inte_d, model_ph_inte_redd_2)
        lrt(model_ph_inte_d, model_ph_inte_redd_3)
        
        model_ph_inte_reds_1 = update(model_ph_inte_s, random = ~ . - vm(Pedigree1, add_inverse_i_sst))
        model_ph_inte_reds_2 = update(model_ph_inte_s, random = ~ . - vm(Pedigree2, add_inverse_i_nst))
        model_ph_inte_reds_3 = update(model_ph_inte_s, random = ~ . - vm(Pedigree, sca_inverse_i))
        
        lrt(model_ph_inte_s, model_ph_inte_reds_1)
        lrt(model_ph_inte_s, model_ph_inte_reds_2)
        lrt(model_ph_inte_s, model_ph_inte_reds_3)
        
        # Ear height:
        model_eh_inte_redd_1 = update(model_eh_inte_d, random = ~ . - vm(Pedigree1, add_inverse_i_sst))
        model_eh_inte_redd_2 = update(model_eh_inte_d, random = ~ . - vm(Pedigree2, add_inverse_i_nst))
        model_eh_inte_redd_3 = update(model_eh_inte_d, random = ~ . - vm(Pedigree, dom_inverse_i))
        
        lrt(model_eh_inte_d, model_eh_inte_redd_1)
        lrt(model_eh_inte_d, model_eh_inte_redd_2)
        lrt(model_eh_inte_d, model_eh_inte_redd_3)
        
        model_eh_inte_reds_1 = update(model_eh_inte_s, random = ~ . - vm(Pedigree1, add_inverse_i_sst))
        model_eh_inte_reds_2 = update(model_eh_inte_s, random = ~ . - vm(Pedigree2, add_inverse_i_nst))
        model_eh_inte_reds_3 = update(model_eh_inte_s, random = ~ . - vm(Pedigree, sca_inverse_i))
        
        lrt(model_eh_inte_s, model_eh_inte_reds_1)
        lrt(model_eh_inte_s, model_eh_inte_reds_2)
        lrt(model_eh_inte_s, model_eh_inte_reds_3)
        
        # Silking:
        model_si_inte_redd_1 = update(model_si_inte_d, random = ~ . - vm(Pedigree1, add_inverse_i_sst))
        model_si_inte_redd_2 = update(model_si_inte_d, random = ~ . - vm(Pedigree2, add_inverse_i_nst))
        model_si_inte_redd_3 = update(model_si_inte_d, random = ~ . - vm(Pedigree, dom_inverse_i))
        
        lrt(model_si_inte_d, model_si_inte_redd_1)
        lrt(model_si_inte_d, model_si_inte_redd_2)
        lrt(model_si_inte_d, model_si_inte_redd_3)
        
        model_si_inte_reds_1 = update(model_si_inte_s, random = ~ . - vm(Pedigree1, add_inverse_i_sst))
        model_si_inte_reds_2 = update(model_si_inte_s, random = ~ . - vm(Pedigree2, add_inverse_i_nst))
        model_si_inte_reds_3 = update(model_si_inte_s, random = ~ . - vm(Pedigree, sca_inverse_i))
        
        lrt(model_si_inte_s, model_si_inte_reds_1)
        lrt(model_si_inte_s, model_si_inte_reds_2)
        lrt(model_si_inte_s, model_si_inte_reds_3)
        
        # Anthesis:
        model_an_inte_redd_1 = update(model_an_inte_d, random = ~ . - vm(Pedigree1, add_inverse_i_sst))
        model_an_inte_redd_2 = update(model_an_inte_d, random = ~ . - vm(Pedigree2, add_inverse_i_nst))
        model_an_inte_redd_3 = update(model_an_inte_d, random = ~ . - vm(Pedigree, dom_inverse_i))
        
        lrt(model_an_inte_d, model_an_inte_redd_1)
        lrt(model_an_inte_d, model_an_inte_redd_2)
        lrt(model_an_inte_d, model_an_inte_redd_3)
        
        model_an_inte_reds_1 = update(model_an_inte_s, random = ~ . - vm(Pedigree1, add_inverse_i_sst))
        model_an_inte_reds_2 = update(model_an_inte_s, random = ~ . - vm(Pedigree2, add_inverse_i_nst))
        model_an_inte_reds_3 = update(model_an_inte_s, random = ~ . - vm(Pedigree, sca_inverse_i))
        
        lrt(model_an_inte_s, model_an_inte_reds_1)
        lrt(model_an_inte_s, model_an_inte_reds_2)
        lrt(model_an_inte_s, model_an_inte_reds_3)
        
        # Late:
        
        # Grain yield:
        model_gy_late_redd_1 = update(model_gy_late_d, random = ~ . - vm(Pedigree1, add_inverse_l_sst))
        model_gy_late_redd_2 = update(model_gy_late_d, random = ~ . - vm(Pedigree2, add_inverse_l_nst))
        model_gy_late_redd_3 = update(model_gy_late_d, random = ~ . - vm(Pedigree, dom_inverse_l))
        
        lrt(model_gy_late_d, model_gy_late_redd_1)
        lrt(model_gy_late_d, model_gy_late_redd_2)
        lrt(model_gy_late_d, model_gy_late_redd_3)
        
        model_gy_late_reds_1 = update(model_gy_late_s, random = ~ . - vm(Pedigree1, add_inverse_l_sst))
        model_gy_late_reds_2 = update(model_gy_late_s, random = ~ . - vm(Pedigree2, add_inverse_l_nst))
        model_gy_late_reds_3 = update(model_gy_late_s, random = ~ . - vm(Pedigree, sca_inverse_l))
        
        lrt(model_gy_late_s, model_gy_late_reds_1)
        lrt(model_gy_late_s, model_gy_late_reds_2)
        lrt(model_gy_late_s, model_gy_late_reds_3)
        
        # Plant height:
        model_ph_late_redd_1 = update(model_ph_late_d, random = ~ . - vm(Pedigree1, add_inverse_l_sst))
        model_ph_late_redd_2 = update(model_ph_late_d, random = ~ . - vm(Pedigree2, add_inverse_l_nst))
        model_ph_late_redd_3 = update(model_ph_late_d, random = ~ . - vm(Pedigree, dom_inverse_l))
        
        lrt(model_ph_late_d, model_ph_late_redd_1)
        lrt(model_ph_late_d, model_ph_late_redd_2)
        lrt(model_ph_late_d, model_ph_late_redd_3)
        
        model_ph_late_reds_1 = update(model_ph_late_s, random = ~ . - vm(Pedigree1, add_inverse_l_sst))
        model_ph_late_reds_2 = update(model_ph_late_s, random = ~ . - vm(Pedigree2, add_inverse_l_nst))
        model_ph_late_reds_3 = update(model_ph_late_s, random = ~ . - vm(Pedigree, sca_inverse_l))
        
        lrt(model_ph_late_s, model_ph_late_reds_1)
        lrt(model_ph_late_s, model_ph_late_reds_2)
        lrt(model_ph_late_s, model_ph_late_reds_3)
        
        # Ear height:
        model_eh_late_redd_1 = update(model_eh_late_d, random = ~ . - vm(Pedigree1, add_inverse_l_sst))
        model_eh_late_redd_2 = update(model_eh_late_d, random = ~ . - vm(Pedigree2, add_inverse_l_nst))
        model_eh_late_redd_3 = update(model_eh_late_d, random = ~ . - vm(Pedigree, dom_inverse_l))
        
        lrt(model_eh_late_d, model_eh_late_redd_1)
        lrt(model_eh_late_d, model_eh_late_redd_2)
        lrt(model_eh_late_d, model_eh_late_redd_3)
        
        model_eh_late_reds_1 = update(model_eh_late_s, random = ~ . - vm(Pedigree1, add_inverse_l_sst))
        model_eh_late_reds_2 = update(model_eh_late_s, random = ~ . - vm(Pedigree2, add_inverse_l_nst))
        model_eh_late_reds_3 = update(model_eh_late_s, random = ~ . - vm(Pedigree, sca_inverse_l))
        
        lrt(model_eh_late_s, model_eh_late_reds_1)
        lrt(model_eh_late_s, model_eh_late_reds_2)
        lrt(model_eh_late_s, model_eh_late_reds_3)
        
        # Silking:
        model_si_late_redd_1 = update(model_si_late_d, random = ~ . - vm(Pedigree1, add_inverse_l_sst))
        model_si_late_redd_2 = update(model_si_late_d, random = ~ . - vm(Pedigree2, add_inverse_l_nst))
        model_si_late_redd_3 = update(model_si_late_d, random = ~ . - vm(Pedigree, dom_inverse_l))
        
        lrt(model_si_late_d, model_si_late_redd_1)
        lrt(model_si_late_d, model_si_late_redd_2)
        lrt(model_si_late_d, model_si_late_redd_3)
        
        model_si_late_reds_1 = update(model_si_late_s, random = ~ . - vm(Pedigree1, add_inverse_l_sst))
        model_si_late_reds_2 = update(model_si_late_s, random = ~ . - vm(Pedigree2, add_inverse_l_nst))
        model_si_late_reds_3 = update(model_si_late_s, random = ~ . - vm(Pedigree, sca_inverse_l))
        
        lrt(model_si_late_s, model_si_late_reds_1)
        lrt(model_si_late_s, model_si_late_reds_2)
        lrt(model_si_late_s, model_si_late_reds_3)
        
        # Anthesis:
        model_an_late_redd_1 = update(model_an_late_d, random = ~ . - vm(Pedigree1, add_inverse_l_sst))
        model_an_late_redd_2 = update(model_an_late_d, random = ~ . - vm(Pedigree2, add_inverse_l_nst))
        model_an_late_redd_3 = update(model_an_late_d, random = ~ . - vm(Pedigree, dom_inverse_l))
        
        lrt(model_an_late_d, model_an_late_redd_1)
        lrt(model_an_late_d, model_an_late_redd_2)
        lrt(model_an_late_d, model_an_late_redd_3)
        
        model_an_late_reds_1 = update(model_an_late_s, random = ~ . - vm(Pedigree1, add_inverse_l_sst))
        model_an_late_reds_2 = update(model_an_late_s, random = ~ . - vm(Pedigree2, add_inverse_l_nst))
        model_an_late_reds_3 = update(model_an_late_s, random = ~ . - vm(Pedigree, sca_inverse_l))
        
        lrt(model_an_late_s, model_an_late_reds_1)
        lrt(model_an_late_s, model_an_late_reds_2)
        lrt(model_an_late_s, model_an_late_reds_3)
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        