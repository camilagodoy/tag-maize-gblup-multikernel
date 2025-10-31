        ##################################################
        #         GBLUP-based multi-kernel models        #   
        #   Estimation of genetic variance components    #
        ##################################################
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Removing variables:
        rm(list=ls())
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Directory, packages, and functions:
        setwd("~/Documents/tag-maize-ss-nss-gblup-multikernel/")
        library("asreml")

        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Getting the blues:
        blues = readRDS("data/single_env_analyses/Final_blues_all_traits.rds")

        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Getting the additive and dominance relationship matrices:
        add_inverse_sst = readRDS("data/genomic_matrices/all_groups/Add_inverse_sst.rds")
        add_inverse_nst = readRDS("data/genomic_matrices/all_groups/Add_inverse_nst.rds")
        dom_inverse = readRDS("data/genomic_matrices/all_groups/Dom_inverse.rds")
        sca_inverse = readRDS("data/genomic_matrices/all_groups/Sca_inverse.rds")
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Fitting the models:
        
        # Grain yield:
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_g.yield_d = asreml(fixed = PV_G.Y ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_sst) + 
                                         vm(Pedigree2, add_inverse_nst) + 
                                         vm(Pedigree, dom_inverse) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +                                
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues,
                                         weights = Weights_G.Y,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        model_g.yield_d = update(model_g.yield_d)

        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_g.yield_s = asreml(fixed = PV_G.Y ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_sst) + 
                                         vm(Pedigree2, add_inverse_nst) + 
                                         vm(Pedigree, sca_inverse) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +                                
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues,
                                         weights = Weights_G.Y,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        model_g.yield_s = update(model_g.yield_s)
       
        # Plant height:
        
        # Considering GEI (identity) and assuming additive and dominance matrices for lines and hybrids:     
        model_p.heigh_d = asreml(fixed = PV_P.H ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_sst) + 
                                         vm(Pedigree2, add_inverse_nst) + 
                                         vm(Pedigree, dom_inverse) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         Year_Location:Pedigree1:Pedigree2,
                                         data = blues,
                                         weights = Weights_P.H,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)

        # Considering GEI (identity) and assuming additive and dominance matrices for lines and hybrids:     
        model_p.heigh_s = asreml(fixed = PV_P.H ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_sst) + 
                                         vm(Pedigree2, add_inverse_nst) + 
                                         vm(Pedigree, sca_inverse) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         Year_Location:Pedigree1:Pedigree2,
                                         data = blues,
                                         weights = Weights_P.H,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        # Ear height:
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_e.heigh_d = asreml(fixed = PV_E.H ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_sst) + 
                                         vm(Pedigree2, add_inverse_nst) + 
                                         vm(Pedigree, dom_inverse) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues,
                                         weights = Weights_E.H,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        model_e.heigh_2 = update(model_e.heigh_d)
        model_e.heigh_2 = update(model_e.heigh_d)
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_e.heigh_s = asreml(fixed = PV_E.H ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_sst) + 
                                         vm(Pedigree2, add_inverse_nst) + 
                                         vm(Pedigree, sca_inverse) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues,
                                         weights = Weights_E.H,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        model_e.heigh_s = update(model_e.heigh_s)

        # Silking:
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_silking_d = asreml(fixed = PV_S.I ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_sst) + 
                                         vm(Pedigree2, add_inverse_nst) + 
                                         vm(Pedigree, dom_inverse) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues,
                                         weights = Weights_S.I,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        model_silking_d = update(model_silking_d)
        
        # Considering GEI (identity + diagonal) and assuming additive and dominance matrices for lines and hybrids: 
        model_silking_s = asreml(fixed = PV_S.I ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_sst) + 
                                         vm(Pedigree2, add_inverse_nst) + 
                                         vm(Pedigree, sca_inverse) + 
                                         Year_Location:Pedigree1 + Year_Location:Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues,
                                         weights = Weights_S.I,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        model_silking_s = update(model_silking_s)

        # Anthesis:
        
        # Considering GEI (diagonal) and assuming additive and dominance matrices for lines and hybrids:     
        model_anthesi_d = asreml(fixed = PV_A.N ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_sst) + 
                                         vm(Pedigree2, add_inverse_nst) + 
                                         vm(Pedigree, dom_inverse) + 
                                         idh(Year_Location):Pedigree1 + idh(Year_Location):Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues,
                                         weights = Weights_A.N,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        model_anthesi_d = update(model_anthesi_d)
        model_anthesi_d = update(model_anthesi_d)

        # Considering GEI (diagonal) and assuming additive and dominance matrices for lines and hybrids:     
        model_anthesi_s = asreml(fixed = PV_A.N ~ Year_Location,
                                         random = ~ vm(Pedigree1, add_inverse_sst) + 
                                         vm(Pedigree2, add_inverse_nst) + 
                                         vm(Pedigree, sca_inverse) + 
                                         idh(Year_Location):Pedigree1 + idh(Year_Location):Pedigree2 +
                                         idh(Year_Location):Pedigree1:Pedigree2,
                                         data = blues,
                                         weights = Weights_A.N,
                                         family = asr_gaussian(dispersion = 1),
                                         na.action = na.method(x = "include", y = "include"), 
                                         workspace = 6e8,
                                         maxit = 100)
        
        model_anthesi_s = update(model_anthesi_s)
        model_anthesi_s = update(model_anthesi_s)
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Extracting the variance components:
        summary(model_g.yield_d)$varcomp
        summary(model_g.yield_s)$varcomp
        summary(model_p.heigh_d)$varcomp
        summary(model_p.heigh_s)$varcomp
        summary(model_e.heigh_d)$varcomp
        summary(model_e.heigh_s)$varcomp
        summary(model_silking_d)$varcomp
        summary(model_silking_s)$varcomp
        summary(model_anthesi_d)$varcomp
        summary(model_anthesi_s)$varcomp
        
        #-------------------------------------------------------------------------------------------------------------------------------------#
        # Checking the significance of each component:
        
        # Grain yield:
        model_gy_redd_1 = update(model_g.yield_d, random = ~ . - vm(Pedigree1, add_inverse_sst))
        model_gy_redd_2 = update(model_g.yield_d, random = ~ . - vm(Pedigree2, add_inverse_nst))
        model_gy_redd_3 = update(model_g.yield_d, random = ~ . - vm(Pedigree, dom_inverse))
        
        lrt(model_g.yield_d, model_gy_redd_1)
        lrt(model_g.yield_d, model_gy_redd_2)
        lrt(model_g.yield_d, model_gy_redd_3)
        
        model_gy_reds_1 = update(model_g.yield_s, random = ~ . - vm(Pedigree1, add_inverse_sst))
        model_gy_reds_2 = update(model_g.yield_s, random = ~ . - vm(Pedigree2, add_inverse_nst))
        model_gy_reds_3 = update(model_g.yield_s, random = ~ . - vm(Pedigree, sca_inverse))
        
        lrt(model_g.yield_s, model_gy_reds_1)
        lrt(model_g.yield_s, model_gy_reds_2)
        lrt(model_g.yield_s, model_gy_reds_3)
        
        # Plant height:
        model_ph_redd_1 = update(model_p.heigh_d, random = ~ . - vm(Pedigree1, add_inverse_sst))
        model_ph_redd_2 = update(model_p.heigh_d, random = ~ . - vm(Pedigree2, add_inverse_nst))
        model_ph_redd_3 = update(model_p.heigh_d, random = ~ . - vm(Pedigree, dom_inverse))
        
        lrt(model_p.heigh_d, model_ph_redd_1)
        lrt(model_p.heigh_d, model_ph_redd_2)
        lrt(model_p.heigh_d, model_ph_redd_3)
        
        model_ph_reds_1 = update(model_p.heigh_s, random = ~ . - vm(Pedigree1, add_inverse_sst))
        model_ph_reds_2 = update(model_p.heigh_s, random = ~ . - vm(Pedigree2, add_inverse_nst))
        model_ph_reds_3 = update(model_p.heigh_s, random = ~ . - vm(Pedigree, sca_inverse))
        
        lrt(model_p.heigh_s, model_ph_reds_1)
        lrt(model_p.heigh_s, model_ph_reds_2)
        lrt(model_p.heigh_s, model_ph_reds_3)
        
        # Ear height:
        model_eh_redd_1 = update(model_e.heigh_d, random = ~ . - vm(Pedigree1, add_inverse_sst))
        model_eh_redd_2 = update(model_e.heigh_d, random = ~ . - vm(Pedigree2, add_inverse_nst))
        model_eh_redd_3 = update(model_e.heigh_d, random = ~ . - vm(Pedigree, dom_inverse))
        
        lrt(model_e.heigh_d, model_eh_redd_1)
        lrt(model_e.heigh_d, model_eh_redd_2)
        lrt(model_e.heigh_d, model_eh_redd_3)
        
        model_eh_reds_1 = update(model_e.heigh_s, random = ~ . - vm(Pedigree1, add_inverse_sst))
        model_eh_reds_2 = update(model_e.heigh_s, random = ~ . - vm(Pedigree2, add_inverse_nst))
        model_eh_reds_3 = update(model_e.heigh_s, random = ~ . - vm(Pedigree, sca_inverse))
        
        lrt(model_e.heigh_s, model_eh_reds_1)
        lrt(model_e.heigh_s, model_eh_reds_2)
        lrt(model_e.heigh_s, model_eh_reds_3)
        
        # Silking:
        model_si_redd_1 = update(model_silking_d, random = ~ . - vm(Pedigree1, add_inverse_sst))
        model_si_redd_2 = update(model_silking_d, random = ~ . - vm(Pedigree2, add_inverse_nst))
        model_si_redd_3 = update(model_silking_d, random = ~ . - vm(Pedigree, dom_inverse))
        
        lrt(model_silking_d, model_si_redd_1)
        lrt(model_silking_d, model_si_redd_2)
        lrt(model_silking_d, model_si_redd_3)
        
        model_si_reds_1 = update(model_silking_s, random = ~ . - vm(Pedigree1, add_inverse_sst))
        model_si_reds_2 = update(model_silking_s, random = ~ . - vm(Pedigree2, add_inverse_nst))
        model_si_reds_3 = update(model_silking_s, random = ~ . - vm(Pedigree, sca_inverse))
        
        lrt(model_silking_s, model_si_reds_1)
        lrt(model_silking_s, model_si_reds_2)
        lrt(model_silking_s, model_si_reds_3)
        
        # Anthesis:
        model_an_redd_1 = update(model_anthesi_d, random = ~ . - vm(Pedigree1, add_inverse_sst))
        model_an_redd_2 = update(model_anthesi_d, random = ~ . - vm(Pedigree2, add_inverse_nst))
        model_an_redd_3 = update(model_anthesi_d, random = ~ . - vm(Pedigree, dom_inverse))
        
        lrt(model_anthesi_d, model_an_redd_1)
        lrt(model_anthesi_d, model_an_redd_2)
        lrt(model_anthesi_d, model_an_redd_3)
        
        model_an_reds_1 = update(model_anthesi_s, random = ~ . - vm(Pedigree1, add_inverse_sst))
        model_an_reds_2 = update(model_anthesi_s, random = ~ . - vm(Pedigree2, add_inverse_nst))
        model_an_reds_3 = update(model_anthesi_s, random = ~ . - vm(Pedigree, sca_inverse))
        
        lrt(model_anthesi_s, model_an_reds_1)
        lrt(model_anthesi_s, model_an_reds_2)
        lrt(model_anthesi_s, model_an_reds_3)
       
        #-------------------------------------------------------------------------------------------------------------------------------------#
        
        
       