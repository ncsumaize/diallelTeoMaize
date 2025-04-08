!WORKSPACE 8 !NOGRAPHICS !RENAME 2 !ARGS  pEarN 4 !OUTFOLDER diallel_out# $A $B #DTA DTS PltHt LfLen LfWth TillN Fleck tEarN HkLfN ShkLen ShkIntLen CobLen KrowN K1rowN EarInternodeLen CobDia pKN tKN tKW weight_avKW50 !OUTFOLDER diallel_out
Title: Landraces_FL1314_all_data_updated_4398taxa_Fest-Multivariate_analysis.
 PlantID  !A      # 13FL990004
 Year     !I      # 2013
 Plot_year  *     # 4
 seed_source !A   # 12FL0166-10
 Expected_Mother  !A      # 166_10
 Mother  !P       # 166_10  #need to declare Mother and Father as pedigree factors, otherwise diallel structure fails since they don't have same levels (some plants are only male or female parents) 
 Father  !P       # 166_10
 Type    !A       # Self
 Column       # 1
 Row         # 4
 Germ            # 1
 Tire_track     # 0
 Shading          # 1507.75
 #Planting_date  !MDY      # 9/27/2013
 #Tassel_date    !MDY      # 12/10/2013
 #Silk_date      !MDY      # 12/26/2013
 !CSKIP 3
 DTA
 DTS
 ASI
 PltHt
 EarHt
 PEHt
 TslLen
 mSpLen
 BZLen
 RSpLenTslLen
 TslpbrCt
 AvBZIntLen
 LfLen
 LfWth
 Rust
 TillN !*10 #increase scale to help convergence
 eRtLdg
 eStkLdg
 RtLdg
 StkLdg
 Fleck
 VnChl
 pEarN
 sEarN
 tEarN
 HkLfN
 ShkLen
 ShkIntLen
 CobLen
 KrowN
 K1rowN
 EarInternodeLen
 CobDia
 HkCol
 pKN
 pKW
 pKW50
 sKN
 sKW
 sKW50
 tKN
 tKW
 weight_avKW50
 F_est
 com !A #cross combination
 dr #dummy variable coding for reciprocal crosses
 TypeYear !A #combination of Type and Year to set residual variance

'Q:/My Drive/Teo and landrace/Landrace/Genetic Architecture Simulation/LR_parent_ped.txt' !SKIP 1 !ALPHA #pedigree
'Q:/My Drive/Teo and landrace/Landrace/Diallel/Landraces_FL1314_all_data_updated_4398taxa_for_diallel.csv' !SKIP 1  !DOPATH $B !MAXITER 75

!PATH 1 #base model
$A ~ mu Year Type Year.Type Shading Year.Tire_track at(Year).pol(Column,-4) at(Year).pol(Row,-4),
       !r at(Type,Outcross).idv(Mother) and(at(Type,Outcross).idv(Father)), #GCA in outcrosses
       at(Type,Self).idv(Mother), #S1 variance
       #!r str(at(Type,Outcross).idv(Mother) and(at(Type,Outcross).idv(Father)) at(Type,Self).idv(Mother) us(2).Mother), #GCA, S1 var, and their cov
       #at(Type,Outcross).Mother.idv(dr) and(at(Type,Outcross).Father.idv(dr)), #reciprocal GCA
       at(Type,Outcross).com , #SCA
       #at(Type,Outcross).com.dr, #reciprocal SCA
       at(Type,Outcross).Mother.idv(Year) and(at(Type,Outcross).Father.idv(Year)) at(Type,Self).Mother.idv(Year),  #GCA*Year S0, S1, not fitting cov bc of convergence problems
       #at(Type,Outcross).Mo.Ye.dr and(at(Type,Outcross).Father.Year.dr), #reciprocal GCA*Year #this term is confounded so leave it out
       at(Type,Outcross).c.Ye #SCA.Year
      # at(Type,Outcross).c.Ye.dr #reciprocal SCA.year
       residual at(Year).id(units) #residual variance same for outbred and S1 individuals each in separate years

#predictions need to be made at specific values of covariates, these are arbitrary
predict Type Tire_track 0 Row 1 Column 1 !AVERAGE Year
predict Type.Mother Tire_track 0 Row 1 Column 1  !AVERAGE Year !PRESENT Type Mother

VPREDICT diallel_out/LR_13_14_diallel.pin !DEFINE
F SCA Type_Outcross.com
F Dominance SCA*4  #Dominance variance
F GCA 5 #at(Type,Outcross).idv(Mother) doesn't match, BE CAREFUL! #V_GCA for moms, V_GCA for dads are equal
F Additive GCA*4 #V(GCA) = V(HS) = 1/4(VA)
F S1_Var 6 #at(Type,Self).nrm(M;nrm(Mother) doesn't match #variance among S0:1 families  = Va + (1/4)Vd + D1  +(1/8)D2
F S1_VarA S1_Var*1.5 - Dominance*0.375 #Additive variance among S1 individuals = 1.5(S1_var - 0.25Dominance) bc S1_var = Va + 0.25Vd + other stuff, S1 individuals have 1.5Va
F Res_13 Residual_1  #within family variance in 2013 pooled over generations
F Res_14 Residual_2  #within family variance in 2014 pooled over generations
F A_Year 8*4 #at(Type,Outcross).nrm(Mo).idv(Year doesn't match, so using the number, BE CAREFUL!
F D_Year Type_Outcross.c.Ye*4 #SCA*Year
F G_out Additive + Dominance #total genetic variance in outcross generation
H D_G_ratio Dominance G_out #ratio of dominance to total variance
F GE_total A_Year + D_Year
F Pheno_out Additive*0.5 + Dominance*0.25 + A_Year*0.25 + D_Year*0.25 + Res_13*0.5 + Res_14*0.5   #Vp = Vgca_m + Vgca_d + Vsca + avg(Vresid_outcross)
F Pheno_S1 S1_Var + 7 + Res_13*0.5 + Res_14*0.5   #Vp for S1 generation can't match at(Type,Self).nrm(Mother).idv(Year  so using NUMBER, BE CAREFUL!
H h2 Additive Pheno_out  #narrow-sense in outbred individuals
H H_out G_out Pheno_out  # broad-sense in outbred
H H_S1 S1_VarA Pheno_S1 #'heritability' of S1 individuals ignoring D1 and D2
H D_P_ratio Dominance Pheno_out #ratio of dominance variance to total PHENO variance
H GE_P_ratio GE_total Pheno_out #proportion of Vp that is Vge
F S1_Exp Additive + Dominance*0.25 #Genetic variance expected among S1 familes based on outbred Va and Vd estimates (assuming D1, D2, etc = 0)
H S1_Obs_Exp S1_Var S1_Exp #Ratio of genetic variance observed among S1 families to expected


!PATH 2 #add heterogeneous residual variances for S0 and S1 generations
$A ~ mu Year Type Year.Type Shading Year.Tire_track at(Year).pol(Column,-4) at(Year).pol(Row,-4),
       !r at(Type,Outcross).idv(Mother) and(at(Type,Outcross).idv(Father)), #GCA in outcrosses
       at(Type,Self).idv(Mother), #S1 variance
       #!r str(at(Type,Outcross).idv(Mother) and(at(Type,Outcross).idv(Father)) at(Type,Self).idv(Mother) us(2).Mother), #GCA, S1 var, and their cov
       #at(Type,Outcross).Mother.idv(dr) and(at(Type,Outcross).Father.idv(dr)), #reciprocal GCA
       at(Type,Outcross).com , #SCA
       #at(Type,Outcross).com.dr, #reciprocal SCA
       at(Type,Outcross).Mother.idv(Year) and(at(Type,Outcross).Father.idv(Year)) at(Type,Self).Mother.idv(Year),  #GCA*Year S0, S1, not fitting cov bc of convergence problems
       #at(Type,Outcross).Mo.Ye.dr and(at(Type,Outcross).Father.Year.dr), #reciprocal GCA*Year #this term is confounded so leave it out
       at(Type,Outcross).c.Ye #SCA.Year
      # at(Type,Outcross).c.Ye.dr #reciprocal SCA.year
       residual at(TypeYear).id(units) #residual variance separate for outbred and S1 individuals each in separate years

#predictions need to be made at specific values of covariates, these are arbitrary
predict Type Tire_track 0 Row 1 Column 1 !AVERAGE Year
predict Type.Mother Tire_track 0 Row 1 Column 1  !AVERAGE Year !PRESENT Type Mother

VPREDICT diallel_out/LR_13_14_diallel.pin !DEFINE
F SCA Type_Outcross.com
F Dominance SCA*4  #Dominance variance
F GCA 7 #at(Type,Outcross).nrm(M;nrm(Mother) not matching #at(Type,Outcross).Mother #V_GCA for moms, V_GCA for dads are equal
F Additive GCA*4 #V(GCA) = V(HS) = 1/4(VA)
F S1_Var 8 #at(Type,Self).nrm(M;nrm(Mother) doesn't match #variance among S0:1 families  = Va + (1/4)Vd + D1  +(1/8)D2
F S1_VarA S1_Var*1.5 - Dominance*0.375 #Additive variance among S1 individuals = 1.5(S1_var - 0.25Dominance) bc S1_var = Va + 0.25Vd + other stuff, S1 individuals have 1.5Va
F ResS0_13 Residual_1  #within FS family variance in 2013(1/2 Va + 3/4 Vd + error) I figured out which residual corresponded to which group only by counting observations and checking against number of effects listed in asreml output
F ResS1_13 Residual_2  #within S0:1 family variance in 2013 (1/2 Va + 1/4 Vd + error)
F ResS0_14 Residual_3  #within FS family variance in 2014 (1/2 Va + 3/4 Vd + error)
F ResS1_14 Residual_4  #within S0:1 family variance in 2014 (1/2 Va + 3/4 Vd + error)
F ErrorS0 ResS0_13*0.5 + ResS0_14*0.5 - Additive*0.5 - Dominance*0.75 #error variance in S0s, removing genetic variances from within-family variance
F ErrorS1 ResS1_13*0.5 + ResS1_14*0.5 - Additive*0.5 - Dominance*0.25 #error variance in S1s, removing genetic variances from within-family variance
F A_Year 10*4 #at(Type,Outcross).nrm(Mo).idv(Year doesn't match, so using the number, BE CAREFUL!
F D_Year Type_Outcross.c.Ye*4 #SCA*Year
F G_out Additive + Dominance #total genetic variance in outcross generation
H D_G_ratio Dominance G_out #ratio of dominance to total variance
F GE_total A_Year + D_Year
F Pheno_out Additive*0.5 + Dominance*0.25 + A_Year*0.25 + D_Year*0.25 + ResS0_13*0.5 + ResS0_14*0.5   #Vp = Vgca_m + Vgca_d + Vsca + avg(Vresid_outcross)
F Pheno_S1 S1_Var + 9 + ResS1_13*0.5 + ResS1_14*0.5   #Vp for S1 generation can't match at(Type,Self).nrm(Mother).idv(Year  so using NUMBER, BE CAREFUL!
H h2 Additive Pheno_out  #narrow-sense in outbred individuals
H H_out G_out Pheno_out  # broad-sense in outbred
H H_S1 S1_VarA Pheno_S1 #'heritability' of S1 individuals ignoring D1 and D2
H D_P_ratio Dominance Pheno_out #ratio of dominance variance to total PHENO variance
H GE_P_ratio GE_total Pheno_out #proportion of Vp that is Vge
F S1_Exp Additive + Dominance*0.25 #Genetic variance expected among S1 familes based on outbred Va and Vd estimates (assuming D1, D2, etc = 0)
H S1_Obs_Exp S1_Var S1_Exp #Ratio of genetic variance observed among S1 families to expected



!PATH 3 #add Covariance between S0 and S1  generations
$A ~ mu Year Type Year.Type Shading Year.Tire_track at(Year).pol(Column,-4) at(Year).pol(Row,-4),
       #!r at(Type,Outcross).Mother and(at(Type,Outcross).Father), #GCA in outcrosses
       #at(Type,Self).Mother, #S1 variance
       !r str(at(Type,Outcross).idv(Mother) and(at(Type,Outcross).idv(Father)) at(Type,Self).idv(Mother) us(2, !GUUU).Mother), #GCA, S1 var, and their cov
       #at(Type,Outcross).Mother.idv(dr) and(at(Type,Outcross).Father.idv(dr)), #reciprocal GCA
       at(Type,Outcross).com , #SCA
       #at(Type,Outcross).com.dr, #reciprocal SCA
       at(Type,Outcross).Mother.idv(Year) and(at(Type,Outcross).Father.idv(Year)) at(Type,Self).Mother.idv(Year),  #GCA*Year S0, S1, not fitting cov bc of convergence problems
       #at(Type,Outcross).Mo.Ye.dr and(at(Type,Outcross).Father.Year.dr), #reciprocal GCA*Year #this term is confounded so leave it out
       at(Type,Outcross).c.Ye #SCA.Year
      # at(Type,Outcross).c.Ye.dr #reciprocal SCA.year
       residual at(TypeYear).id(units) #residual variance separate for outbred and S1 individuals each in separate years

#predictions need to be made at specific values of covariates, these are arbitrary
predict Type Tire_track 0 Row 1 Column 1 !AVERAGE Year
predict Type.Mother Tire_track 0 Row 1 Column 1  !AVERAGE Year !PRESENT Type Mother

VPREDICT diallel_out/LR_13_14_diallel.pin !DEFINE
F SCA Type_Outcross.com
F Dominance SCA*4  #Dominance variance
F GCA us(2).Mother;us(2)[1] #at(Type,Outcross).Mother #V_GCA for moms, V_GCA for dads are equal
F Additive GCA*4 #V(GCA) = V(HS) = 1/4(VA)
F S1_Var us(2).Mother;us(2)[3] #at(Type,Self).Mother  #variance among S0:1 families  = Va + (1/4)Vd + D1  +(1/8)D2
F S1_VarA S1_Var*1.5 - Dominance*0.375 #Additive variance among S1 individuals = 1.5(S1_var - 0.25Dominance) bc S1_var = Va + 0.25Vd + other stuff, S1 individuals have 1.5Va
F CS0S1 us(2).Mother;us(2)[2] #Covariance of additive breeding values in S0,S1 generations
R rS0S1 us(2).Mother;us #form S0/S1 breeding value correlation from covariance and variances
F ResS0_13 Residual_1  #within FS family variance in 2013(1/2 Va + 3/4 Vd + error) I figured out which residual corresponded to which group only by counting observations and checking against number of effects listed in asreml output
F ResS1_13 Residual_2  #within S0:1 family variance in 2013 (1/2 Va + 1/4 Vd + error)
F ResS0_14 Residual_3  #within FS family variance in 2014 (1/2 Va + 3/4 Vd + error)
F ResS1_14 Residual_4  #within S0:1 family variance in 2014 (1/2 Va + 3/4 Vd + error)
F ErrorS0 ResS0_13*0.5 + ResS0_14*0.5 - Additive*0.5 - Dominance*0.75 #error variance in S0s, removing genetic variances from within-family variance
F ErrorS1 ResS1_13*0.5 + ResS1_14*0.5 - Additive*0.5 - Dominance*0.25 #error variance in S1s, removing genetic variances from within-family variance
F GCA_Year 11 #at(Type,Outcross).Mother.Year #GCA *Year
F SCA_Year Type_Outcross.c.Ye #SCA*Year
F G_out Additive + Dominance  #total genetic variance in outcross generation
H D_G_ratio Dominance G_out #ratio of dominance to total variance
F AD  Additive + Dominance #Additive + Dominance variance
H D_AD_ratio Dominance AD #ratio of dominance to total of additive and dominance variance
F GE_total GCA_Year*4 + SCA_Year*4
F Pheno_out Additive*0.5 + Dominance*0.25  + GCA_Year*2 + SCA_Year + ResS0_13*0.5 + ResS0_14*0.5   #Vp = Vgca_m + Vgca_d + Vsca + avg(Vresid_outcross)
F Pheno_S1 S1_Var + 10 + ResS1_13*0.5 + ResS1_14*0.5   #Vp for S1 generation
H h2 Additive Pheno_out  #narrow-sense in outbred individuals
H H_out G_out Pheno_out  # broad-sense in outbred
H H_S1 S1_VarA Pheno_S1 #'heritability' of S1 individuals ignoring D1 and D2
H D_P_ratio Dominance Pheno_out #ratio of dominance variance to total PHENO variance
H GE_P_ratio GE_total Pheno_out #proportion of Vp that is Vge
F S1_Exp Additive + Dominance*0.25 #Genetic variance expected among S1 familes based on outbred Va and Vd estimates (assuming D1, D2, etc = 0)
H S1_Obs_Exp S1_Var S1_Exp #Ratio of genetic variance observed among S1 families to expected
F D1 CS0S1*4 - Additive*2
F D2 S1_Var*8 - Additive*8 - Dominance*2 - D1*8
H D1_G D1 G_out #ratio of D1 to total outcross variance
H D2_G D2 G_out #ratio of D2 to total outcross variance



!PATH 4 #Full model, add reciprocal variances
$A ~ mu Year Type Year.Type Shading Year.Tire_track at(Year).pol(Column,-4) at(Year).pol(Row,-4),
       #!r at(Type,Outcross).Mother and(at(Type,Outcross).Father), #GCA in outcrosses
       #at(Type,Self).Mother, #S1 variance
       !r str(at(Type,Outcross).idv(Mother) and(at(Type,Outcross).idv(Father)) at(Type,Self).idv(Mother) us(2, !GUUU).Mother), #GCA, S1 var, and their cov
       at(Type,Outcross).Mother.idv(dr) and(at(Type,Outcross).Father.idv(dr)), #reciprocal GCA
       at(Type,Outcross).com , #SCA
       at(Type,Outcross).com.dr, #reciprocal SCA
       at(Type,Outcross).Mother.idv(Year) and(at(Type,Outcross).Father.idv(Year)) at(Type,Self).Mother.idv(Year),  #GCA*Year S0, S1, not fitting cov bc of convergence problems
       #at(Type,Outcross).Mo.Ye.dr and(at(Type,Outcross).Father.Year.dr), #reciprocal GCA*Year #this term is confounded so leave it out
       at(Type,Outcross).c.Ye, #SCA.Year
       at(Type,Outcross).c.Ye.dr #reciprocal SCA.year
       residual at(TypeYear).id(units) #residual variance separate for outbred and S1 individuals each in separate years

#predictions need to be made at specific values of covariates, these are arbitrary
predict Type Tire_track 0 Row 1 Column 1 !AVERAGE Year
predict Type.Mother Tire_track 0 Row 1 Column 1  !AVERAGE Year !PRESENT Type Mother

VPREDICT diallel_out/LR_13_14_diallel.pin !DEFINE
F SCA Type_Outcross.com
F Dominance SCA*4  #Dominance variance
F GCA us(2).Mother;us(2)[1] #at(Type,Outcross).Mother #V_GCA for moms, V_GCA for dads are equal
F Additive GCA*4 #V(GCA) = V(HS) = 1/4(VA)
F S1_Var us(2).Mother;us(2)[3] #at(Type,Self).Mother  #variance among S0:1 families  = Va + (1/4)Vd + D1  +(1/8)D2
F S1_VarA S1_Var*1.5 - Dominance*0.375 #Additive variance among S1 individuals = 1.5(S1_var - 0.25Dominance) bc S1_var = Va + 0.25Vd + other stuff, S1 individuals have 1.5Va
F CS0S1 us(2).Mother;us(2)[2] #Covariance of additive breeding values in S0,S1 generations
R rS0S1 us(2).Mother;us #form S0/S1 breeding value correlation from covariance and variances
F ResS0_13 Residual_1  #within FS family variance in 2013(1/2 Va + 3/4 Vd + error) I figured out which residual corresponded to which group only by counting observations and checking against number of effects listed in asreml output
F ResS1_13 Residual_2  #within S0:1 family variance in 2013 (1/2 Va + 1/4 Vd + error)
F ResS0_14 Residual_3  #within FS family variance in 2014 (1/2 Va + 3/4 Vd + error)
F ResS1_14 Residual_4  #within S0:1 family variance in 2014 (1/2 Va + 3/4 Vd + error)
F ErrorS0 ResS0_13*0.5 + ResS0_14*0.5 - Additive*0.5 - Dominance*0.75 #error variance in S0s, removing genetic variances from within-family variance
F ErrorS1 ResS1_13*0.5 + ResS1_14*0.5 - Additive*0.5 - Dominance*0.25 #error variance in S1s, removing genetic variances from within-family variance
F GCA_recip  12   # reciprocal GCA, can't match it any other way, be careful
F SCA_recip Type_Outcross.com.dr  #reciprocal SCA
F RecipA GCA_recip*4 #reciprocal additive
F RecipD SCA_recip*4 #reciprocal dominance
F RecipT RecipA + RecipD #total reciprocal genetic variance
F GCA_Year 14 #at(Type,Outcross).Mother.Year #GCA *Year
F SCA_Year Type_Outcross.c.Ye #SCA*Year
F Recip_Yr Type_Outcross.c.Ye.d*4 #Reciprocal*year variance but we can't estimate GCA_Yr_reciprocal
F G_out Additive + Dominance + RecipT #total genetic variance in outcross generation
H D_G_ratio Dominance G_out #ratio of dominance to total variance
F AD  Additive + Dominance #Additive + Dominance variance
H D_AD_ratio Dominance AD #ratio of dominance to total of additive and dominance variance
H Recip_G_ratio RecipT G_out
F GE_total GCA_Year*4 + SCA_Year*4 + Recip_Yr
F Pheno_out Additive*0.5 + Dominance*0.25 + RecipA*0.5 + RecipD*0.25 + GCA_Year*2 + SCA_Year + Recip_Yr + ResS0_13*0.5 + ResS0_14*0.5   #Vp = Vgca_m + Vgca_d + Vsca + avg(Vresid_outcross)
F Pheno_S1 S1_Var + 13 + ResS1_13*0.5 + ResS1_14*0.5   #Vp for S1 generation
H h2 Additive Pheno_out  #narrow-sense in outbred individuals
H H_out G_out Pheno_out  # broad-sense in outbred
H H_S1 S1_VarA Pheno_S1 #'heritability' of S1 individuals ignoring D1 and D2
H D_P_ratio Dominance Pheno_out #ratio of dominance variance to total PHENO variance
H Recip_P_ratio RecipT Pheno_out  #proportion of Vp that is Vrecip
H GE_P_ratio GE_total Pheno_out #proportion of Vp that is Vge
F S1_Exp Additive + Dominance*0.25 #Genetic variance expected among S1 familes based on outbred Va and Vd estimates (assuming D1, D2, etc = 0)
H S1_Obs_Exp S1_Var S1_Exp #Ratio of genetic variance observed among S1 families to expected
F D1 CS0S1*4 - Additive*2
F D2 S1_Var*8 - Additive*8 - Dominance*2 - D1*8
H D1_G D1 G_out #ratio of D1 to total outcross variance
H D2_G D2 G_out #ratio of D2 to total outcross variance
