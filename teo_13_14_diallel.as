!WORKSPACE 8 !NOGRAPHICS !RENAME 2 !ARGS SeedWeight 3 !OUTFOLDER diallel_out # $A $B DTA DTS PlantHt LeafLen LeafWd TillNum EarsPerSideBran LatBranLen LBIL EarLen SeedsPerRow FC_Width_mm SeedsPerEar EarInternodeLen PredSeedPerPlant PredTotalSeedWeight SeedWeight MalePerc ProporYoked ProporNonShatter ProporPairedSpikelets PredSeedPerPlant PredTotalSeedWeight #
Title: Teosinte 13-14 Diallel analysis sequence of models
#Plot,Year,Type,F_est,Mother,Father,Shading,X_axis,Y_axis,Border,DTA,DTS,PlantHt,LeafLen,LeafWd,TillNum,EarsPerSideBran,LatBranNodeNum,LatBranLen,LBIL,EarLen,SeedsPerRow,FC_Width_mm,SeedsPerEar,EarInternodeLen,PredSeedPerPlant,PredTotalSeedWeight,SeedWeight,CulmDiam,BranNum,MalePerc,Clustering,ASI,FruitcaseWt,ProporYoked,ProporNonShatter,ProporPairedSpikelets,Perc_Germ,FC_Length_mm,FC_LWR,FC_Triangularity,com,dr,TypeYear
#A0003,2014,self,0.5894369,PC_N14_ID2,PC_N14_ID2,2802.29,0,2,E,68,62,129,32,4,10,6,2,480,240,38.5,4,3.907524,8,4.81,3931,88,0.0224,1.38,9,100,0.0125,6,0.0448,0,0.038,0.015,95,5.760111,1.479436,1.397664773,PC_N14_ID2xPC_N14_ID2,1,self2014
#A0004,2014,self,0.4263937,PC_J48_ID2,PC_J48_ID2,2978.09,0,3,E,70,70,86,43,5,6,8,4,200,50,50,4.5,3.455967,9,5.56,3999,85,0.0213,0.8,6,0,0.04,0,0.0425,0,0,0,100,6.257701,1.819066,1.39098873,PC_J48_ID2xPC_J48_ID2,1,self2014
#A0008,2014,outcross,-0.03638047,PC_J13_ID1,PC_K60_ID1,2145.39,0,7,E,64,61,152,40,5.4,13,14,3,520,173.33,46.5,4.5,3.694183,9,5.17,6616,143,0.0216,1.55,10,100,0.0269,3,0.0432,0,0.006,0,95,5.505086,1.498049,1.407056807,PC_J13_ID1xPC_K60_ID1,1,outcross2014
#A0013,2014,self,0.4377412,PC_N10_ID2,PC_N10_ID2,3634.49,0,12,E,70,67,125,41,5.6,8,8,2,310,155,55,4.75,4.134554,9.5,5.79,4580,137,0.0299,1.41,7,100,0.0258,3,0.0598,0,0.012,0.012,95,5.776356,1.40113,1.4310004,PC_N10_ID2xPC_N10_ID2,1,self2014
 Plot  !A      # A0013
 Year  !I      # 2014
 Type  !A  !SORT    # self
 F_est        # 0.4377412
 Mother  !P #!A !SORT #     # PC_N10_ID2
 Father  !P #!A !AS Mother  # PC_N10_ID2 !AS Mother needed to match mom/dad levels if not using pedigree
 Shading        # 3634.49
 X_axis  #!I       # 0
 Y_axis  #!I      # 12
 Border  !A      # E
 DTA        # 70
 DTS        # 67
 PlantHt        # 125
 LeafLen        # 41
 LeafWd        # 5.6
 TillNum         # 8
 EarsPerSideBran         # 8
 LatBranNodeNum         # 2
 LatBranLen        # 310
 LBIL        # 155
 EarLen       # 55
 SeedsPerRow        # 4.75
 FC_Width_mm        # 4.134554
 SeedsPerEar        # 9.5
 EarInternodeLen        # 5.79
 PredSeedPerPlant  !/100 #need to reduce the scale of measurement for convergence, varcomps need to be multiplied by 10000 to be on original scale      # 4580
 PredTotalSeedWeight     # 137
 SeedWeight !*100  #need to increase scale for convergence, varcomps need to be divided by 10000     # 0.0299
 CulmDiam        # 1.41
 BranNum         # 7
 MalePerc        # 100
 Clustering        # 0.0258
 ASI         # 3
 FruitcaseWt        # 0.0598
 ProporYoked         # 0
 ProporNonShatter        # 0.012
 ProporPairedSpikelets        # 0.012
 Perc_Germ    # 95
 FC_Length_mm     # 5.776356
 FC_LWR        # 1.40113
 FC_Triangularity        # 1.4310004
 com  !A  !LL 21   # PC_N10_ID2xPC_N10_ID2
 dr         # 1
 TypeYear  !A      # self2014

Teo_parental_pedigree.csv !ALPHA !SKIP 1   #model the pedigree relationships among the parents, some teo parents are half-sibs
teo_imp_02102017_traits-2_for_diallel.csv  !SKIP 1  !DOPATH $B !MAXITER 75

!PATH 1 #base model, including the relationships among the parents
$A ~ mu Year Type Year.Type Shading Year.Border at(Year).pol(X_axis,-4) at(Year).pol(Y_axis,-4),
       !r at(Type,Outcross).nrm(Mother) and(at(Type,Outcross).nrm(Father)), #GCA in outcrosses
       at(Type,Self).nrm(Mother), #S1 variance
       #!r str(at(Type,Outcross).idv(Mother) and(at(Type,Outcross).idv(Father)) at(Type,Self).idv(Mother) us(2).idv(Mother)), #GCA, S1 var, and their cov
       #at(Type,Outcross).Mother.dr and(at(Type,Outcross).Father.dr), #reciprocal GCA
       at(Type,Outcross).com , #SCA
       #at(Type,Outcross).com.dr, #reciprocal SCA
       at(Type,Outcross).Moth.idv(Year) and(at(Type,Outcross).Fath.idv(Year)),  #GCA*Year
       at(Type,Self).Moth.idv(Year),  #S1*Year
       #at(Type,Outcross).Mo.Ye.dr and(at(Type,Outcross).Father.Year.dr), #reciprocal GCA*Year #this term is confounded so leave it out
       at(Type,Outcross).c.Ye #SCA.Year
       #at(Type,Outcross).c.Ye.dr #reciprocal SCA.year
       residual at(Year).id(units) #residual variance separate for separate years, but common across generations
#predictions need to be made at specific values of covariates, these are arbitrary
predict Type Border 1 X_axis 1 Y_axis 1 !AVERAGE Year
predict Type.Mother Border 1 X_axis 1 Y_axis 1 !AVERAGE Year !PRESENT Type Mother
VPREDICT diallel_out/teo_13-14_diallel.pin !DEFINE
F SCA Type_Outcross.com
F Dominance SCA*4  #Dominance variance
F GCA at(Type,Outcross).nrm(M;nrm(Mother) #at(Type,Outcross).Mother #V_GCA for moms, V_GCA for dads are equal
F Additive GCA*4 #V(GCA) = V(HS) = 1/4(VA)
F S1_Var 6 #at(Type,Self).nrm(Mothe;nrm(Mother) doesn't match #variance among S0:1 families  = Va + (1/4)Vd + D1  +(1/8)D2
F S1_VarA S1_Var*1.5 - Dominance*0.375 #Additive variance among S1 individuals = 1.5(S1_var - 0.25Dominance) bc S1_var = Va + 0.25Vd + other stuff, S1 individuals have 1.5Va
F Res_13 Residual_1  #within family variance in 2013 pooled over generations
F Res_14 Residual_2  #within family variance in 2014 pooled over generations
F A_Year 8*4 #at(Type,Outcross).id(Mo).idv(Year doesn't match, so using the number, BE CAREFUL!
F D_Year Type_Outcross.c.Ye*4 #SCA*Year
F G_out Additive + Dominance #total genetic variance in outcross generation
H D_G_ratio Dominance G_out #ratio of dominance to total variance
F GE_total A_Year + D_Year
F Pheno_out Additive*0.5 + Dominance*0.25 + A_Year*0.25 + D_Year*0.25 + Res_13*0.5 + Res_14*0.5   #Vp = Vgca_m + Vgca_d + Vsca + avg(Vresid_outcross)
F Pheno_S1 S1_Var + 7 + Res_13*0.5 + Res_14*0.5   #Vp for S1 generation can't match at(Type,Self).id(Mother).idv(Year  so using NUMBER, BE CAREFUL!
H h2 Additive Pheno_out  #narrow-sense in outbred individuals
H H_out G_out Pheno_out  # broad-sense in outbred
H H_S1 S1_VarA Pheno_S1 #'heritability' of S1 individuals ignoring D1 and D2
H D_P_ratio Dominance Pheno_out #ratio of dominance variance to total PHENO variance
H GE_P_ratio GE_total Pheno_out #proportion of Vp that is Vge
F S1_Exp Additive + Dominance*0.25 #Genetic variance expected among S1 familes based on outbred Va and Vd estimates (assuming D1, D2, etc = 0)
H S1_Obs_Exp S1_Var S1_Exp #Ratio of genetic variance observed among S1 families to expected

!PATH 2 #fit separate error variances for S0 and S1 generations within each year
$A ~ mu Year Type Year.Type Shading Year.Border at(Year).pol(X_axis,-4) at(Year).pol(Y_axis,-4),
       !r at(Type,Outcross).nrm(Mother) and(at(Type,Outcross).nrm(Father)), #GCA in outcrosses
       at(Type,Self).nrm(Mother), #S1 variance
       #!r str(at(Type,Outcross).id(Mother) and(at(Type,Outcross).id(Father)) at(Type,Self).id(Mother) us(2).id(Mother)), #GCA, S1 var, and their cov
       #at(Type,Outcross).Mother.dr and(at(Type,Outcross).Father.dr), #reciprocal GCA
       at(Type,Outcross).com , #SCA
       #at(Type,Outcross).com.dr, #reciprocal SCA
       at(Type,Outcross).Moth.idv(Year) and(at(Type,Outcross).Fath.idv(Year)) at(Type,Self).Moth.idv(Year),  #GCA*Year S0, S1, not fitting cov bc of convergence problems
       #at(Type,Outcross).Mo.Ye.dr and(at(Type,Outcross).Father.Year.dr), #reciprocal GCA*Year #this term is confounded so leave it out
       at(Type,Outcross).c.Ye #SCA.Year
       #at(Type,Outcross).c.Ye.dr #reciprocal SCA.year
       residual at(TypeYear).id(units) #residual variance separate for outbred and S1 individuals each in separate years
#predictions need to be made at specific values of covariates, these are arbitrary
predict Type Border 1 X_axis 1 Y_axis 1 !AVERAGE Year
predict Type.Mother Border 1 X_axis 1 Y_axis 1 !AVERAGE Year !PRESENT Type Mother
VPREDICT diallel_out/teo_13-14_diallel.pin !DEFINE
F SCA Type_Outcross.com
F Dominance SCA*4  #Dominance variance
F GCA 7 #at(Type,Outcross).nrm(M;nrm(Mother) not matching #V_GCA for moms, V_GCA for dads are equal
F Additive GCA*4 #V(GCA) = V(HS) = 1/4(VA)
F S1_Var 8 #at(Type,Self).nrm(Mothe;nrm(Mother) doesn't match #variance among S0:1 families  = Va + (1/4)Vd + D1  +(1/8)D2
F S1_VarA S1_Var*1.5 - Dominance*0.375 #Additive variance among S1 individuals = 1.5(S1_var - 0.25Dominance) bc S1_var = Va + 0.25Vd + other stuff, S1 individuals have 1.5Va
F ResS0_13 Residual_1  #within FS family variance in 2013(1/2 Va + 3/4 Vd + error) I figured out which residual corresponded to which group only by counting observations and checking against number of effects listed in asreml output
F ResS1_13 Residual_2  #within S0:1 family variance in 2013 (1/2 Va + 1/4 Vd + error)
F ResS0_14 Residual_3  #within FS family variance in 2014 (1/2 Va + 3/4 Vd + error)
F ResS1_14 Residual_4  #within S0:1 family variance in 2014 (1/2 Va + 3/4 Vd + error)
F ErrorS0 ResS0_13*0.5 + ResS0_14*0.5 - Additive*0.5 - Dominance*0.75 #error variance in S0s, removing genetic variances from within-family variance
F ErrorS1 ResS1_13*0.5 + ResS1_14*0.5 - Additive*0.5 - Dominance*0.25 #error variance in S1s, removing genetic variances from within-family variance
F A_Year 10*4 #at(Type,Outcross).id(Mo).idv(Year doesn't match, so using the number, BE CAREFUL!
F D_Year Type_Outcross.c.Ye*4 #SCA*Year
F G_out Additive + Dominance #total genetic variance in outcross generation
H D_G_ratio Dominance G_out #ratio of dominance to total variance
F GE_total A_Year + D_Year
F Pheno_out Additive*0.5 + Dominance*0.25 + A_Year*0.25 + D_Year*0.25 + ResS0_13*0.5 + ResS0_14*0.5   #Vp = Vgca_m + Vgca_d + Vsca + avg(Vresid_outcross)
F Pheno_S1 S1_Var + 9 + ResS1_13*0.5 + ResS1_14*0.5   #Vp for S1 generation can't match at(Type,Self).id(Mother).idv(Year  so using NUMBER, BE CAREFUL!
H h2 Additive Pheno_out  #narrow-sense in outbred individuals
H H_out G_out Pheno_out  # broad-sense in outbred
H H_S1 S1_VarA Pheno_S1 #'heritability' of S1 individuals ignoring D1 and D2
H D_P_ratio Dominance Pheno_out #ratio of dominance variance to total PHENO variance
H GE_P_ratio GE_total Pheno_out #proportion of Vp that is Vge
F S1_Exp Additive + Dominance*0.25 #Genetic variance expected among S1 familes based on outbred Va and Vd estimates (assuming D1, D2, etc = 0)
H S1_Obs_Exp S1_Var S1_Exp #Ratio of genetic variance observed among S1 families to expected

!PATH 3 #add Covariance between S0 and S1  generations
$A ~ mu Year Type Year.Type Shading Year.Border at(Year).pol(X_axis,-4) at(Year).pol(Y_axis,-4),
       #!r at(Type,Outcross).Mother and(at(Type,Outcross).Father), #GCA in outcrosses
       #at(Type,Self).Mother, #S1 variance
       !r str(at(Type,Outcross).nrm(Mother) and(at(Type,Outcross).nrm(Father)) at(Type,Self).nrm(Mother) us(2, !GUUU).nrm(Mother)), #GCA, S1 var, and their cov
       #at(Type,Outcross).Mother.dr and(at(Type,Outcross).Father.dr), #reciprocal GCA
       at(Type,Outcross).com , #SCA
       #at(Type,Outcross).com.dr, #reciprocal SCA
       at(Type,Outcross).Moth.idv(Year) and(at(Type,Outcross).Fath.idv(Year)) at(Type,Self).Moth.idv(Year),  #GCA*Year S0, S1, not fitting cov bc of convergence problems
       #at(Type,Outcross).Mo.Ye.dr and(at(Type,Outcross).Father.Year.dr), #reciprocal GCA*Year #this term is confounded so leave it out
       at(Type,Outcross).c.Ye #SCA.Year
       #at(Type,Outcross).c.Ye.dr #reciprocal SCA.year
       residual at(TypeYear).id(units) #residual variance separate for outbred and S1 individuals each in separate years
#predictions need to be made at specific values of covariates, these are arbitrary
predict Type Border 1 X_axis 1 Y_axis 1 !AVERAGE Year
predict Type.Mother Border 1 X_axis 1 Y_axis 1 !AVERAGE Year !PRESENT Type Mother
#VPREDICT 'Q:/My Drive/Teo and landrace/Landrace/Diallel/diallel_out/LR_13-14_diallel.pin' !DEFINE
VPREDICT diallel_out/teo_13-14_diallel.pin !DEFINE
F SCA Type_Outcross.com
F Dominance SCA*4  #Dominance variance
F GCA us(2).nrm(Mother);us(2)[1] #at(Type,Outcross).Mother #V_GCA for moms, V_GCA for dads are equal
F Additive GCA*4 #V(GCA) = V(HS) = 1/4(VA)
F S1_Var us(2).nrm(Mother);us(2)[3] #at(Type,Self).Mother  #variance among S0:1 families  = Va + (1/4)Vd + D1  +(1/8)D2
F S1_VarA S1_Var*1.5 - Dominance*0.375 #Additive variance among S1 individuals = 1.5(S1_var - 0.25Dominance) bc S1_var = Va + 0.25Vd + other stuff, S1 individuals have 1.5Va
F CS0S1 us(2).nrm(Mother);us(2)[2] #Covariance of additive breeding values in S0,S1 generations
R rS0S1 us(2).nrm(Mother);us(2) #form S0/S1 breeding value correlation from covariance and variances
F ResS0_13 Residual_1  #within FS family variance in 2013(1/2 Va + 3/4 Vd + error) I figured out which residual corresponded to which group only by counting observations and checking against number of effects listed in asreml output
F ResS1_13 Residual_2  #within S0:1 family variance in 2013 (1/2 Va + 1/4 Vd + error)
F ResS0_14 Residual_3  #within FS family variance in 2014 (1/2 Va + 3/4 Vd + error)
F ResS1_14 Residual_4  #within S0:1 family variance in 2014 (1/2 Va + 3/4 Vd + error)
F ErrorS0 ResS0_13*0.5 + ResS0_14*0.5 - Additive*0.5 - Dominance*0.75 #error variance in S0s, removing genetic variances from within-family variance
F ErrorS1 ResS1_13*0.5 + ResS1_14*0.5 - Additive*0.5 - Dominance*0.25 #error variance in S1s, removing genetic variances from within-family variance
F A_Year 11*4 #at(Type,Outcross).id(Mo).idv(Year doesn't match, so using the number, BE CAREFUL!
F D_Year Type_Outcross.c.Ye*4 #SCA*Year
F G_out Additive + Dominance #total genetic variance in outcross generation
H D_G_ratio Dominance G_out #ratio of dominance to total variance
F GE_total A_Year + D_Year
F Pheno_out Additive*0.5 + Dominance*0.25 + A_Year*0.25 + D_Year*0.25 + ResS0_13*0.5 + ResS0_14*0.5   #Vp = Vgca_m + Vgca_d + Vsca + avg(Vresid_outcross)
F Pheno_S1 S1_Var + 10 + ResS1_13*0.5 + ResS1_14*0.5   #Vp for S1 generation can't match at(Type,Self).id(Mother).idv(Year  so using NUMBER, BE CAREFUL!
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
$A ~ mu Year Type Year.Type Shading Year.Border at(Year).pol(X_axis,-4) at(Year).pol(Y_axis,-4),
       #!r at(Type,Outcross).Mother and(at(Type,Outcross).Father), #GCA in outcrosses
       #at(Type,Self).Mother, #S1 variance
       !r str(at(Type,Outcross).nrm(Mother) and(at(Type,Outcross).nrm(Father)) at(Type,Self).nrm(Mother) us(2, !GUUU).nrm(Mother)), #GCA, S1 var, and their cov
       at(Type,Outcross).nrm(Mother).idv(dr) and(at(Type,Outcross).nrm(Father).idv(dr)), #reciprocal GCA
       at(Type,Outcross).com , #SCA
       at(Type,Outcross).com.dr, #reciprocal SCA
       at(Type,Outcross).Moth.idv(Year) and(at(Type,Outcross).Fath.idv(Year)) at(Type,Self).Moth.idv(Year),  #GCA*Year S0, S1
       #at(Type,Outcross).Mo.Ye.dr and(at(Type,Outcross).Father.Year.dr), #reciprocal GCA*Year #this term is confounded so leave it out
       at(Type,Outcross).c.Ye, #SCA.Year
       at(Type,Outcross).c.Ye.dr #reciprocal SCA.year
       residual at(TypeYear).id(units) #residual variance separate for outbred and S1 individuals each in separate years
#predictions need to be made at specific values of covariates, these are arbitrary
predict Type Border 1 X_axis 1 Y_axis 1 !AVERAGE Year
predict Type.Mother Border 1 X_axis 1 Y_axis 1 !AVERAGE Year !PRESENT Type Mother
#VPREDICT 'Q:/My Drive/Teo and landrace/Landrace/Diallel/diallel_out/LR_13-14_diallel.pin' !DEFINE
VPREDICT diallel_out/teo_13-14_diallel.pin !DEFINE
F SCA Type_Outcross.com
F Dominance SCA*4  #Dominance variance
F GCA us(2).nrm(Mother);us(2)[1] #at(Type,Outcross).Mother #V_GCA for moms, V_GCA for dads are equal
F Additive GCA*4 #V(GCA) = V(HS) = 1/4(VA)
F S1_Var us(2).nrm(Mother);us(2)[3] #at(Type,Self).Mother  #variance among S0:1 families  = Va + (1/4)Vd + D1  +(1/8)D2
F S1_VarA S1_Var*1.5 - Dominance*0.375 #Additive variance among S1 individuals = 1.5(S1_var - 0.25Dominance) bc S1_var = Va + 0.25Vd + other stuff, S1 individuals have 1.5Va
F CS0S1 us(2).nrm(Mother);us(2)[2] #Covariance of additive breeding values in S0,S1 generations
R rS0S1 us(2).nrm(Mother);us(2) #form S0/S1 breeding value correlation from covariance and variances
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
