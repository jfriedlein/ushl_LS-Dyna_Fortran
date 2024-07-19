*KEYWORD
*PARAMETER
$
$#   prmr1      val1     prmr2      val2     prmr3      val3     prmr4      val4
R dtmax   1
R dt0     1
$ Section and material ID
$ Use "1" for both to use user-defined shell element
$ Use "2" for both to use LS-Dyna axisymmetric element 15 and Mat24 for elasticity
I secid   1
I MID     1
$
$ Number of quadrature points
I NQP     4
$
$ Elasticity
R BULKMOD  164285.71
R SHEARMOD 80232.558
R YOUNGMOD  207000.0
R POISSON       0.29
R DENSITY    7.83e-9
$
*CONTROL_ACCURACY
$#     osu       inn    pidosu      iacc    
         0         1         0         1
*CONTROL_IMPLICIT_AUTO
$#   iauto    iteopt    itewin     dtmin     dtmax     dtexp     kfail    kcycle
         1        15        10     0.001    &dtmax       0.0         0         0
*CONTROL_IMPLICIT_GENERAL
$#  imflag       dt0    imform      nsbs       igs     cnstn      form    zero_v
         1      &dt0         2         1         1         0         0         0
*CONTROL_IMPLICIT_SOLUTION
$#  nsolvr    ilimit    maxref     dctol     ectol     rctol     lstol    abstol
        12         1        201.00000E+81.00000E+81.00000E-8       0.9-1.0000E-8
$#   dnorm    diverg     istif   nlprint    nlnorm   d3itctl     cpchk     
         1         1         1         2         2         0         0
$#  arcctl    arcdir    arclen    arcmth    arcdmp    arcpsi    arcalf    arctim
         0         0       0.0         1         2         0         0         0
$#   lsmtd     lsdir      irad      srad      awgt      sred    
         4         2       0.0       0.0       0.0       0.0
*CONTROL_SOLUTION
$#    soln       nlq     isnan     lcint     lcacc     ncdcf     
         0         0         0      1000         0         1
*DATABASE_ELOUT
$#      dt    binary      lcur     ioopt   option1   option2   option3   option4
       1.0         3         0         1         9         0         0         0
*DATABASE_NODFOR
$#      dt    binary      lcur     ioopt     
       1.0         3         0         1
*DATABASE_NODAL_FORCE_GROUP
$#    nsid       cid   
         1         0
*DATABASE_NODOUT
$#      dt    binary      lcur     ioopt   option1   option2       
       1.0         3         0         1       0.0         0
*DATABASE_SECFORC
$#      dt    binary      lcur     ioopt     
       1.0         3         0         1
*DATABASE_BINARY_D3PLOT
$#      dt      lcdt      beam     npltc    psetid      
       1.0         0         0         0         0
$#   ioopt      rate    cutoff    window      type      pset    
         0       0.0       0.0       0.0         0         0
$
*INCLUDE
3_ELFORM15_2Dax_elastic.inc
3_UEL_Q1X_ushl_2Dax_linElastic.inc
$
*INCLUDE
3_NoaR_ax_Norris_mX.k
2_NoaR_ax_Norris_m2.inc
*END
