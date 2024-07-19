*KEYWORD
*PARAMETER
$#   prmr1      val1     prmr2      val2     prmr3      val3     prmr4      val4
R dt      1e-6
R dt0     0.1
R dtmax   0.5
$
$ Section and material ID
$ Use "1" for both to use user-defined shell element
$ Use "2" for both to use LS-Dyna axisymmetric element 15 and Mat24 for elasticity
I secid   1
I MID     1
$
$ Number of quadrature points
I NQP     4
I NIP     4
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
         1        20         5  -0.00001    &dtmax       0.0         0         0
*CONTROL_IMPLICIT_GENERAL
$#  imflag       dt0    imform      nsbs       igs     cnstn      form    zero_v
         1      &dt0         2         1         1         0         0         0
*CONTROL_IMPLICIT_SOLVER
$#  lsolvr    lprint     negev     order      drcm    drcprm   autospc   autotol
         2         0         2         0         4       0.0         1       0.0
$#  lcpack    mtxdmp    iparm1    rparm1    rparm2      
         2         0       5001.00000E-9     0.001
$#  emxdmp    rdcmem      
         0      0.85
*CONTROL_SOLUTION
$#    soln       nlq     isnan     lcint     lcacc     ncdcf     
         0         0         0      1000         0         1
*DATABASE_ELOUT
$#      dt    binary      lcur     ioopt   option1   option2   option3   option4
&dt                3         0         1         9         0         0         0
*DATABASE_NODOUT
$#      dt    binary      lcur     ioopt   option1   option2       
&dt                3         0         1       0.0         0
*DATABASE_SECFORC
$#      dt    binary      lcur     ioopt     
&dt                3         0         1
*DATABASE_BINARY_D3PLOT
$#      dt      lcdt      beam     npltc    psetid      
&dt                0         0         0         0
$#   ioopt      rate    cutoff    window      type      pset    
         0       0.0       0.0       0.0         0         0
*DATABASE_EXTENT_BINARY
$#   neiph     neips    maxint    strflg    sigflg    epsflg    rltflg    engflg
        15        15         3         0         1         1         1         1
$#  cmpflg    ieverp    beamip     dcomp      shge     stssz    n3thdt   ialemat
         0         0         0         1         1         1         2         1
$# nintsld   pkp_sen      sclp     hydro     msscl     therm    intout    nodout
         0         0       1.0         0         0         0STRESS
$ HYDRO=1: output the "internal energy density" (incompatible with paraview)
$         0         0       1.0         1         0         0                    
$#    dtdt    resplt     neipb     quadr     cubic     
         0         0         0         0         0
*DATABASE_ELOUT
$#      dt    binary      lcur     ioopt   option1   option2   option3   option4
      0.01         3         0         1         0         0         0         0
$
$ ******************************************************************************
*PARAMETER
$#   prmr1      val1     prmr2      val2     prmr3      val3     prmr4      val4
R DT_OUT  1e-6
$
*PARAMETER
$
$#   prmr1      val1     prmr2      val2     prmr3      val3     prmr4      val4
R MID     1
R HG_onoff0
I NQP     4
$
*INCLUDE
3_ELFORM15_2Dax_elastic.inc
3_UEL_Q1X_ushl_2Dax_linElastic.inc
$
$ ******************************************************************************
$
*INCLUDE
2_HTV_mesh_t15e-1_m3_PID.inc
3_HTV_main.k
*END
