*KEYWORD
*PARAMETER
$
$#   prmr1      val1     prmr2      val2     prmr3      val3     prmr4      val4
R tmax    10
R loadmax 0.1
R dtout   1e-3
$ Section and material ID
$ Use "1" for both to use user-defined shell element
$ Use "2" for both to use LS-Dyna plane strain element 13 and Mat24 for elasticity
I secid   1
I MID     1
$
$ Number of quadrature points
I NQP     4
I NIP     4
I HG_onoff         0
$
*PARAMETER
$#   prmr1      val1     prmr2      val2     prmr3      val3     prmr4      val4
$ Elasticity
R BULKMOD  164285.71
R SHEARMOD 80232.558
R YOUNGMOD  207000.0
R POISSON       0.29
R DENSITY    7.83e-9
$
*INCLUDE
3_UEL_Q1X_ushl_2D_linElastic.inc
3_ELFORM_planeStrain.inc
$
*CONTROL_IMPLICIT_GENERAL
$#  imflag       dt0    imform      nsbs       igs     cnstn      form    zero_v
         1       1.0         2         1         1         0         0         0
$
*CONTROL_TIMESTEP
$#  dtinit    tssfac      isdo    tslimt     dt2ms      lctm     erode     ms1st
       0.0      0.90         0       0.0         0         0         0         0
$#  dt2msf   dt2mslc     imscl    unused    unused     rmscl    unused      ihdo
       0.0         0         0                           0.0                   0
*CONTROL_IMPLICIT_AUTO
$#   iauto    iteopt    itewin     dtmin     dtmax     dtexp     kfail    kcycle
         1        15         5      0.01       1.0       0.0         0         0
*CONTROL_IMPLICIT_SOLUTION
$#  nsolvr    ilimit    maxref     dctol     ectol     rctol     lstol    abstol
        12         4        16     0.001      0.011.00000E10       0.91.0000E-20
$#   dnorm    diverg     istif   nlprint    nlnorm   d3itctl     cpchk     
         1         2         1         2         2         0         0
$#  arcctl    arcdir    arclen    arcmth    arcdmp    arcpsi    arcalf    arctim
         0         0       0.0         1         2         0         0         0
$#   lsmtd     lsdir      irad      srad      awgt      sred    
         2         2       0.0       0.0       0.0       0.0
*CONTROL_IMPLICIT_SOLVER
$#  lsolvr    lprint     negev     order      drcm    drcprm   autospc   autotol
         2         0         2         0         4       0.0         1       0.0
$#  lcpack    mtxdmp    iparm1    rparm1    rparm2      
         2         0       5001.00000E-9     0.001
$#  emxdmp    rdcmem      
         0      0.85
*CONTROL_TERMINATION
$#  endtim    endcyc     dtmin    endeng    endmas     nosol     
&tmax              0       0.0       0.0       0.0         0
*DATABASE_ELOUT
$#      dt    binary      lcur     ioopt   option1   option2   option3   option4
&dtout             3         0         1         9         0         0         0
*DATABASE_NODFOR
$#      dt    binary      lcur     ioopt     
&dtout             3         0         1
*DATABASE_NODOUT
$#      dt    binary      lcur     ioopt   option1   option2       
&dtout             3         0         1       0.0         0
*DATABASE_SECFORC
$#      dt    binary      lcur     ioopt     
&dtout             3         0         1
*DATABASE_BINARY_D3PLOT
$#      dt      lcdt      beam     npltc    psetid      
&dtout             0         0         0         0
$#   ioopt      rate    cutoff    window      type      pset    
         0       0.0       0.0       0.0         0         0
*DATABASE_NODAL_FORCE_GROUP
$#    nsid       cid   
        10         0
$
*BOUNDARY_PRESCRIBED_MOTION_SET
$#    nsid       dof       vad      lcid        sf       vid     death     birth
        10         2         2         1       1.0         01.00000E28       0.0
*CONSTRAINED_GLOBAL
$#      tc        rc       dir         x         y         z       tol   
         1         0         1       0.0       0.0       0.0       0.0
         2         0         2       0.0       0.0       0.0       0.0
$
*PART
$#                                                                         title
shell_4p
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         1    &secid      &mid         0         0         0         2         0
*DEFINE_CURVE_TITLE
Linear_Load
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
         1         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
                 0.0                 0.0
&tmax               &loadmax
100                 &loadmax
$
*INCLUDE
2_OET.inc
*END
