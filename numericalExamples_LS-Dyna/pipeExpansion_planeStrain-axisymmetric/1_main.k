$# LS-DYNA Keyword file created by LS-PrePost(R) V4.7.20 -24Aug2020
$# Created on May-23-2022 (12:53:43)
*KEYWORD
*PARAMETER
$
$#   prmr1      val1     prmr2      val2     prmr3      val3     prmr4      val4
R tmax    10
$ Section and material ID
$ Use "1" for both to use user-defined shell element
$ Use "2" for both to use LS-Dyna plane strain element 13 and Mat24 for elasticity
I secid   2
I MID     2
$
$ Number of quadrature points
I NQP     4
I NIP     4
I HG_onoff         0
$
*PARAMETER
$#   prmr1      val1     prmr2      val2     prmr3      val3     prmr4      val4
$ Elasticity
R BULKMOD  499.66667
R SHEARMOD 1
R YOUNGMOD 2.998
R POISSON  0.499
R DENSITY    7.85e-9
$
*INCLUDE
3_UEL_Q1X_ushl_2D_linElastic.inc
3_ELFORM_planeStrain.inc
$3_ELFORM_Axi.inc
$
*CONSTRAINED_GLOBAL
$#      tc        rc       dir         x         y         z       tol   
         1         0         1       0.0       0.0       0.0       0.0
         2         0         2       0.0       0.0       0.0       0.0
$
*CONTROL_TERMINATION
$#  endtim    endcyc     dtmin    endeng    endmas     nosol     
&tmax              0       0.0       0.01.000000E8         0 
*DEFINE_CURVE_TITLE
Curve_load
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
         1         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
                 0.0                 0.0
&tmax                                 -1
*COMMENT  *LOAD_SEGMENT_SET_ID
$#      id                                                               heading
         0Apply internal pressure
$#    ssid      lcid        sf        at  
         1         1      -1.0       0.0
*USER_LOADING
$     LCID       PID       AXI     UNDEF
         1         2         0         1
*CONTROL_ACCURACY
$#     osu       inn    pidosu      iacc    
         0         1         0         1
*CONTROL_IMPLICIT_AUTO
$#   iauto    iteopt    itewin     dtmin     dtmax     dtexp     kfail    kcycle
         1        15         5       0.1       1.0       0.0         0         0
*CONTROL_IMPLICIT_GENERAL
$#  imflag       dt0    imform      nsbs       igs     cnstn      form    zero_v
         1       1.0         2         1         1         0         0         0
*CONTROL_IMPLICIT_SOLUTION
$#  nsolvr    ilimit    maxref     dctol     ectol     rctol     lstol    abstol
$        12         1        15    0.0001     0.001  1.000E10       0.91.0000E-20
        12         4        15    0.0001     0.001  1.000E10       0.91.0000E-20
$#   dnorm    diverg     istif   nlprint    nlnorm   d3itctl     cpchk     
         1         1         1         2         2         0         0
$#  arcctl    arcdir    arclen    arcmth    arcdmp    arcpsi    arcalf    arctim
         0         0       0.0         1         2         0         0         0
$#   lsmtd     lsdir      irad      srad      awgt      sred    
         4         2       0.0       0.0       0.0       0.0
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
$
*DATABASE_ELOUT
$#      dt    binary      lcur     ioopt   option1   option2   option3   option4
       1.0         3         0         1         9         0         0         0
*DATABASE_NODFOR
$#      dt    binary      lcur     ioopt     
       1.0         3         0         1
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
*DATABASE_EXTENT_BINARY
$#   neiph     neips    maxint    strflg    sigflg    epsflg    rltflg    engflg
        11        11         1         0         1         1         1         1
$#  cmpflg    ieverp    beamip     dcomp      shge     stssz    n3thdt   ialemat
         0         0         0         1         1         1         2         1
$# nintsld   pkp_sen      sclp     hydro     msscl     therm    intout    nodout
         0         0       1.0         1         0         0                    
$#    dtdt    resplt     neipb     quadr     cubic     
         0         0         0         0         0
*PART
$#                                                                         title

$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         1    &secid      &mid         0         0         0         0         0
*PART
$#                                                                         title

$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         2    &secid      &mid         0         0         0         0         0
*include
2_mesh_PID.inc
*END

