$# LS-DYNA Keyword file created by LS-PrePost(R) V4.6.15 - 16Aug2019
$# Created on Feb-18-2021 (16:40:51)
*KEYWORD
*PARAMETER
$
$#   prmr1      val1     prmr2      val2     prmr3      val3     prmr4      val4
R tmax    100
I NIP     2
I secid   2
I MID     2
*CONTROL_TERMINATION
$#  endtim    endcyc     dtmin    endeng    endmas     nosol     
&tmax              0       0.0       0.01.000000E8         0
*DATABASE_CROSS_SECTION_PLANE_ID
$#    csid                                                                 title
         2                                                                      
$#    psid       xct       yct       zct       xch       ych       zch    radius
         0      -1.0      26.5      -1.0      -1.0      28.0      -1.0       0.0
$#    xhev      yhev      zhev      lenl      lenm        id     itype     
      -1.0       0.0       2.0       2.0       8.0         0         0
*BOUNDARY_PRESCRIBED_MOTION_SET_ID
$#      id                                                               heading
         1BC_uy
$#    nsid       dof       vad      lcid        sf       vid     death     birth
         1         2         2         1       1.0         01.00000E28       0.0
*CONSTRAINED_GLOBAL
$#      tc        rc       dir         x         y         z       tol  
         1         0         1       0.0       0.0       0.0       0.0
         2         0         2       0.0       0.0       0.0       0.0
*PART
$#                                                                         title
                                                                                
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         1    &secid      &MID         0         0         0         0         0
$
*DEFINE_CURVE_TITLE
Curve_uy
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
         1         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
                 0.0                 0.0
          &tmax                      8.0
*END
