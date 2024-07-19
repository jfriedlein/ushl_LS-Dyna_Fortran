$# LS-DYNA Keyword file created by LS-PrePost(R) V4.7.20 -24Aug2020
$# Created on May-4-2022 (11:32:16)
*PARAMETER
$
$#   prmr1      val1     prmr2      val2     prmr3      val3     prmr4      val4
R tmax    100
R pmax    1
$
*CONTROL_IMPLICIT_SOLUTION
$#  nsolvr    ilimit    maxref     dctol     ectol     rctol     lstol    abstol
        12         1        20     0.001      0.011.00000E10       0.91.0000E-20
$#   dnorm    diverg     istif   nlprint    nlnorm   d3itctl     cpchk     
         1         1         1         2         2         0         0
$#  arcctl    arcdir    arclen    arcmth    arcdmp    arcpsi    arcalf    arctim
         0         0       0.0         1         2         0         0         0
$  &arcnode         2       0.0         3         1       0.1      -0.1   &arctim
$#   lsmtd     lsdir      irad      srad      awgt      sred    
         4         2       0.0       0.0       0.0       0.0
*CONSTRAINED_GLOBAL
$#      tc        rc       dir         x         y         z       tol   
         1         0         1       0.0       0.0       0.0       0.0
         7         0         2       0.0     -28.0       0.0       0.0
         7,        0,        2,      0.0,<28.+shThck>,   0.0,      0.0
         7         0         1     195.0       0.0       0.0       0.0
$
*DATABASE_CROSS_SECTION_PLANE_ID
$#    csid                                                                 title
         2tool upper
$#    psid       xct       yct       zct       xch       ych       zch    radius
         2       0.0      18.0       0.0       0.0       0.0      0.0     200.0
$#    xhev      yhev      zhev      lenl      lenm        id     itype     
       0.0       0.0       0.0       0.0       0.0         0         0
*DATABASE_CROSS_SECTION_PLANE_ID
$#    csid                                                                 title
         3tool lower
$#    psid       xct       yct       zct       xch       ych       zch    radius
         3       0.0     -12.0       0.0       0.0       0.0       0.0     200.0
$#    xhev      yhev      zhev      lenl      lenm        id     itype     
       0.0       0.0       0.0       0.0       0.0         0         0
*CONTACT_2D_AUTOMATIC_SURFACE_TO_SURFACE_MORTAR_TITLE
$#     cid                                                                 title
         2Matrize_Basisblech
$#    sids      sidm     sfact      freq        fs        fd        dc     membs
         1         2       1.0         1      0.15       0.0       0.0         6
$#  tbirth    tdeath       sos       som       nds       ndm       cof      init
       0.0       0.0       1.0       1.0         0         0         0         0
$#      vc       vdc       ipf     slide    istiff   tiedgap    igapcl    tietyp
       0.0       0.0         0         0         0       0.0         0         0
$#  sldsos    sldsom      
       0.0       0.0
*SET_PART_LIST_TITLE
sheet
$#     sid       da1       da2       da3       da4    solver      
         1       0.0       0.0       0.0       0.0MECH
$#    pid1      pid2      pid3      pid4      pid5      pid6      pid7      pid8
         1         0         0         0         0         0         0         0
*SET_PART_LIST_TITLE
Tool upper
$#     sid       da1       da2       da3       da4    solver      
         2       0.0       0.0       0.0       0.0MECH
$#    pid1      pid2      pid3      pid4      pid5      pid6      pid7      pid8
         2         0         0         0         0         0         0         0
*CONTACT_2D_AUTOMATIC_SURFACE_TO_SURFACE_MORTAR_TITLE
$#     cid                                                                 title
         2Matrize_Basisblech
$#    sids      sidm     sfact      freq        fs        fd        dc     membs
         1         3       1.0         1      0.15       0.0       0.0         6
$#  tbirth    tdeath       sos       som       nds       ndm       cof      init
       0.0       0.0       1.0       1.0         0         0         0         0
$#      vc       vdc       ipf     slide    istiff   tiedgap    igapcl    tietyp
       0.0       0.0         0         0         0       0.0         0         0
$#  sldsos    sldsom      
       0.0       0.0
*SET_PART_LIST_TITLE
Tool lower
$#     sid       da1       da2       da3       da4    solver      
         3       0.0       0.0       0.0       0.0MECH
$#    pid1      pid2      pid3      pid4      pid5      pid6      pid7      pid8
         3         0         0         0         0         0         0         0
$
*CONTROL_TERMINATION
$#  endtim    endcyc     dtmin    endeng    endmas     nosol     
&tmax              0       0.0       0.01.000000E8         0
*DEFINE_CURVE_TITLE
pressure
$#    lcid      sidr       sfa       sfo      offa      offo    dattyp     lcint
         1         0       1.0       1.0       0.0       0.0         0         0
$#                a1                  o1  
                 0.0                 0.0
               &tmax              -&pmax
*COMMENT  *LOAD_SEGMENT_SET
$#    ssid      lcid        sf        at  
         1         1       1.0       0.0
*USER_LOADING
$     LCID       PID       AXI     UNDEF
         2         2         1         0
*MAT_ELASTIC_TITLE
Werkzeug elastisch
$#     mid        ro         e        pr        da        db  not used        
        117.83000E-9  210000.0       0.3       0.0       0.0       0.0
*PART
$#                                                                         title
sheet
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         1    &secid      &MID         0         0         0         0         0
*PART
$#                                                                         title
sheet
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
        10    &secid      &MID         0         0         0         0         0
*PART
$#                                                                         title
tool upper
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         2        11        11         0         0         0         0         0
*PART
$#                                                                         title
tool lower
$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid
         3        11        11         0         0         0         0         0
*SECTION_SHELL_TITLE
rod - axisym
$#   secid    elform      shrf       nip     propt   qr/irid     icomp     setyp
        11        15       1.0         4       3.0         0         0         1
$#      t1        t2        t3        t4      nloc     marea      idof    edgset
       0.0       0.0       0.0       0.0       1.0       0.0       0.0         0
*END
