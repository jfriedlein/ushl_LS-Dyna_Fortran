$# LS-DYNA Keyword file created by LS-PrePost(R) V4.8.30 - 20Apr2022
$# Created on Jul-19-2023 (07:59:00)
*KEYWORD   
*PARAMETER
$
$#   prmr1      val1     prmr2      val2     prmr3      val3     prmr4      val4
I NQP              4
I HG_onoff         0
$
R BULKMOD          0
R SHEARMOD         0
R YOUNGMOD         0
R POISSON          0
R DENSITY          0
$
*MAT_ELASTIC
$#     mid        ro         e        pr        da        db  not used        
         1  &DENSITY &YOUNGMOD  &POISSON       0.0       0.0       0.0
*CONSTRAINED_GLOBAL
$#      tc        rc       dir         x         y         z       tol   
         3         7         3       0.0       0.0       0.0       0.0
*SECTION_SHELL
$#   secid    elform      shrf       nip     propt   qr/irid     icomp     setyp
         1       101       1.0         2       1.0         0         0         1
$#      t1        t2        t3        t4      nloc     marea      idof    edgset
       1.0       1.0       1.0       1.0       0.0       0.0       0.0         0
$#    nipp     nxdof      iunf      ihgf      itaj       lmc      nhsv      iloc
         0         0         0         0         0        16         7         1
$#      p1        p2        p3        p4        p5        p6        p7        p8
$   ELTYPE       DIM       NQP   HGSCALE     DEBUG
         1         2      &NQP &HG_onoff         2
$#      p1        p2        p3        p4        p5        p6        p7        p8
         0  &BULKMOD &SHEARMOD         0         0         0         0         0
*DATABASE_EXTENT_BINARY
$#   neiph     neips    maxint    strflg    sigflg    epsflg    rltflg    engflg
         7         7         3         0         1         1         1         1
$#  cmpflg    ieverp    beamip     dcomp      shge     stssz    n3thdt   ialemat
         0         0         0         1         1         1         2         1
$# nintsld   pkp_sen      sclp     hydro     msscl     therm    intout    nodout
         0         0       1.0         0         0         0STRESS
$ HYDRO=1: output the "internal energy density" (incompatible with paraview)
$         0         0       1.0         1         0         0                    
$#    dtdt    resplt     neipb     quadr     cubic     
         0         0         0         0         0
*END
