*PARAMETER
$
$#   prmr1      val1     prmr2      val2     prmr3      val3     prmr4      val4
I NIP              2
R YOUNGMOD         0
R POISSON          0
R DENSITY          0
*SECTION_SHELL
$#   secid    elform      shrf       nip     propt   qr/irid     icomp     setyp
         2        13       1.0      &NIP       1.0         0         0         1
$#      t1        t2        t3        t4      nloc     marea      idof    edgset
       1.0       1.0       1.0       1.0       0.0       0.0       0.0         0
*MAT_PIECEWISE_LINEAR_PLASTICITY
$#     mid        ro         e        pr      sigy      etan      fail      tdel
         2  &DENSITY &YOUNGMOD  &POISSON      1e10       0.01.00000E21       0.0
$#       c         p      lcss      lcsr        vp  
       0.0       0.0         0         0       0.0
$#    eps1      eps2      eps3      eps4      eps5      eps6      eps7      eps8
       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0
$#     es1       es2       es3       es4       es5       es6       es7       es8
       0.0       0.0       0.0       0.0       0.0       0.0       0.0       0.0
