setup D0RunII p20.17.03
setup lhapdf v5_8_1 -q GCC3_4_4
cp $LHAPDF_DIR/bin/lhapdf-config .
g77 -pthread -UD0  -fexceptions -fPIC -fdollar-ok -fno-automatic \
-fno-second-underscore -ffixed-line-length-132 -fno-globals \
-fno-strict-aliasing -B${LINKTIME_DIR}/usr/X11R6/lib \
-B${RUNTIME_DIR}/usr/lib -B${LINKTIME_DIR}/lib \
-B${LINKTIME_DIR}/usr/lib -Wl,-rpath-link \
-Wl,${LINKTIME_DIR}/usr/X11R6/lib -Wl,-rpath-link \
-Wl,${RUNTIME_DIR}/usr/lib -Wl,-rpath-link -Wl,${LINKTIME_DIR}/lib \
-Wl,-rpath-link -Wl,${LINKTIME_DIR}/usr/lib \
VecBosProd.f \
$PYTHIA_DIR/src/dummies/fh*.f $PYTHIA_DIR/src/dummies/up*.f \
$PYTHIA_DIR/src/dummies/ssmssm.f $PYTHIA_DIR/src/dummies/sugra.f \
$PYTHIA_DIR/src/dummies/visaje.f \
$PYTHIA_DIR/lib/libpythia.a $LHAPDF_DIR/lib/libLHAPDF.a -lstdc++ 
