#!/bin/bash

numnum=81920

ZERORK_DATA_DIR=${ZERORK_DATA_DIR:-../../data}

for rctrs in 128 256 512 1024 2048
do
    nm=h2
    mf=hydrogen/h2_v1b_mech.txt
    tf=hydrogen/h2_v1a_therm.txt
    ./gpuMultiOdeFunc.x ${ZERORK_DATA_DIR}/mechanisms/${mf} \
                  ${ZERORK_DATA_DIR}/mechanisms/${tf} \
                  ./${nm}.cklog \
                  ./data/state10_p20T0900.inp \
                  /dev/null $rctrs $((${numnum}/$rctrs))

    nm=dme
    mf=${nm}/dme_24_mech.txt
    tf=${nm}/dme_24_therm.txt
    ./gpuMultiOdeFunc.x ${ZERORK_DATA_DIR}/mechanisms/${mf} \
                  ${ZERORK_DATA_DIR}/mechanisms/${tf} \
                  ./${nm}.cklog \
                  ./data/state159_p20T0900.inp \
                  /dev/null $rctrs $((${numnum}/$rctrs))

    nm=n-heptane_reduced
    mf=${nm}/heptanesymp159_mec.txt
    tf=${nm}/heptanesymp_therm.txt
    ./gpuMultiOdeFunc.x ${ZERORK_DATA_DIR}/mechanisms/${mf} \
                  ${ZERORK_DATA_DIR}/mechanisms/${tf} \
                  ./${nm}.cklog \
                  ./data/state159_p20T0900.inp \
                  /dev/null $rctrs $((${numnum}/$rctrs))

    nm=iso-octane
    mf=iso-octane/species874/ic8_ver3_mech.txt
    tf=iso-octane/species874/prf_v3_therm_dat.txt
    ./gpuMultiOdeFunc.x ${ZERORK_DATA_DIR}/mechanisms/${mf} \
                  ${ZERORK_DATA_DIR}/mechanisms/${tf} \
                  ./${nm}.cklog \
                  ./data/state874_p20T0900.inp \
                  /dev/null $rctrs $((${numnum}/$rctrs))

    nm=plog
    mf=plog/plog_test.mech
    tf=plog/plog_test.therm
    ./gpuMultiOdeFunc.x ${ZERORK_DATA_DIR}/mechanisms/${mf} \
                  ${ZERORK_DATA_DIR}/mechanisms/${tf} \
                  ./${nm}.cklog \
                  ./data/state10_p20T0900.inp \
                  /dev/null $rctrs $((${numnum}/$rctrs))

    nm=ford
    mf=real_order/ford_propane.mech
    tf=ideal/const_specific_heat.therm
    ./gpuMultiOdeFunc.x ${ZERORK_DATA_DIR}/mechanisms/${mf} \
                  ${ZERORK_DATA_DIR}/mechanisms/${tf} \
                  ./${nm}.cklog \
                  ./data/state10_p20T0900.inp \
                  /dev/null $rctrs $((${numnum}/$rctrs))

    nm=rord
    mf=real_order/rord_propane.mech
    tf=ideal/const_specific_heat.therm
    ./gpuMultiOdeFunc.x ${ZERORK_DATA_DIR}/mechanisms/${mf} \
                  ${ZERORK_DATA_DIR}/mechanisms/${tf} \
                  ./${nm}.cklog \
                  ./data/state10_p20T0900.inp \
                  /dev/null $rctrs $((${numnum}/$rctrs))



done

