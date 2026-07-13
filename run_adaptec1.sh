#!/usr/bin/env bash
set -e

PROJECT_ROOT=/home/djw/Desktop/GlobalPlacement_huawei
LIMBO_ROOT=${PROJECT_ROOT}/thrid_party/Limbo
AUX_PATH=${PROJECT_ROOT}/testbench/ispd2005/adaptec3/adaptec3.aux

cd ${PROJECT_ROOT}

LIBBOOKSHELF=$(find ${LIMBO_ROOT} -name "libbookshelfparser*" | head -n 1 || true)

if [ -z "${LIBBOOKSHELF}" ]; then
    echo "[Info] libbookshelfparser not found. Building Limbo first..."
    cd ${LIMBO_ROOT}
    rm -rf build
    mkdir build
    cd build
    cmake ..
    make -j$(nproc)

    cd ${PROJECT_ROOT}
    LIBBOOKSHELF=$(find ${LIMBO_ROOT} -name "libbookshelfparser*" | head -n 1 || true)
fi

if [ -z "${LIBBOOKSHELF}" ]; then
    echo "[Error] Still cannot find libbookshelfparser after building Limbo."
    exit 1
fi

LIMBO_LIB_DIR=$(dirname "${LIBBOOKSHELF}")

echo "[Info] Using Limbo root: ${LIMBO_ROOT}"
echo "[Info] Using Limbo library directory: ${LIMBO_LIB_DIR}"

rm -rf build
mkdir build
cd build

cmake .. \
  -DLIMBO_ROOT=${LIMBO_ROOT} \
  -DCMAKE_LIBRARY_PATH=${LIMBO_LIB_DIR}

make -j$(nproc)

./global_placer \
  --aux ${AUX_PATH} \

