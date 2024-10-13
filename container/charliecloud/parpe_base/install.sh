#!/bin/bash -e

export DEBIAN_FRONTEND=noninteractive

apt-get clean && apt-get update
apt-get install -y apt-utils

echo "================ Installing locales ======================="
apt-get install -q locales

dpkg-divert --local --rename --add /sbin/initctl
locale-gen en_US en_US.UTF-8
dpkg-reconfigure locales

echo "HOME=$HOME"

echo "================= parPE requirements ============"
# using openmpi coming with libboost-all-dev instead of libmpich-dev
apt-get install -q -y \
  clang \
  cmake \
  curl \
  coinor-libipopt-dev \
  gfortran \
  git \
  hdf5-tools \
  libatlas-base-dev \
  libboost-all-dev \
  libceres-dev \
  libhdf5-dev \
  libomp-dev \
  nano \
  python3-dev \
  python3-numpy \
  python3-pip \
  python3-venv \
  libspdlog-dev \
  swig \
  unzip \
  wget

echo "================= Cleaning package lists ==================="
apt-get clean
apt-get autoclean
apt-get autoremove
