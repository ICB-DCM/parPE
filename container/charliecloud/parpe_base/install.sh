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
cd /u18

echo "================= Updating package lists ==================="
apt-get update

echo "================= Adding some global settings ==================="
mv gbl_env.sh /etc/profile.d/
mkdir -p "$HOME/.ssh/"
mv config "$HOME/.ssh/"
mv 90forceyes /etc/apt/apt.conf.d/
touch "$HOME/.ssh/known_hosts"
mkdir -p /etc/drydock

echo "================= Installing basic packages ==================="
apt-get install -q \
  build-essential \
  curl \
  gcc \
  clang \
  gettext \
  libxml2-dev \
  libxslt1-dev \
  make \
  nano \
  openssh-client \
  openssl \
  software-properties-common \
  sudo  \
  texinfo \
  zip \
  unzip \
  wget \
  rsync \
  psmisc \
  netcat-openbsd \
  vim \
  python-lxml

echo "================= Installing Python packages ==================="
apt-get install -q -y \
  python-pip \
  python-dev \
  python3-pip \
  python3-dev

echo "================= Installing Git ==================="
add-apt-repository ppa:git-core/ppa -y
apt-get update
apt-get install -q -y git


#echo "================= Adding JQ 1.5.1 ==================="
#apt-get install -q jq


#echo "================ Adding ansible 2.4.3.0 ===================="
#sudo pip install -q 'ansible==2.4.3.0'

#echo "================ Adding boto 2.48.0 ======================="
#sudo pip install -q 'boto==2.48.0'

#echo "================ Adding boto3 ======================="
#sudo pip install -q 'boto3==1.6.16'

#export PK_VERSION=1.2.2
#echo "================ Adding packer $PK_VERSION  ===================="
#export PK_FILE=packer_"$PK_VERSION"_linux_amd64.zip

#echo "Fetching packer"
#echo "-----------------------------------"
#rm -rf /tmp/packer
#mkdir -p /tmp/packer
#wget -nv https://releases.hashicorp.com/packer/$PK_VERSION/$PK_FILE
#unzip -o $PK_FILE -d /tmp/packer
#sudo chmod +x /tmp/packer/packer
#mv /tmp/packer/packer /usr/bin/packer

#echo "Added packer successfully"
#echo "-----------------------------------"

#echo "================= Adding awscli 1.14.64 ============"
#sudo pip install -q 'awscli==1.14.64'

#echo "================= parPE requirements ============"
apt-get install gfortran libmpich-dev libatlas-base-dev libboost-all-dev libhdf5-dev cmake libceres-dev coinor-libipopt-dev libcpputest-dev gcovr valgrind swig3.0 python3.6 python3-venv hdf5-tools
# for setuptools to find:
ln -s /usr/bin/swig3.0 /usr/bin/swig
python3 -m pip install --upgrade pip
pip3 install -U setuptools pkgconfig wheel

echo "================= Intalling Shippable CLIs ================="

git clone https://github.com/Shippable/node.git nodeRepo
./nodeRepo/shipctl/x86_64/Ubuntu_16.04/install.sh
rm -rf nodeRepo

echo "Installed Shippable CLIs successfully"
echo "-------------------------------------"

echo "================= Cleaning package lists ==================="
apt-get clean
apt-get autoclean
apt-get autoremove
