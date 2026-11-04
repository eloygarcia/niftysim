#FROM nvidia/cuda:12.6.1-devel-ubuntu24.04
FROM nvidia/cuda:13.4.3-devel-ubuntu24.04
# FROM nvidia/cuda:13.4.3-devel-ubuntu26.04 ## problems with cpp compiler

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC
# RUN apt-get install -y tzdata

RUN set -ex;		 
RUN apt-get update; 	
RUN apt-get install -y build-essential gcc g++ wget git curl
RUN apt-get install -y libx11-dev libxt-dev libxext-dev libxrender-dev
RUN apt-get install -y libgl1-mesa-dev mesa-common-dev libglu1-mesa-dev freeglut3-dev
RUN apt-get install -y libssl-dev openssl;
RUN apt-get install -y tcl tk

RUN apt-get install -y cmake
RUN apt-get install -y qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools
RUN apt-get install -y vtk9 libvtk9-dev

## copy niftysim
COPY ./ /home/niftysim/
RUN mkdir -p /home/niftysim/build/
RUN cd /home/niftysim/build/ 

RUN cmake -S /home/niftysim \
          -B /home/niftysim/build \
          -G "Unix Makefiles"  \
          -DUSE_VIZ=ON   \
          -DUSE_CUDA=ON   \
          -DUSE_NAMESPACE_STD=ON   \
          -DUSE_GPU_GP_CONTACT=OFF   \
          -DCUDA_SDK_COMMON_INCLUDE_DIR="/home/niftysim/cuda samples/Common" 

RUN make

