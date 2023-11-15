
# FROM ubuntu:18.04
FROM nvidia/cuda:11.7.1-devel-ubuntu20.04 
#nvcc -V

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC
# RUN apt-get install -y tzdata

RUN set -ex;		\ 
apt-get update; 	\
apt-get install -y g++ curl libzmq3-dev wget git; \
apt-get install -y libx11-dev libxt-dev xserver-xorg-dev xorg-dev ;
    
# cnake
 RUN wget http://www.cmake.org/files/v3.12/cmake-3.12.0.tar.gz \
    && tar xzf cmake-3.12.0.tar.gz
# old version 2.8.12.2
# RUN wget http://www.cmake.org/files/v2.8/cmake-2.8.12.2.tar.gz \
#    && tar xzf cmake-2.8.12.2.tar.gz
    
RUN cd cmake-3.12.0 \
    && ./configure --prefix=/usr/local \
    && make \
    && make install


# Install OpenGL
# Debian, Ubuntu
# https://en.wikibooks.org/wiki/OpenGL_Programming/Installation/Linux
RUN apt-get update && apt-get install --yes build-essential libgl1-mesa-dev

# cmake 
#RUN apt-get update && apt-get -y install cmake protobuf-compiler cmake-gui

# Download & build Tcls
# https://www.tcl.tk/doc/howto/compile.html#unix
RUN wget http://prdownloads.sourceforge.net/tcl/tcl8.6.6-src.tar.gz && tar -zxvf tcl8.6.6-src.tar.gz
RUN cd tcl8.6.6/unix && ./configure && make && make install
#
# Download & build Tk
# https://www.tcl.tk/doc/howto/compile.html
RUN wget http://prdownloads.sourceforge.net/tcl/tk8.6.6-src.tar.gz && tar -zxvf tk8.6.6-src.tar.gz
RUN cd tk8.6.6/unix && ./configure && make && make install
 
# vtk 
# Download & extract VTK
#RUN mkdir -p
RUN wget http://www.vtk.org/files/release/8.2/VTK-8.2.0.tar.gz && tar -zxvf VTK-8.2.0.tar.gz
#  old version 7.1.1

RUN mkdir -p /home/vtk-build2
RUN cd /home/vtk-build2/ && cmake \
   -DCMAKE_BUILD_TYPE:STRING=Release \
   -DBUILD_TESTING:BOOL=OFF \
   -DVTK_WRAP_PYTHON:BOOL=OFF \
   -DVTK_WRAP_PYTHON_SIP:BOOL=OFF \
   -DVTK_WRAP_TCL:BOOL=ON \
   -DVTK_PYTHON_VERSION:STRING=3 \
   -DVTK_USE_TK:BOOL=ON \
   /VTK-8.2.0/ && make

## copy niftysim
COPY ./ $HOME/home/niftysim/
RUN mkdir -p /home/niftysim/release/
RUN cd /home/niftysim/release/

#COPY ~/Escritorio/Oliver/Libraries/cuda-samples/  $HOME/cuda-samples/
#RUN cd /home/
#RUN git clone git@github.com:NVIDIA/cuda-samples.git

RUN cd /home/niftysim/release/ && cmake -G "Unix Makefiles" -DUSE_VIZ:BOOL=ON \-DUSE_CUDA:BOOL=ON \-DUSE_NAMESPACE_STD:BOOL=ON / -DVTK_DIR:STRING=/home/vtk-build2/ \-DCUDA_SDK_COMMON_INCLUDE_DIR:STRING=/home/niftysim/cuda-samples/Common/ \-DUSE_NAMESPACE_STD:BOOL=ON \-DUSE_GPU_GP_CONTACT:BOOL=ON \/home/niftysim  && make

#RUN make 


# vtk 
#RUN apt-get update && apt-get -y install vtk7

### COPIED
### https://github.com/lukin0110/docker-vtk-python/blob/master/Dockerfile
###################################################################################################################
# Download sources & setup supporting libraries that are needed to build VTK
###################################################################################################################
# Download & extract VTK
# RUN wget http://www.vtk.org/files/release/7.1/VTK-7.1.1.tar.gz && tar -zxvf VTK-7.1.1.tar.gz

# Download, extract & build CMake
# http://www.vtk.org/Wiki/VTK/Configure_and_Build
# RUN wget http://www.cmake.org/files/v2.8/cmake-2.8.12.2.tar.gz \
#    && tar xzf cmake-2.8.12.2.tar.gz
#
#RUN cd cmake-2.8.12.2 \
#    && ./configure --prefix=/usr/local \
#    && make 
 #   && make install



# Download & build Tcls
# https://www.tcl.tk/doc/howto/compile.html#unix
# RUN wget http://prdownloads.sourceforge.net/tcl/tcl8.6.6-src.tar.gz && tar -zxvf tcl8.6.6-src.tar.gz
# RUN cd tcl8.6.6/unix && ./configure && make && make install
#
# Download & build Tk
# https://www.tcl.tk/doc/howto/compile.html
# RUN wget http://prdownloads.sourceforge.net/tcl/tk8.6.6-src.tar.gz && tar -zxvf tk8.6.6-src.tar.gz
# RUN cd tk8.6.6/unix && ./configure && make && make install
###################################################################################################################
# /end setup
###################################################################################################################

###################################################################################################################
# Building VTK with python interfaces
# http://ghoshbishakh.github.io/blog/blogpost/2016/03/05/buid-vtk.html
###################################################################################################################
# RUN mkdir /vtk-build2
# RUN cd /vtk-build2/ && cmake \
#   -DCMAKE_BUILD_TYPE:STRING=Release \
#   -DBUILD_TESTING:BOOL=OFF \
#   -DVTK_WRAP_PYTHON:BOOL=OFF \
#   -DVTK_WRAP_PYTHON_SIP:BOOL=OFF \
#   -DVTK_WRAP_TCL:BOOL=ON \
#   -DVTK_PYTHON_VERSION:STRING=3 \
#   -DVTK_USE_TK:BOOL=ON \
#   /tmpbuild/VTK-7.0.0
# 
# Build VTK
# RUN cd /vtk-build2/ && make
#
# Now install the python bindings
#  RUN cd /vtk-build2/Wrapping/Python && make && make install

# Set environment variable to add the VTK libs to the Shared Libraries
# http://tldp.org/HOWTO/Program-Library-HOWTO/shared-libraries.html
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib:/vtk-build2/lib
# ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/usr/local/lib:/vtk-build2/lib
###################################################################################################################
# /end VTK Build
###################################################################################################################    

        
## installing cuda toolkit ubuntun 22.04 x86_64
# RUN wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-ubuntu2204.pin
# RUN mv cuda-ubuntu2204.pin /etc/apt/preferences.d/cuda-repository-pin-600
#RUN wget https://developer.download.nvidia.com/compute/cuda/12.0.0/local_installers/cuda-repo-ubuntu2204-12-0-local_12.0.0-525.60.13-1_amd64.debs
#RUN dpkg -i cuda-repo-ubuntu2204-12-0-local_12.0.0-525.60.13-1_amd64.deb
#RUN cp /var/cuda-repo-ubuntu2204-12-0-local/cuda-*-keyring.gpg /usr/share/keyrings/
#RUN apt-get update
#RUN apt-get -y install cuda
 

