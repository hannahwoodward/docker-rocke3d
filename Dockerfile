# fedora:37 required for GCC 12.x, so as to prevent fortran compiler errors
FROM fedora:37
LABEL maintainer=woodwardsh

ENV HOME=/home/app
ENV MODELDIR=modelE2_planet_1.0
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
WORKDIR ${HOME}

# --- Install dependencies (layer shared with exocam image) ---
RUN dnf install -y gcc gcc-gfortran git nano wget which xz netcdf.x86_64 netcdf-fortran.x86_64 netcdf-devel.x86_64 netcdf-fortran-devel.x86_64 openmpi.x86_64 openmpi-devel.x86_64 'perl(File::Copy)'

# Install other dependencies
RUN dnf install -y ksh task-spooler

# Ensure mpif90 in path & set LD_LIBRARY_PATH for mk_diags
ENV PATH=$PATH:/usr/lib64/openmpi/bin
ENV LD_LIBRARY_PATH=/usr/lib64/openmpi/lib

# --- Download latest model code, SOCRATES (radiation), spectral files ---
RUN wget https://simplex.giss.nasa.gov/snapshots/modelE2_planet_1.0.latest.tgz
RUN wget https://simplex.giss.nasa.gov/snapshots/socrates_1710.tar.xz
RUN wget https://www.giss.nasa.gov/staff/mway/spectral_files.tgz    # 1.6GB
RUN wget https://www.giss.nasa.gov/staff/mway/stellar_spectra.tgz

# --- Extract model code and install some diagnostic tools ---
# NB adding --no-check-certificate due to certificate errors from portal.nccs.nasa.gov during get_input_data
RUN mkdir $MODELDIR && \
    tar -xf modelE2_planet_1.0.latest.tgz -C $MODELDIR --strip-components=1 && \
    sed -i "s|wget \$DATA_PORTAL_URL|wget --no-check-certificate \$DATA_PORTAL_URL|" $MODELDIR/exec/get_input_data && \
    rm modelE2_planet_1.0.latest.tgz && \
    cd $MODELDIR/decks && \
    make config ModelE_Support=$HOME/ModelE_Support

# --- Adjust .modelErc ---
RUN sed -i "s|#\{0,1\} \{0,1\}SOCRATESPATH=.*|SOCRATESPATH=${HOME}\/socrates|" .modelErc && \
    sed -i "s|#\{0,1\} \{0,1\}NETCDFHOME=.*|NETCDFHOME=/usr|" .modelErc && \
    sed -i "s|#\{0,1\} \{0,1\}MPI=.*|MPI=YES|" .modelErc && \
    sed -i "s|#\{0,1\} \{0,1\}MPIDISTR=.*|MPIDISTR=openmpi|" .modelErc && \
    sed -i "s|#\{0,1\} \{0,1\}COMPILER=.*|COMPILER=gfortran|" .modelErc && \
    sed -i "s|#\{0,1\} \{0,1\}UMASK=.*|UMASK=022|" .modelErc && \
    sed -i "s|#\{0,1\} \{0,1\}VERBOSE_OUTPUT=.*|VERBOSE_OUTPUT=NO|" .modelErc && \
    sed -i "s|#\{0,1\} \{0,1\}MPIDIR=.*|MPIDIR=/usr/lib64/openmpi|" .modelErc && \
    sed -i "s|#\{0,1\} \{0,1\}MPILIBDIR=.*|MPILIBDIR=/usr/lib64/openmpi/lib|" .modelErc && \
    sed -i "s|#\{0,1\} \{0,1\}MPIINCLUDEDIR=.*|MPIINCLUDEDIR=/usr/include/openmpi-x86_64|" .modelErc

# --- Change bottom ocean model latitude from 4 to 1 degree ---
RUN sed -i "s|J1O=4|J1O=1|" $MODELDIR/model/OCNDYN.f

# --- Copy over additional model files ---
COPY src/model/ATMDYN_HIGHLOWTEMP.f $MODELDIR/model/ATMDYN_HIGHLOWTEMP.f
COPY src/model/SEAICE_DRVLOWTEMP.f $MODELDIR/model/SEAICE_DRVLOWTEMP.f

# --- Create required dirs ---
RUN mkdir -p ModelE_Support/{exec,huge_space,prod_decks,prod_input_files,prod_runs}

# --- Unpack SOCRATES spectral_files, stellar spectra and code files ---
RUN mkdir -p socrates/{spectral_files,stellar_spectra} && \
    tar -xf socrates_1710.tar.xz -C socrates --strip-components=1 --no-same-owner && \
    tar -xzf spectral_files.tgz -C socrates/spectral_files --strip-components=1 --no-same-owner && \
    tar -xzf stellar_spectra.tgz -C socrates/stellar_spectra --strip-components=1 --no-same-owner && \
    rm socrates_1710.tar.xz spectral_files.tgz stellar_spectra.tgz

# --- Compile SOCRATES ---
RUN sed -i "s|Mk_cmd|Mk_cmd_gfortran|" socrates/make/Makefile && \
    sed -i "s|INCCDF_PATH.*|INCCDF_PATH     = /usr/lib64/gfortran/modules|" socrates/make/Mk_cmd_gfortran && \
    sed -i "s|LIBCDF_PATH.*|LIBCDF_PATH     = /usr/lib64|" socrates/make/Mk_cmd_gfortran && \
    sed -i "s|LIBCDF_NAME.*|LIBCDF_NAME     = netcdff -lnetcdf|" socrates/make/Mk_cmd_gfortran && \
    cd socrates && \
    ./build_code && \
    sed -i "s|LIBCDF_PATH=.*|LIBCDF_PATH=/usr/lib64|" set_rad_env

# --- Add GISS diagnostic tools ---
RUN cd $MODELDIR/model/mk_diags && \
    make diags

# Add mk_diags to path
ENV PATH=$PATH:$HOME/$MODELDIR/model/mk_diags

# Add custom diagnostic tools
COPY src/bin/ ${HOME}/bin
RUN chmod +x ${HOME}/bin/scaleaccm
ENV PATH=$PATH:$HOME/bin

# --- Copy over example decks ---
RUN mkdir ${HOME}/${MODELDIR}/decks/examples
COPY src/decks/ ${HOME}/${MODELDIR}/decks/examples

# --- Copy over test runs ---
COPY src/test-earth.sh test-earth.sh
COPY src/test-planet.sh test-planet.sh

# --- Add entrypoint.sh ---
COPY src/entrypoint.sh $HOME
RUN chmod +x ${HOME}/entrypoint.sh

ENTRYPOINT ["/home/app/entrypoint.sh"]
CMD ["/bin/bash"]
