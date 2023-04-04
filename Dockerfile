FROM fedora:latest

ENV HOME=/home/app
ENV MODELDIR=modelE2_planet_1.0
WORKDIR ${HOME}

# --- Install dependencies ---
RUN yum -y install gcc gcc-gfortran git nano wget xz netcdf.x86_64 netcdf-fortran.x86_64 netcdf-devel.x86_64 netcdf-fortran-devel.x86_64 openmpi.x86_64 openmpi-devel.x86_64 && \
    echo "module load mpi/openmpi-x86_64" >> /root/.bashrc && \
    dnf install -y 'perl(File::Copy)'

# --- Download latest model code, SOCRATES (radiation), spectral files ---
RUN wget https://simplex.giss.nasa.gov/snapshots/modelE2_planet_1.0.latest.tgz
RUN wget https://simplex.giss.nasa.gov/snapshots/socrates_1710.tar.xz
RUN wget https://www.giss.nasa.gov/staff/mway/spectral_files.tgz    # 1.6GB
RUN wget https://www.giss.nasa.gov/staff/mway/stellar_spectra.tgz

# --- Extract model code ---
RUN mkdir $MODELDIR && \
    tar -xf modelE2_planet_1.0.latest.tgz -C $MODELDIR --strip-components=1 && \
    rm modelE2_planet_1.0.latest.tgz

# --- Prepare computer for installing some diagnostic tools and running the model ---
RUN cd $MODELDIR/decks && \
    make config ModelE_Support=$HOME/ModelE_Support

# --- Adjust .modelErc ---
RUN sed -i "s|#\{0,1\} \{0,1\}SOCRATESPATH=.*|SOCRATESPATH=${HOME}\/ModelE_Support\/socrates|" .modelErc && \
    sed -i "s|#\{0,1\} \{0,1\}NETCDFHOME=.*|NETCDFHOME=/usr|" .modelErc && \
    sed -i "s|#\{0,1\} \{0,1\}MPI=.*|MPI=YES|" .modelErc && \
    sed -i "s|#\{0,1\} \{0,1\}MPIDISTR=.*|MPIDISTR=openmpi|" .modelErc && \
    sed -i "s|#\{0,1\} \{0,1\}COMPILER=.*|COMPILER=gfortran|" .modelErc && \
    sed -i "s|#\{0,1\} \{0,1\}UMASK=.*|UMASK=022|" .modelErc && \
    sed -i "s|#\{0,1\} \{0,1\}VERBOSE_OUTPUT=.*|VERBOSE_OUTPUT=NO|" .modelErc && \
    sed -i "s|#\{0,1\} \{0,1\}MPIDIR=.*|MPIDIR=/usr/lib64/openmpi|" .modelErc && \
    sed -i "s|#\{0,1\} \{0,1\}MPILIBDIR=.*|MPILIBDIR=/usr/lib64/openmpi/lib|" .modelErc && \
    sed -i "s|#\{0,1\} \{0,1\}MPIINCLUDEDIR=.*|MPIINCLUDEDIR=/usr/include/openmpi-x86_64|" .modelErc

# --- Create required dirs ---
RUN mkdir -p ModelE_Support/{exec,huge_space,prod_decks,prod_input_files,prod_runs}

# --- Unpack SOCRATES spectral_files, stellar spectra and code files ---
RUN mkdir -p ModelE_Support/socrates/{spectral_files,stellar_spectra} && \
    tar -xf socrates_1710.tar.xz -C ModelE_Support/socrates --strip-components=1 && \
    tar -xzf spectral_files.tgz -C ModelE_Support/socrates/spectral_files --strip-components=1 && \
    tar -xzf stellar_spectra.tgz -C ModelE_Support/socrates/stellar_spectra --strip-components=1 && \
    rm socrates_1710.tar.xz spectral_files.tgz stellar_spectra.tgz

# --- Create GISS diagnostic tools ---
RUN cd $MODELDIR/model/mk_diags && \
    make diags && \
    echo "export PATH=$PATH:$HOME/${MODELDIR}/model/mk_diags" >> /root/.bashrc
