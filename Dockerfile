FROM neurodebian:stretch
MAINTAINER <switt4@uwo.ca>

ARG DEBIAN_FRONTEND="noninteractive"

ENV LANG="en_US.UTF-8" \
    LC_ALL="en_US.UTF-8" \
    ND_ENTRYPOINT="/neurodocker/startup.sh"
RUN export ND_ENTRYPOINT="/neurodocker/startup.sh" \
    && apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
           apt-utils \
           bzip2 \
           ca-certificates \
           curl \
           locales \
           unzip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen \
    && dpkg-reconfigure --frontend=noninteractive locales \
    && update-locale LANG="en_US.UTF-8" \
    && chmod 777 /opt && chmod a+s /opt \
    && mkdir -p /neurodocker \
    && if [ ! -f "$ND_ENTRYPOINT" ]; then \
         echo '#!/usr/bin/env bash' >> "$ND_ENTRYPOINT" \
    &&   echo 'set -e' >> "$ND_ENTRYPOINT" \
    &&   echo 'export USER="${USER:=`whoami`}"' >> "$ND_ENTRYPOINT" \
    &&   echo 'if [ -n "$1" ]; then "$@"; else /usr/bin/env bash; fi' >> "$ND_ENTRYPOINT"; \
    fi \
    && chmod -R 777 /neurodocker && chmod a+s /neurodocker

ENTRYPOINT ["/neurodocker/startup.sh"]

ENV FSLDIR="/opt/fsl-5.0.11" \
    PATH="/opt/fsl-5.0.11/bin:$PATH"
ENV FSLOUTPUTTYPE NIFTI_GZ
RUN apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
           bc \
           dc \
           file \
           libfontconfig1 \
           libfreetype6 \
           libgl1-mesa-dev \
           libgl1-mesa-dri \
           libglu1-mesa-dev \
           libgomp1 \
           libice6 \
           libxcursor1 \
           libxft2 \
           libxinerama1 \
           libxrandr2 \
           libxrender1 \
           libxt6 \
           sudo \
           wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && echo "Downloading FSL ..." \
    && mkdir -p /opt/fsl-5.0.11 \
    && curl -fsSL --retry 5 https://fsl.fmrib.ox.ac.uk/fsldownloads/fsl-5.0.11-centos6_64.tar.gz \
    | tar -xz -C /opt/fsl-5.0.11 --strip-components 1 \
    && sed -i '$iecho Some packages in this Docker container are non-free' $ND_ENTRYPOINT \
    && sed -i '$iecho If you are considering commercial use of this container, please consult the relevant license:' $ND_ENTRYPOINT \
    && sed -i '$iecho https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Licence' $ND_ENTRYPOINT \
    && sed -i '$isource $FSLDIR/etc/fslconf/fsl.sh' $ND_ENTRYPOINT \
    && echo "Installing FSL conda environment ..." \
    && bash /opt/fsl-5.0.11/etc/fslconf/fslpython_install.sh -f /opt/fsl-5.0.11

ENV CONDA_DIR="/opt/miniconda-latest" \
    PATH="/opt/miniconda-latest/bin:$PATH"
RUN export PATH="/opt/miniconda-latest/bin:$PATH" \
    && echo "Downloading Miniconda installer ..." \
    && conda_installer="/tmp/miniconda.sh" \
    && curl -fsSL --retry 5 -o "$conda_installer" https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash "$conda_installer" -b -p /opt/miniconda-latest \
    && rm -f "$conda_installer" \
    && conda update -yq -nbase conda \
    && conda config --system --prepend channels conda-forge \
    && conda config --system --set auto_update_conda false \
    && conda config --system --set show_channel_urls true \
    && sync && conda clean --all && sync \
    && conda install -y -q --name base \
           "python=3.7.1" \
           "numpy=1.15.4" \
           "pandas=0.23.4" \
           "scipy=1.1.0" \
           "nilearn=0.5.2" \
           "matplotlib=2.2.2" \
           "scikit-learn=0.19.1" \
           "traits=4.6.0" \
	   "fslpy=2.6.0" \
    && sync && conda clean --all && sync

RUN conda install -c conda-forge fslpy
RUN pip install fslpy

#RUN apt-get update -qq \
#    && apt-get install -y -q --no-install-recommends \
#           dirmngr \
#           gnupg2 \
#    && apt-get clean \
#    && rm -rf /var/lib/apt/lists/* \
#    && curl -fsSL --retry 5 http://neuro.debian.net/lists/stretch.us-nh.full \
#       > /etc/apt/sources.list.d/neurodebian.sources.list \
#    && curl -sSL https://dl.dropbox.com/s/zxs209o955q6vkg/neurodebian.gpg | apt-key add - \
#    && (apt-key adv --refresh-keys --keyserver hkp://pool.sks-keyservers.net:80 0xA5D32F012649A5A9 || true) \
#    && apt-get -qq update

#RUN apt-get update -qq \
#    && apt-get install -y -q --no-install-recommends \
#           python3-fsl \
#    && apt-get clean \
#    && rm -rf /var/lib/apt/lists/*

RUN apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
           git \
           gcc \
           pigz \
           liblzma-dev \
           libc-dev \
           git-annex-standalone \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN python -c "from matplotlib import font_manager" && \
    sed -i 's/\(backend *: \).*$/\1Agg/g' $( python -c "import matplotlib; print(matplotlib.matplotlib_fname())" )

RUN mkdir -p /opt/templates && cd /opt/templates && wget https://www.dropbox.com/s/k8k5a5ba94p65ki/tpl-MNI152NLin2009cAsym_res-02_atlas-HOCPA_desc-th0_dseg.nii.gz
ENV TEMPLATEDIR="/opt/templates" \
    PATH="/opt/templates:$PATH"

ADD confounder.py /opt/confounder.py

ENTRYPOINT ["/opt/confounder.py"]