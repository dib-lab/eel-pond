5. Annotating your transcriptome
================================

dammit!

::

    sudo apt-get update
    sudo apt-get install git ruby hmmer unzip build-essential \
        infernal ncbi-blast+ liburi-escape-xs-perl emboss liburi-perl \
        libsm6 libxrender1 libfontconfig1 parallel transdecoder


    sudo gem install crb-blast

    cd
    curl -LO http://busco.ezlab.org/v1/files/BUSCO_v1.22.tar.gz
    tar -xvzf BUSCO_v1.22.tar.gz
    chmod +x BUSCO_v1.22/*.py
    export PATH=$HOME/BUSCO_v1.22:$PATH

    # @CTB put path in env

    pip install -U setuptools
    pip install numpy
    pip install dammit

    cd
    curl -LO http://last.cbrc.jp/last-658.zip
    unzip last-658.zip
    cd last-658
    make
    export PATH=$HOME/last-658/src:$PATH
    export PATH=$HOME/last-658/scripts:$PATH

    libfreetype6-dev?
