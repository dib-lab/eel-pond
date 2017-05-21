=================================================================
Installation and configuration: Docker for Continuous Integration
=================================================================

(This page will some day include links to the various services and
repos we are using to manage the CI service for eel-pond.)

.. shell start

.. ::

   set -x
   set -e

Install software
----------------

On the new machine, run the following commands to update the base
software:
::

    apt-get update
    apt-get -y install screen git curl gcc make g++ python-dev unzip \
            default-jre pkg-config libncurses5-dev r-base-core r-cran-gplots \
            python-matplotlib python-pip python-virtualenv sysstat fastqc \
            trimmomatic bowtie samtools blast2 wget bowtie2 openjdk-8-jre
    apt-get -y install sudo
.. ::

Install `khmer <http://khmer.readthedocs.org>`__ from its source code.
::

    cd ~/
    python2.7 -m virtualenv pondenv
    source pondenv/bin/activate
    cd pondenv
    pip install -U setuptools
    git clone https://github.com/dib-lab/khmer.git
    cd khmer
    make install

Installing Trinity
~~~~~~~~~~~~~~~~~~

To install Trinity:
::

    cd ${HOME}

    wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.3.2.tar.gz \
     -O trinity.tar.gz
    tar xzf trinity.tar.gz
    cd trinityrnaseq*/
    make |& tee trinity-build.log

Assuming it succeeds, modify the path appropriately in your virtualenv
activation setup:
::

    echo export PATH=$PATH:$(pwd) >> ~/pondenv/bin/activate

You will also need to set the default Java version to 1.8
::

   sudo update-alternatives --set java /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java

Next: :doc:`1-quality`
