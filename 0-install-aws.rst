===================================================
Installation and configuration: Amazon Web Services
===================================================

Boot up an m4.xlarge machine from Amazon Web Services running Ubuntu
15.10 LTS (e.g. us-west AMI ami-05384865); increase the root volume
size to ~100 GB.  The m4.xlarge machines have 16 GB of RAM, and 4
CPUs, and will be enough to complete the assembly of the Nematostella
data set. If you are using your own data, be aware of your space
requirements and obtain an appropriately sized machine ("instance")
and storage ("volume").

.. shell start

.. ::

   set -x
   set -e

Install software
----------------

On the new machine, run the following commands to update the base
software:
::

   sudo apt-get update && \
   sudo apt-get -y install screen git curl gcc make g++ python-dev unzip \
            default-jre pkg-config libncurses5-dev r-base-core r-cran-gplots \
            python-matplotlib python-pip python-virtualenv sysstat fastqc \
            trimmomatic bowtie samtools blast2 wget bowtie2 openjdk-8-jre \
            hmmer

Install `khmer <http://khmer.readthedocs.org>`__ from its source code.
::

   cd ~/
   python2.7 -m virtualenv pondenv
   source pondenv/bin/activate
   cd pondenv
   pip install -U setuptools
   git clone --branch v2.0 https://github.com/dib-lab/khmer.git
   cd khmer
   make install

The use of ``virtualenv`` allows us to install Python software without having
root access. If you come back to this protocol in a different terminal session
you will need to run::

        source ~/pondenv/bin/activate

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

You will also need to set the default Java version to 1.8::

  sudo update-alternatives --set java /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java


Install transrate
-----------------

We use `transrate <http://hibberdlab.com/transrate/getting_started.html>`__
to evaluate assemblies.  Install!
::

  cd
  curl -LO https://bintray.com/artifact/download/blahah/generic/transrate-1.0.3-linux-x86_64.tar.gz
  tar -zxf transrate-1.0.3-linux-x86_64.tar.gz
  echo 'export PATH=$PATH:"$HOME/transrate-1.0.3-linux-x86_64"' >> ~/pondenv/bin/activate
  source ~/.profile
  curl -LO ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.3.0/ncbi-blast-2.3.0+-x64-linux.tar.gz
  tar -zxf ncbi-blast-2.3.0+-x64-linux.tar.gz
  echo 'export PATH="$HOME/ncbi-blast-2.3.0+/bin:$PATH"' >> ~/pondenv/bin/activate
  source ~/.profile

Install busco
-------------

Install stuff:

::

  cd
  git clone https://gitlab.com/ezlab/busco.git
  cd busco
  echo "export PATH=$PATH:$(pwd)" >> ~/pondenv/bin/activate
  curl -OL http://busco.ezlab.org/datasets/metazoa_odb9.tar.gz
  curl -OL http://busco.ezlab.org/datasets/eukaryota_odb9.tar.gz
  tar -xzvf metazoa_odb9.tar.gz 
  tar -xzvf eukaryota_odb9.tar.gz

Load your data onto /mnt/data
-----------------------------

Load your data into ``/mnt/work/data``.  You may need to make the
``/mnt/`` directory writeable by doing
::

   sudo chmod a+rwxt /mnt

first, and then creating the subdirectories
::

   cd /mnt
   mkdir -p work work/data
   cd /mnt/work/data

.. ::


   cd /mnt/work
   curl -O https://s3.amazonaws.com/public.ged.msu.edu/mrnaseq-subset.tar
   cd data
   tar xvf ../mrnaseq-subset.tar

Define your $PROJECT variable to be the location of your work
directory; in this case, it will be ``/mnt/work``::

  export PROJECT=/mnt/work

Now load your data in!

.. note::

   If you want to try things out with a small test data set, you can use
   a subset of the Nematostella data from Tulin et al. (2013)::

      cd /mnt/work
      curl -O https://s3.amazonaws.com/public.ged.msu.edu/mrnaseq-subset.tar
      cd data
      tar xvf ../mrnaseq-subset.tar

Check that your data is where it should be
------------------------------------------

Check::

   ls $PROJECT/data

If you see all the files you think you should, good!  Otherwise, debug.

If you're using the Tulin et al. data provided in the snapshot above,
you should see a bunch of files like::

   0Hour_ATCACG_L002_R1_001.fastq.gz

Next: :doc:`1-quality`
