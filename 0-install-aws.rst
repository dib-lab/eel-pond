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
            trimmomatic bowtie samtools blast2 wget bowtie2
.. ::

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

.. ::


   cd /mnt/work
   curl -O https://s3.amazonaws.com/public.ged.msu.edu/mrnaseq-subset.tar
   cd data
   tar xvf ../mrnaseq-subset.tar

Define your $PROJECT variable to be the location of your work
directory; in this case, it will be ``/mnt/work``::

  export PROJECT=/mnt/work

Check::

   ls $PROJECT/data

If you see all the files you think you should, good!  Otherwise, debug.

If you're using the Tulin et al. data provided in the snapshot above,
you should see a bunch of files like::

   0Hour_ATCACG_L002_R1_001.fastq.gz

Next: :doc:`1-quality`
