===================================================
Installation and configuration: UC Davis Farm Cluster
===================================================

Many of us work with data on high performance computing clusters supported by academic institutions. This set of installation instructions is specific for the UC Davis Farm cluster*, which is currently running Ubuntu 14.04.5 LTS. We will make use of some available pre-installed modules. Other software you will need to install in your home directory without root privaleges.

This tutorial will be using a subset of `Nematostella vectensis <https://en.wikipedia.org/wiki/Starlet_sea_anemone>`__ mRNAseq data from `Tulin et al (2013) <http://evodevojournal.biomedcentral.com/articles/10.1186/2041-9139-4-16>`__. The total size of the compressed subset is ~120MB. If you are using your own data, be aware that the space and memory resource requirements will be substantially higher, so please plan accordingly.

Don't run programs on the head node. This slows down the operating system for all the other users and someone will get mad. We don't want that to happen. 

So, let's log in to a node that will be dedicated for you.

::

    srun -p high -t 24:00:00 --mem=40000 --pty bash

This login will expire after 24 hrs, so if you forget, not a big deal. When you're done, it's generally good etiquette to type ``exit`` to free up the resources for others.

Load modules
----------------

We will need several modules for this tutorial:

::

    module load bio/1.0 rsem/1.2.23 trinity/2.2.0
    
.. ::

To see what additional modules are available on the cluster:

::

    module avail

Use this to see what programs were loaded in the bio/1.0 module (which appears to be a py2.7 `conda <http://conda.pydata.org/docs/using/using.html>`__ environment):

::

    conda list

There is some software that we will need that is not installed as a module. We can do that in your home directory, but will have to be careful since you do not have permission to modify system files (no root privaleges). 

Let's download miniconda so we can set up python3 software environment:

::

    cd
    curl -LO https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b

Then source your ``.bashrc`` to put it into your `$PATH <http://unix.stackexchange.com/questions/26047/how-to-correctly-add-a-path-to-path>`__:

::

    source ~/.bashrc

Create a python 3 environment:

::
    conda create -y -n eel-pond anaconda python=3.5
    source activate eel-pond
    # If you want to exit out of the environment, type:
    # source deactivate

You are now in a python 3 environment and can install programs.

Install transrate
-----------------

We use `transrate <http://hibberdlab.com/transrate/getting_started.html>`__
to evaluate assemblies.  Install!
::

  cd
  curl -LO https://bintray.com/artifact/download/blahah/generic/transrate-1.0.3-linux-x86_64.tar.gz
  tar -zxf transrate-1.0.3-linux-x86_64.tar.gz
  echo 'export PATH=$PATH:"$HOME/transrate-1.0.3-linux-x86_64"' >> ~/pondenv/bin/activate
  curl -LO ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.3.0/ncbi-blast-2.3.0+-x64-linux.tar.gz
  tar -zxf ncbi-blast-2.3.0+-x64-linux.tar.gz
  echo 'export PATH="$HOME/ncbi-blast-2.3.0+/bin:$PATH"' >> ~/pondenv/bin/activate
  source ~/pondenv/bin/activate

Install busco

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
  source ~/pondenv/bin/activate
  
  Install `shmlast <https://github.com/camillescott/shmlast>`__

::

    conda install -y --file <(curl https://raw.githubusercontent.com/camillescott/shmlast/master/environment.txt)
    pip install --upgrade pip
    pip install shmlast

Install last

::

    cd
    curl -LO http://last.cbrc.jp/last-658.zip
    unzip last-658.zip
    pushd last-658 && make && make install prefix=~ && popd

Install the proper version of GNU parallel:

::

    cd 
    (wget -O - pi.dk/3 || curl pi.dk/3/ || fetch -o - http://pi.dk/3) | bash
    sudo cp /home/ubuntu/bin/parallel /usr/bin/parallel

Transdecoder

::

    cd
    curl -LO https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz
    tar -xvzf 2.0.1.tar.gz
    cd TransDecoder-2.0.1; make
    
BUSCO

::

    cd
    git clone https://gitlab.com/ezlab/busco.git

Put everything in the path:

::

    echo export PATH=$HOME/last-658/src:$PATH >> /home/ubuntu/miniconda3/bin/activate
    echo export PATH=$HOME/last-658/scripts:$PATH >> /home/ubuntu/miniconda3/bin/activate
    echo export PATH=$HOME/busco:$PATH >> /home/ubuntu/miniconda3/bin/activate
    echo export PATH=$HOME/TransDecoder-2.0.1:$PATH >> /home/ubuntu/miniconda3/bin/activate

Install the proper version of matplotlib

::

    pip install https://pypi.python.org/packages/source/m/matplotlib/matplotlib-1.5.1.tar.gz

Finally, install dammit from the refactor/1.0 branch

::

    pip install https://github.com/camillescott/dammit/archive/refactor/1.0.zip
    
Install databases (this step alone takes ~15-20 min)
# Is there a faster install?
# Don't need everything?

::

    dammit databases --install

By default, the metazoan busco group will be installed. For the eukaryota database, use this:

::

    dammit databases --install --busco-group eukaryota


Get the data
-----------------------------

First, create a working directory and subdirectories:

::

    cd
    mkdir -p work work/data
    cd ~/work
    curl -O https://s3.amazonaws.com/public.ged.msu.edu/mrnaseq-subset.tar
    cd data
    tar xvf ../mrnaseq-subset.tar

Define your $PROJECT variable to be the location of your work
directory; in this case, it will be ``~/work``:
::

    export PROJECT=~/work

Check that your data is where it should be
------------------------------------------

Check::

   ls $PROJECT/data

If you see all the files you think you should, good!  Otherwise, debug.

If you're using the Tulin et al. data provided in the snapshot above,
you should see a bunch of files like::

   0Hour_ATCACG_L002_R1_001.fastq.gz
   
To analyze the entire `Tulin et al. (2013) <http://evodevojournal.biomedcentral.com/articles/10.1186/2041-9139-4-16>`__ data set (if you're feeling ambitious), the files are located in my home directory on the farm cluster here:
 
::
 
    ls /home/ljcohen/Nematostella

Since they are located in my home directory, and thus read only to you, you will need to copy them to your own directory

::

    cp /home/ljcohen/Nematostella/*.gz ~/work/data/

Farm uses the `slurm workload management scheduling system <https://slurm.schedmd.com/sbatch.html>`__.  After you run through this tutorial and become familiar with how the programs run and the expected output, you can write scripts and submit these commands as slurm jobs so that they will run while you can walk away from the computer. The scrolling output you would normally see on the screen will be automatically saved to slurm output files for you to review later.

Example script, requesting 32 GB RAM on 1 node with 16 processors for 4 hrs at high priority: 

::

        #!/bin/bash -l
        #SBATCH -D /home/ljcohen/osmotic_salmon/sbatch_files/
        #SBATCH -J salmon
        #SBATCH -t 4:00:00
        #SBATCH -N 1
        #SBATCH -n 1
        #SBATHC -p high
        #SBATCH -c 16
        #SBATCH --mem=32000
        
        module load <blah>
        
        <command>
        <command>

To run this script, save as (for example) ``salmon.sh`` then submit:

::

       sbatch salmon.sh
       
After the job finished, it will produce an output file named with the job ID, e.g. ``slurm-10654264.out``. To inspect the status of the job, type this:

::

        watch squeue -u ljcohen

References
-------------
* https://wiki.cse.ucdavis.edu/support/systems/farm
* https://github.com/WhiteheadLab/Lab_Wiki/wiki/Using-the-farm-cluster
* https://github.com/RILAB/lab-docs/wiki/Using-Farm

Disclaimer*
-------------

While this set of instructions is moderately relevant to other cluster hpc systems, you will likely need to make modifications. We encourage you to contact your hpc administrators for assistance if you have questions. They are generally friendly people and like to hear from users. :) They will be able to provide helpful suggestions for how to get software running on your hpc system.




Next: :doc:`1-quality`
