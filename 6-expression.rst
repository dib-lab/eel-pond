==========================================
6. Analyzing RNAseq expression with Salmon
==========================================

Salmon!

We will use `salmon <http://salmon.readthedocs.org/en/latest/>`__ to quantify differential expression. `Salmon <https://github.com/COMBINE-lab/salmon>`__ is a new breed of software for quantifying RNAseq reads that is both really fast and takes transcript length into consideration (`Patro et al. 2015 <http://biorxiv.org/content/early/2015/06/27/021592>`__).

Install
=======

::
  
  cd
  curl -LO https://github.com/COMBINE-lab/salmon/releases/download/v0.7.2/Salmon-0.7.2_linux_x86_64.tar.gz
  tar -xvzf Salmon-0.7.2_linux_x86_64.tar.gz
  cd Salmon*/bin
  echo export PATH=$PATH:$(pwd) >> ~/pondenv/bin/activate
  source ~/pondenv/bin/activate
  






http://angus.readthedocs.io/en/2016/rob_quant/tut.html

Run Salmon
==========

First, build the index

::

  cd /mnt/work/quant/
  ln -s /mnt/work/assembly/trinity_out_dir/Trinity.fasta .
  salmon index --index nema --transcripts Trinity.fasta --type quasi

Then, run the salmon command:

::
  
  for R1 in $(ls *R1*.qc.fq.gz)
  do
    sample=$(basename $R1 .extract.extract.qc.fq.gz)
    echo $sample
    echo $R1
    R2=${R1/R1/R2}
    echo $R2
    salmon quant -i nema -p 2 -l IU -1 <(gunzip -c $R1) -2 <(gunzip -c $R2) -o nema_quants/${sample}
  done

(Note that --libType must come before the read files!)

This will create a bunch of directories named something like ``0Hour_ATCACG_L002001.quant``, containing a bunch of files. Take a look at what files there are:

::
  
    find 0Hour_ATCACG_L002_R1_001 -type f

The two most interesting files are ``salmon_quant.log`` and ``quant.sf``. The latter contains the counts; the former contains the log information from running things.
