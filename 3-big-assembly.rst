==============================
3. Running the Actual Assembly
==============================

.. docker::

   RUN apt-get -y install wget
   RUN pip install -U setuptools && pip install -U khmer==2.0

.. shell start

Make sure you've got the PROJECT location defined, and your data is there:
::

   set -u
   printf "\nMy diginormed files are in $PROJECT/diginorm/, and consist of $(ls -1 ${PROJECT}/diginorm/*.keep.abundfilt.fq.gz | wc -l) files\n\n"
   set +u

**Important:** If you get an error above or the count of files is
wrong...  STOP!! Revisit the `installation instructions
<install.html>`__ for your compute platform!

Also, be sure you have loaded the right Python packages
::

   source ~/pondenv/bin/activate
   

Build the files to assemble
---------------------------

Let's make another working directory for the assembly
::

   cd ${PROJECT}
   mkdir -p assembly
   cd assembly

For paired-end data, Trinity expects two files, 'left' and 'right';
there can be orphan sequences present, however.  So, below, we split
all of our interleaved pair files in two, and then add the single-ended
seqs to one of 'em. :
::

    for file in ../diginorm/*.pe.qc.keep.abundfilt.fq.gz
    do
       split-paired-reads.py ${file}
    done
   
    cat *.1 > left.fq
    cat *.2 > right.fq
   
    gunzip -c ../diginorm/orphans.keep.abundfilt.fq.gz >> left.fq

Assembling with Trinity
-----------------------

Run the assembler!
::

    Trinity --left left.fq \
     --right right.fq --seqType fq --max_memory 14G \
     --CPU 2

Note that these last two parts (``--max_memory 14G --CPU 2``)
configure the maximum amount of memory and CPUs to
use.  You can increase (or decrease) them based on what machines you
are running on.

Once this completes, you'll have an assembled transcriptome in
``${PROJECT}/assembly/trinity_out_dir/Trinity.fasta``.

.. shell stop

Next: :doc:`4-evaluating-assembly`
