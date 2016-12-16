=================================
2. Applying Digital Normalization
=================================

In this section, we'll apply `digital normalization
<http://arxiv.org/abs/1203.4802>`__ and `variable-coverage k-mer
abundance trimming <https://peerj.com/preprints/890/>`__ to the reads
prior to assembly.  This has the effect of reducing the computational
cost of assembly `without negatively affecting the quality of the
assembly <https://peerj.com/preprints/505/>`__.

.. shell start

.. ::

   set -x
   set -e
   source /home/ubuntu/work/bin/activate

.. docker::

   RUN pip install -U setuptools && pip install -U khmer==2.0

.. note::

   You'll need ~15 GB of RAM for this, or more if you have a LOT of data.

Make sure you've got the PROJECT location defined, and your data is there:
::

   set -u
   printf "\nMy QC-trimmed files are in $PROJECT/quality/, and consist of $(ls -1 ${PROJECT}/quality/*.qc.fq.gz | wc -l) files\n\n"
   set +u

**Important:** If you get an error above or the count of files is
wrong...  STOP!! Revisit the `installation instructions
<install.html>`__ for your compute platform!

Also, be sure you have loaded the right Python packages::

  source ~/pondenv/bin/activate

Run digital normalization
-------------------------

Make a new working directory for digital normalization and link in the files:
::
   
   cd ${PROJECT}
   mkdir -p diginorm
   cd diginorm
   ln -s ../quality/*.qc.fq.gz .
   
Apply digital normalization to the paired-end reads
::

   normalize-by-median.py -p -k 20 -C 20 -M 4e9 \
     --savegraph normC20k20.ct -u orphans.qc.fq.gz \
     *.pe.qc.fq.gz

Note the ``-p`` in the normalize-by-median command -- when run on
PE data, that ensures that no paired ends are orphaned.  The ``-u`` tells
noralize-by-median that the following filename is unpaired.

Also note the ``-M`` parameter.  This specifies how much memory diginorm
should use, and should be less than the total memory on the computer
you're using. (See `choosing hash
sizes for khmer
<http://khmer.readthedocs.org/en/latest/choosing-hash-sizes.html>`__
for more information.)

Trim off likely erroneous k-mers
--------------------------------

.. ::

   echo 2-diginorm filter-abund `date` >> ${HOME}/times.out

Now, run through all the reads and trim off low-abundance parts of
high-coverage reads
::

   filter-abund.py -V -Z 18 normC20k20.ct *.keep && \
      rm *.keep normC20k20.ct

This will turn some reads into orphans when their partner read is
removed by the trimming.

Rename files
~~~~~~~~~~~~

You'll have a bunch of ``keep.abundfilt`` files -- let's make things prettier.

.. ::
   
   echo 2-diginorm extract `date` >> ${HOME}/times.out

First, let's break out the orphaned and still-paired reads
::

   for file in *.pe.*.abundfilt
   do 
      extract-paired-reads.py ${file} && \
            rm ${file}
   done

We can combine all of the orphaned reads into a single file
::

   gzip -9c orphans.qc.fq.gz.keep.abundfilt > orphans.keep.abundfilt.fq.gz && \
       rm orphans.qc.fq.gz.keep.abundfilt
   for file in *.pe.*.abundfilt.se
   do
      gzip -9c ${file} >> orphans.keep.abundfilt.fq.gz && \
           rm ${file}
   done

We can also rename the remaining PE reads & compress those files
::

   for file in *.abundfilt.pe
   do
      newfile=${file%%.fq.gz.keep.abundfilt.pe}.keep.abundfilt.fq
      mv ${file} ${newfile}
      gzip ${newfile}
   done

This leaves you with a bunch of files named ``*.keep.abundfilt.fq.gz``,
which represent the paired-end/interleaved reads that remain after
both digital normalization and error trimming, together with
``orphans.keep.abundfilt.fq.gz``

.. ::

   echo 2-diginorm DONE `date` >> ${HOME}/times.out

.. shell stop

Next: :doc:`3-big-assembly`.
