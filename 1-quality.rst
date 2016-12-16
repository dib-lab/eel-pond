================================================
1. Quality Trimming and Filtering Your Sequences
================================================

.. shell start

Make sure you've got the PROJECT location defined, and your data is there:
::

   set -u
   printf "\nMy raw data is in $PROJECT/data/, and consists of $(ls -1 ${PROJECT}/data/*.fastq.gz | wc -l) files\n\n"
   set +u

**Important:** If you get an error above or the count of files is
wrong...  STOP!! Revisit the `installation instructions
<install.html>`__ for your compute platform!

Also, be sure you have loaded the right Python packages::

  source ~/pondenv/bin/activate

Link your data into your working directory
------------------------------------------

Change into your project directory and make a workspace for quality trimming:
::
  
   cd ${PROJECT}
   mkdir -p quality
   cd quality

Now, link the data files into your new workspace
::

   ln -s ../data/*.fastq.gz .

(Linking with ``ln`` avoids making a copy of the files.)

Check to make sure it worked::

   printf "I see $(ls -1 *.fastq.gz | wc -l) files here.\n"

You can also do an ``ls`` to list the files.

If you see only one entry, ``*.fastq.gz``, then the ln command above
didn't work properly.  One possibility is that your files aren't in
your data directory; another is that their names don't end with
``.fastq.gz``.

.. note::

   This protocol takes many hours (days!) to run, so you might not want
   to run it on all the data the first time.  If you're using the
   example data, you can work with a subset of it by running this command
   instead of the ``ln -s`` command above::

      cd ${PROJECT}/data
      mkdir -p extract
      for file in *.fastq.gz
      do
          gunzip -c ${file} | head -400000 | gzip \
              > extract/${file%%.fastq.gz}.extract.fastq.gz
      done

   This will pull out the first 100,000 reads of each file (4 lines per record)
   and put them in the new ``data/extract`` directory.  Then, do::

     cd ../quality/
     ln -s ../data/extract/*.fastq.gz .

   to work with the subset data.

Run FastQC on all your files
----------------------------

.. note::

   We can use FastQC to look at the quality of
   your sequences::

      fastqc *.fastq.gz

Find the right Illumina adapters
--------------------------------

You'll need to know which Illumina sequencing adapters were used for
your library in order to trim them off. Below, we will use the TruSeq3-PE.fa
adapters
::

   wget https://anonscm.debian.org/cgit/debian-med/trimmomatic.git/plain/adapters/TruSeq3-PE.fa

.. note::

   You'll need to make sure these are the right adapters for your
   data.  If they are the right adapters, you should see that some of
   the reads are trimmed; if they're not, you won't see anything
   get trimmed.
   

Adapter trim each pair of files
-------------------------------

(From this point on, you may want to be running things inside of
screen, so that you can leave it running while you go do something
else.)

.. @CTB using screen

Run:
::

   rm -f orphans.qc.fq.gz

   for filename in *_R1_*.fastq.gz
   do
        # first, make the base by removing fastq.gz
        base=$(basename $filename .fastq.gz)
        echo $base
        
        # now, construct the R2 filename by replacing R1 with R2
        baseR2=${base/_R1_/_R2_}
        echo $baseR2
        
        # finally, run Trimmomatic
        TrimmomaticPE ${base}.fastq.gz ${baseR2}.fastq.gz \
           ${base}.qc.fq.gz s1_se \
           ${baseR2}.qc.fq.gz s2_se \
           ILLUMINACLIP:TruSeq3-PE.fa:2:40:15 \
           LEADING:2 TRAILING:2 \
           SLIDINGWINDOW:4:2 \
           MINLEN:25
        
        # save the orphans
        gzip -9c s1_se s2_se >> orphans.qc.fq.gz
        rm -f s1_se s2_se
   done


The paired sequences output by this set of commands will be in the
files ending in ``.qc.fq.gz``, with any orphaned sequences all together
in ``orphans.qc.fq.gz``.

Interleave the sequences
------------------------

Next, we need to take these R1 and R2 sequences and convert them into
interleaved form, for the next step.  To do this, we'll use scripts
from the `khmer package <http://khmer.readthedocs.org>`__, which we
installed above.

Now let's use a for loop again - you might notice this is only a minor
modification of the previous for loop...
::

   for filename in *_R1_*.qc.fq.gz
   do
        # first, make the base by removing .extract.fastq.gz
        base=$(basename $filename .qc.fq.gz)
        echo $base

        # now, construct the R2 filename by replacing R1 with R2
        baseR2=${base/_R1_/_R2_}
        echo $baseR2

        # construct the output filename
        output=${base/_R1_/}.pe.qc.fq.gz

        (interleave-reads.py ${base}.qc.fq.gz ${baseR2}.qc.fq.gz | \
            gzip > $output) && rm ${base}.qc.fq.gz ${baseR2}.qc.fq.gz
   done

.. ::

   echo 1-quality DONE `date` >> ${HOME}/times.out

The final product of this is now a set of files named
``*.pe.qc.fq.gz`` that are paired-end / interleaved and quality
filtered sequences, together with the file ``orphans.qc.fq.gz`` that
contains orphaned sequences.

Finishing up
------------

Make the end product files read-only::

   chmod u-w *.pe.qc.fq.gz orphans.qc.fq.gz

to make sure you don't accidentally delete them.

Since you linked your original data files into the ``quality`` directory, you
can now do ::

   rm *.fastq.gz

to remove them from this location; you don't need them for any future steps.

Things to think about
~~~~~~~~~~~~~~~~~~~~~

Note that the filenames, while ugly, are conveniently structured with the
history of what you've done to them.  This is a good strategy to keep
in mind.

Evaluate the quality of your files with FastQC again
----------------------------------------------------

.. note::

   We can once again use FastQC to look at the
   quality of your newly-trimmed sequences::

     fastqc *.pe.qc.fq.gz

Next step: :doc:`2-diginorm`.
