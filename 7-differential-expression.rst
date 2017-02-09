7. Extracting differentially expressed genes with edgeR
=======================================================

Run this `edgeR <https://bioconductor.org/packages/release/bioc/html/edgeR.html>`__ script (`nema.salmon.R
<https://raw.githubusercontent.com/dib-lab/eel-pond/DE/nema.salmon.R>`__)
that loads all this. For more information, see the `edgeR manual <http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf>`__

::

   curl -O -L https://raw.githubusercontent.com/dib-lab/eel-pond/DE/nema.salmon.R
   Rscript nema.salmon.R

These will produce two plots, nema-edgeR-MDS.pdf and nema-edgeR-MA-plot.pdf.

----

You can see the plot outputs for the whole data set (all the reads) here:

* `nema-edgeR-MDS.pdf <https://github.com/dib-lab/eel-pond/blob/DE/edgeR_output/nema-edgeR-MDS.pdf>`__
* `nema-edgeR-MA-plot.pdf <https://github.com/dib-lab/eel-pond/blob/DE/edgeR_output/nema-edgeR-MA-plot.pdf>`__ (0 vs 6 hour)
