encode3-curation
================

This is my collection of IPython notebooks to document some of the work I've done
to try and curate the Wold Lab's data submissions to the ENCODE-DCC.

Originally I just started sharing it so I could show my work to our super data wrangler.

However I then started writing some notebooks that might be considered useful documentation.

 * <a href="http://nbviewer.ipython.org/github/detrout/encode3-curation/blob/master/download-encodedcc.ipynb">Batch downlader</a>
 * <a href="http://nbviewer.ipython.org/github/detrout/encode3-curation/blob/master/ENCODESparql.ipynb">Using SPARQL to search metadata</a>

Git Annex
=========

There were a few pdfs that I wanted to track with in this repository, but
partially because some were still embargoed and because committing binary files to git
has issues I decided to use git-annex to track them.

Homepage: <a href="http://git-annex.branchable.com/">git-annex</a> 

git-annex is available in Debian, Ubuntu, and OS X Homebrew.

to download the pdfs you would do:

    git clone https://github.com/detrout/encode3-curation.git
    cd encode3-curiation
    git annex init
    git annex get <pdf name>

