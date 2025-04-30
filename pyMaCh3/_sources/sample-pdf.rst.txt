Sample PDF
==========

This module deals with sampling from the posterior density function of your 
particular experimental model at different points, given your data. 
        
In order to do this, you will generally need to create a SamplePDF object derived from :py:class:`pyMaCh3._pyMaCh3.fitter.SamplePDFFDBase` 
for each sample of events for your experiment. For some more details on this you can see `the wiki page <https://github.com/mach3-software/MaCh3/wiki/04.-Making-a-samplePDF-experiment-class>`_ on this. 
The code examples there are written using c++ however the general ideas are the same. 
Happy sampling!

.. automodapi:: pyMaCh3._pyMaCh3.sample_pdf
   :members:
   :undoc-members:
   :show-inheritance:
