# Statistical Raman Analysis of Carbon Nanomaterials

Raman spectroscopy is a common tool for the analysis of carbon nanomaterials, however like many micro-scale measurements it is not automatically clear how spectra collected from single points correspond to the entire bulk material.
By collecting many points this problem can be addressed. This was the subject of a research project with full details prublished [here](https://pubs.acs.org/doi/10.1021/acsanm.0c02361).
Please reference this publication if the code or data are used.

A code has been developed to fit many independent Raman spectra from carbon nanomaterials, returning key peak parameters from the D, G and 2D peaks. This Python code is available for download and use.
Requires Python3.7 and modules:
 - [lmfit](https://lmfit.github.io/lmfit-py/)
 - [matplotlib](https://matplotlib.org/)
 - [tkinter](https://docs.python.org/3/library/tkinter.html)

To download the complete program that requires no installation of Python or modules click on the 'New Releases' tab to the right to download a .zip folder. After extraction use `FittingRamanMap3.4.exe` to open the analysis program. It is recommended to create a shortcut to this file, as moving the `.exe` file from the parent folder will cause the program to crash.

The most up to-date version, 4.1, is only available as a python `.py` file at this time. This version includes new flexibility to open Renishaw map data without any formatting or data modifications.
