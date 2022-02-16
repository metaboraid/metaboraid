# MetaboRAID

MetaboRAID is an R package enabling rapid metabolite identification using mass spectrometry MS1 and MS2 data. MetaboRAID includes functions for estimating adducts, neutral masses, and merging of spectra.

In order to install MetaboRAID, first make sure you have devtools installed in R:

```r
install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("CAMERA")
```

Install the latest stable version of gfortran for your macOS from:

```r
https://gcc.gnu.org/wiki/GFortranBinaries
```

Then run:

```r
library(devtools)
library(CAMERA)
```

MetaboRAID requires certain tools that need to be installed before running the analysis. You need to have [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) available on the system. If you already have conda on your system, make sure to update it before the following steps:
```bash
conda update conda
```


Now, Switch to R and run:

```r
install_github("metaboraid/metaboraid")
library(metaboraid)
metaboraid::install_tools()
```

An environment will be created and all the required tools will be installed in that environment.

you are now ready to start the analysis.

Alternatively you can do this manually by first creating an environment called `metaboraid_package`:

```bash
conda create -n  metaboraid_package
```
and

```bash
conda activate metaboraid_package
```

You can then install the tools using:

```bash
conda install metaboraid -c payamemami -c bioconda -c conda-forge -c anaconda -c defaults
```

For installation of CFM-ID (only supported on Linux and IOS), use:

```bash
conda create -n  metaboraid_package_cfm

conda activate metaboraid_package_cfm

conda install cfm -c bioconda -c conda-forge -c anaconda -c defaults
```

If you are using Windows, please install [CFM](http://cfmid.wishartlab.com/) manually and when running `run_cfm`, provide the absolute path to `cfm-id` file.
