# Access and process MODIS data with MODIS and MODISExtra 

The aim of `MODISExtra` was to modify a few lines of code in the [`MODIS`](https://r-forge.r-project.org/R/?group_id=1252) package, to address some research needs. There are also a few new functions to perform e.g. temporal resampling (`interpolate_raster`), convert MODIS DN to physical values (`convert_dn_modis`) or MODIS Quality Flags (`convert_qf_modis`)
__ The package is essentially for personal use and is at its early development stage. Use with care!  

## To install MODISExtra, use:

```r
devtools::install_github("antoinestevens/MODISExtra")
```

The last version of the `MODIS` package should be installed with:

```r
install.packages("MODIS", repos="http://R-Forge.R-project.org")
```
