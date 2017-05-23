#' @title World map
#' @description A dataset of class SpatialPolygonsDataFrame of the World map
#' @name wrld
#' @docType data
#' @author M. Iturbide 
#' @keywords map
NULL

#' @title Bioclim RasterStack for baseline and future data
#' @description A \code{\link[raster]{stack}} dataset containing bioclimatic variables 
#' (bio1, bio4, bio5, bio9, bio15, bio18, bio19) for the reference period 1971-2000 (E-OBS dataset, Haylock et al., 2008, v12)
#' and future period 2071-2100, given by 7 different regional climate models from the ENSEMBLES project 
#' (van der Linden & Mitchell, 2009).
#' 
#' 
#' @section Bioclimatic Variables:
#' \tabular{rr}{
#' \strong{ID} \tab \strong{Variable definition}\cr
#' BIO1 \tab Annual Mean Temperature\cr
#' BIO4 \tab Temperature Seasonality\cr
#' BIO5 \tab Max Temperature of Warmest Month\cr
#' BIO9 \tab Mean Temperature of Driest Quarter\cr
#' BIO15 \tab Precipitation Seasonality\cr
#' BIO18 \tab Precipitation of Warmest Quarter\cr
#' BIO19 \tab Precipitation of Coldest Quarter\cr
#' \cr
#' \cr
#' }
#' 
#' 
#' @section Regional Climate Models: 
#' \tabular{rrrr}{
#' \strong{Acronym}\tab \strong{RCM} \tab \strong{Driving GCM} \tab \strong{Reference}\cr
#' CNRM \tab	 ALADIN \tab  ARPEGE \tab	 Radu et al. (2008) \cr
#' DMI \tab HIRHAM \tab	 ARPEGE \tab Christensen et al. (2008b) \cr
#' ETHZ \tab CLM \tab HadCM3Q0 \tab Jaeger et al. (2008)\cr
#' HC \tab HadRM3Q0 \tab HadCM3Q0 \tab	Haugen & Haakensatd (2005)\cr
#' ICTP \tab RegCM3 \tab ECHAM5-r3 \tab Pal et al. (2007)\cr
#' MPI \tab M-REMO \tab ECHAM5-r3 \tab	Jacob (2001)\cr
#' SMHI-BCM \tab	 RCA \tab BCM \tab Samuelsson et al. (2011)\cr
#' }
#' @name biostack
#' @docType data
#' @author M. Iturbide 
#' @keywords variables
#' @references 
#' Haylock, M.R., Hofstra, N., Klein Tank, A.M.G., Klok, E.J., Jones, P.D. & New, M.
#' (2008) A European daily high-resolution gridded data set of surface temperature and
#' precipitation for 1950-2006. Journal of Geophysical Research 113.
#' 
#' van der Linden, P. & Mitchell, J. (2009) ENSEMBLES: Climate Change and its Impacts:
#' Summary of research and results from the ENSEMBLES project - European Environment Agency 
#' (EEA). Tech. rep., Met Office Hadley Centre, FitzRoy Road, Exeter EX1 3PB, UK.
#' 
#' Radu, R., Deue, M. & Somot, S. (2008) Spectral nudging in a spectral regional climate
#' model. Tellus A 60, 898-910.
#' 
#' Christensen, O.B., Drews, M., Christensen, J.H., Dethloff, K., Ketelsen, K., Hebestadt,
#' I. & Rinke, A. (2008b) The HIRHAM Regional Climate Model. Version 5 (beta). Tech.
#' Rep. 06-17, Danish Meteorological Institute (DMI).
#' 
#' Jaeger, E.B., Anders, I., Luthi, D., Rockel, B., Schar, C. & Seneviratne, S.I. (2008) 
#' Analysis of ERA40-driven CLM simulations for Europe. Meteorologische Zeitschrift pp. 349-367.
#' 
#' Haugen, J.E. & Haakensatd, H. (2005) Validation of hirham version 2 with 50 km and 
#' 25 km resolution. Tech. Rep. 9, Regional Climate Development Under Global Warming (RegClim).
#' 
#' Pal, J., Giorgi, F., Bi, X., Elguindi, N., Solmon, F., Rauscher, S., Gao, X., Francisco,
#' R., Zakey, A., Winter, J., Ashfaq, M., Syed, F., Sloan, L., Bell, J., Diffenbaugh, N.,
#' Karmacharya, J., Konare, A., Martinez, D., da Rocha, R. & Steiner, A. (2007) Regiona
#' lClimate Modeling for the Developing World: The ICTP RegCM3 and RegCNET. Bulletin of 
#' the American Meteorological Society 88, 1395-1409.
#' 
#' Jacob, D. (2001) A note to the simulation of the annual and inter-annual variability of the
#' water budget over the Baltic Sea drainage basin. Meteorology and Atmospheric Physics 77, 61-73.
#' 
#' Samuelsson, P., Jones, C.G., Willen, U., Ullerstig, A., Gollvik, S., Hansson, U.,
#' Jansson ,C., Kjellstrom, E., Nikulin, G. & Wyser, K. (2011) The Rossby Centre Regional Cliate model RCA3: model
#' description and performance. Tellus A 63, 4-23.
NULL

#' @title Oak distribution
#' @description A dataset consisting of a list with two
#' data frames with xy coordinates corresponding to the distribution of two 
#' phylogenetic groups of oaks (H11 and H1)
#' @name Oak_phylo2
#' @docType data
#' @keywords xy records
#' @source  Oak_phylo2 is a modified subset of the \strong{Quercus sp Europe Petit 2002} 
#' database (Petit et al., 2002b), which is available in the \strong{Georeferenced Database 
#' of Genetic Diversity} or \strong{(GD)^2}. 
#' @references Petit, R. J. \emph{et al} 2002. Chloroplast DNA variation in european 
#' white oaks: Phylogeography and patterns of diversity based on data from over 2600 populations. 
#' Forest Ecology and Management 156 (1-3), 5-26.
#' 
#' Evolution of Trees and Forest Communities \cr 
#' Ten years of the EVOLTREE network \cr
#' Evoltree E-Lab - An information system for forest genetics \cr
#' ISBN: 978-2-9519296-3-9
NULL

#' @title Quercus pubsencens distribution
#' @description A data frame with xy coordinates of Quercus pubesnces distribution
#' @name Q_pubescens
#' @docType data
#' @keywords xy records
#' @source  Q_pubescens is a modified subset of occurrences obtained from \strong{GBIF.org} 
#' @references GBIF.org (14th March 2017) GBIF Occurrence Download http://doi.org/10.15468/dl.4ss6vr
NULL

