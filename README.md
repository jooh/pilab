This code provides some general tools for univariate and multivariate
modelling of data in Matlab. The code is angled toward fMRI
datasets but is sufficiently general to support other use cases too.

# organisation
The core components of pilab are a set of container classes for storing and
manipulating data and a set of GLM-based classes for data analysis. These
classes generally implement low-level, basic functionality and are used by
more high-level data-to-results functions for representational similarity
analysis and decoding (see subfolders). 

The basic design logic is that new features can be supported by
sub-classing existing containers or GLM classes, while new specific
applications can be supported by adding new analysis functions. Finally,
batch processing is supported by a large number of modules for Automatic
Analysis (see my clone of this repo on github).

## GLM-based classes
The key heavy lifting in pilab is done by a set of GLM-based classes. A
number of key analysis methods (including representational similarity
analysis) are implemented as special (sub) cases of the general linear
model.

## containers
Data and metadata are stored together in BaseVolum and MriVolume instances.
The basic representation of data is 2D (samples by features), and meta data
is similarly supported for either dimension. A number of flexible methods
are available for indexing data according to meta data, and (in the case of
MriVolume) for converting samples (1 by nfeatures vectors) to volumetric 3D
matrix form, or for writing out to disk as nifti.

All this flexibility and metadata support comes at a substantial
performance cost. Most analysis functions therefore convert the container
data to a raw matrix form or to a lighter GLM class for
performance-dependent applications.

## analysis functions

### representational similarity analysis
The bulk of analysis functionality was developed in support of
representational similarity analysis. There are functions to compute
representational dissimilarity matrices from activation patterns of raw
timecourses (rsa/roidata2rdmvol, rsa/roidata2rdmvol_lindisc) and for comparing data
and model RDMs (rsa/roidata_rsa).

### decoding
Information-based decoding using linear discriminant analysis is supported
by decoding/roidata_lindisc.

### random effects analysis
Results from the rsa and decoding functions are saved in a common format,
which enables a standardised solution for random effects group analysis
(glm/roidata_rfx).
