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
Analysis (AA, see [my clone](https://github.com/jooh/automaticanalysis) of
this repo).

## GLM-based classes
The key heavy lifting in pilab is done by a set of GLM-based classes. A
number of key analysis methods (including representational similarity
analysis) are implemented as special (sub) cases of the general linear
model. Convenient class methods enable bootstrapping and permutation
testing.

## containers
Data and metadata are stored together in Volume classes, which provide a
convenient, pandas-like (but not as performant) interface. The basic
representation of data is 2D (samples by features), and metadata is similarly
supported for either dimension. A number of flexible methods are available for
indexing data according to metadata, and (in the case of MriVolume) for
converting samples (1 by nfeatures vectors) to volumetric 3D matrix form, or for
writing out to disk in nifti format.

This flexibility and metadata support comes at a substantial performance
cost. Most analysis functions therefore convert the container to a raw
matrix or to a lighter GLM class for performance-dependent applications.

## analysis functions
A few high level functions are used to provide a simple data-to-results
interface. 

### representational similarity analysis
The bulk of pilab was developed for the purposes of representational
similarity analysis (Kriegeskorte et al., 2008). Accordingly, there are
functions to compute representational dissimilarity matrices from
activation patterns of raw timecourses (rsa/roidata2rdmvol,
rsa/roidata2rdmvol_lindisc) and for comparing data and model RDMs
(rsa/roidata_rsa).

### decoding
Information-based decoding using linear discriminant analysis is supported
by rsa/roidata_lindisc.

### random effects analysis
Results from the rsa and decoding functions are saved in a common format,
which enables a standardised solution for random effects group analysis
(glm/roidata_rfx).

# about pattern information analysis and this repository
pilab is being developed by [Johan Carlin](mailto:j.carlin@ucl.ac.uk). Any
comments and queries are welcome. 

This code is an implementation of analysis methods originally proposed by [Niko
Kriegeskorte](http://www.mrc-cbu.cam.ac.uk/people/nikolaus.kriegeskorte) and
collaborators. The reference section below lists relevant academic publications,
which should be cited in any publications resulting from the use of this code.

Note that pilab is pre-alpha code undergoing active development. Absolutely
no guarantees are provided concerning the accuracy of results at this
stage. Use at your own risk, and please validate results extensively.

# references
Nili, H., Wingfield, C., Su, L., Walther, A., & Kriegeskorte, N. (2014). [A
toolbox for representational similarity
analysis](http://dx.doi.org/10.1371/journal.pcbi.1003553). PLoS Computational
Biology, 10, e1003553.

Kriegeskorte, N., Goebel, R., & Bandettini, P. A. (2006).
[Information-based functional brain
mapping](http://dx.doi.org/10.1073/pnas.0600244103). Proceedings of the
National Academy of Sciences, 103, 3863–3868.

Kriegeskorte, N., & Bandettini, P. A. (2007). [Analyzing for information,
not activation, to exploit high-resolution
fMRI](http://dx.doi.org/10.1016/j.neuroimage.2007.02.022). NeuroImage, 38,
649–662.

Kriegeskorte, N., Mur, M., & Bandettini, P. A. (2008). [Representational
similarity analysis - connecting the branches of systems
neuroscience](http://dx.doi.org/10.3389/neuro.06.004.2008).  Frontiers in
Systems Neuroscience, 2, 1–28.
