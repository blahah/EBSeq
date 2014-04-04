EBSeq (experimental)
=====

Clone of EBSeq from [bioconductor](http://www.bioconductor.org/packages/devel/bioc/html/EBSeq.html) with **experimental features**. Please use the original EBSeq unless you specifically need one of the experimental features being developed here.

## Experimental features

1. Allow hyperparameters (`alpha`, `beta` and `P`) to be fixed, rather than optimised by maximum likelihood with expectation maximisation. This drastically reduces runtime for complex experiments (e.g for 6 conditions \* 3 samples, runtime reduces from 3 days to 35 seconds). **BEWARE**: there may be a corresponding decrease in the accuracy of the fitted model.

## Credits

Most of the code is (c) [Ning Leng](http://www.biostat.wisc.edu/~ningleng/), with modifications by [Richard Smith-Unna](https://github.com/Blahah).

EBSeq is *free software* under the [Artistic License 2.0](http://opensource.org/licenses/Artistic-2.0)
