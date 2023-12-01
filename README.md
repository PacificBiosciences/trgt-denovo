<h1 align="center"></></h1>

<h1 align="center">TRGT-denovo</h1>

<h3 align="center">Calling <em>de novo</em> tandem repeats in trios</h3>

***

TRGT-denovo is a companion tool of [TRGT (Tandem Repeat Genotyper)](https://github.com/PacificBiosciences/trgt) that does targeted *de novo* calling of tandem repeat mutations from PacBio Hifi Data in trios. It uses the output generated by TRGT to do so. 

Authors: [Tom Mokveld](https://github.com/tmokveld), [Egor Dolzhenko](https://github.com/egor-dolzhenko)

## Early version warning

Please note that TRGT-denovo is still in early development and is subject to signficant changes that can affect anything from input / output file formats to program behavior.

## Availability

* [Latest release with binary](https://github.com/PacificBiosciences/trgt-denovo/releases/latest)

To build TRGT-denovo you need a working C compiler. It was tested on Linux with Clang 13.0.0 & GCC 11.3.0 and on Mac OSX (M1) with Clang 15.0.7 & GCC 14.0.0.

## Documentation

* [Example run](docs/example.md)
* [Output](docs/output.md)
* [Command-line interface](docs/cli.md)

## Need help?
If you notice any missing features, bugs, or need assistance with analyzing the output of TRGT-denovo, 
please don't hesitate to open a GitHub issue.

## Support information
TRGT-denovo is a pre-release software intended for research use only and not for use in diagnostic procedures. 
While efforts have been made to ensure that TRGT-denovo lives up to the quality that PacBio strives for, we make no warranty regarding this software.

As TRGT-denovo is not covered by any service level agreement or the like, please do not contact a PacBio Field Applications Scientists or PacBio Customer Service for assistance with any TRGT-denovo release. 
Please report all issues through GitHub instead. 
We make no warranty that any such issue will be addressed, to any extent or within any time frame.

### DISCLAIMER
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
