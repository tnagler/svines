Fixes issue with false positive `stringop-overread` by enforcing CXX11 standard 
until this is resolved upstream in the Eigen C++ library.

## Test environments
* ubuntu 20.04 (release, devel)
* Windows Server 2019 (release)
* win-builder (devel)
* macOS (release)

## R CMD check results

0 errors | 0 warnings | 2 notes
