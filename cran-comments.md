# cran-comments.md

## Test environments
* Windows 10 Home x64, R 3.6.1
* ubuntu 14.05.5 (on travis-ci), R (oldrel, release and devel)
* Mac OS X 10.13.3 (on travis-ci), R (oldrel and release)
* win-builder, R (oldrelease, release and devel)
* macOS Mojave 10.14.6, R 3.5.1

### R CMD check

-- R CMD check results -------------------------------- fdaACF 0.1.0 ----
Duration: 25.2s

0 errors v | 0 warnings v | 0 notes v



## Travis CI
- linux
    - oldrel pass
    - release pass
    - devel pass
- osx
    - oldrel pass
    - release pass

* devel/osx has not been tested due to an issue with the current r-devel in osx (seen on https://travis-ci.community/t/r-devel-el-capitan-signed-pkg-fails-to-download-when-testing-on-r-devel/672/2). As CRAN doesn’t publish OSX binary packages for r-devel, I have decided to not test it.
  
## win-builder
- oldrelease pass
- release pass
- devel pass