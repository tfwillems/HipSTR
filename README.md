# HipSTR
**H**aplotype-based **i**mputation, **p**hasing and genotyping of **STR**s

#### Author: Thomas Willems <twillems@mit.edu>
#### License: GNU v2

## Introduction

## Installation
HipSTR requires a standard c++ compiler as well as Java version 1.7 or later.
To obtain HipSTR and all of its associated  submodules, use:

    % git clone --recursive https://github.com/ekg/vcflib.git

To build, use Make:

    % cd vcflib
    % make

On Mac, before running Make, change the line in *vcflib/smithwaterman/Makefile* from

    % LDFLAGS=-Wl,-s
to

    % LDFLAGS=-Wl


## Usage

