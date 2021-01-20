
import Pkg
Pkg.add("OffsetArrays")
Pkg.add("Plots")
Pkg.add("ArgParse")


#!/usr/local/bin/julia

# arrays with arbitrary starting index
using OffsetArrays
# statistics package
using Statistics
# plotting
using Plots
# delimited fields file output
using DelimitedFiles
# command line argument parsing
using ArgParse
# LAPACK
using LinearAlgebra


