module SINTRA_CHEOPS

using Dates, Statistics, LinearAlgebra, Distributed
using PythonCall, FITSIO, CSV, JSON, Tar, DataFrames
using SatelliteToolbox, ProgressMeter
using ImageContrastAdjustment, ImageBinarization, ImageFiltering



export download_tles, get_cheops_visits, download_images


# TODO: replace this with an environment variable or something
"Directory with subdirectories for FITS and TLEs"
const datadir = "c:/Users/danie/code/sure/cheops/data"

include("downloading.jl")
include("image_processing.jl")
include("ephemeris_calculation.jl")
include("summary_generation.jl")

end # module SINTRA_CHEOPS
