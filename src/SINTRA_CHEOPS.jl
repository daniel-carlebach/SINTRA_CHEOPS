module SINTRA_CHEOPS

using Dates, Statistics, LinearAlgebra, Distributed
using PythonCall, FITSIO, CSV, JSON, Tar, DataFrames
using SatelliteToolbox, ProgressMeter
using ImageContrastAdjustment, ImageBinarization, ImageFiltering

export download_tles, download_satcat, download_cheops_visits, download_images, read_fits

"Path of file containing information about CHEOPS visits"
cheops_visits_path = ""

"Path of file containing catalog of tracked resident space objects"
satcat_path = ""

"Directory containing FITS files with CHEOPS images"
fits_dir = ""

"Directory containing files with TLE tracking data during CHEOPS observations"
tle_dir = ""

include("downloading.jl")
include("image_processing.jl")
include("ephemeris_calculation.jl")
include("summary_generation.jl")

end # module SINTRA_CHEOPS
