spacetrack_email = ""
spacetrack_password = ""
cookie_path = ""

"""
Make cookies for space-track.org at a given `path`.
You must set `SINTRA_CHEOPS.spacetrack_email` and `SINTRA_CHEOPS.spacetrack_password`
before using this function.
"""
function make_cookies(path)
    @assert !isempty(spacetrack_email)
    @assert !isempty(spacetrack_password)

    url = "https://www.space-track.org/ajaxauth/login"
    run(`curl -c $path -b $path $url -d "identity=$spacetrack_email&password=$spacetrack_password"`)

    global cookie_path = path
end

"""
Download TLEs for all tracked objects with ephocs between the start of the day 
on `start_date` and the end of the day on `end_date`.
You must set `SINTRA_CHEOPS.spacetrack_email` and `SINTRA_CHEOPS.spacetrack_password`
or have already made cookies before using this function.
"""
function download_tles(start_date::DateTime, end_date::DateTime, filename)
    start_str = floor(start_date, Day) |> Date |> string
    end_str = ceil(end_date, Day) |> Date |> string

    if isempty(cookie_path)
        make_cookies(tempname())
    end

    url = "https://www.space-track.org/basicspacedata/query/class/gp_history/EPOCH/$start_str--$end_str/format/3le/predicates/OBJECT_NAME,TLE_LINE1,TLE_LINE2/"

    run(`curl --cookie $cookie_path $url -o $filename`)
end

"""
Download TLEs for all tracked objects with ephocs within 1 day of the start and end of
the CHEOPS observation identified by `file_key`. Results are saved to 
`[SINTRA_CHEOPS.tle_dir]/[file_key].txt`. If `SINTRA_CHEOPS.tle_dir` has not been set yet,
it is set to a temporary directory that will be deleted next time Julia is closed.
You must set `SINTRA_CHEOPS.spacetrack_email` and `SINTRA_CHEOPS.spacetrack_password`
or have already made cookies before using this function.
"""
function download_tles(file_key)
    fits = read_fits(file_key)
    start_date = minimum(fits.time) - Day(1)
    end_date = maximum(fits.time) + Day(1)

    if isempty(tle_dir)
        global tle_dir = mktempdir()
    end

    download_tles(start_date, end_date, joinpath(tle_dir, "$file_key.txt"))
end

"""
Download JSON containing all cataloged resident space objects to 
`SINTRA_CHEOPS.satcat_path`. If `SINTRA_CHEOPS.satcat_path` has not been set yet,
it is set to a temporary path that will be deleted next time Julia is closed.
You must set `SINTRA_CHEOPS.spacetrack_email` and `SINTRA_CHEOPS.spacetrack_password`
or have already made cookies before using this function.
"""
function download_satcat()
    if isempty(cookie_path)
        make_cookies(tempname())
    end
    if isempty(satcat_path)
        global satcat_path = tempname()
    end

    url = "https://www.space-track.org/basicspacedata/query/class/satcat"

    run(`curl --cookie $cookie_path $url -o $satcat_path`)
end

"""
Convert string containing right ascension in hours, minutes, seconds and
declination in degrees, arcminutes, arcseconds to a 2-tuple containing
both in radians.
"""
function hms_dms_to_rad_rad(hms_dms)
    key = r"(\d{2}):(\d{2}):(\d{2}\.*\d*) \/ ([+-]\d{2}):(\d{2}):(\d{2}\.*\d*)"
    captures = match(key, hms_dms).captures
    components = parse.(Float64, captures)

    δsign = sign(components[4])
    components[4] = abs(components[4])

    α = components[1:3] ⋅ [π / 12, π / 720, π / 43200]
    δ = components[4:6] ⋅ [π / 180, π / 10800, π / 648000]

    return (α, δ * δsign)
end

"""
Get information for all CHEOPS observations and save the result to 
`SINTRA_CHEOPS.cheops_visits_path`. If `SINTRA_CHEOPS.cheops_visits_path` has 
not been set yet, it is set to a temporary path that will be deleted next time 
Julia is closed.
"""
function download_cheops_visits()
    # Query all CHEOPS visits with current data revision
    cheops_visits = pyimport("dace_query.cheops" => "Cheops").query_database(
                        filters=pydict(Dict("data_arch_rev" => pydict(Dict("equal" => pylist([3]))))),
                        output_format="pandas"
                    ) |> PyTable |> DataFrame # leaders and best

    sort!(cheops_visits, :date_mjd_start)

    # Calculate UTC start time from MJD start times
    cheops_visits.utc_start = cheops_visits.date_mjd_start .+ 2400000.5 .|> julian2datetime
    cheops_visits.utc_end = cheops_visits.date_mjd_end .+ 2400000.5 .|> julian2datetime

    αδ = map(hms_dms_to_rad_rad, cheops_visits.obj_pos_coordinates_hms_dms)
    cheops_visits.α = first.(αδ)
    cheops_visits.δ = last.(αδ)

    if isempty(cheops_visits_path)
        global cheops_visits_path = tempname()
    end

    CSV.write(cheops_visits_path, cheops_visits)
end

"""
Download corrected subarray FITS file for a CHEOPS observation with a given `file_key`
and extract it to `[SINTRA_CHEOPS.fits_dir]/[file_key].fits`. If `SINTRA_CHEOPS.fits_dir` 
has not been set yet, it is set to a temporary directory that will be deleted
next time Julia is closed.
"""
function download_images(file_key)
    if isempty(cheops_visits_path)
        download_cheops_visits()
    end

    cheops_visits = CSV.read(cheops_visits_path, DataFrame)

    rootpath = cheops_visits[cheops_visits.file_key.==file_key, :file_rootpath][1]
    visit_products = pyimport("dace_query.cheops" => "Cheops").list_data_product(
                         rootpath, output_format="pandas") |> PyTable |> DataFrame

    files = [product for product in visit_products.file if contains(product, "COR_SubArray_V0300")]

    tar_path = tempname()
    pyimport("dace_query.cheops" => "Cheops").download_files(files=files, file_type="files", output_filename=tar_path)
    extracted_tar = Tar.extract(tar_path)

    if isempty(fits_dir)
        global fits_dir = mktempdir()
    end

    for (root, dirs, fits_files) in walkdir(extracted_tar)
        for file in fits_files
            mv(joinpath(root, file), joinpath(fits_dir, "$file_key.fits"))
        end
    end
end

"""
Read `.fits` file stored at `[SINTRA_CHEOPS.fits_dir]/[file_key].fits`.
`SINTRA_CHEOPS.fits_dir` must be set before calling this function.
Images have their values mapped from 0 to 1.
"""
function read_fits(file_key)
    # Get a list of the entries in the FITS file
    @assert !isempty(fits_dir)
    f = FITS(joinpath(fits_dir, "$file_key.fits"))

    # Images are always in the second element
    imgs = read(f[2])
    # Metadata is either in the 3rd or 10th element
    df = (typeof(f[3]) == TableHDU ? f[3] : f[10]) |> DataFrame

    # Parse time strings into DateTime objects
    # DateTime objects only support millisecond precision, not microsecond
    df.time = DateTime.(map(str -> str[1:end-3], df.UTC_TIME))

    # Add images to DataFrame with 
    df.image = eachslice(imgs, dims=3) .|> adjust_image

    return df
end
