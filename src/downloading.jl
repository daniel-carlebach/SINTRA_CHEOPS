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
