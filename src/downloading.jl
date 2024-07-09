# query_database() = pyimport("dace_query.cheops" => "Cheops").query_database()
# list_data_product() = pyimport("dace_query.cheops" => "Cheops").list_data_product()
# download_files() = pyimport("dace_query.cheops" => "Cheops").download_files()

# TODO: make this use environment variables for username and password
"""
Make cookies for space-track.org at a given `path`.
This has my real login info hard coded into it, so please don't publically share it.
"""
function make_cookies(path)
    run(`curl -c $path -b $path https://www.space-track.org/ajaxauth/login -d "identity=dcarleba@umich.edu&password="`)
end

"""
Download TLEs for all tracked objects with ephocs between the start of the day 
on `start_date` and the end of the day on `end_date`.
Cookies are automatically generated and deleted each time this function is called.
"""
function download_tles(start_date::DateTime, end_date::DateTime, filename)
    start_str = floor(start_date, Day) |> Date |> string
    end_str = ceil(end_date, Day) |> Date |> string

    cookie_path = tempname()
    make_cookies(cookie_path)

    url = "https://www.space-track.org/basicspacedata/query/class/gp_history/EPOCH/$start_str--$end_str/format/3le/predicates/OBJECT_NAME,TLE_LINE1,TLE_LINE2/"

    run(`curl --cookie $cookie_path $url -o $filename`)
end

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

"Get information for all CHEOPS observations and return result as DataFrame"
function get_cheops_visits()
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

    return cheops_visits
end

function download_images(file_key, output_dir)
    cheops_visits = CSV.read(
        joinpath(dirname(pathof(SINTRA_CHEOPS)), "..", "data", "cheops_visits.csv"),
        DataFrame
    )
    rootpath = cheops_visits[cheops_visits.file_key.==file_key, :file_rootpath][1]
    visit_products = pyimport("dace_query.cheops" => "Cheops").list_data_product(
                         rootpath, output_format="pandas") |> PyTable |> DataFrame
    files = [product for product in visit_products.file if contains(product, "COR_SubArray_V0300")]
    tar_path = tempname()
    pyimport("dace_query.cheops" => "Cheops").download_files(files=files, file_type="files", output_filename=tar_path)
    extracted_tar = Tar.extract(tar_path)
    for (root, dirs, fits_files) in walkdir(extracted_tar)
        for file in fits_files
            mv(joinpath(root, file), joinpath(output_dir, file))
        end
    end
end