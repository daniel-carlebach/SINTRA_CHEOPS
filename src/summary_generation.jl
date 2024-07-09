"""
Make two tables `radar_perspective` and `cheops_perspective` synthesizing data
from TLE analysis and streak detetection. `radar_perspective` contains one
row for each timestep where an RSO is expected to be in the field of view.
`cheops_perspective` contains one row for each image CHEOPS produces.

- `filekey` is a unique identifier for a CHEOPS observation period. It will match the
regex `CH_PR\\d{6}_TG\\d{6}_V0300` if you are following the DACE CHEOPS database.
- `tle_analysis_result` should be filename of a CSV containing the output of 
[`full_tle_analysis`](@ref).
- `streak_frame_idx` should be a vector containing the indicies of the frames in the
current CHEOPS observation that contain streaks. The first frame is index 1.
"""
function make_tables(filekey, tle_analysis_result, streak_frame_idx)
    fits = read_fits(joinpath(datadir, "fits", filekey * ".fits"))
    fits.visible_streak .= false
    fits[streak_frame_idx, :visible_streak] .= true
    select!(fits, Not(:image))

    cheops_tle, tles = read_tle_separate_cheops(joinpath(datadir, "tle", "$filekey.txt"))

    cheops_visits = CSV.read(joinpath(dirname(pathof(SINTRA_CHEOPS)), "..", "data", "cheops_visits.csv"), DataFrame)
    tle_analysis_result.right_ascension .= rad2deg(cheops_visits[cheops_visits.file_key.==filekey, :α][1])
    tle_analysis_result.declination .= rad2deg(cheops_visits[cheops_visits.file_key.==filekey, :δ][1])

    tle_analysis_result.angle = rad2deg.(tle_analysis_result.angle)

    tle_analysis_result.inclination = missings(Float64, nrow(tle_analysis_result))
    tle_analysis_result.raan = missings(Float64, nrow(tle_analysis_result))
    tle_analysis_result.eccentricity = missings(Float64, nrow(tle_analysis_result))
    tle_analysis_result.argument_of_perigee = missings(Float64, nrow(tle_analysis_result))
    tle_analysis_result.mean_anomaly = missings(Float64, nrow(tle_analysis_result))
    tle_analysis_result.mean_motion = missings(Float64, nrow(tle_analysis_result))
    tle_analysis_result.cheops_sun_angle = missings(Float64, nrow(tle_analysis_result))

    for row in eachrow(tle_analysis_result)
        tle = tles[getproperty.(tles, :satellite_number).==row.sat_number][end]
        row.inclination = tle.inclination
        row.raan = tle.raan
        row.eccentricity = tle.eccentricity
        row.argument_of_perigee = tle.argument_of_perigee
        row.mean_anomaly = tle.mean_anomaly
        row.mean_motion = tle.mean_motion
        row.cheops_sun_angle = angle_cheops_to_sun(tle, cheops_tle, row.time) |> rad2deg
    end

    satcat = JSON.parsefile(joinpath(datadir, "satcat.json")) |> DataFrame
    max_id = maximum(parse.(Int, satcat.NORAD_CAT_ID))
    tle_analysis_result = tle_analysis_result[tle_analysis_result.sat_number.<max_id, :]

    exp_time_float = cheops_visits[cheops_visits.file_key.==filekey, :obs_exptime][1]
    exp_time = round(Int64, exp_time_float * 1000) |> Millisecond

    tle_analysis_result.exp_time .= exp_time_float

    tle_analysis_result.sat_name = replace(map(
            id -> first(satcat[satcat.NORAD_CAT_ID.==id, :OBJECT_NAME]),
            string.(tle_analysis_result.sat_number)
        ), nothing => missing)

    tle_analysis_result.rcs_label = replace(map(
            id -> first(satcat[satcat.NORAD_CAT_ID.==id, :RCS_SIZE]),
            string.(tle_analysis_result.sat_number)
        ), nothing => missing)

    time_diff = tle_analysis_result.time .- reshape(fits.time, 1, :)
    transit_during_exposure = time_diff .> Minute(0) .&& time_diff .< exp_time
    cheops_frame = findfirst.(eachrow(transit_during_exposure))
    tle_analysis_result.cheops_frame = replace(cheops_frame, nothing => missing)

    transit_streak_visible = transit_during_exposure .&& (fits.visible_streak)'

    tle_analysis_result.visible_by_cheops = any(transit_streak_visible, dims=2) |> vec

    dropmissing!(tle_analysis_result, :cheops_frame)

    fits.transit_expected = any(transit_during_exposure, dims=1) |> vec

    fits.cheops_frame = 1:size(fits)[1]

    function first_row_with_duration(group)
        duration_typed = maximum(group.time) - minimum(group.time)
        return DataFrame(transit_duration=round(duration_typed, Millisecond).value / 1000)
    end

    transit_duration_df = combine(groupby(tle_analysis_result, :cheops_frame), group -> size(group)[1] < 1 ? missing : first_row_with_duration(group))

    tle_analysis_per_frame = combine(groupby(tle_analysis_result, :cheops_frame), group -> size(group)[1] < 1 ? missing : group[1, :])
    if "time" in names(tle_analysis_per_frame)
        rename!(tle_analysis_per_frame, :time => :transit_start_time)
    end
    if "x1" in names(tle_analysis_per_frame)
        rename!(tle_analysis_per_frame, :x1 => :transit_start_time)
    end
    leftjoin!(tle_analysis_per_frame, transit_duration_df, on=:cheops_frame)
    leftjoin!(fits, tle_analysis_per_frame, on=:cheops_frame)

    leftjoin!(tle_analysis_result, fits[!, [:cheops_frame, :LOS_TO_SUN_ANGLE, :LOS_TO_MOON_ANGLE, :LOS_TO_EARTH_ANGLE, :LATITUDE, :LONGITUDE]], on=:cheops_frame)

    radar_perspective = tle_analysis_result
    cheops_perspective = fits

    return (radar_perspective, cheops_perspective)
end