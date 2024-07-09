"112.5 arcsec or 0.03125°"
const fov_half_angle = π / 5760

function read_tle_separate_cheops(filename)
    all_tles = read_tles_from_file(filename)
    cheops_tle = all_tles[getproperty.(all_tles, :satellite_number).==44874][end]
    tles = all_tles[getproperty.(all_tles, :satellite_number).!=44874]
    return (cheops_tle, tles)
end

"Calculate the position of the object described by a TLE in the J2000 reference frame"
function j2000_position(tle::TLE, start_time::DateTime, end_time::DateTime; timestep=Millisecond(100))
    sat_start = Nanosecond(start_time - tle_epoch(DateTime, tle)).value * 1e-9
    sat_end = Nanosecond(end_time - tle_epoch(DateTime, tle)).value * 1e-9
    timestep_val = Nanosecond(timestep).value * 1e-9

    orbp = Propagators.init(Val(:SGP4), tle)

    teme = Propagators.propagate!.(orbp, sat_start:timestep_val:sat_end) .|> first

    rot_matrix = r_eci_to_eci(TEME(), J2000(), date_to_jd(start_time))

    return map(pos -> rot_matrix * pos, teme)
end

"""
Calculate the angle between CHEOPS' line of sight, as defined by `α` and `δ`, and
an object whose position is described by `sat_tle`.
"""
function angle_to_los(sat_tle, cheops_tle, α, δ, start_time, end_time; timestep=Millisecond(100))
    sat_position = j2000_position(sat_tle, start_time, end_time; timestep=timestep)
    cheops_position = j2000_position(cheops_tle, start_time, end_time; timestep=timestep)

    sat_displacement = map((u, v) -> u - v, sat_position, cheops_position)
    sat_direction = sat_displacement ./ norm.(sat_displacement)

    star_direction = [cos(δ) * cos(α), cos(δ) * sin(α), sin(δ)]

    return map(direction -> acos(star_direction ⋅ direction), sat_direction)
end

"""
Calculate the angle between CHEOPS' line of sight, as defined by `α` and `δ`, and
an object whose position is described by `sat_tle`. Also includes the distance from
CHEOPS to the other object.
"""
function angle_distance(sat_tle, cheops_tle, α, δ, start_time, end_time; timestep=Millisecond(100))
    sat_position = j2000_position(sat_tle, start_time, end_time; timestep=timestep)
    cheops_position = j2000_position(cheops_tle, start_time, end_time; timestep=timestep)

    sat_displacement = map((u, v) -> u - v, sat_position, cheops_position)
    sat_distance = norm.(sat_displacement)
    sat_direction = sat_displacement ./ sat_distance

    star_direction = [cos(δ) * cos(α), cos(δ) * sin(α), sin(δ)]

    return (map(direction -> acos(star_direction ⋅ direction), sat_direction), sat_distance)
end

"""
Find times between `start_time` and `end_time` when the angle between CHEOPS' 
line of sight and the direction from CHEOPS to the other RSO is less than [`fov_half_angle`](@ref).
CHEOPS' line of sight is defined by the equatorial coordinates of the star it is looking at:
`α` and `δ`.
"""
function find_transits(sat_tle, cheops_tle, α, δ, start_time, end_time; timestep=Millisecond(100))
    try
        angle, dist = angle_distance(sat_tle, cheops_tle, α, δ, start_time, end_time; timestep=timestep)
        transit_mask = angle .< fov_half_angle
        times = start_time .+ (timestep * findall(transit_mask))
        return DataFrame(sat_name=sat_tle.name, sat_number=sat_tle.satellite_number, time=times, angle=angle[transit_mask], distance=dist[transit_mask])
    catch
        @error "$(sat_tle.satellite_number) produced an error between $start_time and $end_time"
        return DataFrame()
    end
end

"""
For a CHEOPS observation period defined by `filekey`, find all of the tracked RSOs
that are expected to transit CHEOPS' field of view. Requires that `[datadir]/fits/[filekey].fits`, 
`[datadir]/tle/[filekey].txt`, and `[datadir]/cheops_visits.csv` all exist and contain the 
correct information. Supprots progress logging and multiprocessing.

`filekey` should match the regex `CH_PR\\d{6}_TG\\d{6}_V0300` if you are
following the DACE CHEOPS database.
"""
function full_tle_analysis(filekey; showprogress=true, multiprocess=true)
    fits = read_fits(joinpath(datadir, "fits", "$filekey.fits"))
    cheops_tle, tles = read_tle_separate_cheops(joinpath(datadir, "tle", "$filekey.txt"))

    start_time, end_time = extrema(fits.time)

    cheops_visits = CSV.read(joinpath(dirname(pathof(SINTRA_CHEOPS)), "..", "data", "cheops_visits.csv"), DataFrame)
    α = cheops_visits[cheops_visits.file_key.==filekey, :α][1]
    δ = cheops_visits[cheops_visits.file_key.==filekey, :δ][1]

    if showprogress && multiprocess
        result_vec = @showprogress desc = "Filekey $filekey" pmap(
            tle -> find_transits(tle, cheops_tle, α, δ, start_time, end_time; timestep=Millisecond(20)),
            tles
        )
    elseif showprogress && !multiprocess
        result_vec = @showprogress desc = "Filekey $filekey" map(
            tle -> find_transits(tle, cheops_tle, α, δ, start_time, end_time; timestep=Millisecond(20)),
            tles
        )
    elseif !showprogress && multiprocess
        result_vec = pmap(
            tle -> find_transits(tle, cheops_tle, α, δ, start_time, end_time; timestep=Millisecond(20)),
            tles
        )
    else
        result_vec = map(
            tle -> find_transits(tle, cheops_tle, α, δ, start_time, end_time; timestep=Millisecond(20)),
            tles
        )
    end

    result = vcat(result_vec...)

    result_unique = combine(groupby(result, [:sat_number, :time]), group -> group[1, :])
    sort!(result_unique, :time)

    return result_unique
end

"""
From the perspective of an RSO described by `obj_tle`, find the angle in
radians between CHEOPS, described by `cheops_tle`, and the sun.
"""
function angle_cheops_to_sun(obj_tle::TLE, cheops_tle::TLE, time::DateTime)
    obj_position = CHEOPSandTLE.j2000_position(obj_tle, time, time)[1]
    cheops_position = CHEOPSandTLE.j2000_position(cheops_tle, time, time)[1]

    jd_time = datetime2julian(time)
    sun_position = r_eci_to_eci(MOD(), J2000(), jd_time) * sun_position_mod(time)

    cheops_direction = cheops_position - obj_position
    cheops_direction /= norm(cheops_direction)

    sun_direction = sun_position - obj_position
    sun_direction /= norm(sun_direction)

    return acos(cheops_direction ⋅ sun_direction)
end