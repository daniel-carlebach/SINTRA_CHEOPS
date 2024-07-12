"Binary mask of circular cropped image that CHEOPS should produce"
const image_mask = BitArray((x - 100)^2 + (y - 100)^2 <= 100^2 for x in 1:200, y in 1:200)

const nθ = 200
const nρ = 200
const θ = range(0, π, nθ)
const cosθ = cos.(θ)
const sinθ = sin.(θ)

const psf_radius = 16
const psf_offset = CartesianIndex.(
    ceil.(Int, [psf_radius * cosθ; -psf_radius * cosθ]),
    ceil.(Int, [psf_radius * sinθ; -psf_radius * sinθ])
) |> unique

"""
Adjust and image so all of the values are between 0 and 1.
All of the values that should be cropped out because CHEOPS images
are circular are set to 0.
"""
function adjust_image(image)
    @assert size(image) == (200, 200)

    mn, mx = extrema(image[image_mask])
    img_scaled = (image .- mn) / (mx - mn)
    img_scaled[.!image_mask] .= 0.0
    return img_scaled
end

"""
Increase the contrast of a CHEOPS image by first applying a median filter
then histogram equalization.
"""
function contrast_image(image)
    @assert size(image) == (200, 200)

    filtered = mapwindow(median, image, (3, 3)) |> adjust_image

    contrasted = zeros(Float64, (200, 200))
    contrasted[image_mask] = adjust_histogram(filtered[image_mask], Equalization(nbins=512))

    return contrasted
end

"""
Create a binary mask of the foreground in a CHEOPS image. [`contrast_image`](@ref)
is applied within this function, so you can pass the original image to this.
"""
function binarize_image(image)
    @assert size(image) == (200, 200)

    contrasted = contrast_image(image)

    thresh = find_threshold(contrasted[image_mask], Entropy(), nbins=512)

    @assert thresh > 0.0

    binarized = contrasted .> thresh
    binarized[.!image_mask] .= false
    return binarized
end

"""
Compute the Hough transform of a binarized CHEOPS image. The Hough transform is
centered around the index (100, 100) within `image`.
"""
function line_hough_transform(image)
    @assert size(image) == (200, 200)
    @assert !any(image[.!image_mask])

    accumulator = zeros(Int, (nθ, nρ))

    for index in findall(image)
        x, y = Tuple(index) .- 100
        for i in 1:nθ
            ρ = clamp(ceil.(Int, x * cosθ[i] + y * sinθ[i]) .+ 100, 1, 200)
            accumulator[i, ρ] += 1
        end
    end

    return accumulator
end

"Compute the circle Hough transform of a binarized CHEOPS image"
function circle_hough_transform(image)
    accumulator = zeros(Int, size(image))

    for index in findall(image)
        circle_indices = index .+ psf_offset
        for circle_index in circle_indices
            tpl = Tuple(circle_index)
            if all(tpl .> 0 .&& tpl .<= 200)
                accumulator[circle_index] += 1
            end
        end
    end

    return accumulator
end

"Find local maxima in a Hough domain with correspond to likely lines or circles"
function find_local_maxima(hough_domain; vote_threshold=60, max_count=5, window=21)
    hough_dodge = hough_domain .+ rand(-0.1:1e-3:0.1, size(hough_domain))
    max_hough = mapwindow(maximum, hough_dodge, (window, window))
    local_maxima = findall((hough_domain .>= vote_threshold) .&& (hough_dodge .== max_hough))
    sort!(local_maxima, lt=(a, b) -> hough_domain[a] > hough_domain[b], rev=true)
    return local_maxima[1:min(max_count, end)]
end

"Find all straight line streaks in a CHEOPS image and return their radius and angle"
function find_streaks(image; remove_circles=false)
    binarized = binarize_image(image)

    if remove_circles
        circle_hough = circle_hough_transform(binarized)
        circle_maxima = find_local_maxima(circle_hough, vote_threshold=70)
        circle_mask = any(
            [(x - tpl[1])^2 + (y - tpl[2])^2 <= psf_radius^2
             for x in 1:200, y in 1:200, tpl in Tuple.(circle_maxima)],
            dims=3)[:, :, 1]
        binarized[circle_mask] .= false
    end

    hough = line_hough_transform(binarized)
    maxima = find_local_maxima(hough, vote_threshold=60, window=35)

    θhat = θ[[Tuple(idx)[1] for idx in maxima]]
    rhat = [Tuple(idx)[2] - 100 for idx in maxima]

    return [(θp, rp) for (θp, rp) in zip(θhat, rhat)]
end