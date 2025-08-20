"""
    CoordsKDTree{TC <: AbstractSkyCoords, K <: NearestNeighbors.KDTree}
A wrapper for a NearestNeighbors.jl `KDTree` that also includes the type of coordinate that it was constructed from as the parameter `TC <: AbstractSkyCoords`. The only field is `tree`, which provides access to the underlying `KDTree`. An internal Cartesian coordinate representation is used, so nearest neighbor queries on the `KDTree` must first be converted into the coordinate `TC`, then to Cartesian coordinates. This conversion is applied automatically in the following methods, which extend those from `NearestNeighbors.jl`.
 - `nn(tree::CoordsKDTree, coord::AbstractSkyCoords)` queries the `tree` for the entry nearest the provided `coord`.
 - `nn(tree::CoordsKDTree, coords::AbstractArray{<:AbstractSkyCoords})` queries the `tree` for the entry nearest each coordinate in `coords`.
"""
struct CoordsKDTree{TC <: AbstractSkyCoords, K <: KDTree}
    tree::K
end
CoordsKDTree{TC}(tree::K) where {TC, K} = CoordsKDTree{TC, K}(tree)
# If the element type of `data` is AbstractSkyCoords, then it contains mixed
# coordinates that need to be converted to be uniform
function CoordsKDTree(data::AbstractArray{AbstractSkyCoords}; kws...)
    TC = typeof(first(data))
    return CoordsKDTree([convert(TC, d) for d in data]; kws...)
end
# For uniform, concrete `eltype(data)`
function CoordsKDTree(data::AbstractArray{TC}; kws...) where {TC <: AbstractSkyCoords}
    data_array = [coords2cart(lon(data[i]), lat(data[i])) for i in eachindex(data)]
    tree = KDTree(data_array, Euclidean(); kws...)
    return CoordsKDTree{TC}(tree)
end

"""
    nn(tree::CoordsKDTree, coord::AbstractSkyCoords)
Queries the `tree` for the nearest entry to the provided `coord`. Returns the index into the tree of the nearest entry and the angular separation between the two coordinates, in radians.

    nn(tree::CoordsKDTree, coords::AbstractArray{<:AbstractSkyCoords})
Returns arrays `(id, sep)` containing the indices and angular separations (in radians) of the closest entries in `tree` for each coordinate in `coords`.
"""
nn(tree::CoordsKDTree{TC}, coord::T) where {TC, T <: AbstractSkyCoords} = nn(tree, convert(TC, coord))
function nn(tree::CoordsKDTree{T}, coord::T) where {T <: AbstractSkyCoords}
    id, sep = nn(tree.tree, coords2cart(lon(coord), lat(coord)))
    sep = 2 * asin(sep / 2) # Convert from cartesian separation to radians
    return id, sep
end
nn(tree::CoordsKDTree{TC}, coords::AbstractArray{T}) where {TC, T <: AbstractSkyCoords} = nn(tree, [convert(TC, coord) for coord in coords])
function nn(tree::CoordsKDTree{T}, coords::AbstractArray{T}) where {T <: AbstractSkyCoords}
    id, sep = nn(tree.tree, [coords2cart(lon(coord), lat(coord)) for coord in coords])
    sep = @. 2 * asin(sep / 2) # Convert from cartesian separation to radians
    return id, sep
end

"""
    match_coords(refcoords::AbstractArray{<:AbstractSkyCoords}, 
                 matchcoords::AbstractArray{<:AbstractSkyCoords})
Finds the nearest entries in `refcoords` to the coordinates contained in `matchcoords`. Returns `(id, sep)`, where
 - `id` is an array giving the indices into `refcoords` that are closest to the elements of `matchcoords`, and 
 - `sep` is an array giving the angular separation between the elements of `matchcoords` and the above matches.

Note that this method creates a `CoordsKDTree` from `refcoords` and then calls the method below. If you plan to use the same `refcoords` to match to many different `matchcoords`s, then you should directly construct `CoordsKDTree(refcoords)` and call the method below.

    match_coords(tree::CoordsKDTree, matchcoords::AbstractArray{<:AbstractSkyCoords})
As above, but uses a pre-constructed `tree::CoordsKDTree` rather than creating one from a reference catalog of coordinates.
"""
match_coords(tree::CoordsKDTree, matchcoords::AbstractArray{<:AbstractSkyCoords}) = nn(tree, matchcoords)
function match_coords(refcoords::AbstractArray{<:AbstractSkyCoords}, matchcoords::AbstractArray{<:AbstractSkyCoords})
    return match_coords(CoordsKDTree(refcoords), matchcoords)
end