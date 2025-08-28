module NearestNeighborsExt

using NearestNeighbors: Euclidean
import NearestNeighbors: KDTree, nn, knn, inrange
using SkyCoords: AbstractSkyCoords, ICRSCoords, CartesianCoords

# Convert arbitrary coordinates into reference coordinate system
# before being passed into KDTree
_convert(x) = vec(convert(CartesianCoords{ICRSCoords}, x))

"""
    KDTree(data::AbstractArray{<:AbstractSkyCoords}; kws...)
*Requires Julia ≥ 1.9 and NearestNeighbors.jl to be loaded (e.g., `using NearestNeighbors`).*

Construct a `KDTree` from NearestNeighbors.jl. The provided `data` are used to construct the tree, and the `kws...` are passed to `NearestNeighbors.KDTree`. An internal Cartesian coordinate representation is used, with a standard coordinate representation of `CartesianCoords{ICRSCoords}`. Coordinate conversions are applied automatically in the following methods, which extend those from `NearestNeighbors.jl`.
 - `nn(tree::KDTree, coord::AbstractSkyCoords)` queries the `tree` for the entry nearest the provided `coord`.
 - `nn(tree::KDTree, coords::AbstractArray{<:AbstractSkyCoords})` queries the `tree` for the entry nearest each coordinate in `coords`.
 - `knn(tree::KDTree, coord::AbstractSkyCoords, k::Int)` queries the `tree` for the `k` coordinates nearest the provided `coord`.
 - `knn(tree::KDTree, coords::AbstractArray{<:AbstractSkyCoords}, k::Int)` queries the `tree` for the `k` coordinates nearest each coordinate in `coords`.
"""
function KDTree(data::AbstractArray{<:AbstractSkyCoords}; kws...)
    if isempty(data)
        throw(ArgumentError("`data` provided to `KDTree` cannot be empty."))
    end
    return KDTree(map(_convert, vec(data)), Euclidean(); kws...)
end

"""
    nn(tree::KDTree, coord::AbstractSkyCoords)
*Requires Julia ≥ 1.9 and NearestNeighbors.jl to be loaded (e.g., `using NearestNeighbors`).*

Queries the `tree` for the nearest entry to the provided `coord`. Returns the index into the tree of the nearest entry and the angular separation between the two coordinates, in radians.

    nn(tree::KDTree, coords::AbstractArray{<:AbstractSkyCoords})
Returns arrays `(id, sep)` containing the indices and angular separations (in radians) of the closest entries in `tree` for each coordinate in `coords`.
"""
nn(tree::KDTree, coord::AbstractSkyCoords) = nn(tree, _convert(coord))
function nn(tree::KDTree, coords::AbstractArray{<:AbstractSkyCoords})
    if isempty(coords)
        throw(ArgumentError("`coords` provided to `nn` cannot be empty."))
    end
    id, sep = nn(tree, map(_convert, vec(coords)))
    sep = @. 2 * asin(sep / 2) # Convert from cartesian separation to radians
    return reshape(id, size(coords)), reshape(sep, size(coords))
end

"""
    knn(tree::KDTree, coord::AbstractSkyCoords, 
        k::Int, sortres::Bool = false)
*Requires Julia ≥ 1.9 and NearestNeighbors.jl to be loaded (e.g., `using NearestNeighbors`).*

Queries the `tree` for the `k` nearest entries to the provided `coord`. Returns vectors `(id, sep)`, which, respectively, contain the indices into the tree of the `k` nearest entries, and the angular separations between `coord` and the `k` nearest entries, in radians. If `sortres` is `true`, the returned neighbors are sorted by separation.

    knn(tree::KDTree, coords::AbstractArray{<:AbstractSkyCoords}, k::Int, sortres::Bool = false)
Returns arrays `(id, sep)` containing vectors of indices and angular separations (in radians) of the `k` closest entries in `tree` for each coordinate in `coords`. If `sortres` is `true`, the returned neighbors are sorted by separation.
"""
function knn(tree::KDTree, coord::AbstractSkyCoords, k::Int, sortres::Bool = false)
    id, sep = knn(tree, _convert(coord), k, sortres)
    @. sep = 2 * asin(sep / 2) # Convert from cartesian separation to radians
    return id, sep
end
function knn(tree::KDTree, coords::AbstractArray{<:AbstractSkyCoords}, k::Int, sortres::Bool = false)
    if isempty(coords)
        throw(ArgumentError("`coords` provided to `knn` cannot be empty."))
    end
    id, sep = knn(tree, map(_convert, vec(coords)), k, sortres)
    for i in eachindex(sep)
        @. sep[i] = 2 * asin(sep[i] / 2) # Convert from cartesian separation to radians
    end
    return reshape(id, size(coords)), reshape(sep, size(coords))
end

"""
    inrange(tree::KDTree, coord::AbstractSkyCoords, seplim::Number)
*Requires Julia ≥ 1.9 and NearestNeighbors.jl to be loaded (e.g., `using NearestNeighbors`).*

Searches for coordinates in the `tree` with angular separations from `coord` less than `seplim`, which must be given in radians. If `tree = KDTree(data)`, returns indices into `data`.

    inrange(tree::KDTree, coords::AbstractArray{<:AbstractSkyCoords}, seplim::Number)
For each coordinate in `coords`, finds coordinates in `tree` that lie within `seplim` radians.
"""
function inrange(tree::KDTree, coord::AbstractSkyCoords, 
                 seplim::Number)
    cart_sep = 2 * sin(seplim / 2) # Cartesian distance = seplim
    return inrange(tree, _convert(coord), cart_sep)
end
function inrange(tree::KDTree, coords::AbstractArray{<:AbstractSkyCoords}, 
                 seplim::Number)
    cart_sep = 2 * sin(seplim / 2) # Cartesian distance = seplim
    return inrange(tree, map(_convert, vec(coords)), cart_sep)
end

"""
    match(refcoords::AbstractArray{<:AbstractSkyCoords}, 
          matchcoords::AbstractArray{<:AbstractSkyCoords};
          nthneighbor::Int = 1)
*Requires Julia ≥ 1.9 and NearestNeighbors.jl to be loaded (e.g., `using NearestNeighbors`).*

Finds the nearest entries in `refcoords` to the coordinates contained in `matchcoords`. The keyword argument `nthneighbor` determines *which* nearest neighbor to search for; typically this should be `1` when matching one set of coordinates to another. Another common use case is setting `nthneighbor = 2` when matching a catalog against itself to find the nearest neighbor of each coordinate in the same catalog. 

Returns `(id, sep)`, where
 - `id` is an array containing indices of the coordinates in `refcoords` that matched with the elements of `matchcoords`, and 
 - `sep` is an array giving the angular separation between the elements of `matchcoords` and the above matches.

Note that this method creates a [`KDTree`](@ref) from `refcoords` and then calls the method below. If you plan to use the same `refcoords` to match to many different `matchcoords`, then you should directly construct `KDTree(refcoords)` and call the method below.

    match(tree::KDTree, matchcoords::AbstractArray{<:AbstractSkyCoords};
          nthneighbor::Int = 1)
As above, but uses a pre-constructed `tree::KDTree` rather than creating one from a reference catalog of coordinates.
"""
function Base.match(tree::KDTree, matchcoords::AbstractArray{<:AbstractSkyCoords};
                    nthneighbor::Int = 1)
    if isempty(matchcoords)
        throw(ArgumentError("`matchcoords` provided to `match` cannot be empty."))
    end
    if nthneighbor == 1
        return nn(tree, matchcoords)
    else
        id, sep = knn(tree, matchcoords, nthneighbor, true)
        return getindex.(id, nthneighbor), getindex.(sep, nthneighbor)
    end
end

function Base.match(refcoords::AbstractArray{<:AbstractSkyCoords}, 
                    matchcoords::AbstractArray{<:AbstractSkyCoords}; kws...)
    if isempty(refcoords)
        throw(ArgumentError("`refcoords` provided to `match` cannot be empty."))
    end
    return match(KDTree(refcoords), matchcoords; kws...)
end

end # module