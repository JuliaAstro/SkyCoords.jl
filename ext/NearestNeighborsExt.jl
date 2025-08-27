module NearestNeighborsExt

using NearestNeighbors: Euclidean
import NearestNeighbors: KDTree, nn, knn

using SkyCoords: AbstractSkyCoords, coords2cart
import SkyCoords: match_coords, CoordsKDTree

export match_coords, CoordsKDTree, nn, knn

CoordsKDTree{TC}(tree::K) where {TC, K} = CoordsKDTree{TC, K}(tree)
# If the element type of `data` is AbstractSkyCoords, then it contains mixed
# coordinates that need to be converted to be uniform
function CoordsKDTree(data::AbstractArray{AbstractSkyCoords}; kws...)
    TC = typeof(first(data))
    return CoordsKDTree(convert.(TC, data); kws...)
end
# For uniform, concrete TC
function CoordsKDTree(data::AbstractArray{TC}; kws...) where {TC <: AbstractSkyCoords}
    if isempty(data)
        throw(ArgumentError("`data` provided to `CoordsKDTree` cannot be empty."))
    end
    data_array = coords2cart.(vec(data))
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
    id, sep = nn(tree.tree, coords2cart(coord))
    sep = 2 * asin(sep / 2) # Convert from cartesian separation to radians
    return id, sep
end
nn(tree::CoordsKDTree{TC}, coords::AbstractArray{T}) where {TC, T <: AbstractSkyCoords} = nn(tree, [convert(TC, coord) for coord in coords])
function nn(tree::CoordsKDTree{T}, coords::AbstractArray{T}) where {T <: AbstractSkyCoords}
    if isempty(coords)
        throw(ArgumentError("`coords` provided to `nn` cannot be empty."))
    end
    id, sep = nn(tree.tree, vec(coords2cart.(coords)))
    sep = @. 2 * asin(sep / 2) # Convert from cartesian separation to radians
    return reshape(id, size(coords)), reshape(sep, size(coords))
end

"""
    knn(tree::CoordsKDTree, coord::AbstractSkyCoords, k::Int, sortres::Bool = false)
Queries the `tree` for the `k` nearest entries to the provided `coord`. Returns vectors `(id, sep)`, which, respectively, contain the indices into the tree of the `k` nearest entries, and the angular separations between `coord` and the `k` nearest entries, in radians. If `sortres` is `true`, the returned neighbors are sorted by separation.

    knn(tree::CoordsKDTree, coords::AbstractArray{<:AbstractSkyCoords}, k::Int, sortres::Bool = false)
Returns arrays `(id, sep)` containing vectors of indices and angular separations (in radians) of the `k` closest entries in `tree` for each coordinate in `coords`. If `sortres` is `true`, the returned neighbors are sorted by separation.
"""
knn(tree::CoordsKDTree{TC}, coord::T, k::Int, sortres::Bool = false) where {TC, T <: AbstractSkyCoords} = knn(tree, convert(TC, coord), k, sortres)
function knn(tree::CoordsKDTree{T}, coord::T, k::Int, sortres::Bool = false) where {T <: AbstractSkyCoords}
    id, sep = knn(tree.tree, coords2cart(coord), k, sortres)
    @. sep = 2 * asin(sep / 2) # Convert from cartesian separation to radians
    return id, sep
end
function knn(tree::CoordsKDTree{TC}, coords::AbstractArray{T}, k::Int, sortres::Bool = false) where {TC, T <: AbstractSkyCoords}
    return knn(tree, convert.(TC, coords), k, sortres)
end
function knn(tree::CoordsKDTree{T}, coords::AbstractArray{T}, k::Int, sortres::Bool = false) where {T <: AbstractSkyCoords}
    if isempty(coords)
        throw(ArgumentError("`coords` provided to `knn` cannot be empty."))
    end
    id, sep = knn(tree.tree, vec(coords2cart.(coords)), k, sortres)
    for i in eachindex(sep)
        @. sep[i] = 2 * asin(sep[i] / 2) # Convert from cartesian separation to radians
    end
    return reshape(id, size(coords)), reshape(sep, size(coords))
end

function match_coords(tree::CoordsKDTree, matchcoords::AbstractArray{<:AbstractSkyCoords};
                      nthneighbor::Int = 1)
    if isempty(matchcoords)
        throw(ArgumentError("`matchcoords` provided to `match_coords` cannot be empty."))
    end
    if nthneighbor == 1
        return nn(tree, matchcoords)
    else
        id, sep = knn(tree, matchcoords, nthneighbor, true)
        return getindex.(id, nthneighbor), getindex.(sep, nthneighbor)
    end
end

function match_coords(refcoords::AbstractArray{<:AbstractSkyCoords}, matchcoords::AbstractArray{<:AbstractSkyCoords}; kws...)
    if isempty(refcoords)
        throw(ArgumentError("`refcoords` provided to `match_coords` cannot be empty."))
    end
    return match_coords(CoordsKDTree(refcoords), matchcoords; kws...)
end

end # module