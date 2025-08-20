"""
    CoordsKDTree{TC <: AbstractSkyCoords, K <: NearestNeighbors.KDTree}
A wrapper for a NearestNeighbors.jl `KDTree` that also includes the type of coordinate that it was constructed from as the parameter `TC <: AbstractSkyCoords`. The only field is `tree`, which provides access to the underlying `KDTree`. An internal Cartesian coordinate representation is used, so nearest neighbor queries on the `KDTree` must first be converted into the coordinate `TC`, then to Cartesian coordinates. This conversion is applied automatically in the following methods, which extend those from `NearestNeighbors.jl`.
 - `nn(tree::CoordsKDTree, coord::AbstractSkyCoords)` queries the `tree` for the entry nearest the provided `coord`.
 - `nn(tree::CoordsKDTree, coords::AbstractArray{<:AbstractSkyCoords})` queries the `tree` for the entry nearest each coordinate in `coords`.
 - `knn(tree::CoordsKDTree, coord::AbstractSkyCoords, k::Int)` queries the `tree` for the `k` coordinates nearest the provided `coord`.
 - `knn(tree::CoordsKDTree, coords::AbstractArray{<:AbstractSkyCoords}, k::Int)` queries the `tree` for the `k` coordinates nearest each coordinate in `coords`.
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
    if isempty(data)
        throw(ArgumentError("`data` provided to `CoordsKDTree` cannot be empty."))
    end
    data_array = vec([coords2cart(lon(data[i]), lat(data[i])) for i in eachindex(data)])
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
    if isempty(coords)
        throw(ArgumentError("`coords` provided to `nn` cannot be empty."))
    end
    id, sep = nn(tree.tree, vec([coords2cart(lon(coord), lat(coord)) for coord in coords]))
    sep = @. 2 * asin(sep / 2) # Convert from cartesian separation to radians
    return reshape(id, size(coords)), reshape(sep, size(coords))
end

"""
    knn(tree::CoordsKDTree, coord::AbstractSkyCoords, k::Int)

Queries the `tree` for the `k` nearest entries to the provided `coord`. Returns vectors `(id, sep)`, which, respectively, contain the indices into the tree of the `k` nearest entries, and the angular separations between `coord` and the `k` nearest entries, in radians.

    knn(tree::CoordsKDTree, coords::AbstractArray{<:AbstractSkyCoords}, k::Int)

Returns arrays `(id, sep)` containing vectors of indices and angular separations (in radians) of the `k` closest entries in `tree` for each coordinate in `coords`.
"""
knn(tree::CoordsKDTree{TC}, coord::T, k::Int) where {TC, T <: AbstractSkyCoords} = knn(tree, convert(TC, coord), k)
function knn(tree::CoordsKDTree{T}, coord::T, k::Int) where {T <: AbstractSkyCoords}
    id, sep = knn(tree.tree, coords2cart(lon(coord), lat(coord)), k)
    @. sep = 2 * asin(sep / 2) # Convert from cartesian separation to radians
    return id, sep
end
knn(tree::CoordsKDTree{TC}, coords::AbstractArray{T}, k::Int) where {TC, T <: AbstractSkyCoords} = knn(tree, [convert(TC, coord) for coord in coords], k)
function knn(tree::CoordsKDTree{T}, coords::AbstractArray{T}, k::Int) where {T <: AbstractSkyCoords}
    if isempty(coords)
        throw(ArgumentError("`coords` provided to `knn` cannot be empty."))
    end
    id, sep = knn(tree.tree, vec([coords2cart(lon(coord), lat(coord)) for coord in coords]), k)
    for i in eachindex(sep)
        @. sep[i] = 2 * asin(sep[i] / 2) # Convert from cartesian separation to radians
    end
    return reshape(id, size(coords)), reshape(sep, size(coords))
end

"""
    match_coords(refcoords::AbstractArray{<:AbstractSkyCoords}, 
                 matchcoords::AbstractArray{<:AbstractSkyCoords};
                 nthneighbor::Int = 1)
Finds the nearest entries in `refcoords` to the coordinates contained in `matchcoords`. The keyword argument `nthneighbor` determines *which* nearest neighbor to search for; typically this should be `1` when matching one set of coordinates to another. Another common use case is setting `nthneighbor = 2` when matching a catalog against itself to find the nearest neighbor of each coordinate in the same catalog. 

Returns `(id, sep)`, where
 - `id` is an array containing indices of the coordinates in `refcoords` that matched with the elements of `matchcoords`, and 
 - `sep` is an array giving the angular separation between the elements of `matchcoords` and the above matches.

Note that this method creates a [`CoordsKDTree`](@ref) from `refcoords` and then calls the method below. If you plan to use the same `refcoords` to match to many different `matchcoords`, then you should directly construct `CoordsKDTree(refcoords)` and call the method below.

    match_coords(tree::CoordsKDTree, matchcoords::AbstractArray{<:AbstractSkyCoords};
                 nthneighbor::Int = 1)
As above, but uses a pre-constructed `tree::CoordsKDTree` rather than creating one from a reference catalog of coordinates.
"""
function match_coords(tree::CoordsKDTree, matchcoords::AbstractArray{<:AbstractSkyCoords};
                      nthneighbor::Int = 1)
    if isempty(matchcoords)
        throw(ArgumentError("`matchcoords` provided to `match_coords` cannot be empty."))
    end
    if nthneighbor == 1
        nn(tree, matchcoords)
    else
        # Get all neighbors out to nthneighbor and return the one with largest separation
        id, sep = knn(tree, matchcoords, nthneighbor)
        id_result = Array{eltype(first(id))}(undef, size(matchcoords))
        sep_result = Array{eltype(first(sep))}(undef, size(matchcoords))
        for i in eachindex(id, sep)
            sval, a = findmax(sep[i])
            id_result[i] = id[i][a]
            sep_result[i] = sval
        end
        return id_result, sep_result
    end
end
function match_coords(refcoords::AbstractArray{<:AbstractSkyCoords}, matchcoords::AbstractArray{<:AbstractSkyCoords}; kws...)
    if isempty(refcoords)
        throw(ArgumentError("`refcoords` provided to `match_coords` cannot be empty."))
    end
    return match_coords(CoordsKDTree(refcoords), matchcoords; kws...)
end