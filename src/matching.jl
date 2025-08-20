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

nn(tree::CoordsKDTree{TC}, point::T) where {TC, T <: AbstractSkyCoords} = nn(tree, convert(TC, point))
function nn(tree::CoordsKDTree{T}, point::T) where {T <: AbstractSkyCoords}
    id, sep = nn(tree.tree, coords2cart(lon(point), lat(point)))
    sep = 2 * asin(sep / 2) # Convert from cartesian separation to radians
    return id, sep
end
nn(tree::CoordsKDTree{TC}, points::AbstractArray{T}) where {TC, T <: AbstractSkyCoords} = nn(tree, [convert(TC, point) for point in points])
function nn(tree::CoordsKDTree{T}, points::AbstractArray{T}) where {T <: AbstractSkyCoords}
    id, sep = nn(tree.tree, [coords2cart(lon(point), lat(point)) for point in points])
    sep = @. 2 * asin(sep / 2) # Convert from cartesian separation to radians
    return id, sep
end

match_coords(tree::CoordsKDTree, matchcat::AbstractArray{<:AbstractSkyCoords}) = nn(tree, matchcat)
function match_coords(refcat::AbstractArray{<:AbstractSkyCoords}, matchcat::AbstractArray{<:AbstractSkyCoords})
    return match_coords(CoordsKDTree(refcat), matchcat)
end