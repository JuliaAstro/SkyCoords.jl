struct CoordKDTree{TC <: AbstractSkyCoords, K <: KDTree}
    tree::K
end
CoordKDTree{TC}(tree::K) where {TC, K} = CoordKDTree{TC, K}(tree)
# If the element type of `data` is AbstractSkyCoords, then it contains mixed
# coordinates that need to be converted to be uniform
function CoordKDTree(data::AbstractArray{AbstractSkyCoords}; kws...)
    TC = typeof(first(data))
    return CoordKDTree([convert(TC, d) for d in data]; kws...)
end
# For uniform, concrete `eltype(data)`
function CoordKDTree(data::AbstractArray{TC}; kws...) where {TC <: AbstractSkyCoords}
    data_array = [coords2cart(lon(data[i]), lat(data[i])) for i in eachindex(data)]
    tree = KDTree(data_array, Euclidean(); kws...)
    return CoordKDTree{TC}(tree)
end

nn(tree::CoordKDTree{TC}, point::T) where {TC, T <: AbstractSkyCoords} = nn(tree, convert(TC, point))
function nn(tree::CoordKDTree{T}, point::T) where {T <: AbstractSkyCoords}
    id, sep = nn(tree.tree, coords2cart(lon(point), lat(point)))
    sep = 2 * asin(sep / 2) # Convert from cartesian separation to radians
    return id, sep
end
nn(tree::CoordKDTree{TC}, points::AbstractArray{T}) where {TC, T <: AbstractSkyCoords} = nn(tree, [convert(TC, point) for point in points])
function nn(tree::CoordKDTree{T}, points::AbstractArray{T}) where {T <: AbstractSkyCoords}
    id, sep = nn(tree.tree, [coords2cart(lon(point), lat(point)) for point in points])
    sep = @. 2 * asin(sep / 2) # Convert from cartesian separation to radians
    return id, sep
end

match_coordinates(tree::CoordKDTree, matchcat::AbstractArray{<:AbstractSkyCoords}) = nn(tree, matchcat)
function match_coordinates(refcat::AbstractArray{<:AbstractSkyCoords}, matchcat::AbstractArray{<:AbstractSkyCoords})
    return match_coords(CoordKDTree(refcat), matchcat)
end