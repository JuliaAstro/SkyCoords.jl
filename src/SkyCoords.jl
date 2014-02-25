module SkyCoords

export ICRSCoords,
       FK5J2000Coords,
       FK5Coords,
       GalCoords,
       to_icrs,
       to_fk5j2000,
       to_gal

# Helper functions -----------------------------------------------------------

# Create rotation matrix about a given axis (x, y, z)
function xrotmat(angle)
    s = sin(angle)
    c = cos(angle)
    [1.  0.  0.;
     0.  c   s ;
     0. -s   c ]
end

function yrotmat(angle)
    s = sin(angle)
    c = cos(angle)
    [c   0. -s ;
     0.  1.  0.;
     s   0.  c ]
end

function zrotmat(angle)
    s = sin(angle)
    c = cos(angle)
    [ c   s   0.;
     -s   c   0.;
      0.  0.  1.]
end

# (lon, lat) -> (x, y, z) unit vector
coords2cart(lon, lat) = [cos(lat)*cos(lon); cos(lat)*sin(lon); sin(lat)]

# (x,y,z) unit vector -> (lon, lat)
function cart2coords(r)
    [atan2(r[2], r[1]); atan2(r[3], sqrt(r[1]*r[1] + r[2]*r[2]))]
end

# Computes the precession matrix from J2000 to the given Julian Epoch.
# Expression from from Capitaine et al. 2003 as expressed in the USNO
# Circular 179.  This should match the IAU 2006 standard from SOFA.
const pzeta =  [ 2.650545; 2306.083227;  0.2988499;  0.01801828; -0.000005971;
                -0.0000003173]
const pz =     [-2.650545; 2306.077181;  1.0927348;  0.01826837; -0.000028596;
                -0.0000002904]
const ptheta = [      0.0; 2004.191903; -0.4294934; -0.04182264; -0.000007089;
                -0.0000001274]
function precess_from_j2000_capitaine(epoch)
    t = (epoch - 2000.0) / 100.0
    tn = 1.0
    zeta = pzeta[1]
    z = pz[1]
    theta = ptheta[1]
    for i = 2:6
        tn *= t
        zeta += pzeta[i] * tn
        z += pz[i] * tn
        theta += ptheta[i] * tn
    end
    zeta /= 3600.0
    z /= 3600.0
    theta /= 3600.0
    zrotmat(-z) * yrotmat(theta) * zrotmat(-zeta)
end

# Constants --------------------------------------------------------------

# TODO: can remove this in julia 0.3 (degrees2radians has become deg2rad)
deg2rad = degrees2radians

# ICRS --> FK5 at J2000 (See USNO Circular 179, section 3.5)
eta0 = deg2rad(-19.9 / 3600000.)
xi0 = deg2rad(9.1 / 3600000.)
da0 = deg2rad(-22.9 / 3600000.)
const icrs_to_fk5j2000 = xrotmat(-eta0) * yrotmat(xi0) * zrotmat(da0)
const fk5j2000_to_icrs = icrs_to_fk5j2000'

# FK5 --> Gal
# North galactic pole and zeropoint of l in FK4/FK5 coordinates. Needed for
# transformations to/from FK4/5
# These are from Reid & Brunthaler 2004
ngp_fk5j2000_ra = deg2rad(192.859508)
ngp_fk5j2000_dec = deg2rad(27.128336)
lon0_fk5j2000 = deg2rad(122.932)
const fk5j2000_to_gal = (zrotmat(pi - lon0_fk5j2000) *
                         yrotmat(pi/2. - ngp_fk5j2000_dec) *
                         zrotmat(ngp_fk5j2000_ra))
const gal_to_fk5j2000 = fk5j2000_to_gal'

# Gal -> ICRS
const gal_to_icrs = gal_to_fk5j2000 * fk5j2000_to_icrs
const icrs_to_gal = gal_to_icrs'

# Types ----------------------------------------------------------------------

# epoch-independent systems
abstract Coords
for syms in ((:ICRSCoords, :ra, :dec),
             (:FK5J2000Coords, :ra, :dec),
             (:GalCoords, :l, :b))
    t, lon, lat = syms
    @eval begin
        immutable ($t) <: Coords
            ($lon)::Float64
            ($lat)::Float64

            function ($t)(($lon)::Real, ($lat)::Real)
                new(mod(float64($lon), 2pi), float64($lat))
            end
        end
    end
end

for t = (:FK5Coords,)
    @eval begin
        immutable ($t) <: Coords
            ra::Float64
            dec::Float64
            epoch::Float64
        
            function ($t)(ra::Real, dec::Real, epoch::Real)
                new(mod(float64(ra), 2pi), float64(dec), float64(epoch))
            end
        end
    end
end



# Functions ---------------------------------------------------------------

# To FK5 at epoch J2000
function to_fk5j2000(c::ICRSCoords)
    r = icrs_to_fk5j2000 * coords2cart(c.ra, c.dec)
    lonlat = cart2coords(r)
    FK5J2000Coords(lonlat[1], lonlat[2])
end
function to_fk5j2000(c::GalCoords)
    r = gal_to_fk5j2000 * coords2cart(c.l, c.b)
    lonlat = cart2coords(r)
    FK5J2000Coords(lonlat[1], lonlat[2])
end 
to_fk5j2000(c::FK5J2000Coords) = c

# To FK5 with general epoch
function to_fk5(c::ICRSCoords, epoch::Real)
    pmat = precess_from_j2000_capitaine(epoch)
    r = icrs_to_fk5j2000 * pmat * coords2cart(c.ra, c.dec)
    lonlat = cart2coords(r)
    FK5Coords(lonlat[1], lonlat[2], epoch)
end

# To Galactic
function to_gal(c::ICRSCoords)
    r = icrs_to_gal * coords2cart(c.ra, c.dec)
    lonlat = cart2coords(r)
    GalCoords(lonlat[1], lonlat[2])
end
function to_gal(c::FK5J2000Coords)
    r = fk5j2000_to_gal * coords2cart(c.ra, c.dec)
    lonlat = cart2coords(r)
    GalCoords(lonlat[1], lonlat[2])
end    
to_gal(c::GalCoords) = c

# To ICRS
function to_icrs(c::GalCoords)
    r = gal_to_icrs * coords2cart(c.l, c.b)
    lonlat = cart2coords(r)
    ICRSCoords(lonlat[1], lonlat[2])
end
function to_icrs(c::FK5J2000Coords)
    r = fk5j2000_to_icrs * coords2cart(c.ra, c.dec)
    lonlat = cart2coords(r)
    ICRSCoords(lonlat[1], lonlat[2])
end
function to_icrs(c::FK5Coords)
    pmat = precess_from_j2000_capitaine(c.epoch)'
    r = pmat * fk5j2000_to_icrs * coords2cart(c.ra, c.dec)
    lonlat = cart2coords(r)
    ICRSCoords(lonlat[1], lonlat[2])
end
to_icrs(c::ICRSCoords) = c

# Vectorize all functions
for t = (:FK5J2000Coords, :ICRSCoords, :GalCoords)
    @eval @vectorize_1arg $t to_gal
    @eval @vectorize_1arg $t to_fk5j2000
    @eval @vectorize_1arg $t to_icrs
end

end # module
