r"""
rho(R) = (R/R0)^-\alpha
M = rho v
R1 = 2 R0
R2 = 2 R1 = 4 R0
M0 = M(R0) = 4/3 pi rho_0 R0^3
rho_0 = rho(r_0)
M(R) = 4/3 pi rho_0 (4-alpha)^-1 (R)^{4-\alpha} R0^\alpha
"""

# # (M(R1)-M(R0)) / M(R0)
# # alpha = 0
# C = 4/3 pi rho_0
# mr0 = C 1/4
# mr1 = C 16/4
# mr1 - mr0 = C 15/4
# (mr1 - mr0)/mr0 = 15
# # alpha = 1
# C = 4/3 pi rho_0 R0^1
# mr0 = C 1/3 R0^4
# mr1 = C 8/3 R0^4
# mr1 - mr0 = C 7/3 R0^4
# (mr1 - mr0)/mr0 = 7
# # alpha = 2
# C = 4/3 pi rho_0 R0^2
# mr0 = C 1/2
# mr1 = C 4/2
# mr1 - mr0 = C 3/2
# (mr1 - mr0)/mr0 = 3
# 
# # (M(R2)-M(R1)) / M(R1)
# # alpha = 0
# mr0 = C 1/4
# mr1 = C 16/4
# mr2 = C 256/4
# mr2 - mr1 = 240/4
# mr1 - mr0 = 15/4
# (mr2-mr1)/(mr1-mr0) = 16
# # alpha = 1
# mr0 = C 1/3
# mr1 = C 8/3
# mr2 = C 64/3
# mr2 - mr1 = 56/3
# mr1 - mr0 = 7/3
# (mr2-mr1)/(mr1-mr0) = 8
# # alpha = 2
# mr0 = C 1/2
# mr1 = C 4/2
# mr2 = C 16/2
# mr2 - mr1 = 16/2-4/2 = 8-2 = 6
# mr1 - mr0 = 3/2
# (mr2-mr1)/(mr1-mr0) = 4

import numpy as np

def m_of_r(R, alpha, r0=1):
    return 4/3. * np.pi*(4-alpha)**-1 * R**(4-alpha) * r0**alpha

r2 = 3
r1 = 2
r0 = 1


# alpha=0
# print("r2=2r1=4r0.  alpha={alpha}.  M(R2) = {0}, M(R1) = {1}, M(R0) = {2}  "
#       "(M(R2)-M(R1))/(M(R1)-M(R0)) = {3}  "
#       "(M(R1)-M(R0))/(M(R0)) = {4}  "
#       .format(m_of_r(r2, alpha, r0),
#               m_of_r(r1, alpha, r0),
#               m_of_r(r0, alpha, r0),
#               (m_of_r(r2, alpha, r0)-m_of_r(r1, alpha, r0))/(m_of_r(r1, alpha, r0)-m_of_r(r0, alpha, r0)),
#               (m_of_r(r1, alpha, r0)-m_of_r(r0, alpha, r0))/(m_of_r(r0, alpha, r0)),
#               alpha=alpha
#              )
#      )
# 
# alpha=1
# print("r2=2r1=4r0.  alpha={alpha}.  M(R2) = {0}, M(R1) = {1}, M(R0) = {2}  "
#       "(M(R2)-M(R1))/(M(R1)-M(R0)) = {3}  "
#       "(M(R1)-M(R0))/(M(R0)) = {4}  "
#       .format(m_of_r(r2, alpha, r0),
#               m_of_r(r1, alpha, r0),
#               m_of_r(r0, alpha, r0),
#               (m_of_r(r2, alpha, r0)-m_of_r(r1, alpha, r0))/(m_of_r(r1, alpha, r0)-m_of_r(r0, alpha, r0)),
#               (m_of_r(r1, alpha, r0)-m_of_r(r0, alpha, r0))/(m_of_r(r0, alpha, r0)),
#               alpha=alpha
#              )
#      )
# 
# alpha=2
# print("r2=2r1=4r0.  alpha={alpha}.  M(R2) = {0}, M(R1) = {1}, M(R0) = {2}  "
#       "(M(R2)-M(R1))/(M(R1)-M(R0)) = {3}  "
#       "(M(R1)-M(R0))/(M(R0)) = {4}  "
#       .format(m_of_r(r2, alpha, r0),
#               m_of_r(r1, alpha, r0),
#               m_of_r(r0, alpha, r0),
#               (m_of_r(r2, alpha, r0)-m_of_r(r1, alpha, r0))/(m_of_r(r1, alpha, r0)-m_of_r(r0, alpha, r0)),
#               (m_of_r(r1, alpha, r0)-m_of_r(r0, alpha, r0))/(m_of_r(r0, alpha, r0)),
#               alpha=alpha
#              )
#      )
# 
# alpha=3
# print("r2=2r1=4r0.  alpha={alpha}.  M(R2) = {0}, M(R1) = {1}, M(R0) = {2}  "
#       "(M(R2)-M(R1))/(M(R1)-M(R0)) = {3}  "
#       "(M(R1)-M(R0))/(M(R0)) = {4}  "
#       .format(m_of_r(r2, alpha, r0),
#               m_of_r(r1, alpha, r0),
#               m_of_r(r0, alpha, r0),
#               (m_of_r(r2, alpha, r0)-m_of_r(r1, alpha, r0))/(m_of_r(r1, alpha, r0)-m_of_r(r0, alpha, r0)),
#               (m_of_r(r1, alpha, r0)-m_of_r(r0, alpha, r0))/(m_of_r(r0, alpha, r0)),
#               alpha=alpha
#              )
#      )

"""
Do some correct integrals to determine how much of the outer envelope flux is
included in the inner aperture.  These should effectively be infinite integrals
over a cylinder, so there's not really much point in doing the spherical caps
(oops).
"""
def rho_sphere(z, theta, radius, alpha=1, r0=1):
    if radius < r0:
        return 1
    else:
        return (radius/r0)**-alpha
     
def rho_cyl(z, theta, radius, alpha=1, r0=1):
    r_eff = (radius**2+z**2)**0.5
    if r_eff < r0:
        return 1
    else:
        return (r_eff/r0)**-alpha


def cyl_integral(function, cyl_radius, height, args):
    """
    Try to do a cylindrical integral of a function
    """
    from scipy.integrate import tplquad

    # z bounds
    qfun = lambda y,z: -height #0
    rfun = lambda y,z:  height #2*np.pi

    # y bounds
    gfun = lambda y: 0
    hfun = lambda y: 2*np.pi

    # x bounds
    llim_a = 0
    ulim_b = cyl_radius

    return tplquad(function, a=llim_a, b=ulim_b, gfun=gfun, hfun=hfun,
                   qfun=qfun, rfun=rfun, args=args)

def ratio_outer_inner(inner_radii, outer_radii, alpha=1, r0=1):
    """
    Determine the flux ratio between an aperture integrated in the inner & outer radius
    """
    inner_1 = cyl_integral(rho_cyl, inner_radii[1], inner_radii[1]*2, (alpha, r0))[0]
    inner_0 = cyl_integral(rho_cyl, inner_radii[0], inner_radii[1]*2, (alpha, r0))[0]
    inner = inner_1 - inner_0
    outer_1 = cyl_integral(rho_cyl, outer_radii[1], outer_radii[1]*2, (alpha, r0))[0]
    outer_0 = cyl_integral(rho_cyl, outer_radii[0], outer_radii[1]*2, (alpha, r0))[0]
    outer = outer_1 - outer_0
    return outer/inner, outer, inner

def sph_integral(function, sphere_radius, args):
    """
    Try to do a spherical integral of a function
    """
    from scipy.integrate import tplquad

    # x bounds
    qfun = lambda y,z: 0
    rfun = lambda y,z: sphere_radius

    # y bounds
    gfun = lambda y: 0
    hfun = lambda y: np.pi

    # z bounds
    llim_a = 0
    ulim_b = 2*np.pi

    def f(r,th,ph):
        return function(r,th,ph,*args) * r**2 * np.sin(th)

    return tplquad(f, a=llim_a, b=ulim_b, gfun=gfun, hfun=hfun,
                   qfun=qfun, rfun=rfun)

def sphcapint(function, sphere_radius, cyl_radius, args):

    from scipy.integrate import quad

    def f(x, args):
        return (sphere_radius**2-x**2)*function(x, 0, 0, *args)

    H = (sphere_radius**2 - cyl_radius**2)**0.5
    #h = sphere_radius - H
    assert H > 1,"why integrate if density constant?"
    assert f(0, args) > 0
    assert f(H, args) > 0
    #assert f(sphere_radius, args) > f(H, args)

    return np.pi*quad(f, H, sphere_radius, args=(args,))[0]


def trivial_integrals(r_core, alpha, gridsize=100, plummer=False):

    zz,yy,xx = np.indices([gridsize]*3, dtype='float')
    center = gridsize/2.
    rr = np.sum([(ii-center)**2 for ii in (xx,yy,zz)], axis=0)**0.5

    if plummer:
        dens = (1+rr**2/r_core**2)**-2.5
    else:
        dens = (rr >= r_core)*(rr/r_core)**-alpha
        dens[(rr < r_core)] = 1.0

    img = dens.sum(axis=0)

    return img

def ratio_outer_inner(inner_radii, outer_radii, alpha=1, r0=1):
    """
    Determine the flux ratio between an aperture integrated in the inner & outer radius
    """
    gridsize = 500
    r_core = 10.
    img = trivial_integrals(r_core=r_core, alpha=alpha, gridsize=gridsize)
    yy,xx = np.indices([gridsize]*2, dtype='float')
    center = gridsize/2.
    rr = np.sum([(ii-center)**2 for ii in (xx,yy,)], axis=0)**0.5

    inner = img[(rr>=inner_radii[0]*r_core) & (rr<inner_radii[1]*r_core)].sum()
    outer = img[(rr>=outer_radii[0]*r_core) & (rr<outer_radii[1]*r_core)].sum()
    return outer/inner, outer, inner

def plummer(r, mtot, a, zh2=2.8, mh=1.66053886e-24):
    """
    Return the density given a Plummer profile
    """
    rho = 3 * mtot / (4*np.pi*a**3) * (1. + r**2/a**2)**(-2.5)
    nh2 = rho / (mh*zh2)
    return nh2


if 'mass_scalings' not in locals():
    from astropy.utils.console import ProgressBar
    mass_scalings = {'2-1to1-0':
                     {alpha: ratio_outer_inner((0,1),(1,2), alpha=alpha)
                      for alpha in ProgressBar((0,1,2,3))},
                     '3-2to2-1':
                     {alpha: ratio_outer_inner((1,2),(2,3), alpha=alpha)
                      for alpha in ProgressBar((0,1,2,3))},
                    }

"""
old, integral version
{'2-1to1-0': {0: (2.0, 50.26548245743669, 25.132741228718345),
  1: (1.065947700104067, 21.708906086812526, 20.36582665799937),
  2: (0.6127552207030073, 10.679812076833347, 17.429165376315385),
  3: (0.3810163398660765, 5.928303633156673, 15.559184772076739)},
 '3-2to2-1': {0: (1.5, 75.39822368615503, 50.26548245743669),
  1: (0.9357720656646151, 20.314587892175695, 21.708906086812526),
  2: (0.5636783713680823, 6.019979077986598, 10.679812076833347),
  3: (0.32721823156524066, 1.9398490310233178, 5.928303633156673)}}
"""
