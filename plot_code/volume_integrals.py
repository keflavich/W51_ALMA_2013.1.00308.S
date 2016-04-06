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


alpha=0
print("r2=2r1=4r0.  alpha={alpha}.  M(R2) = {0}, M(R1) = {1}, M(R0) = {2}  "
      "(M(R2)-M(R1))/(M(R1)-M(R0)) = {3}  "
      "(M(R1)-M(R0))/(M(R0)) = {4}  "
      .format(m_of_r(r2, alpha, r0),
              m_of_r(r1, alpha, r0),
              m_of_r(r0, alpha, r0),
              (m_of_r(r2, alpha, r0)-m_of_r(r1, alpha, r0))/(m_of_r(r1, alpha, r0)-m_of_r(r0, alpha, r0)),
              (m_of_r(r1, alpha, r0)-m_of_r(r0, alpha, r0))/(m_of_r(r0, alpha, r0)),
              alpha=alpha
             )
     )

alpha=1
print("r2=2r1=4r0.  alpha={alpha}.  M(R2) = {0}, M(R1) = {1}, M(R0) = {2}  "
      "(M(R2)-M(R1))/(M(R1)-M(R0)) = {3}  "
      "(M(R1)-M(R0))/(M(R0)) = {4}  "
      .format(m_of_r(r2, alpha, r0),
              m_of_r(r1, alpha, r0),
              m_of_r(r0, alpha, r0),
              (m_of_r(r2, alpha, r0)-m_of_r(r1, alpha, r0))/(m_of_r(r1, alpha, r0)-m_of_r(r0, alpha, r0)),
              (m_of_r(r1, alpha, r0)-m_of_r(r0, alpha, r0))/(m_of_r(r0, alpha, r0)),
              alpha=alpha
             )
     )

alpha=2
print("r2=2r1=4r0.  alpha={alpha}.  M(R2) = {0}, M(R1) = {1}, M(R0) = {2}  "
      "(M(R2)-M(R1))/(M(R1)-M(R0)) = {3}  "
      "(M(R1)-M(R0))/(M(R0)) = {4}  "
      .format(m_of_r(r2, alpha, r0),
              m_of_r(r1, alpha, r0),
              m_of_r(r0, alpha, r0),
              (m_of_r(r2, alpha, r0)-m_of_r(r1, alpha, r0))/(m_of_r(r1, alpha, r0)-m_of_r(r0, alpha, r0)),
              (m_of_r(r1, alpha, r0)-m_of_r(r0, alpha, r0))/(m_of_r(r0, alpha, r0)),
              alpha=alpha
             )
     )

alpha=3
print("r2=2r1=4r0.  alpha={alpha}.  M(R2) = {0}, M(R1) = {1}, M(R0) = {2}  "
      "(M(R2)-M(R1))/(M(R1)-M(R0)) = {3}  "
      "(M(R1)-M(R0))/(M(R0)) = {4}  "
      .format(m_of_r(r2, alpha, r0),
              m_of_r(r1, alpha, r0),
              m_of_r(r0, alpha, r0),
              (m_of_r(r2, alpha, r0)-m_of_r(r1, alpha, r0))/(m_of_r(r1, alpha, r0)-m_of_r(r0, alpha, r0)),
              (m_of_r(r1, alpha, r0)-m_of_r(r0, alpha, r0))/(m_of_r(r0, alpha, r0)),
              alpha=alpha
             )
     )

"""
Do some correct integrals to determine how much of the outer envelope flux is
included in the inner aperture.  These should effectively be infinite integrals
over a cylinder, so there's not really much point in doing the spherical caps
(oops).
"""

def cyl_integral(function, sphere_radius, cyl_radius, args):
    """
    Try to do a cylindrical integral of a function
    """
    from scipy.integrate import tplquad

    assert sphere_radius > cyl_radius

    # x bounds
    qfun = lambda y,z: 0
    rfun = lambda y,z: 2*np.pi

    # y bounds
    gfun = lambda y: 0
    hfun = lambda y: cyl_radius

    # z bounds
    if cyl_radius > sphere_radius:
        # in the case that we want to integrate out to infinity....
        H = cyl_radius
    else:
        H = (sphere_radius**2 - cyl_radius**2)**0.5
        assert H > 0

    llim_a = -H/2
    ulim_b = H/2

    return tplquad(function, a=llim_a, b=ulim_b, gfun=gfun, hfun=hfun,
                   qfun=qfun, rfun=rfun, args=args)

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

def rho(radius, theta, z, alpha=1, r0=1):
    if radius < r0:
        return 1
    else:
        return (radius/r0)**-alpha

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
