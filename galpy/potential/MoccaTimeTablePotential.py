from .SphericalPotential import SphericalPotential
from astropy import units
import numpy as np
import numba
from scipy.interpolate import interp1d, RegularGridInterpolator
import matplotlib.pyplot as plt

class MoccaTimeTablePotential(SphericalPotential):
    def __init__(self, potTableList, times, kepAmp, amp=1.0, ro=None, vo=None, num_radius_points=1000000):
        SphericalPotential.__init__(self, amp=amp, ro=ro, vo=vo)
        self.hasC=True

        self.kepAmp = kepAmp
        self.eps = 0.01

        # Regular radius grid spanning all potential tables
        all_radii = np.hstack([table[:,0] for table in potTableList])  # Combine all radii
        print("All Radii:", all_radii)
        self.radii = np.linspace(np.min(all_radii), np.max(all_radii), num_radius_points)
        print("Stored Radii", self.radii)
        self.timePoints = times

        # Create a 3D array for potential values
        self.potential_3d = np.zeros((len(self.timePoints), len(self.radii)))
        print(kepAmp)
        self.massInterp = interp1d(times, kepAmp)
        # print("Mass Interp at t=0 (10000):", self.massInterp(0))
        # print("Mass Interp at t=500 (10500):", self.massInterp(500))
        # print("Mass Interp at t=1000 (11000):", self.massInterp(1000))

        # print("Mass diff, interp vs snapshot t=0:", abs(self.massInterp(0)-kepAmp[0]))
        # print("Mass diff, interp vs snapshot t=500:", abs(self.massInterp(500)-kepAmp[1]))
        # print("Mass diff, interp vs snapshot t=1000:", abs(self.massInterp(1000)-kepAmp[2]))

        rmax = []
        
        for i, table in enumerate(potTableList):
            # Split the table into radii and potential columns
            radii, potentials = table[:, 0], table[:, 1]
            rmax.append(sorted(radii)[-2])
            interpolator = interp1d(radii, potentials, bounds_error=False, fill_value="extrapolate")
            self.potential_3d[i] = interpolator(self.radii)
        self.maxRad = min(rmax)
        print("maxRad:", self.maxRad)
        print("All max radii:", rmax)

        # Expose required attributes for C integration
        self.radius_grid = self.radii  # Radius grid
        self.time_grid = self.timePoints  # Time grid
        print("Time Grid:", self.time_grid)
        self.kepAmp_grid = self.kepAmp
        self.potentials = self.potential_3d.flatten()  # Flattened potential array
        print("self.potentials:", self.potentials)
        self.nr = len(self.radius_grid)  # Number of radius points
        self.nt = len(self.time_grid)  # Number of time points
        self.max_radius = self.maxRad  # Maximum radius
        self._rmin = self.radius_grid[0]  # Minimum radius
        self._rmax = self.radius_grid[-1]  # Maximum radius
        self._tmin = self.time_grid[0]
        self._tmax = self.time_grid[-1]
        self._kepmin = self.kepAmp_grid[0]
        self._kepmax = self.kepAmp_grid[-1]
        self._amp = amp
        # self._total_mass = self.kepAmp  # Keplerian amplitude (mass)



        # Create a RegularGridInterpolator for fast lookup
        self.interpolator = RegularGridInterpolator((self.timePoints, self.radii), self.potential_3d, bounds_error=False, fill_value=None)




    def _revaluate(self, r, t=0.):
        r = np.atleast_1d(r)
        query_points = np.array([[t, ri] for ri in r])
        interpolated_values = self.interpolator(query_points)

        mass = self.massInterp(t)
        keplerian_values = (-mass/ r)
        result = np.where(r <= self.maxRad, interpolated_values, keplerian_values)

        # print("Eval Time:", t)
        # print("Eval Radius:", r)
        # print("Eval Query Points:", query_points)
        # print("Eval Interp Points:", interpolated_values)
        # print("Eval Kepler Result")
        # print("Eval Kep Amp Mass:", mass)
        # print("Eval Result:", result)
        # print()

        return result[0] if result.size == 1 else result

    def _rforce(self, r, t=0.):
        r = np.atleast_1d(r)
        delta_r = np.maximum(self.eps * r, 1e-5)
        dphi_dr_within = (self._revaluate(r + delta_r, t) - self._revaluate(r - delta_r, t)) / (2 * delta_r)
        mass = self.massInterp(t)
        dphi_dr_keplerian = mass / (r ** 2)
        result = -np.where(r <= self.radii[-1], dphi_dr_within, dphi_dr_keplerian)
        # print("Force Calc:", result[0] if result.size == 1 else result)
        return result[0] if result.size == 1 else result

    def _r2deriv(self, r, t=0.):
        r = np.atleast_1d(r)
        delta_r = np.maximum(self.eps * r, 1e-5)
        d2phi_dr2_within = (
            self._revaluate(r + delta_r, t) - 2 * self._revaluate(r, t) + self._revaluate(r - delta_r, t)
        ) / (delta_r ** 2)
        mass = self.massInterp(t)
        d2phi_dr2_keplerian = 2 * mass / (r ** 3)
        result = np.where(r <= self.radii[-1], d2phi_dr2_within, d2phi_dr2_keplerian)
        return result[0] if result.size == 1 else result