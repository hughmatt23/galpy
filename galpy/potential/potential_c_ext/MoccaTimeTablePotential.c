#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp2d.h>
#include <galpy_potentials.h>
#include <interp_2d.h>

#include <stdatomic.h>

// MoccaTimeTablePotential: Arguments passed via args array
// args layout: [amp, nr, nt, r_min, r_max, t_min, t_max, flatten(radii), flatten(times), flatten(potential)]
double MoccaTimeTablePotentialrevaluate(double r, double t,
                                        struct potentialArg * potentialArgs) {
    double * args = potentialArgs->args;
    // Get bounds and interpolators
    double amp = *(args++);
    double rmin = *(args);
    double rmax = *(args + 1);
    int ntime = (int)*(args + 2);  // Explicitly cast to int

    double *times = args + 3;          // Start of time array
    double *masses = times + ntime;    // Start of mass array
    double tmin = times[0];
    double tmax = times[ntime-1];

    // Allocate a separate array (optional if you need isolation)
    // double *timesArray = (double *)malloc(ntime * sizeof(double));
    // memcpy(timesArray, times, ntime * sizeof(double)); // Copy only ntime elements
    double timesArray[ntime];
    double massArray[ntime];
    for (int i=0; i<ntime; i++){
        timesArray[i] = times[i];
        massArray[i] = masses[i];
    }

    // Interpolate the mass (kepAmp) at the given time
    // Extrapolate mass using slope at boundaries
    double mass;
    if (t < tmin) {
        // Extrapolate mass at t < tmin using the slope at tmin
        double mass_slope = gsl_interp_eval_deriv(*potentialArgs->interp1d, timesArray, massArray, tmin, *potentialArgs->acc1d);
        mass = gsl_interp_eval(*potentialArgs->interp1d, timesArray, massArray, tmin, *potentialArgs->acc1d) + mass_slope * (t - tmin);
    } else if (t > tmax) {
        // Extrapolate mass at t > tmax using the slope at tmax
        double mass_slope = gsl_interp_eval_deriv(*potentialArgs->interp1d, timesArray, massArray, tmax, *potentialArgs->acc1d);
        mass = gsl_interp_eval(*potentialArgs->interp1d, timesArray, massArray, tmax, *potentialArgs->acc1d) + mass_slope * (t - tmax);
    } else {
        // Interpolate mass normally
        mass = gsl_interp_eval(*potentialArgs->interp1d, timesArray, massArray, t, *potentialArgs->acc1d);
    }


    // Handle cases where r or t is out of bounds
    if (r < rmin) {
        return 0.0;  // Return zero potential for r < rmin or t < tmin
    } else if (r > rmax) {
        // If r > rmax, use Keplerian approximation
        return -mass / r;
    }

    // Use interp_2d to interpolate the potential at (r, t)
    return interp_2d_eval(potentialArgs->i2d, r, t, potentialArgs->accx, potentialArgs->accy);
}

atomic_bool initialized = false;  // Global flag
double MoccaTimeTablePotentialrforce(double r, double t,
                                     struct potentialArg * potentialArgs) {
    double * args = potentialArgs->args;

    // Get bounds and interpolators
    double amp = *(args++);
    double rmin = *(args);
    double rmax = *(args + 1);
    int ntime = (int)*(args + 2);  // Explicitly cast to int

    double *times = args + 3;          // Start of time array
    double *masses = times + ntime;    // Start of mass array
    double tmin = times[0];
    double tmax = times[ntime-1];


    // Allocate a separate array (optional if you need isolation)
    // double *timesArray = (double *)malloc(ntime * sizeof(double));
    // memcpy(timesArray, times, ntime * sizeof(double)); // Copy only ntime elements
    double timesArray[ntime];
    double massArray[ntime];
    for (int i=0; i<ntime; i++){
        timesArray[i] = times[i];
        massArray[i] = masses[i];
    }

    // Interpolate the mass (kepAmp) at the given time
    // Extrapolate mass using slope at boundaries
    double mass;
    if (t < tmin) {
        // Extrapolate mass at t < tmin using the slope at tmin
        double mass_slope = gsl_interp_eval_deriv(*potentialArgs->interp1d, timesArray, massArray, tmin, *potentialArgs->acc1d);
        mass = gsl_interp_eval(*potentialArgs->interp1d, timesArray, massArray, tmin, *potentialArgs->acc1d) + mass_slope * (t - tmin);
    } else if (t > tmax) {
        // Extrapolate mass at t > tmax using the slope at tmax
        double mass_slope = gsl_interp_eval_deriv(*potentialArgs->interp1d, timesArray, massArray, tmax, *potentialArgs->acc1d);
        mass = gsl_interp_eval(*potentialArgs->interp1d, timesArray, massArray, tmax, *potentialArgs->acc1d) + mass_slope * (t - tmax);
    } else {
        // Interpolate mass normally
        mass = gsl_interp_eval(*potentialArgs->interp1d, timesArray, massArray, t, *potentialArgs->acc1d);
    }

    // Use numerical differentiation for the radial force
    // At small radii the step size is small to match, but as the stars escape and r increases, the max step size is set by rmax
    double dr = fmin(0.003 * r, 0.003 * rmax);;  // Small step for numerical differentiation

    double fplus = interp_2d_eval(potentialArgs->i2d, t, r+dr, potentialArgs->accx, potentialArgs->accy);
    double fminus = interp_2d_eval(potentialArgs->i2d, t, r-dr, potentialArgs->accx, potentialArgs->accy);
    double fplus2 = interp_2d_eval(potentialArgs->i2d, t, r+2*dr, potentialArgs->accx, potentialArgs->accy);
    double fminus2 = interp_2d_eval(potentialArgs->i2d, t, r-2*dr, potentialArgs->accx, potentialArgs->accy);
    double force = -(-fplus2 + 8*fplus - 8*fminus + fminus2)/(12*dr);

    // Handle cases where r or t is out of bounds
    if (r < rmin) {
        return 0.0;
    } else if (r > rmax) {
        return - (mass) / (r * r);
    }
    return force;
}

double MoccaTimeTablePotentialr2deriv(double r, double t,
                                      struct potentialArg * potentialArgs) {
    double * args = potentialArgs->args;
    // Get bounds and interpolators
    double amp = *(args++);
    double rmin = *(args);
    double rmax = *(args + 1);
    int ntime = (int)*(args + 2);  // Explicitly cast to int

    double *times = args + 3;          // Start of time array
    double *masses = times + ntime;    // Start of mass array
    double tmin = times[0];
    double tmax = times[ntime-1];

    // Allocate a separate array (optional if you need isolation)
    // double *timesArray = (double *)malloc(ntime * sizeof(double));
    // memcpy(timesArray, times, ntime * sizeof(double)); // Copy only ntime elements
    double timesArray[ntime];
    double massArray[ntime];
    for (int i=0; i<ntime; i++){
        timesArray[i] = times[i];
        massArray[i] = masses[i];
    }

  
    // Interpolate the mass (kepAmp) at the given time
    // Extrapolate mass using slope at boundaries
    double mass;
    if (t < tmin) {
        // Extrapolate mass at t < tmin using the slope at tmin
        double mass_slope = gsl_interp_eval_deriv(*potentialArgs->interp1d, timesArray, massArray, tmin, *potentialArgs->acc1d);
        mass = gsl_interp_eval(*potentialArgs->interp1d, timesArray, massArray, tmin, *potentialArgs->acc1d) + mass_slope * (t - tmin);
    } else if (t > tmax) {
        // Extrapolate mass at t > tmax using the slope at tmax
        double mass_slope = gsl_interp_eval_deriv(*potentialArgs->interp1d, timesArray, massArray, tmax, *potentialArgs->acc1d);
        mass = gsl_interp_eval(*potentialArgs->interp1d, timesArray, massArray, tmax, *potentialArgs->acc1d) + mass_slope * (t - tmax);
    } else {
        // Interpolate mass normally
        mass = gsl_interp_eval(*potentialArgs->interp1d, timesArray, massArray, t, *potentialArgs->acc1d);
    }

    // Handle cases where r or t is out of bounds
    if (r < rmin) {
        return 0.0;
    } else if (r > rmax) {
        return 2.0 * (mass) / (r * r * r);
    }

    // Use numerical differentiation for the second derivative
    double dr = 1e-5;  // Small step for numerical differentiation
    double fplus = interp_2d_eval(potentialArgs->i2d, r+dr, t, potentialArgs->accx, potentialArgs->accy);
    double fminus = interp_2d_eval(potentialArgs->i2d, r-dr, t, potentialArgs->accx, potentialArgs->accy);
    double f0 = interp_2d_eval(potentialArgs->i2d, r, t, potentialArgs->accx, potentialArgs->accy);

    return (fplus - 2 * f0 + fminus) / (dr * dr);
}

double MoccaTimeTablePotentialrdens(double r, double t,
                                    struct potentialArg * potentialArgs) {
    return M_1_PI / 4.0 * (
        MoccaTimeTablePotentialr2deriv(r, t, potentialArgs) -
        2.0 * MoccaTimeTablePotentialrforce(r, t, potentialArgs) / r
    );
}