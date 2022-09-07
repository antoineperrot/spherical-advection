# spherical-advection

This repo is work done for fun on trying to apply the parametric Kalman Filter (Pannekoucke et al. 2016) on the sphere for a transport equation.


## [Notebook 1](https://github.com/antoineperrot/spherical-advection/blob/master/1.%20Semi%20lagrangian%20method%20-%20transport%20on%20R%C2%B2.ipynb): Semi-Lagrangian method for transport on $\mathbb{R}^2$

This notebook presents a numerical experiment where a semi-lagrangian method is employed to simulate the transport by the wind of an initial state.

<p align="center">
  <img src="https://github.com/antoineperrot/spherical-advection/blob/master/figures/2D_Advection_SemiLagrangian.png" width="800" />
</p>

__observation:__ the wind field (black arrows) transports correctly the initial state, and forms some spirals around critical points.

## [Notebook 2](https://github.com/antoineperrot/spherical-advection/blob/master/2.%20Semi%20lagrangian%20method%20-%20transport%20on%20S%C2%B2.ipynb): Semi-Lagrangian method for transport on $\mathbb{S}^2$

This time, the semi-lagrangian method is applied for a spherical domain. Unfortunately, it is a non-conservative method, meaning that additional treatment must be applied to conserve the total mass along one timestep.

<p align="center">
  <img src="https://github.com/antoineperrot/spherical-advection/blob/master/figures/semilag_transport_on_s2.png" width="800" />
</p>

__observation:__ the initial state (a high concentration on the north pole) has been distorded by the flow.

## [Notebook 3](https://github.com/antoineperrot/spherical-advection/blob/master/3.%20Isotropic%20GRF%20on%20S%C2%B2%20-%20using%20sqrt%20B.ipynb): Generating isotropic gaussian random fields on the sphere using $\mathbf{B}^{1/2}$

In order to prepare future data assimilation experiments (with an EnKF), it is necessary to know how to generate isotropic gaussian random field (GRF).
In this notebook, we do so by building the covariance matrix of an isotropic GRF (on a lat-lon grid), then computing its Cholesky decomposition to apply it to a centered normal vector.
This method is numerically very expensive and therefore limited to very low-resolution grid.

<p align="center">
  <img src="https://github.com/antoineperrot/spherical-advection/blob/master/figures/grf_sqrt_B.png" width="800" />
</p>

Then, a sanity check indicates that the random fields generated are not perfectly isotropic (longer meridional length-scale)

## [Notebook 4](https://github.com/antoineperrot/spherical-advection/blob/master/4.%20Isotropic%20GRF%20on%20S%C2%B2%20-%20using%20spherical%20harmonics.ipynb): Generating isotropic GRF on the sphere with spherical harmonics

In this notebook we use a lemma that relies on spherical harmonics to generate isotropic GRF at high-resolution.
<p align="center">
  <img src="https://github.com/antoineperrot/spherical-advection/blob/master/figures/grf_spherical_harmonics.png" width="800" />
</p>

Then, a sanity chech confirms the isotropic character of the GRF.


## [Notebook 5](https://github.com/antoineperrot/spherical-advection/blob/master/5.%20PKF%20VLATcov%20model%20for%20S%C2%B2.ipynb): Plotting VLATcov modeled correlation function on the sphere

In this notebook we use the VLATcov covariance model to represent homogeneous correlation functions on the sphere.

<p align="center">
  <img src="https://github.com/antoineperrot/spherical-advection/blob/master/figures/homogeneous_correlation_functions_with_VLATcov_model.png" width="800" />
</p>

## [Notebook 6](https://github.com/antoineperrot/spherical-advection/blob/master/6.%20Enkf%20and%20PKF%20for%20Advection%20on%20the%20sphere%20.ipynb): Enkf and PKF for Advection on the sphere

In this last notebook we compare the error statistics forecasted by an EnKF and the PKF for the transport on the sphere (to be continued...).
<p align="center">
  <img src="https://github.com/antoineperrot/spherical-advection/blob/master/figures/enkf_and_pkf_forecasts_spherical_transport.png" width="800" />
</p>


