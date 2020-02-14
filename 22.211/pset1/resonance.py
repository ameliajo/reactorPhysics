# Resonance
#
# Nuclear resonances in the SLBW model

import scipy
from scipy.integrate import quad
from scipy.special import wofz
from scipy.constants import physical_constants
sqrt = scipy.sqrt

MODEL_CHOICES = {"wide", "narrow", "infinite dilute"}
KB = physical_constants['Boltzmann constant in eV/K'][0]

_N_RESONANCES = 3
_E0S = (6.673491e+0, 2.087152e+1, 3.668212e+1)
_GAMMA_NS = (1.475792e-3, 1.009376e-2, 3.354568e-2)
_GAMMA_YS = (2.300000e-2, 2.286379e-2, 2.300225e-2)


class BaseResonance(object):
	"""A nuclear resonance in the SLBW
	(Single Level Breit-Wigner) model

	Parameters:
	-----------
	e0:         float; energy (eV) at which resonance is centered
	gamma_n:    float; FWHM of scattering resonance
	gamma_y:    float; FWHM of absorption resonance
	potential:  float; potential scattering cross section (barns)
	            [Default: 11.2934 - appropriate for U238]
	a:          float; atomic mass of nuclide
	            [Default: 238 - appropriate for U238]

	Attributes:
	-----------
	[All Parameters]
	gamma:      float; total FWHM (scattering + absorption)
	r0:         float
	"""
	
	def __init__(self, e0, gamma_n, gamma_y, potential=11.2934, a=238):
		self.a = a
		self.e0 = e0
		self.gamma_n = gamma_n
		self.gamma_y = gamma_y
		self.gamma = gamma_n + gamma_y
		self.r0 = 2603911/e0*((a + 1)/a)**2
		self.potential = potential
		self.peak = self.sigma_y(e0)
		self.q0 = 2*sqrt(self.potential*self.r0)
	
	def _assert_x(self, e, x):
		"""A handy function for assertions.

		Parameters:
		-----------
		e:              NoneType or float; energy in eV
		x:              NoneType or float

		Returns:
		--------
		x:              float; distance in resonance widths
		"""
		errstr = "You must provide one of either 'e' in eV, or x."
		assert (e is not None) or (x is not None), errstr
		assert (e is None) or (x is None), errstr
		if x is None:
			x = self.x(e)
		return x
	
	def _assert_model_choices(self, model, sigma_b):
		"""Make sure that the user requests a valid model

		Parameters:
		-----------
		model:          str; name of the model to validate
		sigma_b:        Nonetype or float; background cross section to validate

		"""
		if model not in MODEL_CHOICES:
			raise ValueError("Model must be in " + str(MODEL_CHOICES))
		else:
			if sigma_b is None and model != "infinite dilute":
				raise ValueError("You must specify a background cross section "
				                 "for a finite dilute resonance integral.")
	
	def x(self, e):
		"""Find the resonance variable `x` at some energy.

		Parameter:
		----------
		e:      float; energy in eV

		Returns:
		--------
		x:      float
		"""
		return 2*(e - self.e0)/self.gamma
	
	def psi(self, x, t=0):
		"""Call the appropriate 'psi' function

		Parameter:
		----------
		x:      float; the value of self.x() at the current energy

		"""
		pass
	
	def chi(self, x, t=0):
		"""Call the appropriate 'chi' function

		Parameter:
		----------
		x:      float; the value of self.x() at the current energy

		"""
		pass
	
	def sigma_y(self, e=None, x=None, t=0):
		"""Capture cross section

		Parameters:
		-----------
		e:              NoneType or float; energy in eV
		x:              NoneType or float
		t:              float; temperature in K
		
		Returns:
		--------
		float; Capture cross section in barns
		"""
		x = self._assert_x(e, x)
		sigma = sqrt(self.e0/e)*self.r0* \
		        self.gamma_n*self.gamma_y/self.gamma**2* \
		        self.psi(x, t)
		return sigma
	
	def sigma_n(self, e=None, x=None, t=0):
		"""Scattering cross section

		You may approximate the scattering cross section
		as the potential cross section to simplify some of the math.
		
		Parameters:
		-----------
		e:              NoneType or float; energy in eV
		x:              NoneType or float
		t:              float; temperature in K
		
		Returns:
		--------
		float; Scattering cross section in barns
		"""
		x = self._assert_x(e, x)
		sigma = self.gamma_n/self.gamma* \
		        (self.gamma_n/self.gamma*self.r0*self.psi(x, t=t) +
		         self.q0*self.chi(x,t=t)) \
		        + self.potential
		return sigma
	
	def resonance_integral(self, e1, e2, model, numerical=True, sigma_b=None):
		self._assert_model_choices(model, sigma_b)
		pass


class FadeevaResonance(BaseResonance):
	"""A resonance using the full Fadeeva function"""
	
	def _xi(self, t):
		"""A temperature-dependent constant used for Doppler broadening

		Parameter:
		----------
		t:          float; temperature in Kelvins

		Returns:
		--------
		xi:         float
		"""
		xi = self.gamma/2*sqrt(self.a/(4*KB*t*self.e0))
		return xi
	
	def _faddeeva(self, x, t):
		"""The complex function used as a basis for the correct
		psi(x) and chi(x) functions, excluding factor of sqrt(pi).

		Parameter:
		----------
		x:          float; distance in resonance widths
		t:          float; temperature in Kelvins

		Returns:
		--------
		w:          float; the complex evaluation of the Faddeeva function
					at temperature `t`
		"""
		xi = self._xi(t)
		w = xi*wofz((x + 1j)*xi)
		return w
	
	def psi(self, x, t=0):
		"""The full psi function in the SLBW model

		Parameters:
		-----------
		x:      float; distance in resonance widths at the current energy
		t:      float, optional; temperature in Kelvins
				[Default: 0]
		"""
		if t > 0:
			return sqrt(scipy.pi)*self._faddeeva(x, t).real
		elif t == 0:
			return 1/(1 + x**2)
		else:
			raise ValueError("Temperature cannot be negative!")
	
	def chi(self, x, t=0):
		"""The full, unapproximated chi function in the SLBW model

		Parameters:
		-----------
		x:      float; distance in resonance widths at the current energy
		t:      float, optional; temperature in Kelvins
				[Default: 0]
		"""
		if t > 0:
			return sqrt(scipy.pi)*self._faddeeva(x, t).imag
		elif t == 0:
			return x/(1 + x**2)
		else:
			raise ValueError("Temperature cannot be negative!")
	
	def resonance_integral(self, e1, e2, model, numerical=True, sigma_b=None, t=0):
		"""A resonance integral computed using variable 1/E flux

		Note: Narrow Resonance model still needs to be written.

		Parameters:
		-----------
		e1:         float; lower energy bound for the integral
		e2:         float; upper energy bound for the integral
		model:      str; one of {"wide", "narrow", or "infinite dilute"}
		numerical:  Boolean; whether to use the numerical solution.
					The analytical solution has not been implemented in this code;
					 until it is, this parameter should be left as True.
					[Default: True]
		sigma_b:    float, optional; the scattering cross section in barns
					Required for "wide" and "narrow" models.
		t:          float, optional; temperature in Kelvins
					[Default: 0]

		Returns:
		--------
		ri:         float; analytic resonance integral
		"""
		if not numerical:
			raise NotImplementedError("No analytic solution has been "
			                          "implemeted at this time.")
		super().resonance_integral(e1, e2, model, numerical, sigma_b)
		sig_y = lambda e: self.sigma_y(e, t=t)
		if model == "infinite dilute":
			f = lambda e: sig_y(e)/e
			ri = quad(f, e1, e2)[0]
		elif model == "wide":
			f = lambda e: sigma_b*sig_y(e)/(sig_y(e) + sigma_b)/e
			ri = quad(f, e1, e2)[0]
		elif model == "narrow":
			sig_t = lambda e: sig_y(e) + self.sigma_n(e)
			f = lambda e: (self.potential + sigma_b)*sig_y(e)/ \
			              (sig_t(e) + sigma_b)/e
			ri = quad(f, e1, e2)[0]
		else:
			raise NotImplementedError(model + " has not been implemented yet.")
		return ri
	
	def flux(self, e, model, sigma_b=None, t=0):
		"""The flux at some energy

		Parameters:
		-----------
		e:          float; the energy in eV
		model:      str; one of {"wide", "narrow", or "infinite dilute"}
		sigma_b:    float, optional; the scattering cross section in barns
					Required for "wide" and "narrow" models.
		t:          float, optional; temperature in Kelvins
					[Default: 0]

		Returns:
		--------
		phi:        float; the flux at energy `e` in the current model
		"""
		self._assert_model_choices(model, sigma_b)
		if model == "infinite dilute":
			phi = 1/e
		elif model == "wide":
			phi = sigma_b/(self.sigma_y(e, t=t) + sigma_b)/e
		elif model == "narrow":
			phi = (sigma_b + self.potential)/ \
			      (self.sigma_n(e) + self.sigma_y(e, t=t) + sigma_b)/e
		else:
			raise ValueError(model + "is not one of " + str(MODEL_CHOICES))
		return phi
	
	def group_capture_cross_section(self, e1, e2, model, sigma_b=None, t=0):
		"""Calculate the group capture cross section sigma_y_g
		over some energy range.

		Parameters:
		-----------
		e1:         float; lower energy bound for the energy group
		e2:         float; upper energy bound for the energy group
		model:      str; one of {"wide", "narrow", or "infinite dilute"}
		sigma_b:    float, optional; the scattering cross section in barns
					Required for "wide" and "narrow" models.
		t:          float, optional; temperature in Kelvins
					[Default: 0]

		Returns:
		--------
		sigma_y_g:  float; the group capture cross section over the energy group

		"""
		ri = self.resonance_integral(e1, e2, model, True, sigma_b, t)
		phi = lambda e: self.flux(e, model, sigma_b, t)
		flux_integral = quad(phi, e1, e2)[0]
		sigma_y_g = ri/flux_integral
		return sigma_y_g


# Create the resonances used in pset
all_faddeeva_resonances = scipy.empty(3, dtype=FadeevaResonance)
for i in range(_N_RESONANCES):
	e0 = _E0S[i]
	gamma_n = _GAMMA_NS[i]
	gamma_y = _GAMMA_YS[i]
	all_faddeeva_resonances[i] = FadeevaResonance(e0, gamma_n, gamma_y)
