import pystan

functions {

// Set up constants 

	real m_electron() { return  510998.910;}				 // Electron mass in eV
	real c() { return  299792458.;}			   	   	    	 // Speed of light in m/s
	real hbarc() { return 0.0000197327;}					 // Reduced Planck constant*speed of light in centimeter*eV
	real g_f() { return 2.30158e-28;}					 // GF in centimeter/eV
	real n_avogadro() { return 6.02214e23;}					 // Avogadro's number
	real m_osmium_target() { return 1.77198e-11;}				 // Target mass osmium
	real q_w() { return x;}	 						 // QW

// Target

   real target(real A, real detector_mass) {
   	return n_avogadro * detector_mass / A;
   }

// Flux

   real flux(real radius, real activity) {
   	return activity / (4 * pi() * radius ^ 2);
   }

// Minimum Neutrino Energy

   real min_neutrino_energy(real kinetic_energy) {
   	return kinetic_energy / 2 * (1 + sqrt(1 + 2 * m_osmium_target / kinetic_energy));
   }

// Differential Cross Section

   real differential_cross_section(real neutrino_energy, real kinetic_energy) {
   	real prefactor;
	real fraction;
	prefactor <- g_f ^ 2 * m_osmium_target * q_w ^ 2 / (4 * pi());
	fraction <- 1 - m_osmium_target * kinetic_energy / (2 * neutrino_energy ^ 2)
	if (neutrion_energy < min_neutrino_energy( kinetic_energy)) {
	   return 0;
	}
	return prefactor * fraction;
   }

// Neutrino Oscillations

   real oscillations(real sin2theta_s, real delta_m, real neutrino_energy, real radius) {
   	return 1 - sin2theta_s ^ 2 * (sin(12700 * delta_m ^ 2 * radius / neutrino_energy)) ^ 2;
   }

// Beta Nu Spectrum

   real beta_nu_spectrum(real neutrino_energy, real Q) {
   	real eta;
	real xi;
	real K;
	real p;
	real S;
	real diff_G;
	real y;
   	eta <- neutrino_energy / m_electron();
	xi <- Q / m_electron();
	K <- xi - eta +1;
	p <- sqrt(K ^ 2 -1);
	if (eta > xi) {
	   return 0;
	}
	S <- K * p * eta ^ 2;
	diff_G <- sqrt((xi - y + 1) ^ 2 -1) * (xi - y + 1) * y ^ 2;

	return S / diff_G;


   }

// Event Rate (without oscillations)

   real event_rate(real kinetic_energy, real radius, real neutrino_energy) {
   	real unripe_rate;
	unripe_rate <- target(A, detector_mass) * flux(radius, activity) * differential_cross_section(neutrino_energy, kinetic_energy) * beta_nu_spectrum(neutrino_energy, Q);
	return unripe_rate;
   }

// Event Rate (with oscillations)

   real osc_event_rate(real kinetic_energy, real radius, real neutrino_energy) {
   	real unripe_osc_rate;
	unripe_osc_rate <- target(A, detector_mass) * flux(radius, activity) * differential_cross_section(neutrino_energy, kinetic_energy) * beta_nu_spectrum(neutrino_energy, Q) * oscillations(sin2theta_s, delta_m, neutrino_energy, radius);
	return unripe_osc_rate;
   }
}

parameters {

	real sin2theta_s;
	real delta_m;   

}

transformed parameters {

}

data {

}

model {

}
