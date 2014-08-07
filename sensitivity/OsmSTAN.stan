
functions{

// Set up constants 

	real m_electron() { return  510998.910;}				 // Electron mass in eV
	real c() { return  299792458.;}			   	   	    	 // Speed of light in m/s
	real hbarc() { return 0.0000197327;}					 // Reduced Planck constant*speed of light in centimeter*eV
	real g_f() { return 2.30158e-28;}					 // GF in centimeter/eV
	real n_avogadro() { return 6.02214e23;}					 // Avogadro's number
	
	real q_w() { return 105.885;}	 					 // QW

	real a() { return 190.23;}						 // Atomic Weight of Osmium
	real flux_constant() { return 1 / (4 * pi());}				 // 1 over 4pi

// Target

   real target(real detector_mass) {
   	return n_avogadro() * detector_mass / a();
   }


// Minimum Neutrino Energy

   real min_neutrino_energy(real kinetic_energy) {
   	return kinetic_energy / 2 * (1 + sqrt(1 + 2 * m_osmium_target / kinetic_energy));
   }

// Differential Cross Section

   real differential_cross_section(real neutrino_energy, real kinetic_energy) {
   	real prefactor;
	real fraction;
	prefactor <- square(g_f) * m_osmium_target * square(q_w) / (4 * pi());
	fraction <- 1 - m_osmium_target * kinetic_energy / (2 * square(neutrino_energy))
	if (neutrino_energy < min_neutrino_energy(kinetic_energy)) {
	   return 0;
	}
	return prefactor * fraction;
   }

// Neutrino Oscillations

   real oscillations(real sin2theta_s, real delta_m, real neutrino_energy, real radius) {
   	return 1 - square(sin2theta_s) * square(sin(12700 * square(delta_m) * radius / neutrino_energy));
   }

// Beta Nu Spectrum

   real beta_nu_spectrum(real neutrino_energy, real Q) {
   	real eta;
	real xi;
	real K;
	real p;
	real S;
   	eta <- neutrino_energy / m_electron();
	xi <- Q / m_electron();
	K <- xi - eta +1;
	p <- sqrt(square(K) - 1);
	if (eta > xi) {
	   return 0;
	}
	S <- K * p * square(eta);
	return S;
   }



// Signal

   real signal(real kinetic_energy, real radius, real neutrino_energy, real A, real detector_mass, real Q, real sin2theta_s, real delta_m) {
	return target(A, detector_mass) * flux_constant * pow(radius, -2) * differential_cross_section(neutrino_energy, kinetic_energy) * beta_nu_spectrum(neutrino_energy, Q) * oscillations(sin2theta_s, delta_m, neutrino_energy, radius);
   }

// Maximum Kinetic Energy

   real tmax(real neutrino_energy) {
   	return neutrino_energy / (1 + m_osmium_target / (2 * neutrino_energy));
   }

}

parameters {

	real <upper = 1> sin2theta_s;
	real <lower = 0> delta_m;
	real <lower = 0, upper = Q> neutrino_energy;
	real <lower = threshhold, upper = tmax(neutrino_energy)> kinetic_energy;
	real <lower = min_r, upper = max_r> radius;   

}

transformed parameters {

	real square(sin2theta_s);
	real square(delta_m);
	
}

data {

     real Q;								// in eV	
     real threshhold;							// in eV
     real m_osmium_target;						// in eV
     real min_r;							// in cm
     real max_r;							// in cm
     real detector_mass;						// in g

}

model {

      real signal_strength;
      signal_strength <- signalreal kinetic_energy, real radius, real neutrino_energy, real A, real detector_mass, real Q, real sin2theta_s, real delta_m);
      increment_log_prob(log(signal_strength));

}