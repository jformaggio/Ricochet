
functions{

// Set up constants 

	real m_electron() { return  510998.910;}				 // Electron mass in eV
	real c() { return  299792458.;}			   	   	    	 // Speed of light in m/s
	real hbarc() { return 0.0000197327;}			         // Reduced Planck constant*speed of light in centimeter*eV
	real g_f() { return 2.30158e-28;}					 // GF in centimeter/eV
	real n_avogadro() { return 6.02214e23;}				 // Avogadro's number
	real weakangle() { return 0.22296;}			            	 // weak mixing angle

// Target

   real target(real detector_mass, real A) {
   	return n_avogadro() * detector_mass / A;
   }

// Calculate weak charge

   real q_w(real A, real Z) {
   	return (A/2. - Z*(1.-4.*weakangle()));
   }


// Minimum Neutrino Energy

   real min_neutrino_energy(real kinetic_energy, real m_target) {
   	return kinetic_energy / 2 * (1 + sqrt(1 + 2 * m_target/ kinetic_energy));
   }

// Differential Cross Section

   real differential_cross_section(real neutrino_energy, real kinetic_energy, real m_target, real targetA, real targetZ) {
   	real prefactor;
	real fraction;
	prefactor <- square(g_f()) * m_target * square(q_w(targetA, targetZ)) / (4 * pi());
	fraction <- 1 - m_target * kinetic_energy / (2 * square(neutrino_energy));
	return prefactor * fraction;
   }

// Neutrino Oscillations

   real oscillations(real sin2theta_s, real delta_m2, real neutrino_energy, real radius) {
   	return 1. - sin2theta_s * square(sin(12700 * delta_m2 * radius / neutrino_energy));
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


// Maximum Kinetic Energy

   real tmax(real neutrino_energy, real m_A) 
   {
   	return neutrino_energy / (1 + m_A / (2 * neutrino_energy));
   }

}

data {

     real Q;								// in eV	
     real threshold;							// in eV
     real targetMass;							// in eV
     real min_r;								// in cm
     real max_r;								// in cm
     real detector_mass;						// in g
     real targetA;							// the target mass number
     real targetZ;						       	// the target atomic number

}

transformed data {

	real t_upper_limit;
	real e_lower_limit;

	t_upper_limit <- Q/(1. + targetMass/(2*Q));
	e_lower_limit <- threshold/2.*(1.+sqrt(1.+2.*targetMass/threshold));

}
parameters {

//	real <upper = 1.> sin2theta_s;
//	real <lower = 0.> delta_m2;
	real <lower = e_lower_limit, upper = Q> neutrino_energy;
	real <lower = threshold, upper = t_upper_limit> kinetic_energy;
	real <lower = min_r, upper = max_r> radius;   

}

transformed parameters {

	 real theFlux;
	 real theSpectrum;
	 real theCrossSection;
	 real theOscillation;
	 real theRate;

	 theFlux <- 1. / (4. * pi() * square(radius));

	 theSpectrum <- beta_nu_spectrum(neutrino_energy, Q);

	 theCrossSection <- differential_cross_section(neutrino_energy, kinetic_energy, targetMass, targetA, targetZ);

//	 theOscillation <- oscillations(sin2theta_s, delta_m2, neutrino_energy, radius);

	 theRate <- target(targetMass, targetA) * theFlux * theSpectrum * theCrossSection;

}

model {

        increment_log_prob(log(theRate));

      	neutrino_energy ~ uniform(min_neutrino_energy(threshold, targetMass), Q);

        kinetic_energy ~ uniform(threshold, tmax(neutrino_energy, targetMass));

}
