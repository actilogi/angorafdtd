/* AUTORIGHTS
Copyright (C) 2006-2018  Ilker R. Capoglu and Di Zhang

    This file is part of the Angora package.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

//Defines several routines that initialize some of the global variables used throughout the simulation.

#include "headers.h"

#include "init_disp.h"

#include "material/Cmat_types.h"
//
extern bool check_mode;

extern const int vacuum;

extern Array<double,4> J_p_x,J_p_y,J_p_z;
extern Array<double,4> Jm1_p_x,Jm1_p_y,Jm1_p_z;
extern Array<double,3> Ex_m1,Ey_m1,Ez_m1;
extern Array<update_coeff_type,3> Cc_X,Cc_Y,Cc_Z;
extern Array<update_coeff_type,4> alpha_X,xi_X,gamma_X,
                                  alpha_Y,xi_Y,gamma_Y,
                                  alpha_Z,xi_Z,gamma_Z;
extern Array<update_coeff_type,4> Pa_X,Pb_X,Pa_Y,Pb_Y,Pa_Z,Pb_Z;
extern Array<omega_p_x_type,4> omega_p_x_indices;
extern Array<omega_p_y_type,4> omega_p_y_indices;
extern Array<omega_p_z_type,4> omega_p_z_indices;
extern Array<tau_p_x_type,4> tau_p_x_indices;
extern Array<tau_p_y_type,4> tau_p_y_indices;
extern Array<tau_p_z_type,4> tau_p_z_indices;
extern Array<Omega_p_x_type,4> Omega_p_x_indices;
extern Array<Omega_p_y_type,4> Omega_p_y_indices;
extern Array<Omega_p_z_type,4> Omega_p_z_indices;
extern int pole_dim_size;

extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;
//

void init_disp(const int& num_poles)
{
	//Initialize dispersion arrays
	//The size of the 4th dimension will be the number of poles,
	//given by num_poles

	//TODO (longer term): Accept a Drude pole parameter at a given
	//pole index for file input. Somehow communicate the maximum
	//pole index in the file read options to this routine before
	//actually placing the file into the grid, since placing any material
	//with dispersion requires properly allocated dispersion arrays.

	//TODO (longer term): Distinguish between Drude-only and Lorentz cases
	//and call the Drude-only update method if there are no Lorentz pole pairs

  //assign global variable
  pole_dim_size = num_poles;

  //no need to allocate if there are no poles
  if (num_poles==0)
    return;

	// do not initialize if in check mode
	if (check_mode)
	  return;

//Jm1_p_x,Jm1_p_y,Jm1_p_z;
//Ex_m1,Ey_m1,Ez_m1;
//Cc_X,Cc_Y,Cc_Z;
//alpha_X,xi_X,gamma_X,
//alpha_Y,xi_Y,gamma_Y,
//alpha_Z,xi_Z,gamma_Z;

  Ex_m1.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1));
  Ey_m1.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1));
  Ez_m1.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper));

	Cc_X.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1));
  Cc_Y.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1));
  Cc_Z.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper));

	J_p_x.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1),Range(0,num_poles-1));
  J_p_y.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1),Range(0,num_poles-1));
  J_p_z.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper),Range(0,num_poles-1));
	Jm1_p_x.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1),Range(0,num_poles-1));
  Jm1_p_y.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1),Range(0,num_poles-1));
  Jm1_p_z.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper),Range(0,num_poles-1));

	alpha_X.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1),Range(0,num_poles-1));
  alpha_Y.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1),Range(0,num_poles-1));
  alpha_Z.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper),Range(0,num_poles-1));
	xi_X.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1),Range(0,num_poles-1));
  xi_Y.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1),Range(0,num_poles-1));
  xi_Z.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper),Range(0,num_poles-1));
	gamma_X.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1),Range(0,num_poles-1));
  gamma_Y.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1),Range(0,num_poles-1));
  gamma_Z.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper),Range(0,num_poles-1));
#if (0)
  // These can be used for Drude-only simulations in the future
  Pa_X.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1),Range(0,num_poles-1));
  Pb_X.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1),Range(0,num_poles-1));
  Pa_Y.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1),Range(0,num_poles-1));
  Pb_Y.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1),Range(0,num_poles-1));
  Pa_Z.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper),Range(0,num_poles-1));
  Pb_Z.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper),Range(0,num_poles-1));
#endif
  omega_p_x_indices.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1),Range(0,num_poles-1));
  omega_p_y_indices.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1),Range(0,num_poles-1));
  omega_p_z_indices.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper),Range(0,num_poles-1));
  tau_p_x_indices.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1),Range(0,num_poles-1));
  tau_p_y_indices.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1),Range(0,num_poles-1));
  tau_p_z_indices.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper),Range(0,num_poles-1));
  Omega_p_x_indices.resize(Range(iback,ifront),Range(jleft,jright+1),Range(klower,kupper+1),Range(0,num_poles-1));
  Omega_p_y_indices.resize(Range(iback,ifront+1),Range(jleft,jright),Range(klower,kupper+1),Range(0,num_poles-1));
  Omega_p_z_indices.resize(Range(iback,ifront+1),Range(jleft,jright+1),Range(klower,kupper),Range(0,num_poles-1));

  //Initialize polarization current arrays
  J_p_x=0;
  J_p_y=0;
  J_p_z=0;
  Jm1_p_x=0;
  Jm1_p_y=0;
  Jm1_p_z=0;
  //Initialize E-field history arrays
  Ex_m1=0;
  Ey_m1=0;
  Ez_m1=0;
#if (0)
  // These can be used for Drude-only simulations in the future

  // Polarization-current update coefficients for vacuum
  Pa_X=1.0;
  Pb_X=0.0;
  Pa_Y=1.0;
  Pb_Y=0.0;
  Pa_Z=1.0;
  Pb_Z=0.0;
#endif

  // (Lorentz) polarization-current update coefficients for vacuum
  alpha_X=1.0;
  xi_X=0.0;
  gamma_X=0.0;
  alpha_Y=1.0;
  xi_Y=0.0;
  gamma_Y=0.0;
  alpha_Z=1.0;
  xi_Z=0.0;
  gamma_Z=0.0;

  //Default parameters for vacuum
  omega_p_x_indices=vacuum;
  omega_p_y_indices=vacuum;
  omega_p_z_indices=vacuum;
  tau_p_x_indices=vacuum;
  tau_p_y_indices=vacuum;
  tau_p_z_indices=vacuum;
  Omega_p_x_indices=vacuum;
  Omega_p_y_indices=vacuum;
  Omega_p_z_indices=vacuum;
}
