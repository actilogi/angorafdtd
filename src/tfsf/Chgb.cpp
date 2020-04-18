/* AUTORIGHTS
Copyright (C) 2006-2012  Ilker R. Capoglu

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

//Defines the class "Chgb" for a TF/SF Hermite-Gaussian source in free space

#include "headers.h"

#include "Chgb.h"

#include "Cpw_fs.h"
#include "Cpw_2l.h"
#include "Cpw_ml.h"

//Uses TinyVector operations
#include <blitz/tinyvec2.h>

//definition of Cwf needed
#include "waveforms/Cwf.h"

extern double dx;
extern int NCELLS_X,NCELLS_Y;

extern double epsilon_r_upper,mu_r_upper,epsilon_r_lower,mu_r_lower;

extern int rank;

extern int number_of_layers;

//gauss-legendre quadrature rule generator
extern void gaussquadrule(const int& n, Array<double,1>& x, Array<double,1>& w);

extern double Hermite(const int& N, const double& x);

extern void MPI_exit(const int& exitcode);

//factorial function (from Cgsmb)
extern int factorial (int num);


Chgb::Chgb(const HGBDataType& MyData, const int& Index)
		:Data(MyData), HGPIndex(Index)
{
	/** incidence, rotation and polarization angles of the beam **/
	theta = Data.THETA;
	phi = Data.PHI;
	alpha = Data.alpha;
	psi = Data.PSI;

	if (cos(theta)>0)
	{//incident from the upper half-space
		epsilon_r_i = epsilon_r_upper;
		mu_r_i = mu_r_upper;
	}
	else
	{//incident from the lower half-space
		epsilon_r_i = epsilon_r_lower;
		mu_r_i = mu_r_lower;
	}

	//half width of the angular spectrum is assumed independent of frequency
	//consequently, the spatial beam width scales with the wavelength
	// "beam_half_width" is given for the center wavelength of the excitation
	double k_center = Data.waveform->W_0()/c;  //this is for free space
	//the wavenumber times the beam half-width in the incidence half space
	double k_times_beam_half_width = k_center*sqrt(epsilon_r_i*mu_r_i)*Data.beam_half_width;  //THIS IS ASSUMED TO BE CONSTANT FOR ALL WAVELENGTHS!
	double s_halfwidth = 1/k_times_beam_half_width; //half-width in the angle space

	//approximate maxima of the Hermite polynomials
	double hermite_max_x = sqrt(pow(2.0,Data.x_order)*factorial(Data.x_order));
	double hermite_max_y = sqrt(pow(2.0,Data.y_order)*factorial(Data.y_order));

	double angular_threshold = 1e-3;

	//the maximum direction cosine that needs to be considered (it could be at most 1)
	//do rough, brute-force calculation
	double sx_loop=(4+Data.x_order)*s_halfwidth;  //As the order increases, the decay is slower [4+Data.x_order is still a bit too far but..]
	double sx_step = 0.01*sx_loop; //step backward by %1 of the starting point
	double sy_loop=(4+Data.y_order)*s_halfwidth;  //As the order increases, the decay is slower [4+Data.y_order is still a bit too far but..]
	double sy_step = 0.01*sy_loop; //step backward by %1 of the starting point
	double ratio_to_maximum_x,ratio_to_maximum_y; //ratio of the spectrum at a certain point to the (roughly) maximum value
	do
	{//decrease sx until the threshold is found
		sx_loop-=sx_step;
		ratio_to_maximum_x = abs(Hermite(Data.x_order,sx_loop/s_halfwidth)*exp(-pow2(sx_loop/s_halfwidth)/2)/hermite_max_x);
	} while (ratio_to_maximum_x<angular_threshold);
	do
	{//decrease sy until the threshold is found
		sy_loop-=sy_step;
		ratio_to_maximum_y = abs(Hermite(Data.y_order,sy_loop/s_halfwidth)*exp(-pow2(sy_loop/s_halfwidth)/2)/hermite_max_y);
	} while (ratio_to_maximum_y<angular_threshold);

	double s_max = min(max(sx_loop,sy_loop),1.0);

	if (sx_loop>1.0)
	{
		if (rank==0)
		{
			cout << "Warning: The spatial spectrum of Hermite-Gaussian beam " << HGPIndex << " has non-propagating components (spectrum falls to " << abs(Hermite(Data.x_order,s_max/s_halfwidth)*exp(-pow2(s_max/s_halfwidth)/2)/hermite_max_x)*100 << "% of its maximum at the edge of the propagation circle in the x direction)" << endl << endl;
		}
	}
	if (sy_loop>1.0)
	{
		if (rank==0)
		{
			cout << "Warning: The spatial spectrum of Hermite-Gaussian beam " << HGPIndex << " has non-propagating components (spectrum falls to " << abs(Hermite(Data.y_order,s_max/s_halfwidth)*exp(-pow2(s_max/s_halfwidth)/2)/hermite_max_y)*100 << "% of its maximum at the edge of the propagation circle in the y direction)" << endl << endl;
		}
	}

	//maximum k in the excitation waveform
	double k_max = (Data.waveform->w_max_20())/(c/sqrt(epsilon_r_i*mu_r_i));

	//maximum beam half-width among all the wavelengths (direction-cosine spacing is determined by this width)
	double minimum_necessary_shift_in_beam_half_width = 10;
	double k_times_beam_waist = minimum_necessary_shift_in_beam_half_width*k_times_beam_half_width;

	if (minimum_necessary_shift_in_beam_half_width<8)
		if (rank==0) cout << "Warning: There may be aliasing!!!!" << endl;

	//x and y extents are determined either by the beam waist or the size of the TF/SF box
	double TFSF_x_extent = (NCELLS_X-(Data.HGBMarginFrontX+Data.HGBMarginBackX))*dx;
	double TFSF_y_extent = (NCELLS_Y-(Data.HGBMarginRightY+Data.HGBMarginLeftY))*dx;
	double k_times_W_x = max(k_times_beam_waist,k_max*TFSF_x_extent);
	double k_times_W_y = max(k_times_beam_waist,k_max*TFSF_y_extent);

	//uninitialized in config file, use sampling theorem for default spacing
	dsx = 2*M_PI/k_times_W_x;
	//The range 0->(dsx*N_X1) is divided into N_X1 regions and sx is placed at the midpoint of each region.
	//This provides consistency with previous publications.
	//The point (sx,sy) on the 2D plane may fall out of the sin(th)<1 circle for N_X1=N_X2=2, but neither N_X1 nor N_X2 are supposed to be that small anyway. Note that dsx<=2*pi/k_times_beam_waist=2*pi/(minimum_necessary_shift_in_beam_half_width*k_times_beam_half_width)
	//						   <(lambda_center/beam_half_width)/(minimum_necessary_shift_in_beam_half_width)
	//						   <1/minimum_necessary_shift_in_beam_half_width
	// Therefore N_X1 is larger than (minimum_necessary_shift_in_beam_half_width), which is larger than 2.
	N_X1 = int(2*s_max/dsx)+1;

	//uninitialized in config file, use sampling theorem for default spacing
	dsy = 2*M_PI/k_times_W_y;
	//(see note above for sx)
	N_X2 = int(2*s_max/dsy)+1;

	double sx,sy;
	double pw_factor; //extra amplitude factor applied to each plane wave (angle-dependent)
	//start adding the plane waves
	for (int i=0;i<N_X1;i++)
	{
		sx = dsx*(i-(N_X1-1.0)/2.0);
		for (int j=0;j<N_X2;j++)
		{
			sy = dsy*(j-(N_X2-1.0)/2.0);

			if ((pow2(sx)+pow2(sy))<=1)  //is there a direction corresponding to sx,sy?
			{
				/** determine the incidence and polarization angles of the individual plane wave **/
				//the plane wave is always incident from the upper half space
				double sz = real(sqrt((complex<double>)(1-pow2(sx)-pow2(sy)))); //cast into complex, in case there are roundoff errors for sx^2+sy^2=1

				theta_pw = acos(sz);	//sz = cos(theta_pw), result of acos is always in the [0,pi] range

				if ((abs(sx)<LIBSTD_DBL_EPSILON*100)&&(abs(sy)<LIBSTD_DBL_EPSILON*100)) //is theta_pw=0?
				{
					phi_pw = 0; //just assign an arbitrary value, since phi_pw is undefined here
				}
				else
				{
					phi_pw = atan2(sy,sx);
				}
				//atan2 is between [-pi,pi]
				//if less than 0, make phi_pw positive
				if (phi_pw<0)
				{
					phi_pw += 2*M_PI;
				}

				double polarization = psi;  //psi is defined like phi with respect to the local x axis
				double angle_meridional = phi_pw - polarization; //angle with respect to the meridional plane
				// the angle with the meridional plane translates differently to psi_pw depending on the incidence half space
				//psi_pw is adjusted to give zero cross-pol component
				//This is different from Cfb, where the relationship is much simpler
				if (cos(theta_pw)>=0) //incident from upper half space
				{
					psi_pw = M_PI + atan2(cos(angle_meridional),cos(theta_pw)*sin(angle_meridional)); //output of atan2 is between [-pi,pi]
				}
				else //incident from lower half space
				{
					psi_pw = atan2(cos(angle_meridional),cos(theta_pw)*sin(angle_meridional)); //output of atan2 is between [-pi,pi]
				}
				/** incidence and polarization angles of the individual plane wave are determined. **/
				pw_factor = sqrt(pow2(sin(angle_meridional))+pow2(cos(angle_meridional)/cos(theta_pw)))*//amplitude of pw component such that the x-projection has the desired value, not the total pw amplitude
						dsx*dsy*
						pow2(k_center/(2*M_PI))* //k^2 is to convert the angular spectrum into inverse (spatial) Fourier transform, (2pi)^2 is the 2D inverse-Fourier transform coefficient
						(2*M_PI*pow2(Data.beam_half_width))*  //(pi/c) factor in the angular spectrum [the i^{m+n} factor is taken care of below]
						Hermite(Data.x_order,sx/s_halfwidth)*exp(-pow2(sx/s_halfwidth)/2)*
						Hermite(Data.y_order,sy/s_halfwidth)*exp(-pow2(sy/s_halfwidth)/2);

				/****************************************************************/
				/** Rotate the incidence directions and electric field vectors **/
				/****************************************************************/
				//ccw (w.r.t. z) azimuthal rotation angle from the global x axis to the local x axis of the beam
				double rot_az = ((sin(theta)>=0)?phi+M_PI/2+alpha:phi-M_PI/2+alpha);

				Array<double,2> rot_theta(3,3),rot_alpha(3,3),rot_tot(3,3);
				// Matrix representing counter-clockwise (right-hand) rotation w.r.t. the z axis by rot_az
				//(http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle)
				rot_alpha = cos(rot_az) , -sin(rot_az) , 0,
							sin(rot_az) ,  cos(rot_az) , 0,
									0  ,          0  , 1;

				double ux = -sin(phi);
				double uy = cos(phi);
				// Matrix representing counter-clockwise (right-hand) rotation w.r.t. the axis (ux,uy,0) by theta
				//(http://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle)
				rot_theta = (cos(theta)+pow2(ux)*(1-cos(theta))) , ux*uy*(1-cos(theta))             , uy*sin(theta),
							ux*uy*(1-cos(theta))             , (cos(theta)+pow2(uy)*(1-cos(theta))) , -ux*sin(theta),
							-uy*sin(theta)                   , ux*sin(theta)                    , cos(theta);
				//net rotation matrix
				firstIndex iii;
				secondIndex jjj;
				thirdIndex kkk;
				//azimuthal rotation is done before the theta rotation
				rot_tot = sum(rot_theta(iii,kkk)*rot_alpha(kkk,jjj),kkk);

				TinyVector<double,3> rot_tot_row1,rot_tot_row2,rot_tot_row3;
				rot_tot_row1 = (rot_tot(0,0)), (rot_tot(0,1)), (rot_tot(0,2));
				rot_tot_row2 = (rot_tot(1,0)), (rot_tot(1,1)), (rot_tot(1,2));
				rot_tot_row3 = (rot_tot(2,0)), (rot_tot(2,1)), (rot_tot(2,2));

				TinyVector<double,3> dircos;
				dircos = (sin(theta_pw)*cos(phi_pw)) , (sin(theta_pw)*sin(phi_pw)) , (cos(theta_pw));
				TinyVector<double,3> k_inc = -dircos;
				// Note that sin(theta_pw) is always >=0, since theta_pw is always between 0 and pi
				TinyVector<double,3> k_inc_lateral;
				k_inc_lateral = (-cos(phi_pw)),(-sin(phi_pw)),(0.0);

				TinyVector<double,3> unit_z;
				unit_z = (0.0),(0.0),(1.0);
				TinyVector<double,3> k_E;
				k_E = cos(psi_pw)*cross(k_inc_lateral,unit_z)+sin(psi_pw)*cross(cross(k_inc_lateral,unit_z),k_inc);
				TinyVector<double,3>  E;
				E = (k_E(0)), (k_E(1)), (k_E(2));

				TinyVector<double,3> dircos_rot;
				dircos_rot = (dot(rot_tot_row1,dircos)), (dot(rot_tot_row2,dircos)), (dot(rot_tot_row3,dircos));
				TinyVector<double,3> E_rot;
				E_rot = (dot(rot_tot_row1,E)), (dot(rot_tot_row2,E)), (dot(rot_tot_row3,E));

				double theta_rot = acos(dircos_rot(thirdDim));
				double phi_rot = atan2(dircos_rot(secondDim),dircos_rot(firstDim));
				TinyVector<double,3> k_inc_rot,k_inc_rot_lateral;
				k_inc_rot = (-sin(theta_rot)*cos(phi_rot)) , (-sin(theta_rot)*sin(phi_rot)) , (-cos(theta_rot));
				if (sin(theta_rot)>=0)//theta_rot could be any wild direction
				{
					k_inc_rot_lateral = (-cos(phi_rot)),(-sin(phi_rot)),(0.0);
				}
				else
				{
					k_inc_rot_lateral = (cos(phi_rot)),(sin(phi_rot)),(0.0);
				}
				double psi_rot = atan2(dot(E_rot,cross(cross(k_inc_rot_lateral,unit_z),k_inc_rot)),dot(E_rot,cross(k_inc_rot_lateral,unit_z)));

				/*************************************************************/
				/** Incidence directions and electric field vectors rotated **/
				/*************************************************************/

				PWDataType PWData;
				PWData.THETA = theta_rot;
				PWData.PHI = phi_rot;
				PWData.PSI = psi_rot;

				PWData.PWMarginBackX = Data.HGBMarginBackX;
				PWData.PWMarginFrontX = Data.HGBMarginFrontX;
				PWData.PWMarginLeftY = Data.HGBMarginLeftY;
				PWData.PWMarginRightY = Data.HGBMarginRightY;
				PWData.PWMarginLowerZ = Data.HGBMarginLowerZ;
				PWData.PWMarginUpperZ = Data.HGBMarginUpperZ;
				PWData.PWOriginX = Data.HGBOriginX;
				PWData.PWOriginY = Data.HGBOriginY;
				PWData.PWOriginZ = Data.HGBOriginZ;
				PWData.L_req = Data.L_req;

				PWData.DisplayWarnings = Data.DisplayWarnings;

				PWData.E0 = Data.E0*pw_factor;

				//waveform of the constituent plane waves is the integral of the image-plane waveform
				// every (i) multiplier in the i^{m+n} factor corresponds to a 90deg phase shift (-1 times the Hilbert transform)
				switch((Data.x_order+Data.y_order) % 4) {
				case  0  : PWData.waveform = Data.waveform;
						break;
				case  1  : PWData.waveform = Data.waveform->HilbertTransform();
						PWData.E0 *= -1;
						break;
				case  2  : PWData.waveform = Data.waveform;
						PWData.E0 *= -1;
						break;
				case  3  : PWData.waveform = Data.waveform->HilbertTransform();
						break;
				}

				//Add the plane wave source to the Hermite-Gaussian-beam object
				if (number_of_layers==1)
				{//attach a free-space plane-wave source
					boost::shared_ptr<Cpw> new_pw_ptr(new Cpw_fs(PWData));
					PlaneWaves.push_back(new_pw_ptr);
				}
				else if (number_of_layers==2)
				{//attach a 2-layered-medium plane-wave source
					if (cos(theta)*cos(theta_rot)<ANGORA_CPW_ML_MIN_GRAZING_COS_Z)
					{//grazing angle too low, some plane waves fall within opposite half space
					 //this is not physical, so don't allow it
						throw HGBGrazingIncidenceException();
					}
					else
					{
						boost::shared_ptr<Cpw> new_pw_ptr(new Cpw_2l(PWData));
						PlaneWaves.push_back(new_pw_ptr);
					}
				}
				else
				{//attach a multilayered-medium plane-wave source
					if (cos(theta)*cos(theta_rot)<ANGORA_CPW_ML_MIN_GRAZING_COS_Z)
					{//grazing angle too low, some plane waves fall within opposite half space
					 //this is not physical, so don't allow it
						throw HGBGrazingIncidenceException();
					}
					else
					{
						boost::shared_ptr<Cpw> new_pw_ptr(new Cpw_ml(PWData));
						PlaneWaves.push_back(new_pw_ptr);
					}
				}
			}
		}
	}
}

void Chgb::CorrectE(const int& n)
{//applies the E-field corrections on the TF/SF box due to the Hermite-Gaussian beam.
//Currently plane waves are used, but this can change in the future
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->CorrectE(n);		//apply correction to the main-grid E-field (using Hinc at n+1/2)
	}
}

void Chgb::CorrectH(const int& n)
{//applies the H-field corrections on the TF/SF box due to the Hermite-Gaussian beam.
//Currently plane waves are used, but this can change in the future
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->CorrectH(n);		//apply correction to the main-grid H-field (using Einc at n+1)
	}
}

void Chgb::WriteScatteredPWDirections(Array<double,1>& PW_THETA, Array<double,1>& PW_PHI) const
{//write the scattering angles (THETA and PHI) of the scattered PWs into PW_THETA and PW_PHI
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->WriteScatteredPWDirection(PW_THETA,PW_PHI);
	}
}

void Chgb::WriteScatteredPWDelaysFromOrigin(Array<double,1>& origindelay_array, const double& FFOriginX, const double& FFOriginY, const double& FFOriginZ) const
{//write the delays (from the origin) of the scattered PWs into origindelay_array
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->WriteScatteredPWDelayFromOrigin(origindelay_array,FFOriginX,FFOriginY,FFOriginZ);
	}
}

void Chgb::WriteScatteredPWFieldAmplitudes(Array<double,1>& E_x_array, Array<double,1>& E_y_array) const
{//write the field amplitudes of the scattered PWs into field_array
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->WriteScatteredPWFieldAmplitude(E_x_array,E_y_array);
	}
}
