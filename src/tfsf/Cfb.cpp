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

//Defines the class "Cfb" for a TF/SF focused-beam source in free space

#include "headers.h"

#include "Cfb.h"

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


Cfb::Cfb(const FBDataType& MyData, const int& Index)
		:Data(MyData), FPIndex(Index)
{
	/** incidence, rotation and polarization angles of the beam **/
	theta = Data.THETA;
	phi = Data.PHI;
	psi = Data.PSI;


	if ((cos(Data.ap_half_angle)<0)&&(sin(Data.ap_half_angle)<0))
	{
#ifdef __GNUG__
//GNU C++ compiler is being used, use the nice predefined variables for the function name
//		InvalidNumericArgumentException<double> exc(__PRETTY_FUNCTION__,1);
		string func_name = __FUNCTION__;
#else
		string func_name = "";
#endif
		throw AngoraInvalidArgumentExceptionWithType<double>(func_name,Data.ap_half_angle,
			"(should be between 0 and 90deg)");
	}

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

	//maximum x and y extents of the TF/SF box
	// for lambda_min and lambda_max, -40 dB might be too tight, so we use -20dB
	double lambda_min = c/sqrt(epsilon_r_i*mu_r_i)/(Data.waveform->w_max_20()/2/M_PI);
	double lambda_max = c/sqrt(epsilon_r_i*mu_r_i)/(Data.waveform->w_min_20()/2/M_PI);

	double minimum_aliasing_distance_in_beam_width = 5; //must be larger than 2
	double beam_waist = minimum_aliasing_distance_in_beam_width*lambda_max/abs(sin(Data.ap_half_angle));  //FIXME:should find a better way for this later

	if (minimum_aliasing_distance_in_beam_width<5)
		if (rank==0) cout << "Warning: There may be aliasing!!!!" << endl;

	//x and y extents are determined either by the beam waist or the size of the TF/SF box
	double W_x = max(beam_waist,(NCELLS_X-(Data.FBMarginFrontX+Data.FBMarginBackX))*dx);
	double W_y = max(beam_waist,(NCELLS_Y-(Data.FBMarginRightY+Data.FBMarginLeftY))*dx);

	if (Data.angular_discretization=="cartesian")
	{
		if (Data.n_1==-1)
		{
			//uninitialized in config file, use sampling theorem for default spacing
			dsx = lambda_min/W_x;
			//The range 0->(dsx*N_X1) is divided into N_X1 regions and sx is placed at the midpoint of each region.
			//This provides consistency with previous publications.
			//The point (sx,sy) on the 2D plane may fall out of the sin(th)<1 circle for N_X1=N_X2=2, but neither N_X1 nor N_X2 are supposed to be that small anyway. Note that dsx<=lambda_min/beam_waist=sin(th_max)/(minimum_aliasing_distance_in_beam_width)*(lambda_min/lambda_max)
			//						   <sin(th_max)/(minimum_aliasing_distance_in_beam_width)
			// Therefore N_X1 is larger than (minimum_aliasing_distance_in_beam_width), which is larger than 2.
			N_X1 = int(2*abs(sin(Data.ap_half_angle))/dsx)+1;
		}
		else
		{
			N_X1 = Data.n_1;
			if (N_X1<=0)
			{
				if (rank==0)
				{
					cout << "Error: Invalid number of x-direction cosines (" << N_X1 << ")" << endl;
				}
				exit(-1);
			}
			dsx = 2*abs(sin(Data.ap_half_angle))/N_X1;
		}
		if (Data.n_2==-1)
		{
			//uninitialized in config file, use sampling theorem for default spacing
			dsy = lambda_min/W_y;
			//(see note above for sx)
			N_X2 = int(2*abs(sin(Data.ap_half_angle))/dsy)+1;
		}
		else
		{
			N_X2 = Data.n_2;
			if (N_X2<=0)
			{
				if (rank==0)
				{
					cout << "Error: Invalid number of y-direction cosines (" << N_X2 << ")" << endl;
				}
				exit(-1);
			}
			dsy = 2*abs(sin(Data.ap_half_angle))/N_X2;
		}

		//form the 2D theta and phi arrays for the incidence directions
		angle_within_ill_cone.resize(N_X1,N_X2);
		angle_within_ill_cone = true;
		theta_array.resize(N_X1,N_X2);
		phi_array.resize(N_X1,N_X2);
		pwfactor.resize(N_X1,N_X2);

		for (int sx_index=0; sx_index<N_X1; sx_index++)
		{
			//(see note above for the placement of sx and sy)
			sx = dsx*(sx_index-(N_X1-1.0)/2.0);
			for (int sy_index=0; sy_index<N_X2; sy_index++)
			{
				//(see note above for the placement of sx and sy)
				sy = dsy*(sy_index-(N_X2-1.0)/2.0);
				if (((pow2(sx)+pow2(sy))>1) //dir1^2+dir2^2 = sin^2(theta_pw) should be less than or equal to 1 for a valid sin(theta_pw)
				   ||((pow2(sx)+pow2(sy))>pow2(abs(sin(Data.ap_half_angle)))))	// it should also be less than [maximum allowable direction cosine]^2
				{
					angle_within_ill_cone(sx_index,sy_index) = false;
				}
				else
				{
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

					theta_array(sx_index,sy_index) = theta_pw;
					phi_array(sx_index,sy_index) = phi_pw;

					//amplitude factor for every plane wave
					pwfactor(sx_index,sy_index) = dsx*dsy/sz // = infinitesimal solid angle = sin(theta)*d_theta*d_phi
												 *sqrt(abs(cos(theta_pw)))  //radiometric factor
												 *(Data.f/(2*M_PI*c));  //other factors (from Richards&Wolf)
				}
			}
		}
	}
	else if (Data.angular_discretization=="radial")
	{
		// Data.n_1 and Data.n_2 *must* be initialized: This should be ensured in read_tfsf
		N_X1 = Data.n_1;
		N_X2 = Data.n_2;
		//spacing of phi_pw (between [0,pi])
		d_phi = M_PI/N_X2;

		//Gauss-Legendre parameters for rho quadrature [ rho=sin(theta_pw), -sin(ap_half_angle)<rho<sin(ap_half_angle) ]
		GLx.resize(N_X1); //Gauss-Legendre coordinates
		GLw.resize(N_X1); //Gauss-Legendre weights

		//calculate GL quadrature rule for the sin(theta_pw) values (between [-sin(ap_half_angle),sin(ap_half_angle)])
		gaussquadrule(N_X1,GLx,GLw);
		//rescale the positions and weights to the correct ranges of sin(theta_pw) (between [-sin(ap_half_angle),sin(ap_half_angle)])
		//see http://en.wikipedia.org/wiki/Gaussian_quadrature for the following formulas
		//GLx and GLw are for sin(theta_pw)
		GLx = abs(sin(Data.ap_half_angle))*(-GLx);	//-GLx because it is calculated in reverse direction in gaussquadrule
		GLw = abs(sin(Data.ap_half_angle))*GLw;

		//form the 2D theta_pw and phi arrays for the incidence directions
		angle_within_ill_cone.resize(N_X1,N_X2); //always true for "radial"
		angle_within_ill_cone = true;
		theta_array.resize(N_X1,N_X2);
		phi_array.resize(N_X1,N_X2);
		pwfactor.resize(N_X1,N_X2);

		for(int rho_index=0; rho_index<N_X1; rho_index++)
		{
			for(int phi_index=0; phi_index<N_X2; phi_index++)
			{
				//the plane wave is always incident from the upper half space
				theta_pw = asin(GLx(rho_index)); //GLx and GLw are for sin(theta_pw), result is between [-pi/2,pi/2]
				phi_pw = (phi_index+0.5)*d_phi;

				/** convert to traditional spherical angles **/
				//This is important, because the polarizations of each plane wave below are determined with respect to the "phi" angle,
				//which only works if "theta_pw" is between 0 and pi.
				if (sin(theta_pw)<0)
				{
					theta_pw = 2*M_PI-theta_pw;
					phi_pw += M_PI;
				}
				theta_array(rho_index,phi_index) = theta_pw;
				phi_array(rho_index,phi_index) = phi_pw;
				//amplitude factor for every plane wave
				pwfactor(rho_index,phi_index) = abs(sin(theta_pw))*GLw(rho_index)*d_phi  //could also use GLx(rho_index) for sin(theta_pw)
												*sqrt(abs(cos(theta_pw)))  //radiometric factor
												*(Data.f/(2*M_PI*c));  //other factors (from Richards&Wolf)
			}
		}
	}
	else
	{
		if (rank==0)
		{
			cout << "Developer error: Angular discretization type (angular_discretization) should be either \"cartesian\" or \"radial\"" << endl;
		}
		exit(-1);
	}

	//start adding the plane waves
	for (int i=0;i<N_X1;i++)
	{
		for (int j=0;j<N_X2;j++)
		{
			if (angle_within_ill_cone(i,j))
			{//do not add PW if it is outside illumination cone
				theta_pw = theta_array(i,j);
				phi_pw = phi_array(i,j);

				double polarization = psi;  //psi is defined like phi with respect to the local x axis
				//psi_pw is adjusted to give zero cross-pol component
				double angle_meridional = phi_pw - polarization; //angle with respect to the meridional plane
				// the angle with the meridional plane translates differently to psi_pw depending on the incidence half space
				if (cos(theta_pw)>=0) //incident from upper half space
				{
					psi_pw = 3*M_PI/2 - angle_meridional;
				}
				else //incident from lower half space
				{
					psi_pw = M_PI/2 + angle_meridional; //pw_pol is ccw wrt the +x axis, like the spherical phi_pw angle
				}

				/****************************************************************/
				/** Rotate the incidence directions and electric field vectors **/
				/****************************************************************/
				//ccw (w.r.t. z) azimuthal rotation angle from the global x axis to the local x axis of the beam
				//The focused beam is rotationally symmetric around its propagation axis (except the polarization), therefore the local x axis is assumed parallel to the xy plane without any loss of generality.
				double rot_az = ((sin(theta)>=0)?phi+M_PI/2:phi-M_PI/2);

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

				PWData.PWMarginBackX = Data.FBMarginBackX;
				PWData.PWMarginFrontX = Data.FBMarginFrontX;
				PWData.PWMarginLeftY = Data.FBMarginLeftY;
				PWData.PWMarginRightY = Data.FBMarginRightY;
				PWData.PWMarginLowerZ = Data.FBMarginLowerZ;
				PWData.PWMarginUpperZ = Data.FBMarginUpperZ;
				PWData.PWOriginX = Data.FBOriginX;
				PWData.PWOriginY = Data.FBOriginY;
				PWData.PWOriginZ = Data.FBOriginZ;
				PWData.L_req = Data.L_req;

				PWData.DisplayWarnings = Data.DisplayWarnings;

				PWData.E0 = Data.E0*pwfactor(i,j);
																		//abs(sin(theta_pw)) because sin(theta_pw) can go negative in radial quadrature
				//waveform of the constituent plane waves is the derivative of the incident plane-wave waveform
				PWData.waveform = Data.waveform->Derivative();

				//Add the plane wave source to the focused-beam object
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
						throw FBGrazingIncidenceException();
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
						throw FBGrazingIncidenceException();
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

void Cfb::CorrectE(const int& n)
{//applies the E-field corrections on the TF/SF box due to the focused beam.
//Currently plane waves are used, but this can change in the future
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->CorrectE(n);		//apply correction to the main-grid E-field (using Hinc at n+1/2)
	}
}

void Cfb::CorrectH(const int& n)
{//applies the H-field corrections on the TF/SF box due to the focused beam.
//Currently plane waves are used, but this can change in the future
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->CorrectH(n);		//apply correction to the main-grid H-field (using Einc at n+1)
	}
}

void Cfb::WriteScatteredPWDirections(Array<double,1>& PW_THETA, Array<double,1>& PW_PHI) const
{//write the scattering angles (THETA and PHI) of the scattered PWs into PW_THETA and PW_PHI
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->WriteScatteredPWDirection(PW_THETA,PW_PHI);
	}
}

void Cfb::WriteScatteredPWDelaysFromOrigin(Array<double,1>& origindelay_array, const double& FFOriginX, const double& FFOriginY, const double& FFOriginZ) const
{//write the delays (from the origin) of the scattered PWs into origindelay_array
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->WriteScatteredPWDelayFromOrigin(origindelay_array,FFOriginX,FFOriginY,FFOriginZ);
	}
}

void Cfb::WriteScatteredPWFieldAmplitudes(Array<double,1>& E_x_array, Array<double,1>& E_y_array) const
{//write the field amplitudes of the scattered PWs into field_array
	for (int i=0; i<NumberOfPlaneWaves(); i++)
	{
		PlaneWaves[i]->WriteScatteredPWFieldAmplitude(E_x_array,E_y_array);
	}
}
