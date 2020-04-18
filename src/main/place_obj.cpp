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

#include "headers.h"

#include "place_obj.h"

#include "material/Cmat.h"
#include "shape/Cshape.h"

extern double dt;

extern const int vacuum;

extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;

extern Array<update_coeff_type,3> Ca_X,Cb_X,Cc_X,Ca_Y,Cb_Y,Cc_Y,Ca_Z,Cb_Z,Cc_Z;
extern Array<update_coeff_type,3> Da_X,Db_X,Da_Y,Db_Y,Da_Z,Db_Z;

extern Array<eps_x_type,3> eps_x_indices;
extern Array<eps_y_type,3> eps_y_indices;
extern Array<eps_z_type,3> eps_z_indices;
extern Array<mu_x_type,3> mu_x_indices;
extern Array<mu_y_type,3> mu_y_indices;
extern Array<mu_z_type,3> mu_z_indices;
extern Array<cond_e_x_type,3> cond_e_x_indices;
extern Array<cond_e_y_type,3> cond_e_y_indices;
extern Array<cond_e_z_type,3> cond_e_z_indices;
extern Array<cond_h_x_type,3> cond_h_x_indices;
extern Array<cond_h_y_type,3> cond_h_y_indices;
extern Array<cond_h_z_type,3> cond_h_z_indices;

extern Array<float,1> eps_x,eps_y,eps_z;
extern Array<float,1> cond_e_x,cond_e_y,cond_e_z;
extern Array<float,1> mu_x,mu_y,mu_z;
extern Array<float,1> cond_h_x,cond_h_y,cond_h_z;

extern Array<bool,3> dispersion_exists_at_Ex_position,dispersion_exists_at_Ey_position,dispersion_exists_at_Ez_position;
extern Array<update_coeff_type,4> alpha_X,xi_X,gamma_X,alpha_Y,xi_Y,gamma_Y,alpha_Z,xi_Z,gamma_Z;
extern Array<omega_p_x_type,4> omega_p_x_indices;
extern Array<omega_p_y_type,4> omega_p_y_indices;
extern Array<omega_p_z_type,4> omega_p_z_indices;
extern Array<tau_p_x_type,4> tau_p_x_indices;
extern Array<tau_p_y_type,4> tau_p_y_indices;
extern Array<tau_p_z_type,4> tau_p_z_indices;
extern Array<Omega_p_x_type,4> Omega_p_x_indices;
extern Array<Omega_p_y_type,4> Omega_p_y_indices;
extern Array<Omega_p_z_type,4> Omega_p_z_indices;
extern Array<float,1> omega_p_x,omega_p_y,omega_p_z,
                      tau_p_x,tau_p_y,tau_p_z,
                      Omega_p_x,Omega_p_y,Omega_p_z;



void place_obj(const Cmat& material, const_Cshape_shared_ptr shapeptr)
{
	int i,j,k;

//	double BlockBack,BlockFront,BlockLeft,BlockRight,BlockLower,BlockUpper;

//	double BlockBack = shapeptr->bounding_box.back_limit;
//	double BlockFront = shapeptr->bounding_box.front_limit;
//	double BlockLeft = shapeptr->bounding_box.left_limit;
//	double BlockRight = shapeptr->bounding_box.right_limit;
//	double BlockLower = shapeptr->bounding_box.lower_limit;
//	double BlockUpper = shapeptr->bounding_box.upper_limit;

	//get the bounding box limits for quicker placement
	int BlockBack = shapeptr->bounding_box_back_cell();
	int BlockFront = shapeptr->bounding_box_front_cell();
	int BlockLeft = shapeptr->bounding_box_left_cell();
	int BlockRight = shapeptr->bounding_box_right_cell();
	int BlockLower = shapeptr->bounding_box_lower_cell();
	int BlockUpper = shapeptr->bounding_box_upper_cell();

	/** TODO: When the resizeAndPreserve() overload with Range arguments is developed,
	do the resizing here, instead of in init.cpp or initgeom.cpp **/
//	//resize the polarization current arrays as necessary
//	if ((material.omega_p_exists())&&(material.omega_p_value()!=0))
//	{
//		//x components
//		if (BlockBack<J_p_x.lbound(firstDim))
//		{
//			J_p_x.resize(Range(BlockBack,J_p_x.ubound(firstDim)),
//									Range(J_p_x.lbound(secondDim),J_p_x.ubound(secondDim)),
//									Range(J_p_x.lbound(thirdDim),J_p_x.ubound(thirdDim)));
//		}
//	}

  //total number of poles in the material
  //(This includes vanishing poles: One vanishing pole
  //is supposed to reset all dispersion arrays to zero.)
  //(No poles means transparent material.)
  int tot_num_lrntz_poles(material.tot_num_lrntz_poles());
  //number of effective (nonzero) poles in the material
  int eff_num_lrntz_poles(material.eff_num_lrntz_poles());
  //maximum number of poles in the simulation grid
  int max_num_lrntz_poles(omega_p_x_indices.extent(fourthDim));

	//place the bulk of the object
	for (i=max(iback,BlockBack); i<=min(ifront,BlockFront); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright+1,BlockRight+1); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper+1,BlockUpper+1); k++)
			{
				if (shapeptr->IsInside(i-0.5,j-1,k-1))
				{
				  eps_x_type &eps_x_idx(eps_x_indices(i,j,k));
				  cond_e_x_type &cond_e_x_idx(cond_e_x_indices(i,j,k));
					if (material.eps_x_exists())
					{
						eps_x_idx=material.eps_x_index();
					}
					if (material.cond_e_x_exists())
					{
						cond_e_x_idx=material.cond_e_x_index();
					}
					/********************************************************************************************************/
          if (eff_num_lrntz_poles>0)
            dispersion_exists_at_Ex_position(i,j,k) = true;
          else
            dispersion_exists_at_Ex_position(i,j,k) = false;

					if (max_num_lrntz_poles==0) {
					  //there is NO dispersion
						Ca_X(i,j,k)=(1-dt*cond_e_x(cond_e_x_idx)/(2.0*eps_x(eps_x_idx)*epsilon_0))/
						            (1+dt*cond_e_x(cond_e_x_idx)/(2.0*eps_x(eps_x_idx)*epsilon_0));
						Cb_X(i,j,k)=dt/eps_x(eps_x_idx)/epsilon_0/dx/(1+dt*cond_e_x(cond_e_x_idx)/(2.0*eps_x(eps_x_idx)*epsilon_0));
					}
					else {
						for (int p(0); p<max_num_lrntz_poles; ++p) {
							omega_p_x_type &omega_p_x_idx(omega_p_x_indices(i,j,k,p));
							tau_p_x_type &tau_p_x_idx(tau_p_x_indices(i,j,k,p));
							Omega_p_x_type &Omega_p_x_idx(Omega_p_x_indices(i,j,k,p));
							if (tot_num_lrntz_poles!=0) {//don't change anything if there are no poles (transparent material)
								if (p<tot_num_lrntz_poles) {
								  //place poles present in the material
									omega_p_x_idx=material.omega_p_x_index(p);
									tau_p_x_idx=material.tau_p_x_index(p);
									Omega_p_x_idx=material.Omega_p_x_index(p);
								}
								else {
                  //zero out remaining poles
									omega_p_x_idx=vacuum;
									tau_p_x_idx=vacuum;
									Omega_p_x_idx=vacuum;
								}
							}
							//Assign alpha_p_x, xi_p_x, gamma_p_x
							alpha_X(i,j,k,p) = (2-pow2(omega_p_x(omega_p_x_idx)*dt))/
                                               (1+1/tau_p_x(tau_p_x_idx)*dt);
                            xi_X(i,j,k,p) = (1/tau_p_x(tau_p_x_idx)*dt-1)/
                                            (1/tau_p_x(tau_p_x_idx)*dt+1);
                            gamma_X(i,j,k,p) = (epsilon_0*pow2(Omega_p_x(Omega_p_x_idx)*dt))/
                                               (1/tau_p_x(tau_p_x_idx)*dt+1)/(2*dt);
#if (0)
							Pa_X(i,j,k,p)=(tau_p_x(tau_p_x_idx)-dt/2.0)/
										(tau_p_x(tau_p_x_idx)+dt/2.0);
							Pb_X(i,j,k,p)=0.5*dx*(1+Pa_X(i,j,k,p))*
									 (tau_p_x(tau_p_x_idx)*pow2(Omega_p_x(Omega_p_x_idx))*epsilon_0*dt/2.0)/
									 (tau_p_x(tau_p_x_idx)+dt/2.0);
#endif
						}
						// Compute gamma_p sum
						// (...)
						// Assign Ca_X, Cb_X, Cc_X
						float gamma_p_sum(0.0);
						for (int p(0); p<max_num_lrntz_poles; ++p) {
						     gamma_p_sum += 0.5*gamma_X(i,j,k,p)*(2*dt);
            }
            Ca_X(i,j,k) = (2*epsilon_0*eps_x(eps_x_idx)-cond_e_x(cond_e_x_idx)*dt)/
                (2*epsilon_0*eps_x(eps_x_idx)+gamma_p_sum+cond_e_x(cond_e_x_idx)*dt);
            Cb_X(i,j,k) = 2*dt/dx/
                (2*epsilon_0*eps_x(eps_x_idx)+gamma_p_sum+cond_e_x(cond_e_x_idx)*dt);
            Cc_X(i,j,k) = gamma_p_sum/
                (2*epsilon_0*eps_x(eps_x_idx)+gamma_p_sum+cond_e_x(cond_e_x_idx)*dt);

#if (0)
						float beta_p_sum(0.0);
						for (int p(0); p<max_num_lrntz_poles; ++p) {
							beta_p_sum += (Pa_X(i,j,k,p)==-1.0?
                             0.0:
                             2/(dx*(1+Pa_X(i,j,k,p)))*Pb_X(i,j,k,p));
						}
						Ca_X(i,j,k)=(1-dt*(cond_e_x(cond_e_x_idx)+beta_p_sum)/(2.0*eps_x(eps_x_idx)*epsilon_0))
								 /(1+dt*(cond_e_x(cond_e_x_idx)+beta_p_sum)/(2.0*eps_x(eps_x_idx)*epsilon_0));
						Cb_X(i,j,k)=dt/eps_x(eps_x_idx)/epsilon_0/dx
								 /(1+dt*(cond_e_x(cond_e_x_idx)+beta_p_sum)/(2.0*eps_x(eps_x_idx)*epsilon_0));
#endif
          }
					/********************************************************************************************************/
				}
			}
		}
	}
	for (i=max(iback,BlockBack); i<=min(ifront+1,BlockFront+1); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright,BlockRight); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper+1,BlockUpper+1); k++)
			{
				if (shapeptr->IsInside(i-1,j-0.5,k-1))
				{
				  eps_y_type &eps_y_idx(eps_y_indices(i,j,k));
				  cond_e_y_type &cond_e_y_idx(cond_e_y_indices(i,j,k));
					if (material.eps_y_exists())
					{
						eps_y_idx=material.eps_y_index();
					}
					if (material.cond_e_y_exists())
					{
						cond_e_y_idx=material.cond_e_y_index();
					}
					/********************************************************************************************************/
          if (eff_num_lrntz_poles>0)
            dispersion_exists_at_Ey_position(i,j,k) = true;
          else
            dispersion_exists_at_Ey_position(i,j,k) = false;

					if (max_num_lrntz_poles==0) {
					  //there is NO dispersion
						Ca_Y(i,j,k)=(1-dt*cond_e_y(cond_e_y_idx)/(2.0*eps_y(eps_y_idx)*epsilon_0))/
						            (1+dt*cond_e_y(cond_e_y_idx)/(2.0*eps_y(eps_y_idx)*epsilon_0));
						Cb_Y(i,j,k)=dt/eps_y(eps_y_idx)/epsilon_0/dx/(1+dt*cond_e_y(cond_e_y_idx)/(2.0*eps_y(eps_y_idx)*epsilon_0));
					}
					else {
						for (int p(0); p<max_num_lrntz_poles; ++p) {
							omega_p_y_type &omega_p_y_idx(omega_p_y_indices(i,j,k,p)); //<- unused here
							tau_p_y_type &tau_p_y_idx(tau_p_y_indices(i,j,k,p));
							Omega_p_y_type &Omega_p_y_idx(Omega_p_y_indices(i,j,k,p));
							if (tot_num_lrntz_poles!=0) {//don't change anything if there are no poles (transparent material)
								if (p<tot_num_lrntz_poles) {
								  //place poles present in the material
									omega_p_y_idx=material.omega_p_y_index(p);
									tau_p_y_idx=material.tau_p_y_index(p);
									Omega_p_y_idx=material.Omega_p_y_index(p);
								}
								else {
                  //zero out remaining poles
									omega_p_y_idx=vacuum;
									tau_p_y_idx=vacuum;
									Omega_p_y_idx=vacuum;
								}
							}
              //Assign alpha_p_y, xi_p_y, gamma_p_y
							alpha_Y(i,j,k,p) = (2-pow2(omega_p_y(omega_p_y_idx)*dt))/
                                               (1+1/tau_p_y(tau_p_y_idx)*dt);
                            xi_Y(i,j,k,p) = (1/tau_p_y(tau_p_y_idx)*dt-1)/
                                            (1/tau_p_y(tau_p_y_idx)*dt+1);
                            gamma_Y(i,j,k,p) = (epsilon_0*pow2(Omega_p_y(Omega_p_y_idx)*dt))/
                                               (1/tau_p_y(tau_p_y_idx)*dt+1)/(2*dt);
#if (0)
							Pa_Y(i,j,k,p)=(tau_p_y(tau_p_y_idx)-dt/2.0)/
										(tau_p_y(tau_p_y_idx)+dt/2.0);
							Pb_Y(i,j,k,p)=0.5*dx*(1+Pa_Y(i,j,k,p))*
									 (tau_p_y(tau_p_y_idx)*pow2(Omega_p_y(Omega_p_y_idx))*epsilon_0*dt/2.0)/
									 (tau_p_y(tau_p_y_idx)+dt/2.0);
# endif
						}
            // Compute gamma_p sum
						// (...)
						// Assign Ca_Y, Cb_Y, Cc_Y
						float gamma_p_sum(0.0);
						for (int p(0); p<max_num_lrntz_poles; ++p) {
						     gamma_p_sum += 0.5*gamma_Y(i,j,k,p)*(2*dt);
            }
            Ca_Y(i,j,k) = (2*epsilon_0*eps_y(eps_y_idx)-cond_e_y(cond_e_y_idx)*dt)/
                (2*epsilon_0*eps_y(eps_y_idx)+gamma_p_sum+cond_e_y(cond_e_y_idx)*dt);
            Cb_Y(i,j,k) = 2*dt/dx/
                (2*epsilon_0*eps_y(eps_y_idx)+gamma_p_sum+cond_e_y(cond_e_y_idx)*dt);
            Cc_Y(i,j,k) = gamma_p_sum/
                (2*epsilon_0*eps_y(eps_y_idx)+gamma_p_sum+cond_e_y(cond_e_y_idx)*dt);
#if(0)
						float beta_p_sum(0.0);
						for (int p(0); p<max_num_lrntz_poles; ++p) {
							beta_p_sum += (Pa_Y(i,j,k,p)==-1.0?
                             0.0:
                             2/(dx*(1+Pa_Y(i,j,k,p)))*Pb_Y(i,j,k,p));
						}
						Ca_Y(i,j,k)=(1-dt*(cond_e_y(cond_e_y_idx)+beta_p_sum)/(2.0*eps_y(eps_y_idx)*epsilon_0))
								 /(1+dt*(cond_e_y(cond_e_y_idx)+beta_p_sum)/(2.0*eps_y(eps_y_idx)*epsilon_0));
						Cb_Y(i,j,k)=dt/eps_y(eps_y_idx)/epsilon_0/dx
								 /(1+dt*(cond_e_y(cond_e_y_idx)+beta_p_sum)/(2.0*eps_y(eps_y_idx)*epsilon_0));
#endif
          }
					/********************************************************************************************************/
				}
			}
		}
	}
	for (i=max(iback,BlockBack); i<=min(ifront+1,BlockFront+1); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright+1,BlockRight+1); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper,BlockUpper); k++)
			{
				if (shapeptr->IsInside(i-1,j-1,k-0.5))
				{
				  eps_z_type &eps_z_idx(eps_z_indices(i,j,k));
				  cond_e_z_type &cond_e_z_idx(cond_e_z_indices(i,j,k));
					if (material.eps_z_exists())
					{
						eps_z_idx=material.eps_z_index();
					}
					if (material.cond_e_z_exists())
					{
						cond_e_z_idx=material.cond_e_z_index();
					}
					/********************************************************************************************************/
          if (eff_num_lrntz_poles>0)
            dispersion_exists_at_Ez_position(i,j,k) = true;
          else
            dispersion_exists_at_Ez_position(i,j,k) = false;

					if (max_num_lrntz_poles==0) {
					  //there is NO dispersion
						Ca_Z(i,j,k)=(1-dt*cond_e_z(cond_e_z_idx)/(2.0*eps_z(eps_z_idx)*epsilon_0))/
						            (1+dt*cond_e_z(cond_e_z_idx)/(2.0*eps_z(eps_z_idx)*epsilon_0));
						Cb_Z(i,j,k)=dt/eps_z(eps_z_idx)/epsilon_0/dx/(1+dt*cond_e_z(cond_e_z_idx)/(2.0*eps_z(eps_z_idx)*epsilon_0));
					}
					else {
						for (int p(0); p<max_num_lrntz_poles; ++p) {
							omega_p_z_type &omega_p_z_idx(omega_p_z_indices(i,j,k,p)); //<- unused here
							tau_p_z_type &tau_p_z_idx(tau_p_z_indices(i,j,k,p));
							Omega_p_z_type &Omega_p_z_idx(Omega_p_z_indices(i,j,k,p));
							if (tot_num_lrntz_poles!=0) {//don't change anything if there are no poles (transparent material)
								if (p<tot_num_lrntz_poles) {
								  //place poles present in the material
									omega_p_z_idx=material.omega_p_z_index(p);
									tau_p_z_idx=material.tau_p_z_index(p);
									Omega_p_z_idx=material.Omega_p_z_index(p);
								}
								else {
                  //zero out remaining poles
									omega_p_z_idx=vacuum;
									tau_p_z_idx=vacuum;
									Omega_p_z_idx=vacuum;
								}
							}
              //Assign alpha_p_z, xi_p_z, gamma_p_z
							alpha_Z(i,j,k,p) = (2-pow2(omega_p_z(omega_p_z_idx)*dt))/
                                               (1+1/tau_p_z(tau_p_z_idx)*dt);
                            xi_Z(i,j,k,p) = (1/tau_p_z(tau_p_z_idx)*dt-1)/
                                            (1/tau_p_z(tau_p_z_idx)*dt+1);
                            gamma_Z(i,j,k,p) = (epsilon_0*pow2(Omega_p_z(Omega_p_z_idx)*dt))/
                                               (1/tau_p_z(tau_p_z_idx)*dt+1)/(2*dt);
#if(0)
							Pa_Z(i,j,k,p)=(tau_p_z(tau_p_z_idx)-dt/2.0)/
										(tau_p_z(tau_p_z_idx)+dt/2.0);
							Pb_Z(i,j,k,p)=0.5*dx*(1+Pa_Z(i,j,k,p))*
									 (tau_p_z(tau_p_z_idx)*pow2(Omega_p_z(Omega_p_z_idx))*epsilon_0*dt/2.0)/
									 (tau_p_z(tau_p_z_idx)+dt/2.0);
#endif
						}
            // Compute gamma_p sum
						// (...)
						// Assign Ca_Z, Cb_Z, Cc_Z
						float gamma_p_sum(0.0);
						for (int p(0); p<max_num_lrntz_poles; ++p) {
						     gamma_p_sum += 0.5*gamma_Z(i,j,k,p)*(2*dt);
            }
            Ca_Z(i,j,k) = (2*epsilon_0*eps_z(eps_z_idx)-cond_e_z(cond_e_z_idx)*dt)/
                (2*epsilon_0*eps_z(eps_z_idx)+gamma_p_sum+cond_e_z(cond_e_z_idx)*dt);
            Cb_Z(i,j,k) = 2*dt/dx/
                (2*epsilon_0*eps_z(eps_z_idx)+gamma_p_sum+cond_e_z(cond_e_z_idx)*dt);
            Cc_Z(i,j,k) = gamma_p_sum/
                (2*epsilon_0*eps_z(eps_z_idx)+gamma_p_sum+cond_e_z(cond_e_z_idx)*dt);
#if(0)
						float beta_p_sum(0.0);
						for (int p(0); p<max_num_lrntz_poles; ++p) {
							beta_p_sum += (Pa_Z(i,j,k,p)==-1.0?
                             0.0:
                             2/(dx*(1+Pa_Z(i,j,k,p)))*Pb_Z(i,j,k,p));
						}
						Ca_Z(i,j,k)=(1-dt*(cond_e_z(cond_e_z_idx)+beta_p_sum)/(2.0*eps_z(eps_z_idx)*epsilon_0))
								 /(1+dt*(cond_e_z(cond_e_z_idx)+beta_p_sum)/(2.0*eps_z(eps_z_idx)*epsilon_0));
						Cb_Z(i,j,k)=dt/eps_z(eps_z_idx)/epsilon_0/dx
								 /(1+dt*(cond_e_z(cond_e_z_idx)+beta_p_sum)/(2.0*eps_z(eps_z_idx)*epsilon_0));
#endif
          }
					/********************************************************************************************************/
				}
			}
		}
	}
	for (i=max(iback,BlockBack); i<=min(ifront+1,BlockFront+1); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright,BlockRight); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper,BlockUpper); k++)
			{
				if (shapeptr->IsInside(i-1,j-0.5,k-0.5))
				{
					if (material.mu_x_exists())
					{
						mu_x_indices(i,j,k)=material.mu_x_index();
					}
					if (material.cond_h_x_exists())
					{
						cond_h_x_indices(i,j,k)=material.cond_h_x_index();
					}
					Da_X(i,j,k)=(1-dt*cond_h_x(cond_h_x_indices(i,j,k))/(2.0*mu_x(mu_x_indices(i,j,k))*mu_0))/(1+dt*cond_h_x(cond_h_x_indices(i,j,k))/(2.0*mu_x(mu_x_indices(i,j,k))*mu_0));
					Db_X(i,j,k)=dt/mu_x(mu_x_indices(i,j,k))/mu_0/dx/(1+dt*cond_h_x(cond_h_x_indices(i,j,k))/(2.0*mu_x(mu_x_indices(i,j,k))*mu_0));
				}
			}
		}
	}
	for (i=max(iback,BlockBack); i<=min(ifront,BlockFront); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright+1,BlockRight+1); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper,BlockUpper); k++)
			{
				if (shapeptr->IsInside(i-0.5,j-1,k-0.5))
				{
					if (material.mu_y_exists())
					{
						mu_y_indices(i,j,k)=material.mu_y_index();
					}
					if (material.cond_h_y_exists())
					{
						cond_h_y_indices(i,j,k)=material.cond_h_y_index();
					}
					Da_Y(i,j,k)=(1-dt*cond_h_y(cond_h_y_indices(i,j,k))/(2.0*mu_y(mu_y_indices(i,j,k))*mu_0))/(1+dt*cond_h_y(cond_h_y_indices(i,j,k))/(2.0*mu_y(mu_y_indices(i,j,k))*mu_0));
					Db_Y(i,j,k)=dt/mu_y(mu_y_indices(i,j,k))/mu_0/dx/(1+dt*cond_h_y(cond_h_y_indices(i,j,k))/(2.0*mu_y(mu_y_indices(i,j,k))*mu_0));
				}
			}
		}
	}
	for (i=max(iback,BlockBack); i<=min(ifront,BlockFront); i++)
	{
		for (j=max(jleft,BlockLeft); j<=min(jright,BlockRight); j++)
		{
			for (k=max(klower,BlockLower); k<=min(kupper+1,BlockUpper+1); k++)
			{
				if (shapeptr->IsInside(i-0.5,j-0.5,k-1))
				{
					if (material.mu_z_exists())
					{
						mu_z_indices(i,j,k)=material.mu_z_index();
					}
					if (material.cond_h_z_exists())
					{
						cond_h_z_indices(i,j,k)=material.cond_h_z_index();
					}
					Da_Z(i,j,k)=(1-dt*cond_h_z(cond_h_z_indices(i,j,k))/(2.0*mu_z(mu_z_indices(i,j,k))*mu_0))/(1+dt*cond_h_z(cond_h_z_indices(i,j,k))/(2.0*mu_z(mu_z_indices(i,j,k))*mu_0));
					Db_Z(i,j,k)=dt/mu_z(mu_z_indices(i,j,k))/mu_0/dx/(1+dt*cond_h_z(cond_h_z_indices(i,j,k))/(2.0*mu_z(mu_z_indices(i,j,k))*mu_0));
				}
			}
		}
	}

	/***************************************/
	/** TODO: Interpolation at boundaries **/
	/***************************************/
}
