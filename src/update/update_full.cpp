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

#include "material/Cmat_types.h"

/** REMOVE LATER **/
extern double dx;
/** REMOVE LATER **/

extern Array<double,3> Ex,Ey,Ez,Ex_m1,Ey_m1,Ez_m1;
extern Array<double,3> Hx,Hy,Hz;

extern Array<update_coeff_type,3> Ca_X,Cb_X,Cc_X,Ca_Y,Cb_Y,Cc_Y,Ca_Z,Cb_Z,Cc_Z;
extern Array<update_coeff_type,3> Da_X,Db_X,Da_Y,Db_Y,Da_Z,Db_Z;

extern Array<double,1> inv_kappa_e_x,inv_kappa_e_y,inv_kappa_e_z,inv_kappa_h_x,inv_kappa_h_y,inv_kappa_h_z;

extern Array<bool,3> dispersion_exists_at_Ex_position,dispersion_exists_at_Ey_position,dispersion_exists_at_Ez_position;
extern Array<double,4> J_p_x,J_p_y,J_p_z,Jm1_p_x,Jm1_p_y,Jm1_p_z;
extern Array<update_coeff_type,4> alpha_X,xi_X,gamma_X,alpha_Y,xi_Y,gamma_Y,alpha_Z,xi_Z,gamma_Z;
extern int pole_dim_size;

extern int Ex_min_index_in_x,Ex_max_index_in_x,Ex_min_index_in_y,Ex_max_index_in_y,Ex_min_index_in_z,Ex_max_index_in_z;
extern int Ey_min_index_in_x,Ey_max_index_in_x,Ey_min_index_in_y,Ey_max_index_in_y,Ey_min_index_in_z,Ey_max_index_in_z;
extern int Ez_min_index_in_x,Ez_max_index_in_x,Ez_min_index_in_y,Ez_max_index_in_y,Ez_min_index_in_z,Ez_max_index_in_z;
extern int Hx_min_index_in_x,Hx_max_index_in_x,Hx_min_index_in_y,Hx_max_index_in_y,Hx_min_index_in_z,Hx_max_index_in_z;
extern int Hy_min_index_in_x,Hy_max_index_in_x,Hy_min_index_in_y,Hy_max_index_in_y,Hy_min_index_in_z,Hy_max_index_in_z;
extern int Hz_min_index_in_x,Hz_max_index_in_x,Hz_min_index_in_y,Hz_max_index_in_y,Hz_min_index_in_z,Hz_max_index_in_z;

namespace{
	int i,j,k;
//	bool dispersion_exists_here;
	double E_old;
	double *E_new,*H_new,*Em1;
	double jp_old;
	double *jp,*jm1;
	double jsum;
	update_coeff_type *alpha_p_x,*xi_p_x,*gamma_p_x,*alpha_p_y,*xi_p_y,*gamma_p_y,*alpha_p_z,*xi_p_z,*gamma_p_z;
};


//updates will all material properties taken into account
void updateE_full(const int& n)
{
	//Update Ex
	for (i=Ex_min_index_in_x; i<=Ex_max_index_in_x; i++){
		for (j=Ex_min_index_in_y; j<=Ex_max_index_in_y; j++){
			for (k=Ex_min_index_in_z; k<=Ex_max_index_in_z; k++)
			{
				if (dispersion_exists_at_Ex_position(i,j,k))
				{
					/** 1. COMPUTE THE Jx POLARIZATION CONTRIB **/
					jsum=0;
                    jp = &(J_p_x(i,j,k,0));
                    jm1 = &(Jm1_p_x(i,j,k,0));
                    alpha_p_x = &(alpha_X(i,j,k,0));
                    xi_p_x = &(xi_X(i,j,k,0));
                    for (int p(0); p<pole_dim_size;
                        ++jp, ++jm1, ++alpha_p_x, ++xi_p_x, ++p) {
                            jsum += 0.5*dx*((1+(*alpha_p_x))*(*jp)+(*xi_p_x)*(*jm1));
                    }
					/** 2. DO ELECTRIC FIELD UPDATE **/
					E_new=&(Ex(i,j,k));
					E_old=(*E_new);
                    Em1 = &(Ex_m1(i,j,k));
					(*E_new)=Ca_X(i,j,k)*E_old
						+ Cb_X(i,j,k)*
						((Hz(i,j,k)-Hz(i,j-1,k))*inv_kappa_e_y(j)
						+(Hy(i,j,k-1)-Hy(i,j,k))*inv_kappa_e_z(k)
					/** 3. ADD THE Jx POLARIZATION CONTRIB **/
					    -jsum)
					/** 3.b. ADD THE Ex^(n-1) CONTRIB **/
					    + Cc_X(i,j,k)*(*Em1);
					/** 4. DO THE Jx UPDATE USING E_old and E_new **/
                    jp = &(J_p_x(i,j,k,0));
                    // Assign reference
                    jm1 = &(Jm1_p_x(i,j,k,0));
                    alpha_p_x = &(alpha_X(i,j,k,0));
                    xi_p_x = &(xi_X(i,j,k,0));
                    gamma_p_x = &(gamma_X(i,j,k,0));
                    for (int p(0); p<pole_dim_size;
                        ++jp, ++jm1, ++alpha_p_x, ++xi_p_x, ++gamma_p_x, ++p) {
                        jp_old = (*jp);
                        (*jp)=(*alpha_p_x)*jp_old+(*xi_p_x)*(*jm1)+(*gamma_p_x)*((*E_new)-(*Em1));
                        *jm1 = jp_old;
                    }
                    *Em1 = E_old;
                }
				else
				{
					E_new=&(Ex(i,j,k));
					(*E_new)=Ca_X(i,j,k)*(*E_new)
						+ Cb_X(i,j,k)*
						((Hz(i,j,k)-Hz(i,j-1,k))*inv_kappa_e_y(j)
						+(Hy(i,j,k-1)-Hy(i,j,k))*inv_kappa_e_z(k));
				}
			}
		}
	}
	//Update Ey
	for (i=Ey_min_index_in_x; i<=Ey_max_index_in_x; i++){
		for (j=Ey_min_index_in_y; j<=Ey_max_index_in_y; j++){
			for (k=Ey_min_index_in_z; k<=Ey_max_index_in_z; k++)
			{
				if (dispersion_exists_at_Ey_position(i,j,k))
				{
					/** 1. COMPUTE THE Jx POLARIZATION CONTRIB **/
					jsum=0;
                    jp = &(J_p_y(i,j,k,0));
                    jm1 = &(Jm1_p_y(i,j,k,0));
                    alpha_p_y = &(alpha_Y(i,j,k,0));
                    xi_p_y = &(xi_Y(i,j,k,0));
                    for (int p(0); p<pole_dim_size;
                       ++jp, ++jm1, ++alpha_p_y, ++xi_p_y, ++p) {
                    jsum += 0.5*dx*((1+(*alpha_p_y))*(*jp)+(*xi_p_y)*(*jm1));
                    }
					/** 2. DO ELECTRIC FIELD UPDATE **/
					E_new=&(Ey(i,j,k));
					E_old=(*E_new);
					Em1 = &(Ey_m1(i,j,k));
					(*E_new)=Ca_Y(i,j,k)*E_old
						+ Cb_Y(i,j,k)*
						((Hx(i,j,k)-Hx(i,j,k-1))*inv_kappa_e_z(k)
						+(Hz(i-1,j,k)-Hz(i,j,k))*inv_kappa_e_x(i)
					/** 3. ADD THE Jy POLARIZATION CONTRIB **/
					    -jsum)
          /** 3.b. ADD THE Ex^(n-1) CONTRIB **/
					    +Cc_Y(i,j,k)*(*Em1);
					/** 4. DO THE Jy UPDATE USING E_old and E_new **/
                    jp = &(J_p_y(i,j,k,0));
                    // Assign reference
                    jm1 = &(Jm1_p_y(i,j,k,0));
                    alpha_p_y = &(alpha_Y(i,j,k,0));
                    xi_p_y = &(xi_Y(i,j,k,0));
                    gamma_p_y = &(gamma_Y(i,j,k,0));
                    for (int p(0); p<pole_dim_size;
                        ++jp, ++jm1, ++alpha_p_y, ++xi_p_y, ++gamma_p_y, ++p) {
                    jp_old = (*jp);
                    (*jp)=(*alpha_p_y)*jp_old+(*xi_p_y)*(*jm1)+(*gamma_p_y)*((*E_new)-(*Em1));
                    *jm1 = jp_old;
                    }
                    *Em1 = E_old;
				}
				else
				{
					E_new=&(Ey(i,j,k));
					(*E_new)=Ca_Y(i,j,k)*(*E_new)
						+ Cb_Y(i,j,k)*
						((Hx(i,j,k)-Hx(i,j,k-1))*inv_kappa_e_z(k)
						+(Hz(i-1,j,k)-Hz(i,j,k))*inv_kappa_e_x(i));
				}
			}
		}
	}
	//Update Ez
	for (i=Ez_min_index_in_x; i<=Ez_max_index_in_x; i++){
		for (j=Ez_min_index_in_y; j<=Ez_max_index_in_y; j++){
			for (k=Ez_min_index_in_z; k<=Ez_max_index_in_z; k++)
			{
				if (dispersion_exists_at_Ez_position(i,j,k))
				{
					/** 1. COMPUTE THE Jx POLARIZATION CONTRIB **/
					jsum=0;
                    jp = &(J_p_z(i,j,k,0));
                    jm1 = &(Jm1_p_z(i,j,k,0));
                    alpha_p_z = &(alpha_Z(i,j,k,0));
                    xi_p_z = &(xi_Z(i,j,k,0));
                    for (int p(0); p<pole_dim_size;
                         ++jp, ++jm1, ++alpha_p_z, ++xi_p_z, ++p) {
                    jsum += 0.5*dx*((1+(*alpha_p_z))*(*jp)+(*xi_p_z)*(*jm1));
                    }
					/** 2. DO ELECTRIC FIELD UPDATE **/
					E_new=&(Ez(i,j,k));
					E_old=(*E_new);
					Em1 = &(Ez_m1(i,j,k));
					(*E_new)=Ca_Z(i,j,k)*E_old
						+ Cb_Z(i,j,k)*
						((Hy(i,j,k)-Hy(i-1,j,k))*inv_kappa_e_x(i)
							+(Hx(i,j-1,k)-Hx(i,j,k))*inv_kappa_e_y(j)
					/** 3. ADD THE Jz POLARIZATION CONTRIB **/
					    -jsum)
                    /** 3.b. ADD THE Ex^(n-1) CONTRIB **/
					    + Cc_Y(i,j,k)*(*Em1);
                            /** 4. DO THE Jz UPDATE USING E_old and E_new **/
                    jp = &(J_p_z(i,j,k,0));
                    // Assign reference
                    jm1 = &(Jm1_p_z(i,j,k,0));
                    alpha_p_z = &(alpha_Z(i,j,k,0));
                    xi_p_z = &(xi_Z(i,j,k,0));
                    gamma_p_z = &(gamma_Z(i,j,k,0));
                    for (int p(0); p<pole_dim_size;
                         ++jp, ++jm1, ++alpha_p_z, ++xi_p_z, ++gamma_p_z, ++p) {
                    jp_old = (*jp);
                    (*jp)=(*alpha_p_z)*jp_old+(*xi_p_z)*(*jm1)+(*gamma_p_z)*((*E_new)-(*Em1));
                    *jm1 = jp_old;
                    }
                    *Em1 = E_old;
				}
				else
				{
					E_new=&(Ez(i,j,k));
					(*E_new)=Ca_Z(i,j,k)*(*E_new)
						+ Cb_Z(i,j,k)*
						((Hy(i,j,k)-Hy(i-1,j,k))*inv_kappa_e_x(i)
							+(Hx(i,j-1,k)-Hx(i,j,k))*inv_kappa_e_y(j));
				}
			}
		}
	}
}

void updateH_full(const int& n)
{
	//Update Hx
	for (i=Hx_min_index_in_x; i<=Hx_max_index_in_x; i++){
		for (j=Hx_min_index_in_y; j<=Hx_max_index_in_y; j++){
			for (k=Hx_min_index_in_z; k<=Hx_max_index_in_z; k++)
			{
				H_new=&(Hx(i,j,k));
				(*H_new)=Da_X(i,j,k)*(*H_new)
					+ Db_X(i,j,k)*
					((Ey(i,j,k+1)-Ey(i,j,k))*inv_kappa_h_z(k)
						+(Ez(i,j,k)-Ez(i,j+1,k))*inv_kappa_h_y(j));
			}
		}
	}

	//Update Hy
	for (i=Hy_min_index_in_x; i<=Hy_max_index_in_x; i++){
		for (j=Hy_min_index_in_y; j<=Hy_max_index_in_y; j++){
			for (k=Hy_min_index_in_z; k<=Hy_max_index_in_z; k++)
			{
				H_new=&(Hy(i,j,k));
				(*H_new)=Da_Y(i,j,k)*(*H_new)
					+ Db_Y(i,j,k)*
					((Ez(i+1,j,k)-Ez(i,j,k))*inv_kappa_h_x(i)
						+(Ex(i,j,k)-Ex(i,j,k+1))*inv_kappa_h_z(k));
			}
		}
	}

    //Update Hz
	for (i=Hz_min_index_in_x; i<=Hz_max_index_in_x; i++){
		for (j=Hz_min_index_in_y; j<=Hz_max_index_in_y; j++){
			for (k=Hz_min_index_in_z; k<=Hz_max_index_in_z; k++)
			{
				H_new=&(Hz(i,j,k));
				(*H_new)=Da_Z(i,j,k)*(*H_new)
					+ Db_Z(i,j,k)*
					((Ex(i,j+1,k)-Ex(i,j,k))*inv_kappa_h_y(j)
						+(Ey(i,j,k)-Ey(i+1,j,k))*inv_kappa_h_x(i));
			}
		}
	}
}
