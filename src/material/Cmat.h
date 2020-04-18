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

//Declaration of the Cmat object that represents a material

#ifndef CMAT_H
#define CMAT_H

//for the vector STL class
#include <vector>

#include "Cmat_excp.h"

#include "Cmat_types.h"

//for the shared_ptr smart pointer (from the Boost library)
#include <boost/shared_ptr.hpp>

typedef boost::shared_ptr<eps_x_type> eps_x_type_ptr;
typedef boost::shared_ptr<eps_y_type> eps_y_type_ptr;
typedef boost::shared_ptr<eps_z_type> eps_z_type_ptr;
typedef boost::shared_ptr<mu_x_type> mu_x_type_ptr;
typedef boost::shared_ptr<mu_y_type> mu_y_type_ptr;
typedef boost::shared_ptr<mu_z_type> mu_z_type_ptr;
typedef boost::shared_ptr<cond_e_x_type> cond_e_x_type_ptr;
typedef boost::shared_ptr<cond_e_y_type> cond_e_y_type_ptr;
typedef boost::shared_ptr<cond_e_z_type> cond_e_z_type_ptr;
typedef boost::shared_ptr<cond_h_x_type> cond_h_x_type_ptr;
typedef boost::shared_ptr<cond_h_y_type> cond_h_y_type_ptr;
typedef boost::shared_ptr<cond_h_z_type> cond_h_z_type_ptr;
typedef vector<omega_p_x_type> omega_p_x_vec;
typedef vector<omega_p_y_type> omega_p_y_vec;
typedef vector<omega_p_z_type> omega_p_z_vec;
typedef vector<tau_p_x_type> tau_p_x_vec;
typedef vector<tau_p_y_type> tau_p_y_vec;
typedef vector<tau_p_z_type> tau_p_z_vec;
typedef vector<Omega_p_x_type> Omega_p_x_vec;
typedef vector<Omega_p_y_type> Omega_p_y_vec;
typedef vector<Omega_p_z_type> Omega_p_z_vec;


class Cmat
{
public:
	Cmat();

	void set_eps(const float& epsilon_r);
	void set_eps_x(const float& epsilon_r_x);
	void set_eps_y(const float& epsilon_r_y);
	void set_eps_z(const float& epsilon_r_z);
	void set_mu(const float& mu_r);
	void set_mu_x(const float& mu_r_x);
	void set_mu_y(const float& mu_r_y);
	void set_mu_z(const float& mu_r_z);
	void set_cond_e(const float& cond_e);
	void set_cond_e_x(const float& cond_e_x);
	void set_cond_e_y(const float& cond_e_y);
	void set_cond_e_z(const float& cond_e_z);
	void set_cond_h(const float& cond_h);
	void set_cond_h_x(const float& cond_h_x);
	void set_cond_h_y(const float& cond_h_y);
	void set_cond_h_z(const float& cond_h_z);
	void set_Omega_p(const float& Omega);
    void add_drude_pole(const float& omega_p,const float& tau_p);
    void add_lorentz_pole(const float& omega_p,
                      const float& delta_p,
                      const float& deps_p);

	bool eps_x_exists() const;
	bool eps_y_exists() const;
	bool eps_z_exists() const;
	bool mu_x_exists() const;
	bool mu_y_exists() const;
	bool mu_z_exists() const;
	bool cond_e_x_exists() const;
	bool cond_e_y_exists() const;
	bool cond_e_z_exists() const;
	bool cond_h_x_exists() const;
	bool cond_h_y_exists() const;
	bool cond_h_z_exists() const;

	eps_x_type eps_x_index() const;
	eps_y_type eps_y_index() const;
	eps_z_type eps_z_index() const;
	mu_x_type mu_x_index() const;
	mu_y_type mu_y_index() const;
	mu_z_type mu_z_index() const;
	cond_e_x_type cond_e_x_index() const;
	cond_e_y_type cond_e_y_index() const;
	cond_e_z_type cond_e_z_index() const;
	cond_h_x_type cond_h_x_index() const;
	cond_h_y_type cond_h_y_index() const;
	cond_h_z_type cond_h_z_index() const;
	omega_p_x_type omega_p_x_index(const int& p) const;
	omega_p_y_type omega_p_y_index(const int& p) const;
	omega_p_z_type omega_p_z_index(const int& p) const;
	tau_p_x_type tau_p_x_index(const int& p) const;
	tau_p_y_type tau_p_y_index(const int& p) const;
	tau_p_z_type tau_p_z_index(const int& p) const;
	Omega_p_x_type Omega_p_x_index(const int& p) const;
	Omega_p_y_type Omega_p_y_index(const int& p) const;
	Omega_p_z_type Omega_p_z_index(const int& p) const;

	float eps_x_value() const;
	float eps_y_value() const;
	float eps_z_value() const;
	float mu_x_value() const;
	float mu_y_value() const;
	float mu_z_value() const;
	float cond_e_x_value() const;
	float cond_e_y_value() const;
	float cond_e_z_value() const;
	float cond_h_x_value() const;
	float cond_h_y_value() const;
	float cond_h_z_value() const;
	float omega_p_x_value(const int& p) const;
	float omega_p_y_value(const int& p) const;
	float omega_p_z_value(const int& p) const;
	float tau_p_x_value(const int& p) const;
	float tau_p_y_value(const int& p) const;
	float tau_p_z_value(const int& p) const;
	float Omega_p_x_value(const int& p) const;
	float Omega_p_y_value(const int& p) const;
	float Omega_p_z_value(const int& p) const;

	int tot_num_lrntz_poles() const {
	  return _num_lrntz_poles;
	}

	int eff_num_lrntz_poles() const;

private:
	eps_x_type_ptr _eps_x_ptr;
	eps_y_type_ptr _eps_y_ptr;
	eps_z_type_ptr _eps_z_ptr;
	mu_x_type_ptr _mu_x_ptr;
	mu_y_type_ptr _mu_y_ptr;
	mu_z_type_ptr _mu_z_ptr;
	cond_e_x_type_ptr _cond_e_x_ptr;
	cond_e_y_type_ptr _cond_e_y_ptr;
	cond_e_z_type_ptr _cond_e_z_ptr;
	cond_h_x_type_ptr _cond_h_x_ptr;
	cond_h_y_type_ptr _cond_h_y_ptr;
	cond_h_z_type_ptr _cond_h_z_ptr;

	int _num_lrntz_poles; //changed only by add_**() methods

	omega_p_x_vec _omega_p_x_vec;
	omega_p_y_vec _omega_p_y_vec;
	omega_p_z_vec _omega_p_z_vec;
	tau_p_x_vec _tau_p_x_vec;
	tau_p_y_vec _tau_p_y_vec;
	tau_p_z_vec _tau_p_z_vec;
	Omega_p_x_vec _Omega_p_x_vec;
	Omega_p_y_vec _Omega_p_y_vec;
	Omega_p_z_vec _Omega_p_z_vec;

};


#endif // CMAT_H
