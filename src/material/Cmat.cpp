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

//Definition of the Cmat object that represents a material

#include "Cmat_excp.h"

#include "Cmat.h"

extern const int vacuum;

extern Array<float,1> eps_x,eps_y,eps_z;
extern Array<float,1> mu_x,mu_y,mu_z;
extern Array<float,1> cond_e_x,cond_e_y,cond_e_z;
extern Array<float,1> cond_h_x,cond_h_y,cond_h_z;
extern Array<float,1> omega_p_x,omega_p_y,omega_p_z,tau_p_x,tau_p_y,tau_p_z,
                      Omega_p_x,Omega_p_y,Omega_p_z;


/** MOVE INTO SEPARATE HEADER **/
eps_x_type add_distinct_eps_x(const float& epsilon_r_x)
{
	eps_x_type num_of_distinct_eps_x = eps_x.size()+1;
	eps_x.resizeAndPreserve(num_of_distinct_eps_x);
	eps_x(num_of_distinct_eps_x-1) = epsilon_r_x;
	return num_of_distinct_eps_x-1;
}

eps_y_type add_distinct_eps_y(const float& epsilon_r_y)
{
	eps_y_type num_of_distinct_eps_y = eps_y.size()+1;
	eps_y.resizeAndPreserve(num_of_distinct_eps_y);
	eps_y(num_of_distinct_eps_y-1) = epsilon_r_y;
	return num_of_distinct_eps_y-1;
}

eps_z_type add_distinct_eps_z(const float& epsilon_r_z)
{
	eps_z_type num_of_distinct_eps_z = eps_z.size()+1;
	eps_z.resizeAndPreserve(num_of_distinct_eps_z);
	eps_z(num_of_distinct_eps_z-1) = epsilon_r_z;
	return num_of_distinct_eps_z-1;
}

mu_x_type add_distinct_mu_x(const float& mu_r_x)
{
	mu_x_type num_of_distinct_mu_x = mu_x.size()+1;
	mu_x.resizeAndPreserve(num_of_distinct_mu_x);
	mu_x(num_of_distinct_mu_x-1) = mu_r_x;
	return num_of_distinct_mu_x-1;
}

mu_y_type add_distinct_mu_y(const float& mu_r_y)
{
	mu_y_type num_of_distinct_mu_y = mu_y.size()+1;
	mu_y.resizeAndPreserve(num_of_distinct_mu_y);
	mu_y(num_of_distinct_mu_y-1) = mu_r_y;
	return num_of_distinct_mu_y-1;
}

mu_z_type add_distinct_mu_z(const float& mu_r_z)
{
	mu_z_type num_of_distinct_mu_z = mu_z.size()+1;
	mu_z.resizeAndPreserve(num_of_distinct_mu_z);
	mu_z(num_of_distinct_mu_z-1) = mu_r_z;
	return num_of_distinct_mu_z-1;
}

cond_e_x_type add_distinct_cond_e_x(const float& sigma_e_x)
{
	cond_e_x_type num_of_distinct_cond_e_x = cond_e_x.size()+1;
	cond_e_x.resizeAndPreserve(num_of_distinct_cond_e_x);
	cond_e_x(num_of_distinct_cond_e_x-1) = sigma_e_x;
	return num_of_distinct_cond_e_x-1;
}

cond_e_y_type add_distinct_cond_e_y(const float& sigma_e_y)
{
	cond_e_y_type num_of_distinct_cond_e_y = cond_e_y.size()+1;
	cond_e_y.resizeAndPreserve(num_of_distinct_cond_e_y);
	cond_e_y(num_of_distinct_cond_e_y-1) = sigma_e_y;
	return num_of_distinct_cond_e_y-1;
}

cond_e_z_type add_distinct_cond_e_z(const float& sigma_e_z)
{
	cond_e_z_type num_of_distinct_cond_e_z = cond_e_z.size()+1;
	cond_e_z.resizeAndPreserve(num_of_distinct_cond_e_z);
	cond_e_z(num_of_distinct_cond_e_z-1) = sigma_e_z;
	return num_of_distinct_cond_e_z-1;
}

cond_h_x_type add_distinct_cond_h_x(const float& sigma_h_x)
{
	cond_h_x_type num_of_distinct_cond_h_x = cond_h_x.size()+1;
	cond_h_x.resizeAndPreserve(num_of_distinct_cond_h_x);
	cond_h_x(num_of_distinct_cond_h_x-1) = sigma_h_x;
	return num_of_distinct_cond_h_x-1;
}

cond_h_y_type add_distinct_cond_h_y(const float& sigma_h_y)
{
	cond_h_y_type num_of_distinct_cond_h_y = cond_h_y.size()+1;
	cond_h_y.resizeAndPreserve(num_of_distinct_cond_h_y);
	cond_h_y(num_of_distinct_cond_h_y-1) = sigma_h_y;
	return num_of_distinct_cond_h_y-1;
}

cond_h_z_type add_distinct_cond_h_z(const float& sigma_h_z)
{
	cond_h_z_type num_of_distinct_cond_h_z = cond_h_z.size()+1;
	cond_h_z.resizeAndPreserve(num_of_distinct_cond_h_z);
	cond_h_z(num_of_distinct_cond_h_z-1) = sigma_h_z;
	return num_of_distinct_cond_h_z-1;
}

omega_p_x_type add_distinct_omega_p_x(const float& omega)
{
	omega_p_x_type num_of_distinct_omega_p_x = omega_p_x.size()+1;
	omega_p_x.resizeAndPreserve(num_of_distinct_omega_p_x);
	omega_p_x(num_of_distinct_omega_p_x-1) = omega;
	return num_of_distinct_omega_p_x-1;
}

omega_p_y_type add_distinct_omega_p_y(const float& omega)
{
	omega_p_y_type num_of_distinct_omega_p_y = omega_p_y.size()+1;
	omega_p_y.resizeAndPreserve(num_of_distinct_omega_p_y);
	omega_p_y(num_of_distinct_omega_p_y-1) = omega;
	return num_of_distinct_omega_p_y-1;
}

omega_p_z_type add_distinct_omega_p_z(const float& omega)
{
	omega_p_z_type num_of_distinct_omega_p_z = omega_p_z.size()+1;
	omega_p_z.resizeAndPreserve(num_of_distinct_omega_p_z);
	omega_p_z(num_of_distinct_omega_p_z-1) = omega;
	return num_of_distinct_omega_p_z-1;
}

tau_p_x_type add_distinct_tau_p_x(const float& tau)
{
	tau_p_x_type num_of_distinct_tau_p_x = tau_p_x.size()+1;
	tau_p_x.resizeAndPreserve(num_of_distinct_tau_p_x);
	tau_p_x(num_of_distinct_tau_p_x-1) = tau;
	return num_of_distinct_tau_p_x-1;
}

tau_p_y_type add_distinct_tau_p_y(const float& tau)
{
	tau_p_y_type num_of_distinct_tau_p_y = tau_p_y.size()+1;
	tau_p_y.resizeAndPreserve(num_of_distinct_tau_p_y);
	tau_p_y(num_of_distinct_tau_p_y-1) = tau;
	return num_of_distinct_tau_p_y-1;
}

tau_p_z_type add_distinct_tau_p_z(const float& tau)
{
	tau_p_z_type num_of_distinct_tau_p_z = tau_p_z.size()+1;
	tau_p_z.resizeAndPreserve(num_of_distinct_tau_p_z);
	tau_p_z(num_of_distinct_tau_p_z-1) = tau;
	return num_of_distinct_tau_p_z-1;
}

Omega_p_x_type add_distinct_Omega_p_x(const float& Omega)
{
	Omega_p_x_type num_of_distinct_Omega_p_x = Omega_p_x.size()+1;
	Omega_p_x.resizeAndPreserve(num_of_distinct_Omega_p_x);
	Omega_p_x(num_of_distinct_Omega_p_x-1) = Omega;
	return num_of_distinct_Omega_p_x-1;
}

Omega_p_y_type add_distinct_Omega_p_y(const float& Omega)
{
	Omega_p_y_type num_of_distinct_Omega_p_y = Omega_p_y.size()+1;
	Omega_p_y.resizeAndPreserve(num_of_distinct_Omega_p_y);
	Omega_p_y(num_of_distinct_Omega_p_y-1) = Omega;
	return num_of_distinct_Omega_p_y-1;
}

Omega_p_z_type add_distinct_Omega_p_z(const float& Omega)
{
	Omega_p_z_type num_of_distinct_Omega_p_z = Omega_p_z.size()+1;
	Omega_p_z.resizeAndPreserve(num_of_distinct_Omega_p_z);
	Omega_p_z(num_of_distinct_Omega_p_z-1) = Omega;
	return num_of_distinct_Omega_p_z-1;
}
/** MOVE INTO SEPARATE HEADER **/

Cmat::Cmat()
 : _num_lrntz_poles(0) {
}

void Cmat::set_eps(const float& epsilon_r)
{
	set_eps_x(epsilon_r);
	set_eps_y(epsilon_r);
	set_eps_z(epsilon_r);
}
void Cmat::set_eps_x(const float& epsilon_r)
{
	_eps_x_ptr.reset(new eps_x_type(add_distinct_eps_x(epsilon_r)));
}
void Cmat::set_eps_y(const float& epsilon_r)
{
	_eps_y_ptr.reset(new eps_y_type(add_distinct_eps_y(epsilon_r)));
}
void Cmat::set_eps_z(const float& epsilon_r)
{
	_eps_z_ptr.reset(new eps_z_type(add_distinct_eps_z(epsilon_r)));
}

void Cmat::set_mu(const float& mu_r)
{
	set_mu_x(mu_r);
	set_mu_y(mu_r);
	set_mu_z(mu_r);
}
void Cmat::set_mu_x(const float& mu_r)
{
	_mu_x_ptr.reset(new mu_x_type(add_distinct_mu_x(mu_r)));
}
void Cmat::set_mu_y(const float& mu_r)
{
	_mu_y_ptr.reset(new mu_y_type(add_distinct_mu_y(mu_r)));
}
void Cmat::set_mu_z(const float& mu_r)
{
	_mu_z_ptr.reset(new mu_z_type(add_distinct_mu_z(mu_r)));
}

void Cmat::set_cond_e(const float& sigma_e)
{
	set_cond_e_x(sigma_e);
	set_cond_e_y(sigma_e);
	set_cond_e_z(sigma_e);
}
void Cmat::set_cond_e_x(const float& sigma_e)
{
	_cond_e_x_ptr.reset(new cond_e_x_type(add_distinct_cond_e_x(sigma_e)));
}
void Cmat::set_cond_e_y(const float& sigma_e)
{
	_cond_e_y_ptr.reset(new cond_e_y_type(add_distinct_cond_e_y(sigma_e)));
}
void Cmat::set_cond_e_z(const float& sigma_e)
{
	_cond_e_z_ptr.reset(new cond_e_z_type(add_distinct_cond_e_z(sigma_e)));
}

void Cmat::set_cond_h(const float& sigma_h)
{
	set_cond_h_x(sigma_h);
	set_cond_h_y(sigma_h);
	set_cond_h_z(sigma_h);
}
void Cmat::set_cond_h_x(const float& sigma_h)
{
	_cond_h_x_ptr.reset(new cond_h_x_type(add_distinct_cond_h_x(sigma_h)));
}
void Cmat::set_cond_h_y(const float& sigma_h)
{
	_cond_h_y_ptr.reset(new cond_h_y_type(add_distinct_cond_h_y(sigma_h)));
}
void Cmat::set_cond_h_z(const float& sigma_h)
{
	_cond_h_z_ptr.reset(new cond_h_z_type(add_distinct_cond_h_z(sigma_h)));
}

/* In test: add Omega_p materials */
void Cmat::set_Omega_p(const float& Omega)
{
  _Omega_p_x_vec.push_back(add_distinct_Omega_p_x(Omega));
  _Omega_p_y_vec.push_back(add_distinct_Omega_p_y(Omega));
  _Omega_p_z_vec.push_back(add_distinct_Omega_p_z(Omega));
}

void Cmat::add_drude_pole(const float& omega_p,const float& tau_p) {
  //increment number of poles
  _num_lrntz_poles++;
  //assign pole parameters
  _omega_p_x_vec.push_back(vacuum); //<-not the same as Lorentz's omega!!
  _omega_p_y_vec.push_back(vacuum);
  _omega_p_z_vec.push_back(vacuum);
  _tau_p_x_vec.push_back(add_distinct_tau_p_x(2*tau_p));
  _tau_p_y_vec.push_back(add_distinct_tau_p_y(2*tau_p));
  _tau_p_z_vec.push_back(add_distinct_tau_p_z(2*tau_p));
  _Omega_p_x_vec.push_back(add_distinct_Omega_p_x(omega_p));
  _Omega_p_y_vec.push_back(add_distinct_Omega_p_y(omega_p));
  _Omega_p_z_vec.push_back(add_distinct_Omega_p_z(omega_p));
}

void Cmat::add_lorentz_pole(const float& omega_p,
                            const float& delta_p,
                            const float& deps_p) {
  //increment number of poles
  _num_lrntz_poles++;
  //assign pole parameters
  _omega_p_x_vec.push_back(add_distinct_omega_p_x(omega_p));
  _omega_p_y_vec.push_back(add_distinct_omega_p_y(omega_p));
  _omega_p_z_vec.push_back(add_distinct_omega_p_z(omega_p));
  _tau_p_x_vec.push_back(add_distinct_tau_p_x(1./delta_p));
  _tau_p_y_vec.push_back(add_distinct_tau_p_y(1./delta_p));
  _tau_p_z_vec.push_back(add_distinct_tau_p_z(1./delta_p));
  _Omega_p_x_vec.push_back(add_distinct_Omega_p_x(sqrt(deps_p)*omega_p));
  _Omega_p_y_vec.push_back(add_distinct_Omega_p_y(sqrt(deps_p)*omega_p));
  _Omega_p_z_vec.push_back(add_distinct_Omega_p_z(sqrt(deps_p)*omega_p));
}

bool Cmat::eps_x_exists() const
{
	return _eps_x_ptr;
}

bool Cmat::eps_y_exists() const
{
	return _eps_y_ptr;
}

bool Cmat::eps_z_exists() const
{
	return _eps_z_ptr;
}

bool Cmat::mu_x_exists() const
{
	return _mu_x_ptr;
}

bool Cmat::mu_y_exists() const
{
	return _mu_y_ptr;
}

bool Cmat::mu_z_exists() const
{
	return _mu_z_ptr;
}

bool Cmat::cond_e_x_exists() const
{
	return _cond_e_x_ptr;
}

bool Cmat::cond_e_y_exists() const
{
	return _cond_e_y_ptr;
}

bool Cmat::cond_e_z_exists() const
{
	return _cond_e_z_ptr;
}

bool Cmat::cond_h_x_exists() const
{
	return _cond_h_x_ptr;
}

bool Cmat::cond_h_y_exists() const
{
	return _cond_h_y_ptr;
}

bool Cmat::cond_h_z_exists() const
{
	return _cond_h_z_ptr;
}

eps_x_type Cmat::eps_x_index() const
{
	if (eps_x_exists())
	{
		return *_eps_x_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("eps_x");
	}
}

eps_y_type Cmat::eps_y_index() const
{
	if (eps_y_exists())
	{
		return *_eps_y_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("eps_y");
	}
}

eps_z_type Cmat::eps_z_index() const
{
	if (eps_z_exists())
	{
		return *_eps_z_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("eps_z");
	}
}

mu_x_type Cmat::mu_x_index() const
{
	if (mu_x_exists())
	{
		return *_mu_x_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("mu_x");
	}
}

mu_y_type Cmat::mu_y_index() const
{
	if (mu_y_exists())
	{
		return *_mu_y_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("mu_y");
	}
}

mu_z_type Cmat::mu_z_index() const
{
	if (mu_z_exists())
	{
		return *_mu_z_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("mu_z");
	}
}

cond_e_x_type Cmat::cond_e_x_index() const
{
	if (cond_e_x_exists())
	{
		return *_cond_e_x_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_e_x");
	}
}

cond_e_y_type Cmat::cond_e_y_index() const
{
	if (cond_e_y_exists())
	{
		return *_cond_e_y_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_e_y");
	}
}

cond_e_z_type Cmat::cond_e_z_index() const
{
	if (cond_e_z_exists())
	{
		return *_cond_e_z_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_e_z");
	}
}

cond_h_x_type Cmat::cond_h_x_index() const
{
	if (cond_h_x_exists())
	{
		return *_cond_h_x_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_h_x");
	}
}

cond_h_y_type Cmat::cond_h_y_index() const
{
	if (cond_h_y_exists())
	{
		return *_cond_h_y_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_h_y");
	}
}

cond_h_z_type Cmat::cond_h_z_index() const
{
	if (cond_h_z_exists())
	{
		return *_cond_h_z_ptr;
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_h_z");
	}
}

omega_p_x_type Cmat::omega_p_x_index(const int& p) const
{
  return _omega_p_x_vec[p];
  //may throw exception if pole doesn't exist
}

omega_p_y_type Cmat::omega_p_y_index(const int& p) const
{
  return _omega_p_y_vec[p];
  //may throw exception if pole doesn't exist
}

omega_p_z_type Cmat::omega_p_z_index(const int& p) const
{
  return _omega_p_z_vec[p];
  //may throw exception if pole doesn't exist
}

tau_p_x_type Cmat::tau_p_x_index(const int& p) const
{
  return _tau_p_x_vec[p];
  //may throw exception if pole doesn't exist
}

tau_p_y_type Cmat::tau_p_y_index(const int& p) const
{
  return _tau_p_y_vec[p];
  //may throw exception if pole doesn't exist
}

tau_p_z_type Cmat::tau_p_z_index(const int& p) const
{
  return _tau_p_z_vec[p];
  //may throw exception if pole doesn't exist
}

Omega_p_x_type Cmat::Omega_p_x_index(const int& p) const
{
  return _Omega_p_x_vec[p];
  //may throw exception if pole doesn't exist
}

Omega_p_y_type Cmat::Omega_p_y_index(const int& p) const
{
  return _Omega_p_y_vec[p];
  //may throw exception if pole doesn't exist
}

Omega_p_z_type Cmat::Omega_p_z_index(const int& p) const
{
  return _Omega_p_z_vec[p];
  //may throw exception if pole doesn't exist
}

float Cmat::eps_x_value() const
{
	if (eps_x_exists())
	{
		return eps_x(*_eps_x_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("eps_x");
	}
}

float Cmat::eps_y_value() const
{
	if (eps_y_exists())
	{
		return eps_y(*_eps_y_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("eps_y");
	}
}

float Cmat::eps_z_value() const
{
	if (eps_z_exists())
	{
		return eps_z(*_eps_z_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("eps_z");
	}
}

float Cmat::mu_x_value() const
{
	if (mu_x_exists())
	{
		return mu_x(*_mu_x_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("mu_x");
	}
}

float Cmat::mu_y_value() const
{
	if (mu_y_exists())
	{
		return mu_y(*_mu_y_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("mu_y");
	}
}

float Cmat::mu_z_value() const
{
	if (mu_z_exists())
	{
		return mu_z(*_mu_z_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("mu_z");
	}
}

float Cmat::cond_e_x_value() const
{
	if (cond_e_x_exists())
	{
		return cond_e_x(*_cond_e_x_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_e_x");
	}
}

float Cmat::cond_e_y_value() const
{
	if (cond_e_y_exists())
	{
		return cond_e_y(*_cond_e_y_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_e_y");
	}
}

float Cmat::cond_e_z_value() const
{
	if (cond_e_z_exists())
	{
		return cond_e_z(*_cond_e_z_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_e_z");
	}
}

float Cmat::cond_h_x_value() const
{
	if (cond_h_x_exists())
	{
		return cond_h_x(*_cond_h_x_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_h_x");
	}
}

float Cmat::cond_h_y_value() const
{
	if (cond_h_y_exists())
	{
		return cond_h_y(*_cond_h_y_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_h_y");
	}
}

float Cmat::cond_h_z_value() const
{
	if (cond_h_z_exists())
	{
		return cond_h_z(*_cond_h_z_ptr);
	}
	else
	{
		throw MaterialPropertyDoesNotExist("cond_h_z");
	}
}

float Cmat::omega_p_x_value(const int& p) const
{
  return omega_p_x(_omega_p_x_vec[p]);
  //may throw exception if pole doesn't exist
}

float Cmat::omega_p_y_value(const int& p) const
{
  return omega_p_y(_omega_p_y_vec[p]);
  //may throw exception if pole doesn't exist
}

float Cmat::omega_p_z_value(const int& p) const
{
  return omega_p_z(_omega_p_z_vec[p]);
  //may throw exception if pole doesn't exist
}

float Cmat::tau_p_x_value(const int& p) const
{
  return tau_p_x(_tau_p_x_vec[p]);
  //may throw exception if pole doesn't exist
}

float Cmat::tau_p_y_value(const int& p) const
{
  return tau_p_y(_tau_p_y_vec[p]);
  //may throw exception if pole doesn't exist
}

float Cmat::tau_p_z_value(const int& p) const
{
  return tau_p_z(_tau_p_z_vec[p]);
  //may throw exception if pole doesn't exist
}

float Cmat::Omega_p_x_value(const int& p) const
{
  return Omega_p_x(_Omega_p_x_vec[p]);
  //may throw exception if pole doesn't exist
}

float Cmat::Omega_p_y_value(const int& p) const
{
  return Omega_p_y(_Omega_p_y_vec[p]);
  //may throw exception if pole doesn't exist
}

float Cmat::Omega_p_z_value(const int& p) const
{
  return Omega_p_z(_Omega_p_z_vec[p]);
  //may throw exception if pole doesn't exist
}

//Returns the number of "effective" poles, meaning that
//the ones with zero pole frequency don't count
int Cmat::eff_num_lrntz_poles() const {
  int num_poles(0);
  for (int p(0); p<_num_lrntz_poles;++p) {
    //Omega_p=(sqrt(delta_eps) x pole frequency) is in the numerator,
    //so should be nonzero for the pole to contribute.
    //Same goes for tau_p (pole relaxation time).
    if ((Omega_p_x(_Omega_p_x_vec[p])!=0.0)&&(tau_p_x(_tau_p_x_vec[p])!=0.0)) {
      num_poles++;
    }
  }
  return num_poles;
}
