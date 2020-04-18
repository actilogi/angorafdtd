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

#ifndef CMATS_H
#define CMATS_H

//Declaration of the class "Cmats" that holds the named-material information

#include "headers.h"

#include "Cmats_excp.h"

//for the C++ STL map class
#include <map>

#include "Cmat.h"


class Cmats
{
 public:
	 void CreateMaterial(const Cmat& material, const string& material_tag);

	 bool lookupMaterialWithTag(const string& material_tag, Cmat& material) const;  //copies the material identifier that corresponds to tag into material, returns false if tag is not found

	 const Cmat operator[] (const string& material_tag) const; //returns the material ID corresponding to the material tag, throws exception if not found

	 int NumberOfMaterials() const
	 {
		 return NamedMaterials.size();
	 }

	 //maximum effective number of Lorentz poles among the materials
	 int max_eff_num_lrntz_poles() const;

private:
	 //C++ STL map object that holds the named-material information
	 map<string,Cmat> NamedMaterials;

	 bool MaterialTagExists(const string& material_tag);
};

#endif // CMATS_H
