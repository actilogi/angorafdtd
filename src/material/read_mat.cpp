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

// reading the config options that specify the material definitions

#include "headers.h"

#include "read_mat.h"

#include "Cmats.h"


void read_mat(Cmats &Materials, const Config& fdtdconfig, const Config& validsettings)
{
	string material_setting_path = "Materials";
	if (fdtdconfig.exists(material_setting_path))
	{
		Setting& Materiallistsettings = read_list_from_group(fdtdconfig.getRoot(),material_setting_path);

		int num_of_newmaterials = Materiallistsettings.getLength();
		for (int newmaterialindex=0; newmaterialindex<num_of_newmaterials; newmaterialindex++)
		{
			Setting& Newmaterialsettings = Materiallistsettings[newmaterialindex];	//go to the newmaterialindex'th new material setting
			//check group for invalid settings
			CheckAngoraGroupSetting(Newmaterialsettings,validsettings);

			if (SettingEnabledForGrid(Newmaterialsettings))		//apply only if enabled for this grid
			{
				//create the material object
				Cmat NewMaterial;

				//is the material transparent? (are the unspecified parameters assigned default values, or left untouched)
				bool transparent;
				if (!read_optional_value_from_group<bool>(Newmaterialsettings,"transparent",transparent))
				{
					transparent = false; //by default, unspecified parameters are assigned default values
				}

				float rel_permittivity;
				if (!read_optional_value_from_group<float>(Newmaterialsettings,"rel_permittivity",rel_permittivity))
				{
					if (!transparent)
					{
						NewMaterial.set_eps(1.0); //default: eps_r = 1
					}// if transparent, don't set the permittivity
				}
				else
				{
					NewMaterial.set_eps(rel_permittivity);
				}

				float rel_permeability;
				if (!read_optional_value_from_group<float>(Newmaterialsettings,"rel_permeability",rel_permeability))
				{
					if (!transparent)
					{
						NewMaterial.set_mu(1.0); //default: mu_r = 1
					}
				}
				else
				{
					NewMaterial.set_mu(rel_permeability);
				}

				float electric_conductivity;
				if (!read_optional_value_from_group<float>(Newmaterialsettings,"electric_conductivity",electric_conductivity))
				{
					if (!transparent)
					{
						NewMaterial.set_cond_e(0); //default: sigma = 0
					}
				}
				else
				{
					NewMaterial.set_cond_e(electric_conductivity);
				}

				float magnetic_conductivity;
				if (!read_optional_value_from_group<float>(Newmaterialsettings,"magnetic_conductivity",magnetic_conductivity))
				{
					if (!transparent)
					{
						NewMaterial.set_cond_h(0); //default: sigma = 0
					}
				}
				else
				{
					NewMaterial.set_cond_h(magnetic_conductivity);
				}

        /********************************/
        /***     READ DRUDE POLES     ***/
        /********************************/
        bool dpf_exists(Newmaterialsettings.exists("drude_pole_frequency")),
             drt_exists(Newmaterialsettings.exists("drude_pole_relaxation_time"));
        //either all should exist or neither
        if (dpf_exists^drt_exists) {
          throw AngoraSettingException(Newmaterialsettings,
            "drude_pole_frequency and drude_pole_relaxation_time should both be specified");
        }
        if (dpf_exists&&drt_exists) {
          Setting& dpf = Newmaterialsettings["drude_pole_frequency"];
          Setting& drt = Newmaterialsettings["drude_pole_relaxation_time"];
          if (dpf.isNumber()&&drt.isNumber()) {
            float drude_pole_frequency(dpf),drude_pole_relaxation_time(drt);
            NewMaterial.add_drude_pole(drude_pole_frequency,drude_pole_relaxation_time);
          }
          else if (dpf.isArray()&&drt.isArray()) {
            int numpoles(dpf.getLength());
            if (numpoles!=drt.getLength()) {
              throw AngoraSettingException(Newmaterialsettings,
                    "drude_pole_frequency and drude_pole_relaxation_time should have "
                    "the same number of elements");
            }
            else {
              if ((!dpf[0].isNumber())||(!drt[0].isNumber())) {
                throw AngoraSettingException(Newmaterialsettings,
                "drude_pole_frequency and drude_pole_relaxation_time should be "
                "arrays of numbers");
              }
              else {
                for (int p(0); p<numpoles; ++p) {
                  float drude_pole_frequency(dpf[p]),drude_pole_relaxation_time(drt[p]);
                  NewMaterial.add_drude_pole(drude_pole_frequency,drude_pole_relaxation_time);
                }
              }
            }
          }
          else {
            throw AngoraSettingException(Newmaterialsettings,
            "drude_pole_frequency and drude_pole_relaxation_time should both be "
            "numbers or arrays of numbers");
          }
        }
        else {
          if (!transparent) {
            NewMaterial.add_drude_pole(0.0,0.0);
          }
        }

        /********************************/
        /***     READ LORENTZ POLES     ***/
        /********************************/
        bool lpf_exists(Newmaterialsettings.exists("lorentz_pole_frequency")),
             ldf_exists(Newmaterialsettings.exists("lorentz_pole_damping_factor")),
             lde_exists(Newmaterialsettings.exists("lorentz_delta_epsilon"));
        //either all should exist or neither
        if (!((lpf_exists&&ldf_exists&&lde_exists)||(!lpf_exists&&!ldf_exists&&!lde_exists))) {
          throw AngoraSettingException(Newmaterialsettings,
            "lorentz_pole_frequency, lorentz_pole_damping_factor, and lorentz_delta_epsilon "
            "should all be specified");
        }
        if (lpf_exists&&ldf_exists&&lde_exists) {
          Setting& lpf = Newmaterialsettings["lorentz_pole_frequency"];
          Setting& ldf = Newmaterialsettings["lorentz_pole_damping_factor"];
          Setting& lde = Newmaterialsettings["lorentz_delta_epsilon"];
          if (lpf.isNumber()&&ldf.isNumber()&&lde.isNumber()) {
            float lorentz_pole_frequency(lpf),lorentz_pole_damping_factor(ldf),
                  lorentz_delta_epsilon(lde);
            if (lorentz_pole_damping_factor==0.0)
              throw AngoraSettingException(Newmaterialsettings,
                    "lorentz_pole_damping_factor cannot be zero");
            NewMaterial.add_lorentz_pole(lorentz_pole_frequency,
                                        lorentz_pole_damping_factor,
                                        lorentz_delta_epsilon);
          }
          else if (lpf.isArray()&&ldf.isArray()&&lde.isArray()) {
            int numpoles(lpf.getLength());
            if (numpoles!=ldf.getLength()||numpoles!=lde.getLength()) {
              throw AngoraSettingException(Newmaterialsettings,
                    "lorentz_pole_frequency, lorentz_pole_damping_factor, and "
                    "lorentz_delta_epsilon should have "
                    "the same number of elements");
            }
            else {
              if ((!lpf[0].isNumber())||(!ldf[0].isNumber())||(!lde[0].isNumber())) {
                throw AngoraSettingException(Newmaterialsettings,
                "lorentz_pole_frequency, lorentz_pole_damping_factor, and "
                "lorentz_delta_epsilon should be arrays of numbers");
              }
              else {
                for (int p(0); p<numpoles; ++p) {
                  float lorentz_pole_frequency(lpf[p]),lorentz_pole_damping_factor(ldf[p]),
                        lorentz_delta_epsilon(lde[p]);
                  if (lorentz_pole_damping_factor==0.0)
                    throw AngoraSettingException(Newmaterialsettings,
                          "lorentz_pole_damping_factor cannot be zero");
                  NewMaterial.add_lorentz_pole(lorentz_pole_frequency,
                                               lorentz_pole_damping_factor,
                                               lorentz_delta_epsilon);
                }
              }
            }
          }
          else {
            throw AngoraSettingException(Newmaterialsettings,
            "lorentz_pole_frequency, lorentz_pole_damping_factor, and "
            "lorentz_delta_epsilon should all be "
            "numbers or arrays of numbers");
          }
        }
        else {
          if (!transparent) {
            NewMaterial.add_lorentz_pole(0.0,0.0,0.0);
          }
        }

				/********************************/

				string material_tag;
				read_value_from_group<string>(Newmaterialsettings,"material_tag",material_tag);

//				AddIsotropicMaterial(NewMaterialId,rel_permittivity,rel_permeability,electric_conductivity,magnetic_conductivity);
				//Add the new material to the material-collector object (even in check mode, since other objects refer to these materials)
				Materials.CreateMaterial(NewMaterial,material_tag);
			}
		}
	}//if exists
}
