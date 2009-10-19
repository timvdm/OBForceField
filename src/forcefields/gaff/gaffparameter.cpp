/*********************************************************************
GAFF force field

Copyright (C) 2006-2009 by Tim Vandermeersch
Copyright (C) 2009 by Frank Peters
 
This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <openbabel/locale.h>
#include <openbabel/tokenst.h>
#include <openbabel/oberror.h>
#include <openbabel/mol.h>
#include <fstream>

#include "gafftype.h"
#include "gaffparameter.h"

using namespace std;

namespace OpenBabel {
  namespace OBFFs {
  
    GAFFParameterDB::GAFFParameterDB(const std::string &filename) : _initialized(false), _filename(filename)
    {
      OBFFParameterDB("GAFF");
      _initialized = ParseParamFile();
    }

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    //
    //  Parse parameter files
    //
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
 
    bool GAFFParameterDB::ParseParamFile()
    {
      // Set the locale for number parsing to avoid locale issues: PR#1785463
      obLocale.SetLocale();

      string s, sbuffer;
      vector<string> vs;
      char buffer[256];
      bool valid_line=true;
      unsigned int i;
      std::vector<OBVariant> parameter;
      OBFFTable * p;
      std::vector<std::string> header;

      // open data/_parFile
      ifstream ifs;
      if (OpenDatafile(ifs, _filename).length() == 0) {
	obErrorLog.ThrowError(__FUNCTION__, "Cannot open parameter file", obError);
	return false;
      }
      header.reserve(8);
        
      valid_line= ifs.getline(buffer, BUFF_SIZE); //first line title block

      header.clear();
      header.push_back("atom type");
      header.push_back("mass");
      header.push_back("atpol");
      p = AddTable("Atom Properties", header);
      valid_line= ifs.getline(buffer, BUFF_SIZE); //next line
      while(valid_line && strlen(buffer)!=0) //read block 2 (atom mass and polarizability)
	{
	  tokenize(vs, buffer," \n\t");
	  parameter.clear();
	  parameter.push_back(OBVariant(vs[0], "type")); //KNDSYM
	  parameter.push_back(OBVariant(strtod(vs[1].c_str(),0),"mass")); // AMASS
	  parameter.push_back(OBVariant(strtod(vs[2].c_str(),0),"atpol")); // ATPOL [A^3]
	  p->AddRow(parameter);
	  valid_line= ifs.getline(buffer, BUFF_SIZE);
	}
      valid_line= ifs.getline(buffer, BUFF_SIZE); //next line (block 3: line contains JSOLTY)

      header.clear();
      header.push_back("bond type");
      header.push_back("type1");
      header.push_back("type2");
      header.push_back("K");
      header.push_back("r0");
      p = AddTable("Bond Harmonic", header); //Table header needs to correspond with header used in query!!
      valid_line= ifs.getline(buffer, BUFF_SIZE); //next line
      while(valid_line && strlen(buffer)!=0) //read block 4 (bonds)
	{
	  sbuffer=buffer;
	  tokenize(vs, sbuffer," -",2);
	  parameter.clear();
	  s = GAFFType::MakeBondName(vs[0],vs[1]);
	  parameter.push_back(OBVariant(s, "bh")); //  the bond stretch potential is bond-harmonic
	  parameter.push_back(OBVariant(vs[0], "type1"));
	  parameter.push_back(OBVariant(vs[1], "type2"));
	  sbuffer = vs[2];
	  vs.clear();
	  tokenize(vs, sbuffer," \n\t",2);
	  parameter.push_back(OBVariant(strtod(vs[0].c_str(),0), "K")); // RK [kcal/mol/(A^2)]
	  parameter.push_back(OBVariant(strtod(vs[1].c_str(),0), "r0")); // REQ [A]
	  p->AddRow(parameter);
	  valid_line= ifs.getline(buffer, BUFF_SIZE); //next line
	}
      // add default value
      parameter.clear();
      s = GAFFType::MakeBondName("X","X");
      parameter.push_back(OBVariant(s, "bh"));
      parameter.push_back(OBVariant("X", "type1"));
      parameter.push_back(OBVariant("X", "type2"));
      parameter.push_back(OBVariant(500.0, "K")); // RK [kcal/mol/(A^2)]
      parameter.push_back(OBVariant(1.1, "r0")); // REQ [A]
      p->AddRow(parameter);

      header.clear();
      header.push_back("angle type");
      header.push_back("type1");
      header.push_back("type2");
      header.push_back("type3");
      header.push_back("ka");
      header.push_back("theta0");
      p = AddTable("Angle Harmonic", header);
      valid_line= ifs.getline(buffer, BUFF_SIZE); //next line
      while(valid_line && strlen(buffer)!=0) //read block 5 (angles)
	{
	  sbuffer=buffer;
	  tokenize(vs, sbuffer," -",3);
	  parameter.clear();
	  s = GAFFType::MakeAngleName(vs[0],vs[1],vs[2]);
	  parameter.push_back(OBVariant(s, "ah")); //potential is angle harmonic
	  parameter.push_back(OBVariant(vs[0], "type1")); //ITT
	  parameter.push_back(OBVariant(vs[1], "type2")); //JTT
	  parameter.push_back(OBVariant(vs[2], "type3")); //KTT
	  sbuffer = vs[3];
	  vs.clear();
	  tokenize(vs, sbuffer," \n\t",2);
	  parameter.push_back(OBVariant(strtod(vs[0].c_str(),0), "K")); // TK [kcal/mol/(rad**2)]
	  parameter.push_back(OBVariant(strtod(vs[1].c_str(),0), "theta0")); // TEQ [degrees]
	  p->AddRow(parameter);
	  valid_line= ifs.getline(buffer, BUFF_SIZE);
	}
      // add default value
      parameter.clear();
      s = GAFFType::MakeAngleName("X","X","X");
      parameter.push_back(OBVariant(s, "ah")); //potential is angle harmonic
      parameter.push_back(OBVariant("X" , "type1")); //ITT
      parameter.push_back(OBVariant("X", "type2")); //JTT
      parameter.push_back(OBVariant("X", "type3")); //KTT
      parameter.push_back(OBVariant(0.0, "K")); // TK [kcal/mol/(rad**2)]
      parameter.push_back(OBVariant(180, "theta0")); // TEQ [degrees]
      p->AddRow(parameter);

      header.clear();
      header.push_back("torsion type");
      header.push_back("type1");
      header.push_back("type2");
      header.push_back("type3");
      header.push_back("type4");
      header.push_back("K");
      header.push_back("d");
      header.push_back("n");
      p = AddTable("Torsion Harmonic", header);
      valid_line= ifs.getline(buffer, BUFF_SIZE);
      while(valid_line && strlen(buffer)!=0) //read block 6 (dihedrals)
	{
	  // Gaff torsional potential (PK/IDIVF) * (1 + cos(PN*phi - GAMMA))
          // Harmonic torsion potetial is:  K * (1 + d * cos( n * phi)) 
	  sbuffer=buffer;
	  tokenize(vs, sbuffer," -",4);
	  parameter.clear();
	  s = GAFFType::MakeTorsionName(vs[0],vs[1],vs[2],vs[3]);
	  parameter.push_back(OBVariant(s, "th")); //potential is torsion harmonic
	  parameter.push_back(OBVariant(vs[0], "type1")); //IPT
	  parameter.push_back(OBVariant(vs[1], "type2")); //JPT
	  parameter.push_back(OBVariant(vs[2], "type3")); //KPT
	  parameter.push_back(OBVariant(vs[3], "type4")); //LPT
	  sbuffer = vs[4];
	  vs.clear();
	  tokenize(vs, sbuffer," \n\t",4);
          // vs[0] = IDIVF, PK = vs[1], GAMMA = vs[2] (always 0 or 180), PN = vs[3]
	  parameter.push_back(OBVariant(strtod(vs[1].c_str(),0)/strtod(vs[0].c_str(),0),"K")); // K
	  if (strtod(vs[2].c_str(),0)==0)
	    parameter.push_back(OBVariant(1.0,"d")); //d
	  else
	    parameter.push_back(OBVariant(-1.0,"d")); //d
	  parameter.push_back(OBVariant(strtod(vs[3].c_str(),0),"n")); // N
	  p->AddRow(parameter);
	  valid_line= ifs.getline(buffer, BUFF_SIZE);
	}

      header.clear();
      header.push_back("oop type");
      header.push_back("type1");
      header.push_back("type2");
      header.push_back("type3");
      header.push_back("type4");
      header.push_back("K");
      header.push_back("d");
      header.push_back("n");
      p = AddTable("Torsion Harmonic OOP" ,header);
      valid_line= ifs.getline(buffer, BUFF_SIZE); //next line 
      while(valid_line && strlen(buffer)!=0) //read block 7(impropers)
	{
	  //         The input is the same as in for the dihedrals except that
	  // 	   the torsional barrier height is NOT divided by the factor
	  //         idivf.  The improper torsions are defined between any four
	  //         atoms not bonded (in a successive fashion) with each other
	  //         as in the case of "regular" or "proper" dihedrals.  Improper
	  //         dihedrals are used to keep certain groups planar and to
	  //         prevent the racemization of certain centers in the united
	  //         atom model.  Consult the above reference for details.
	  //         From Townhee help:
	  // Improper torsions are not automatically generated by the Towhee code as the rules for determining where they are applied are not always straight-forward. Amber param96 exclusively uses the Stereocenter version of the improper torsions, and they are typically centered on an sp2 atom in order to enforce planarity with its three neighbors. Only one improper torsion allowed to be centered on any atom. These torsions are listed in the Amber literature as i-j-k-l where the angle is the dihedral between i-k-l and j-k-l, and the bonding pattern is i, j, and l are all bonded to atom k, and are also not bonded to each other. In the towhee_input file this stereo improper torsion is listed only for atom k, and the atom order there is l, i, j. Remember that you can set the improper type to 0 to have the code automatically determine the improper type (so long as inpstyle is 2).
	  
	  // Gaff oop-torsional potential PK * (1 + cos(PN*phi - GAMMA))
          // Harmonic torsion potetial is:  K * (1 + d * cos( n * phi)) 
	  sbuffer=buffer;
	  tokenize(vs, sbuffer," -",4);
	  parameter.clear();
	  s = GAFFType::MakeOOPName(vs[0],vs[1],vs[2],vs[4]);
	  parameter.push_back(OBVariant(s, "th")); //potential torsion harmonic (same as for proper dihedrals torsions)
	  parameter.push_back(OBVariant(vs[0], "type1")); //IPT
	  parameter.push_back(OBVariant(vs[1], "type2")); //JPT
	  parameter.push_back(OBVariant(vs[2], "type3")); //KPT
	  parameter.push_back(OBVariant(vs[3], "type4")); //LPT
	  sbuffer = vs[4];
	  vs.clear();
	  tokenize(vs, sbuffer," \n\t",3);
	  parameter.push_back(OBVariant(strtod(vs[0].c_str(),0),"K")); // K
	  if (strtod(vs[2].c_str(),0)==0)
	    parameter.push_back(OBVariant(1.0,"d")); //d
	  else
	    parameter.push_back(OBVariant(-1.0,"d")); //d
	  parameter.push_back(OBVariant(strtod(vs[2].c_str(),0),"n")); // n
	  p->AddRow(parameter);
	  valid_line= ifs.getline(buffer, BUFF_SIZE);
	}

      header.clear();
      header.push_back("bond type");
      header.push_back("type1");
      header.push_back("type2");
      header.push_back("A");
      header.push_back("B");
      p = AddTable("Hydrogen Bond", header);
      valid_line= ifs.getline(buffer, BUFF_SIZE); //next line 
      while(valid_line && strlen(buffer)!=0) //read block 8 (10-12 H-BOND)
	{
	  sbuffer=buffer;
	  tokenize(vs, sbuffer," -",2);
	  parameter.clear();
	  s = GAFFType::MakeBondName(vs[0],vs[1]);
	  parameter.push_back(OBVariant(s, "lj10-12")); //  In AMBER forcefields 10-12 Lennard-Jones potentials are used for Hydrogen Bonds (not in GAFF)
	  parameter.push_back(OBVariant(vs[0], "type1")); //KT1
	  parameter.push_back(OBVariant(vs[1], "type2")); //KT2
	  sbuffer = vs[2];
	  vs.clear();
	  tokenize(vs, sbuffer," \n\t",2);
	  parameter.push_back(OBVariant(strtod(vs[0].c_str(),0), "A"));
	  parameter.push_back(OBVariant(strtod(vs[1].c_str(),0), "B"));
	  p->AddRow(parameter);
	  valid_line= ifs.getline(buffer, BUFF_SIZE);
	}

      valid_line= ifs.getline(buffer, BUFF_SIZE); //next line 
      while(valid_line && strlen(buffer)!=0) //read block 9 (equivalent atoms for vdW potential)
	{
	  // not used. We assume the types are all listed individually in block 10
	  valid_line= ifs.getline(buffer, BUFF_SIZE);
	}
      valid_line= ifs.getline(buffer, BUFF_SIZE); //next line 

      const double two_power_minus_one_sixth = pow(2.,-1./6.);
      header.clear();
      header.push_back("atom type");
      header.push_back("type1");
      header.push_back("sigma");
      header.push_back("epsilon");
      p = AddTable("LJ6_12", header);
      valid_line= ifs.getline(buffer, BUFF_SIZE); //'RE' means vdW radius and well-depth are read
      while(valid_line && strlen(buffer)!=0) //read block 10 (vdWaals)
	{
	  // van der Waals potental GAFF: EDEP * ( (R/r)^12 - 2*(R/r)^6 ), Geometric mixing rules are used, 1-4 interactions reduced by factor 1/2
	  // LJ6-12 potential: 4*epsilon * ( (sigma/r)^12 - (sigma/r)^6 )

	  tokenize(vs, buffer," \n\t");
	  parameter.clear();
	  parameter.push_back(OBVariant(vs[0], "type"));
	  parameter.push_back(OBVariant(two_power_minus_one_sixth*strtod(vs[1].c_str(),0), "sigma")); // sigma = 2^(-1/6) * R (Angstrom)
	  parameter.push_back(OBVariant(strtod(vs[2].c_str(),0), "epsilon")); // epsilon = EDEP (kcal/mol)
	  p->AddRow(parameter);
	  valid_line= ifs.getline(buffer, BUFF_SIZE);
	}
      // add default value
      parameter.clear();
      parameter.push_back(OBVariant("X", "type"));
      parameter.push_back(OBVariant(1.4870, "sigma")); // R
      parameter.push_back(OBVariant(0.0157, "epsilon")); // EDEP (kcal/mol)
      p->AddRow(parameter);
      
      if (ifs)
	ifs.close();
      
      // return the locale to the original one
      obLocale.RestoreLocale();
      return true;
    }
  } // end namespace OpenBabel
}
//! \brief GAFF force field

