/*********************************************************************
forcefieldmmff94.cpp - MMFF94 force field

Copyright (C) 2006-2008 by Tim Vandermeersch <tim.vandermeersch@gmail.com>
 
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


#include "common.h"
#include "parameter.h"

#include <iomanip>
#include <openbabel/mol.h>

using namespace std;

namespace OpenBabel {
namespace OBFFs {

  MMFF94Common::MMFF94Common()
  {
    m_database = static_cast<OBParameterDB*>(new MMFF94SimpleParameterDB("/home/timvdm/OBForceField/git_repo/data/mmff94.ff"));
  }
  
  MMFF94Common::~MMFF94Common()
  {
  
  }

  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Setup Functions
  //
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
 
  // The MMFF94 article doesn't seem to include information about how 
  // aromaticity is perceived. This function was written by studying the 
  // MMFF_opti.log file, trail-and-error and using the MMFF94 validation
  // set to check the results (If all atom types are assigned correctly, 
  // aromatic rings are probably detected correctly)
  bool MMFF94Common::PerceiveAromaticity(/*const*/ OBMol &mol)
  {
    bool done = false; // not done actually....
    OBAtom *ringatom;
    OBBond *ringbond;
    vector<OBRing*> vr;
    vr = mol.GetSSSR();
    
    vector<OBRing*>::iterator ri;
    vector<int>::iterator rj;
    int n, index, ringsize, first_rj, prev_rj, pi_electrons;
    first_rj = prev_rj = index = 0;
    for (ri = vr.begin();ri != vr.end();++ri) { // for each ring
      ringsize = (*ri)->Size();
      
      n = 1;
      pi_electrons = 0;
      for(rj = (*ri)->_path.begin();rj != (*ri)->_path.end();rj++) { // for each ring atom
        index = *rj;
        ringatom = mol.GetAtom(index);
        
        // is the bond to the previous ring atom double?
        if (n > 1) {
          ringbond = mol.GetBond(prev_rj, index);
          if (!ringbond) {
            prev_rj = index;
            continue;
          }
          if (ringbond->GetBondOrder() == 2) {
            pi_electrons += 2;
            prev_rj = index;
            n++;
            continue;
          }
          prev_rj = index;
        } else {
          prev_rj = index;
          first_rj = index;
        }
	
        // does the current ring atom have a exocyclic double bond?
        FOR_NBORS_OF_ATOM (nbr, ringatom) {
          if ((*ri)->IsInRing(nbr->GetIdx()))
            continue;

          if (!IsAromatic(&*nbr))
            continue;

          ringbond = mol.GetBond(nbr->GetIdx(), index);
          if (!ringbond) {
            continue;
          }
          if (ringbond->GetBondOrder() == 2)
            pi_electrons++;
        }

        // is the atom N, O or S in 5 rings
        if (ringsize == 5 &&
            ringatom->GetIdx() == (*ri)->GetRootAtom()) 
          pi_electrons += 2;

        n++;
      
      } // for each ring atom
      
      // is the bond from the first to the last atom double?
      ringbond = mol.GetBond(first_rj, index);
      if (ringbond) {
        if (ringbond->GetBondOrder() == 2) 
          pi_electrons += 2;
      }

      if ((pi_electrons == 6) && ((ringsize == 5) || (ringsize == 6))) {
        // mark ring atoms as aromatic
        for(rj = (*ri)->_path.begin();rj != (*ri)->_path.end();rj++) {
          if (!IsAromatic(mol.GetAtom(*rj)))
            done = true;
          m_aromAtoms[mol.GetAtom(*rj)->GetIdx()-1] = 1;
        }
        // mark all ring bonds as aromatic
        FOR_BONDS_OF_MOL (bond, mol) 
          if((*ri)->IsMember(&*bond))
            m_aromBonds[bond->GetIdx()] = 1;
      }
    }

    return done;
  }
  
  // Symbolic atom typing is skipped
  // 
  // atom typing is based on:
  //   MMFF94 I - Table III
  //   MMFF94 V - Table I
  //
  int MMFF94Common::GetType(OBAtom *atom)
  {
    OBMol *mol = atom->GetParent();
    OBBond *bond;
    int oxygenCount, nitrogenCount, sulphurCount, doubleBondTo;
    ////////////////////////////////
    // Aromatic Atoms
    ////////////////////////////////
    if (IsAromatic(atom)) {
      if (atom->IsInRingSize(5)) {
        bool isAromatic = false;
        vector<OBAtom*> alphaPos, betaPos;
        vector<OBAtom*> alphaAtoms, betaAtoms;

        if (atom->IsSulfur()) {
          return 44; // Aromatic 5-ring sulfur with pi lone pair (STHI)
        }
        if (atom->IsOxygen()) {
          return 59; // Aromatic 5-ring oxygen with pi lone pair (OFUR)
        }
        if (atom->IsNitrogen()) {
          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (nbr->IsOxygen() && (nbr->GetValence() == 1)) {
              return 82; // N-oxide nitrogen in 5-ring alpha position, 
              // N-oxide nitrogen in 5-ring beta position, 
              // N-oxide nitrogen in other 5-ring  position, 
              // (N5AX, N5BX, N5OX) 
            }
          }
        }
        FOR_NBORS_OF_ATOM (nbr, atom) {
          if (!IsAromatic(mol->GetBond(atom, &*nbr)) || !nbr->IsInRingSize(5))
            continue;
   
          if (IsInSameRing(atom, &*nbr)) {
            alphaPos.push_back(&*nbr);
          }
          
          FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
            if (nbrNbr->GetIdx() == atom->GetIdx())
              continue;
            if (!IsAromatic(mol->GetBond(&*nbr, &*nbrNbr)) || !nbrNbr->IsInRingSize(5))
              continue;
             
            isAromatic = true;
	    
            if (IsInSameRing(atom, &*nbrNbr)) {
              betaPos.push_back(&*nbrNbr);
            }
          }
        }
	
        if (isAromatic) {
          
	  
          for (unsigned int i = 0; i < alphaPos.size(); i++) {
            if (alphaPos[i]->IsSulfur()) {
              alphaAtoms.push_back(alphaPos[i]);
            } else if (alphaPos[i]->IsOxygen()) {
              alphaAtoms.push_back(alphaPos[i]);
            } else if (alphaPos[i]->IsNitrogen() && (alphaPos[i]->GetValence() == 3)) {
              bool IsNOxide = false;
              FOR_NBORS_OF_ATOM (nbr, alphaPos[i]) {
                if (nbr->IsOxygen() && (nbr->GetValence() == 1)) {
                  IsNOxide = true;
                }
              }

              if (!IsNOxide) {
                alphaAtoms.push_back(alphaPos[i]);
              }
            }
          }
          for (unsigned int i = 0; i < betaPos.size(); i++) {
            if (betaPos[i]->IsSulfur()) {
              betaAtoms.push_back(betaPos[i]);
            } else if (betaPos[i]->IsOxygen()) {
              betaAtoms.push_back(betaPos[i]);
            } else if (betaPos[i]->IsNitrogen() && (betaPos[i]->GetValence() == 3)) {
              bool IsNOxide = false;
              FOR_NBORS_OF_ATOM (nbr, betaPos[i]) {
                if (nbr->IsOxygen() && (nbr->GetValence() == 1)) {
                  IsNOxide = true;
                }
              }

              if (!IsNOxide) {
                betaAtoms.push_back(betaPos[i]);
              }
            }
          }
          if (!betaAtoms.size()) {
            nitrogenCount = 0;
            FOR_NBORS_OF_ATOM (nbr, atom) {
              if (nbr->IsNitrogen() && (nbr->GetValence() == 3)) {
                if ((nbr->BOSum() == 4) && IsAromatic(&*nbr)) {
                  nitrogenCount++;
                } else if ((nbr->BOSum() == 3) && !IsAromatic(&*nbr)) {
                  nitrogenCount++;
                }
              }
            }
            if (nitrogenCount >= 2) {
              return 80; // Aromatic carbon between N's in imidazolium (CIM+)
            }
          }
          if (!alphaAtoms.size() && !betaAtoms.size()) {
            if (atom->IsCarbon()) {
              // there is no S:, O:, or N:
              // this is the case for anions with only carbon and nitrogen in the ring
              return 78; // General carbon in 5-membered aromatic ring (C5)
            } else if (atom->IsNitrogen()) {
              if (atom->GetValence() == 3) {
                // this is the N: atom
                return 39; // Aromatic 5 ring nitrogen with pi lone pair (NPYL)
              } else {
                // again, no S:, O:, or N:
                return 76; // Nitrogen in 5-ring aromatic anion (N5M)
              }
            }
          }
          if (alphaAtoms.size() == 2) {
            if (atom->IsCarbon() && IsInSameRing(alphaAtoms[0], alphaAtoms[1])) {
              if (alphaAtoms[0]->IsNitrogen() && alphaAtoms[1]->IsNitrogen()) {
                if ((alphaAtoms[0]->GetValence() == 3) && (alphaAtoms[1]->GetValence() == 3)) {
                  return 80; // Aromatic carbon between N's in imidazolium (CIM+)
                }
              }
            }
          }
          if (alphaAtoms.size() && !betaAtoms.size()) {
            if (atom->IsCarbon()) {
              return 63; // Aromatic 5-ring C, alpha to N:, O:, or S: (C5A)
            } else if (atom->IsNitrogen()) {
              if (atom->GetValence() == 3) {
                return 81; // Posivite nitrogen in 5-ring alpha position (N5A+)
              } else {
                return 65; // Aromatic 5-ring N, alpha to N:, O:, or S: (N5A)
              }
            }
          }
          if (!alphaAtoms.size() && betaAtoms.size()) {
            if (atom->IsCarbon()) {
              return 64; // Aromatic 5-ring C, beta to N:, O:, or S: (C5B)
            } else if (atom->IsNitrogen()) {
              if (atom->GetValence() == 3) {
                return 81; // Posivite nitrogen in 5-ring beta position (N5B+)
              } else {
                return 66; // Aromatic 5-ring N, beta to N:, O:, or S: (N5B)
              }
            }
          }
          if (alphaAtoms.size() && betaAtoms.size()) {
            for (unsigned int i = 0; i < alphaAtoms.size(); i++) {
              for (unsigned int j = 0; j < betaAtoms.size(); j++) {
                if (!IsInSameRing(alphaAtoms[i], betaAtoms[j])) {
                  if (atom->IsCarbon()) {
                    return 78; // General carbon in 5-membered aromatic ring (C5)
                  } else if (atom->IsNitrogen()) {
                    return 79; // General nitrogen in 5-membered aromatic ring (N5)
                  }
                }
              }
            }
            for (unsigned int i = 0; i < alphaAtoms.size(); i++) {
              if (alphaAtoms[i]->IsSulfur() || alphaAtoms[i]->IsOxygen()) {
                if (atom->IsCarbon()) {
                  return 63; // Aromatic 5-ring C, alpha to N:, O:, or S: (C5A)
                } else if (atom->IsNitrogen()) {
                  return 65; // Aromatic 5-ring N, alpha to N:, O:, or S: (N5A)
                }
              }
            }
            for (unsigned int i = 0; i < betaAtoms.size(); i++) {
              if (betaAtoms[i]->IsSulfur() || betaAtoms[i]->IsOxygen()) {
                if (atom->IsCarbon()) {
                  return 64; // Aromatic 5-ring C, beta to N:, O:, or S: (C5B)
                } else if (atom->IsNitrogen()) {
                  return 66; // Aromatic 5-ring N, beta to N:, O:, or S: (N5B)
                }
              }
            }
	    
            if (atom->IsCarbon()) {
              return 78; // General carbon in 5-membered aromatic ring (C5)
            } else if (atom->IsNitrogen()) {
              return 79; // General nitrogen in 5-membered aromatic ring (N5)
            }
          }
        }
      }
    
      if (atom->IsInRingSize(6)) {
	
        if (atom->IsCarbon()) {
          return 37; // Aromatic carbon, e.g., in benzene (CB)
        } else if (atom->IsNitrogen()) {
          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (nbr->IsOxygen() && (nbr->GetValence() == 1)) {
              return 69; // Pyridinium N-oxide nitrogen (NPOX)
            }
          }
	  
          if (atom->GetValence() == 3) {
            return 58; // Aromatic nitrogen in pyridinium (NPD+)
          } else {
            return 38; // Aromatic nitrogen with sigma lone pair (NPYD)
          }
        }
      }
    }
    
    ////////////////////////////////
    // Hydrogen
    ////////////////////////////////
    if (atom->GetAtomicNum() == 1) {
      FOR_NBORS_OF_ATOM (nbr, atom) {
        if (nbr->IsCarbon()) {
          return 5; // Hydrogen attatched to carbon (HC)
        }
        if (nbr->GetAtomicNum() == 14) {
          return 5; // Hydrogen attatched to silicon (HSI)
        }
        if (nbr->IsOxygen()) {
          if (nbr->BOSum() == 3) {
            if (nbr->GetValence() == 3) {
              return 50; // Hydrogen on oxonium oxygen (HO+)
            } else {
              return 52; // Hydrogen on oxenium oxygen (HO=+)
            }
          }
	  
          int hydrogenCount = 0;
          FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
            if (nbrNbr->IsHydrogen()) {
              hydrogenCount++;
              continue;
            }
            if (nbrNbr->IsCarbon()) {
              if (IsAromatic(&*nbrNbr)) {
                return 29; // phenol
              }
                     
              FOR_NBORS_OF_ATOM (nbrNbrNbr, &*nbrNbr) {
                if (nbrNbrNbr->GetIdx() == nbr->GetIdx())
                  continue;

                bond = mol->GetBond(&*nbrNbr, &*nbrNbrNbr);
                if (bond->IsDouble()) {
                  if (nbrNbrNbr->IsOxygen()) {
                    return 24; // Hydroxyl hydrogen in carboxylic acids (HOCO)
                  }
                  if (nbrNbrNbr->IsCarbon() || nbrNbrNbr->IsNitrogen()) {
                    return 29; // Enolic or phenolic hydroxyl hydrogen,
                    // Hydroxyl hydrogen in HO-C=N moiety (HOCC, HOCN)
                  }
                }
              }
            }
            if (nbrNbr->IsPhosphorus()) {
              return 24; // Hydroxyl hydrogen in H-O-P moiety (HOP)
            }
            if (nbrNbr->IsSulfur()) {
              return 33; // Hydrogen on oxygen attached to sulfur (HOS)
            }
	  
          }
          if (hydrogenCount == 2) {
            return 31; // Hydroxyl hydrogen in water (HOH)
          }

          return 21; // Hydroxyl hydrogen in alcohols, Generic hydroxyl hydrogen (HOR, HO)
        }
        if (nbr->IsNitrogen()) {
          switch (GetType(&*nbr)) {
          case 81:
            return 36; // Hydrogen on imidazolium nitrogen (HIM+)
          case 68:
            return 23; // Hydrogen on N in N-oxide (HNOX)
          case 67:
            return 23; // Hydrogen on N in N-oxide (HNOX)
          case 62:
            return 23; // Generic hydrogen on sp3 nitrogen, e.g., in amines (HNR)
          case 56:
            return 36; // Hydrogen on guanimdinium nitrogen (HGD+)
          case 55:
            return 36; // Hydrogen on amidinium nitrogen (HNN+)
          case 43:
            return 28; // Hydrogen on NSO, NSO2, or NSO3 nitrogen, Hydrogen on N triply bonded to C (HNSO, HNC%)
          case 39:
            return 23; // Hydrogen on nitrogen in pyrrole (HPYL)
          case 8:
            return 23; // Generic hydrogen on sp3 nitrogen, e.g., in amines, Hydrogen on nitrogen in ammonia (HNR, H3N)
          }

          if (nbr->BOSum() == 4) {
            if (nbr->GetValence() == 2) {
              return 28; // Hydrogen on N triply bonded to C (HNC%)
            } else {
              return 36; // Hydrogen on pyridinium nitrogen, Hydrogen on protonated imine nitrogen (HPD+, HNC+)
            }
          }

          if (nbr->GetValence() == 2) {
            FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
              if (nbrNbr->IsHydrogen())
                continue;
              
              bond = mol->GetBond(&*nbr, &*nbrNbr);
              if (bond->IsDouble()) {
                if (nbrNbr->IsCarbon() || nbrNbr->IsNitrogen()) {
                  return 27; // Hydrogen on imine nitrogen, Hydrogen on azo nitrogen (HN=C, HN=N) 
                }

                return 28; // Generic hydrogen on sp2 nitrogen (HSP2)
              }
            }
          }
	  
          FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
            if (nbrNbr->IsHydrogen())
              continue;
	    
            if (nbrNbr->IsCarbon()) {
              if (IsAromatic(&*nbrNbr)) {
                return 28; // deloc. lp pair
              }
	      
              FOR_NBORS_OF_ATOM (nbrNbrNbr, &*nbrNbr) {
                if (nbrNbrNbr->GetIdx() == nbr->GetIdx())
                  continue;
              
                bond = mol->GetBond(&*nbrNbr, &*nbrNbrNbr);
                if (bond->IsDouble()) {
                  if (nbrNbrNbr->IsCarbon() || nbrNbrNbr->IsNitrogen() || nbrNbrNbr->IsOxygen() || nbrNbrNbr->IsSulfur()) {
                    return 28; // Hydrogen on amide nitrogen, Hydrogen on thioamide nitrogen,
                    // Hydrogen on enamine nitrogen, Hydrogen in H-N-C=N moiety (HNCO, HNCS, HNCC, HNCN)
                  }
                }
              }
            }
            if (nbrNbr->IsNitrogen()) {
              FOR_NBORS_OF_ATOM (nbrNbrNbr, &*nbrNbr) {
                if (nbrNbrNbr->GetIdx() == nbr->GetIdx())
                  continue;
              
                bond = mol->GetBond(&*nbrNbr, &*nbrNbrNbr);
                if (bond->IsDouble()) {
                  if (nbrNbrNbr->IsCarbon() || nbrNbrNbr->IsNitrogen()) {
                    return 28; // Hydrogen in H-N-N=C moiety, Hydrogen in H-N-N=N moiety (HNNC, HNNN)
                  }
                }
              }
            }
            if (nbrNbr->IsSulfur()) {
              FOR_NBORS_OF_ATOM (nbrNbrNbr, &*nbrNbr) {
                if (nbrNbrNbr->GetIdx() == nbr->GetIdx())
                  continue;
              
                if (nbrNbrNbr->IsOxygen() || (nbrNbrNbr->GetValence() == 1)) {
                  return 28; // Hydrogen on NSO, NSO2 or NSO3 nitrogen (HNSO)
                }
              }
            }
          }
              
          return 23; // Generic hydrogen on sp3 nitrogen e.g., in amines,
          // Hydrogen on nitrogen in pyrrole, Hydrogen in ammonia,
          // Hydrogen on N in N-oxide (HNR, HPYL, H3N, HNOX)
        }
        if (nbr->IsSulfur() || nbr->IsPhosphorus()) {
          return 71; // Hydrogen attached to sulfur, Hydrogen attached to >S= sulfur doubly bonded to N,
          // Hydrogen attached to phosphorus (HS, HS=N, HP)
        }
      }
    }

    ////////////////////////////////
    // Lithium
    ////////////////////////////////
    if (atom->GetAtomicNum() == 3) {
      // 0 neighbours
      if (atom->GetValence() == 0) {
        return 92; // Lithium cation (LI+)
      }
    }
 
    ////////////////////////////////
    // Carbon
    ////////////////////////////////
    if (atom->GetAtomicNum() == 6) {
      // 4 neighbours
      if (atom->GetValence() == 4) {
        if (atom->IsInRingSize(3)) {
          return 22; // Aliphatic carbon in 3-membered ring (CR3R)
        } 
	
        if (atom->IsInRingSize(4)) {
          return 20; // Aliphatic carbon in 4-membered ring (CR4R)
        }
        
        return 1; // Alkyl carbon (CR)
      }
      // 3 neighbours
      if (atom->GetValence() == 3) {
        int N2count = 0;
        int N3count = 0;
        oxygenCount = sulphurCount = doubleBondTo = 0;

        FOR_NBORS_OF_ATOM (nbr, atom) {
          bond = mol->GetBond(&*nbr, atom);
          if (bond->IsDouble()) {
            doubleBondTo = nbr->GetAtomicNum();
          }

          if (nbr->GetValence() == 1) {
            if (nbr->IsOxygen()) {
              oxygenCount++;
            } else if (nbr->IsSulfur()) {
              sulphurCount++;
            }
          } else if (nbr->GetValence() == 3) {
            if (nbr->IsNitrogen()) {
              N3count++;
            }
          } else if ((nbr->GetValence() == 2) && bond->IsDouble()) {
            if (nbr->IsNitrogen()) {
              N2count++;
            }
          }
        }
        if ((N3count >= 2) && (doubleBondTo == 7) && !N2count) {
          // N3==C--N3
          return 57; // Guanidinium carbon, Carbon in +N=C-N: resonance structures (CGD+, CNN+)
        }
        if ((oxygenCount == 2) || (sulphurCount == 2)) {
          // O1-?-C-?-O1 or S1-?-C-?-S1
          return 41; // Carbon in carboxylate anion, Carbon in thiocarboxylate anion (CO2M, CS2M)
        }
        if (atom->IsInRingSize(4) && (doubleBondTo == 6)) {
	        return 30; // Olefinic carbon in 4-membered ring (CR4E)
        }
        if ((doubleBondTo ==  7) || (doubleBondTo ==  8) || 
            (doubleBondTo == 15) || (doubleBondTo == 16)) {
          // C==N, C==O, C==P, C==S
          return 3; // Generic carbonyl carbon, Imine-type carbon, Guanidine carbon,
          // Ketone or aldehyde carbonyl carbon, Amide carbonyl carbon,
          // Carboxylic acid or ester carbonyl carbon, Carbamate carbonyl carbon,
          // Carbonic acid or ester carbonyl carbon, Thioester carbonyl (double
          // bonded to O or S), Thioamide carbon (double bonded to S), Carbon
          // in >C=SO2, Sulfinyl carbon in >C=S=O, Thiocarboxylic acid or ester 
          // carbon, Carbon doubly bonded to P (C=O, C=N, CGD, C=OR, C=ON, COO,
          // COON, COOO, C=OS, C=S, C=SN, CSO2, CS=O, CSS, C=P)
        }
	
        return 2; // Vinylic Carbon, Generic sp2 carbon (C=C, CSP2)
	
      }
      // 2 neighbours
      if (atom->GetValence() == 2) {
        return 4; // Acetylenic carbon, Allenic caron (CSP, =C=)
      }
      // 1 neighbours
      if (atom->GetValence() == 1) {
        return 60; // Isonitrile carbon (C%-)
      }
    }

    ////////////////////////////////
    // Nitrogen
    ////////////////////////////////
    if (atom->GetAtomicNum() == 7) {
      // 4 neighbours
      if (atom->GetValence() == 4) {
        FOR_NBORS_OF_ATOM (nbr, atom) {
          if (nbr->IsOxygen() && (nbr->GetValence() == 1)) {
            return 68; // sp3-hybridized N-oxide nitrogen (N3OX)
          }
        }

        return 34; // Quaternary nitrogen (NR+)
      }
      // 3 neighbours
      if (atom->GetValence() == 3) {
        if (atom->BOSum() >= 4) { // -N(=O)(=O) is a valid nitro group
          oxygenCount = nitrogenCount = doubleBondTo = 0;
	  
          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (nbr->IsOxygen() && (nbr->GetValence() == 1)) {
              oxygenCount++;
            }
            if (nbr->IsNitrogen()) {
              bond = mol->GetBond(&*nbr, atom);
              if (bond->IsDouble()) {
                doubleBondTo = 7;
              }
            }
            if (nbr->IsCarbon()) {
              bond = mol->GetBond(&*nbr, atom);
              if (bond->IsDouble()) {
                FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
                  if (nbrNbr->IsNitrogen() && (nbrNbr->GetValence() == 3)) {
                    nitrogenCount++;
                  }
                }
              }
            }
          }

          if (oxygenCount == 1) {
            return 67; // sp2-hybridized N-oxide nitrogen (N2OX)
          }
          if (oxygenCount >= 2) {
            return 45; // Nitrogen in nitro group, Nitrogen in nitrate group (NO2, NO3)
          }

          if (nitrogenCount == 1) {
            return 54; // Iminium nitrogen (N+=C)
          }
          if (nitrogenCount == 2) {
            return 55; // Either nitrogen in N+=C-N: (NCN+)
          }
          if (nitrogenCount == 3) {
            return 56; // Guanidinium nitrogen (NGD+)
          }
	  
          if (doubleBondTo == 7) {
            return 54; // Positivly charged nitrogen doubly bonded to nitrogen (N+=N)
          }
        }
	
        if (atom->BOSum() == 3) {
          bool IsAmide = false;
          bool IsSulfonAmide = false;
          bool IsNNNorNNC = false;
          int tripleBondTo = 0;
          doubleBondTo = 0;
	  
          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (nbr->IsSulfur() || nbr->IsPhosphorus()) {
              oxygenCount = 0;
	      
              FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
                if (nbrNbr->IsOxygen() && (nbrNbr->GetValence() == 1)) {
                  oxygenCount++;
                }
              }
              if (oxygenCount >= 2) {
                IsSulfonAmide = true;
                //return 43; // Sulfonamide nitrogen (NSO2, NSO3)
              }
            }
          }
	
          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (nbr->IsCarbon()) {
              FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
                bond = mol->GetBond(&*nbr, &*nbrNbr);
                if (bond->IsDouble() && (nbrNbr->IsOxygen() || nbrNbr->IsSulfur())) {
                  IsAmide = true;
                  //return 10; // Amide nitrogen, Thioamide nitrogen (NC=O, NC=S)
                }
              }
            }
          }
	  
          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (nbr->IsCarbon()) {
              int N2count = 0;
              int N3count = 0;
              oxygenCount = sulphurCount = 0;

              FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
                bond = mol->GetBond(&*nbr, &*nbrNbr);
                if (bond->IsDouble()) {
                  doubleBondTo = nbrNbr->GetAtomicNum();
                }
                if (IsAromatic(bond)) {
                  if ((nbrNbr->GetAtomicNum() == 7) || (nbrNbr->GetAtomicNum() == 6)) {
                    doubleBondTo = nbrNbr->GetAtomicNum();
                  }
                }
                if (bond->IsTriple()) {
                  tripleBondTo = nbrNbr->GetAtomicNum();
                }
                if (nbrNbr->IsNitrogen() && (nbrNbr->GetValence() == 3)) {
                  int nbrOxygen = 0;
                  FOR_NBORS_OF_ATOM (nbrNbrNbr, &*nbrNbr) {
                    if (nbrNbrNbr->IsOxygen()) {
                      nbrOxygen++;
                    }
                  }
                  if (nbrOxygen < 2) {
                    N3count++;
                  }
                }
                if (nbrNbr->IsNitrogen() && (nbrNbr->GetValence() == 2) && (bond->IsDouble() || IsAromatic(bond))) {
                  N2count++;
                }
                if (IsAromatic(&*nbrNbr)) {
                  if (nbrNbr->IsOxygen()) {
                    oxygenCount++;
                  }
                  if (nbrNbr->IsSulfur()) {
                    sulphurCount++;
                  }
                }
              }
              if (N3count == 3) {
                return 56; // Guanidinium nitrogen (NGD+)
              }
	
              if (!IsAmide && !IsSulfonAmide && !oxygenCount && !sulphurCount && IsAromatic(&*nbr)) {
                return 40;
              }

              if ((N3count == 2) && (doubleBondTo == 7) && !N2count) {
                return 55; // Either nitrogen in N+=C-N: (NCN+)
              }
            }

            if (nbr->IsNitrogen()) {
              nitrogenCount = 0;
              FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
                bond = mol->GetBond(&*nbr, &*nbrNbr);
                if (bond->IsDouble()) {
                  if (nbrNbr->IsCarbon()) {
                    oxygenCount = sulphurCount = 0;
                    FOR_NBORS_OF_ATOM (nbrNbrNbr, &*nbrNbr) {
                      if (nbrNbrNbr->IsOxygen()) {
                        oxygenCount++;
                      }
                      if (nbrNbrNbr->IsSulfur()) {
                        sulphurCount++;
                      }
                      if (nbrNbrNbr->IsSulfur()) {
                        nitrogenCount++;
                      }
                    }
                    if (!oxygenCount && !sulphurCount && (nitrogenCount == 1)) {
                      bool bondToAromC = false;
                      FOR_NBORS_OF_ATOM (nbr2, atom) {
                        if (IsAromatic(&*nbr2) && nbr2->IsCarbon() && nbr2->IsInRingSize(6)) {
                          bondToAromC = true;
                        }
                      }
                      if (!bondToAromC) {
                        IsNNNorNNC = true;
                      }
                    }
                  }
                  if (nbrNbr->IsNitrogen()) {
                    bool bondToAromC = false;
                    FOR_NBORS_OF_ATOM (nbr2, atom) {
                      if (IsAromatic(&*nbr2) && nbr2->IsCarbon() && nbr2->IsInRingSize(6)) {
                        bondToAromC = true;
                      }
                    }
                    if (!bondToAromC) {
                      IsNNNorNNC = true;
                    }
                  }
                }
              }
            }	    
          }
          
          if (IsSulfonAmide) {
            return 43; // Sulfonamide nitrogen (NSO2, NSO3)
          }
          if (IsAmide) {
            return 10; // Amide nitrogen, Thioamide nitrogen (NC=O, NC=S)
          }
 
          if ((doubleBondTo ==  6) || (doubleBondTo == 7) ||(doubleBondTo == 15) || (tripleBondTo == 6)) {
            return 40; // Enamine or aniline nitrogen (deloc. lp), Nitrogen in N-C=N with deloc. lp,
            // Nitrogen in N-C=N with deloc. lp, Nitrogen attached to C-C triple bond
            // (NC=C, NC=N, NC=P, NC%C)
          }
          if (tripleBondTo == 7) {
            return 43; // Nitrogen attached to cyano group (NC%N)
          }
          if (IsNNNorNNC) {
            return 10; // Nitrogen in N-N=C moiety with deloc. lp
            // Nitrogen in N-N=N moiety with deloc. lp (NN=C, NN=N)
          }
	
          return 8; // Amine nitrogen (NR)
        }
      }
      // 2 neighbours
      if (atom->GetValence() == 2) {
        if (atom->BOSum() == 4) {
          FOR_NBORS_OF_ATOM (nbr, atom) {
            bond = mol->GetBond(&*nbr, atom);
            if (bond->IsTriple()) {
              return 61; // Isonitrile nitrogen (NR%)
            }
          }

          return 53; // Central nitrogen in C=N=N or N=N=N (=N=)
        } 
	
        if (atom->BOSum() == 3) {
          doubleBondTo = 0;

          FOR_NBORS_OF_ATOM (nbr, atom) {
            bond = mol->GetBond(&*nbr, atom);
            if (nbr->IsOxygen() && bond->IsDouble() && (nbr->GetValence() == 1)) {
              return 46; // Nitrogen in nitroso group (N=O)
            }
            if ((nbr->IsCarbon() || nbr->IsNitrogen()) && bond->IsDouble()) {
              return 9; // Iminie nitrogen, Azo-group nitrogen (N=C, N=N)
            }
          }
          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (nbr->IsSulfur()) {
              oxygenCount = 0;
	      
              FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
                if (nbrNbr->IsOxygen() && (nbrNbr->GetValence() == 1)) {
                  oxygenCount++;
                }
              }
              if (oxygenCount >= 2) {
                return 43; // Sulfonamide nitrogen (NSO2, NSO3)
              }
            }
          }	
        } 
	
        if (atom->BOSum() == 2) {
          oxygenCount = sulphurCount = 0;

          FOR_NBORS_OF_ATOM (nbr, atom) {
            if (nbr->IsSulfur()) {
              FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
                if (nbrNbr->IsOxygen() && (nbrNbr->GetValence() == 1)) {
                  oxygenCount++;
                }
              }
              if (oxygenCount == 1) {
                return 48; // Divalent nitrogen replacing monovalent O in SO2 group (NSO)
              }
            }
          }

          return 62; // Anionic divalent nitrogen (NM)
        } 
      }
      // 1 neighbours
      if (atom->GetValence() == 1) {
       	FOR_NBORS_OF_ATOM (nbr, atom) {
          bond = mol->GetBond(&*nbr, atom);
          if (bond->IsTriple()) {
            return 42; // Triply bonded nitrogen (NSP)
          }
          if (nbr->IsNitrogen() && (nbr->GetValence() == 2)) {
            return 47; // Terminal nitrogen in azido or diazo group (NAZT)
          }
        }
      }
    }

    ////////////////////////////////
    // Oxygen
    ////////////////////////////////
    if (atom->GetAtomicNum() == 8) {
      // 3 neighbours
      if (atom->GetValence() == 3) {
        return 49; // Oxonium oxygen (O+)
      }
      // 2 neighbours
      if (atom->GetValence() == 2) {
        int hydrogenCount = 0;
        FOR_NBORS_OF_ATOM (nbr, atom) {
          if (nbr->IsHydrogen()) {
            hydrogenCount++;
          }
        }

        if (hydrogenCount == 2) {
          // H--O--H
          return 70; // Oxygen in water (OH2)
        }
        if (atom->BOSum() == 3) {
          return 51; // Oxenium oxygen (O=+)
        }
        
        return 6; // Generic divalent oxygen, Ether oxygen, Carboxylic acid or ester oxygen,
        // Enolic or phenolic oxygen, Oxygen in -O-C=N- moiety, Divalent oxygen in
        // thioacid or ester, Divalent nitrate "ether" oxygen, Divalent oxygen in
        // sulfate group, Divalent oxygen in sulfite group, One of two divalent
        // oxygens attached to sulfur, Divalent oxygen in R(RO)S=O, Other divalent
        // oxygen attached to sulfur, Divalent oxygen in phosphate group, Divalent
        // oxygen in phosphite group, Divalent oxygen (one of two oxygens attached
        // to P), Other divalent oxygen (-O-, OR, OC=O, OC=C, OC=N, OC=S, ONO2, 
        // ON=O, OSO3, OSO2, OSO, OS=O, -OS, OPO3, OPO2, OPO, -OP)

        // 59 ar
      }
      // 1 neighbour
      if (atom->GetValence() == 1) {
        oxygenCount = sulphurCount = 0;
        
        FOR_NBORS_OF_ATOM (nbr, atom) {
          bond = mol->GetBond(&*nbr, atom);

          if (nbr->IsCarbon() || nbr->IsNitrogen()) {
            FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
              if (nbrNbr->IsOxygen() && (nbrNbr->GetValence() == 1)) {
                oxygenCount++;
              }
              if (nbrNbr->IsSulfur() && (nbrNbr->GetValence() == 1)) {
                sulphurCount++;
              }
            }
          }
          // O---H
          if (nbr->IsHydrogen()) {
            return 35;
          }
          // O-?-C
          if (nbr->IsCarbon()) {
            if (oxygenCount == 2) {
              // O-?-C-?-O
              return 32; // Oxygen in carboxylate group (O2CM)
            }
            if (bond->IsSingle()) { 
              // O--C
              return 35; // Oxide oxygen on sp3 carbon, Oxide oxygen on sp2 carbon (OM, OM2)
            } else { 
              // O==C
              return 7; // Generic carbonyl oxygen, Carbonyl oxygen in amides,
              // Carbonyl oxygen in aldehydes and ketones, Carbonyl
              // oxygen in acids or esters (O=C, O=CN, O=CR, O=CO)
            }
          }
          // O-?-N
          if (nbr->IsNitrogen()) {
            if (oxygenCount >= 2) { 
              // O-?-N-?-O
              return 32; // Oxygen in nitro group, Nitro-group oxygen in nitrate,
              // Nitrate anion oxygen (O2N, O2NO, O3N)
            }
            if (bond->IsSingle()) { 
              // O--N
              return 32; // Oxygen in N-oxides (ONX)
            } else { 
              // O==N
              return 7; // Nitroso oxygen (O=N)
            }
          }
          // O-?-S
          if (nbr->IsSulfur()) {
            if (sulphurCount == 1) { 
              // O1-?-S-?-S1
              return 32; // Terminal oxygen in thiosulfinate anion (OSMS)
            }
            if (bond->IsSingle()) { 
              // O--S
              return 32; // Single terminal oxygen on sulfur, One of 2 terminal O's on sulfur, 
              // One of 3 terminal O's on sulfur, Terminal O in sulfate anion, 
              // (O-S, O2S, O3S, O4S)
            } else { 
              // O==S

              // are all sulfur nbr atoms carbon?
              bool isSulfoxide = true;
              FOR_NBORS_OF_ATOM (nbr2, &*nbr) {
                if (atom == &*nbr2)
                  continue;

                if (nbr->GetBond(&*nbr2)->IsDouble() && !nbr2->IsOxygen())
                  return 7; // O=S on sulfur doubly bonded to, e.g., C (O=S=)

                if (nbr2->IsOxygen() && nbr2->GetValence() == 1)
                  isSulfoxide = false;
              }

              if (isSulfoxide)
                return 7; // Doubly bonded sulfoxide oxygen (O=S)
              else
                return 32; // (O2S, O3S, O4S)
            }
          }

          return 32; // Oxygen in phosphine oxide, One of 2 terminal O's on sulfur, 
          // One of 3 terminal O's on sulfur, One of 4 terminal O's on sulfur, 
          // Oxygen in perchlorate anion (OP, O2P, O3P, O4P, O4Cl)
        }
      }
    }
    
    ////////////////////////////////
    // Flourine
    ////////////////////////////////
    if (atom->GetAtomicNum() == 9) {
      // 1 neighbour
      if (atom->GetValence() == 1) {
        return 11; // Fluorine (F)
      }
      // 0 neighbours
      if (atom->GetValence() == 0) {
        return 89; // Fluoride anion (F-)
      }
    }
    
    ////////////////////////////////
    // Sodium
    ////////////////////////////////
    if (atom->GetAtomicNum() == 11) {
      return 93; // Sodium cation (NA+)
    }
    
    ////////////////////////////////
    // Magnesium
    ////////////////////////////////
    if (atom->GetAtomicNum() == 12) {
      return 99; // Dipositive magnesium cation (MG+2)
    }
 
    ////////////////////////////////
    // Silicon
    ////////////////////////////////
    if (atom->GetAtomicNum() == 14) {
      return 19; // Silicon (SI)
    }
 
    ////////////////////////////////
    // Phosphorus
    ////////////////////////////////
    if (atom->GetAtomicNum() == 15) {
      if (atom->GetValence() == 4) {
        return 25; // Phosphate group phosphorus, Phosphorus with 3 attached oxygens,
        // Phosphorus with 2 attached oxygens, Phosphine oxide phosphorus,
        // General tetracoordinate phosphorus (PO4, PO3, PO2, PO, PTET)
      }
      if (atom->GetValence() == 3) {
        return 26; // Phosphorus in phosphines (P)
      }
      if (atom->GetValence() == 2) {
        return 75; // Phosphorus doubly bonded to C (-P=C)
      }
    }
    
    ////////////////////////////////
    // Sulfur
    ////////////////////////////////
    if (atom->GetAtomicNum() == 16) {
      // 4 neighbours
      if (atom->GetValence() == 4) {
        return 18; // Sulfone sulfur, Sulfonamide sulfur, Sulfonate group sulfur,
        // Sulfate group sulfur, Sulfur in nitrogen analog of sulfone 
        // (SO2, SO2N, SO3, SO4, SNO)
      }
      // 3 neighbours
      if (atom->GetValence() == 3) {
        oxygenCount = sulphurCount = doubleBondTo = 0;

        FOR_NBORS_OF_ATOM (nbr, atom) {
          bond = mol->GetBond(&*nbr, atom);
          if (bond->IsDouble()) {
            doubleBondTo = nbr->GetAtomicNum();
          }

          if (nbr->GetValence() == 1) {
            if (nbr->IsOxygen()) {
              oxygenCount++;
            } else if (nbr->IsSulfur()) {
              sulphurCount++;
            }
          } 
        }

        if (oxygenCount == 2) {
          if (doubleBondTo == 6) {
            return 18; // Sulfone sulfur, doubly bonded to carbon (=SO2)
          }
          return 73; // Sulfur in anionic sulfinate group (SO2M)
        }
        if (oxygenCount && sulphurCount)
          return 73; // Tricoordinate sulfur in anionic thiosulfinate group (SSOM)
  
        //if ((doubleBondTo == 6) || (doubleBondTo == 8))
        return 17; // Sulfur doubly bonded to carbon, Sulfoxide sulfur (S=C, S=O)
      }
      // 2 neighbours
      if (atom->GetValence() == 2) {
        doubleBondTo = 0;

        FOR_NBORS_OF_ATOM (nbr, atom) {
          if (nbr->IsOxygen()) {
            bond = mol->GetBond(&*nbr, atom);
            if (bond->IsDouble()) {
              doubleBondTo = 8;
            }
          }
        }

        if (doubleBondTo == 8)
          return 74; // Sulfinyl sulfur, e.g., in C=S=O (=S=O)
	
        return 15; // Thiol, sulfide, or disulfide sulfor (S)
      }
      // 1 neighbour
      if (atom->GetValence() == 1) {
        sulphurCount = doubleBondTo = 0;

        FOR_NBORS_OF_ATOM (nbr, atom) {
          FOR_NBORS_OF_ATOM (nbrNbr, &*nbr) {
            if (nbrNbr->IsSulfur() && (nbrNbr->GetValence() == 1)) {
              sulphurCount++;
            }
          }
          bond = mol->GetBond(&*nbr, atom);
          if (bond->IsDouble()) {
            doubleBondTo = nbr->GetAtomicNum();
          }
        }

        if ((doubleBondTo == 6) && (sulphurCount != 2)) {
          return 16; // Sulfur doubly bonded to carbon (S=C)
        }

        return 72; // Terminal sulfur bonded to P, Anionic terminal sulfur,
        // Terminal sulfur in thiosulfinate group (S-P, SM, SSMO)
      }

      // 44 ar
    }
    
    ////////////////////////////////
    // Clorine
    ////////////////////////////////
    if (atom->GetAtomicNum() == 17) {
      // 4 neighbour
      if (atom->GetValence() == 4) {
        oxygenCount = 0;
        
        FOR_NBORS_OF_ATOM (nbr, atom) {
          if (nbr->IsOxygen()) {
            oxygenCount++;
          }
        }
        if (oxygenCount == 4)
          return 77; // Perchlorate anion chlorine (CLO4)
      }
      // 1 neighbour
      if (atom->GetValence() == 1) {
        return 12; // Chlorine (CL)
      }
      // 0 neighbours
      if (atom->GetValence() == 0) {
        return 90; // Chloride anion (CL-)
      }
    }
    
    ////////////////////////////////
    // Potasium
    ////////////////////////////////
    if (atom->GetAtomicNum() == 19) {
      return 94; // Potasium cation (K+)
    }
    
    ////////////////////////////////
    // Calcium
    ////////////////////////////////
    if (atom->GetAtomicNum() == 20) {
      // 0 neighbours
      if (atom->GetValence() == 0) {
        return 96; // Dipositive calcium cation (CA+2)
      }
    }
 
    ////////////////////////////////
    // Iron
    ////////////////////////////////
    if (atom->GetAtomicNum() == 26) {
      return 87; // Dipositive iron (FE+2)
      return 88; // Tripositive iron (FE+3)
    }
    
    ////////////////////////////////
    // Copper
    ////////////////////////////////
    if (atom->GetAtomicNum() == 29) {
      return 97; // Monopositive copper cation (CU+1)
      return 98; // Dipositive copper cation (CU+2)
    }
    
    ////////////////////////////////
    // Zinc
    ////////////////////////////////
    if (atom->GetAtomicNum() == 30) {
      return 95; // Dipositive zinc cation (ZN+2)
    }
 
    ////////////////////////////////
    // Bromine
    ////////////////////////////////
    if (atom->GetAtomicNum() == 35) {
      // 1 neighbour
      if (atom->GetValence() == 1) {
        return 13; // Bromine (BR)
      }
      // 0 neighbours
      if (atom->GetValence() == 0) {
        return 91; // Bromide anion (BR-)
      }
    }
 
    ////////////////////////////////
    // Iodine
    ////////////////////////////////
    if (atom->GetAtomicNum() == 53) {
      // 1 neighbour
      if (atom->GetValence() == 1) {
        return 14; // Iodine (I)
      }
    }
 


    return 0;
  }

  bool MMFF94Common::SetTypes(/*const*/ OBMol &mol)
  {
    // mark all atoms and bonds as non-aromatic 
    m_aromAtoms.resize(mol.NumAtoms(), 0);
    m_aromBonds.resize(mol.NumBonds(), 0);

    // It might be needed to run this function more than once...
    bool done = true;
    while (done) {
      done = PerceiveAromaticity(mol);
    }
   
    OBAtomIterator iter;
    for (OBAtom *atom = mol.BeginAtom(iter); atom; atom = mol.NextAtom(iter))
      m_types.push_back(GetType(atom));
    
    return true;
  }
  
  bool MMFF94Common::SetFormalCharges(/*const*/ OBMol &mol)
  {
    m_fCharges.resize(mol.NumAtoms(), 0.0);

    FOR_ATOMS_OF_MOL (atom, mol) {
      int type = m_types.at(atom->GetIdx()-1);
       
      bool done = false;
      switch (type) {
      case 34:
      case 49:
      case 51:
      case 54:
      case 58:
      case 92:
      case 93:
      case 94:
      case 97:
        m_fCharges[atom->GetIdx()-1] = 1.0;
        done = true;
        break;
      case 35:
      case 62:
      case 89:
      case 90:
      case 91:
        m_fCharges[atom->GetIdx()-1] = -1.0;
        done = true;
        break;
      case 55:
        m_fCharges[atom->GetIdx()-1] = 0.5;
        done = true;
        break;
      case 56:
        m_fCharges[atom->GetIdx()-1] = 1.0 / 3.0;
        done = true;
        break;
      case 87:
      case 95:
      case 96:
      case 98:
      case 99:
        m_fCharges[atom->GetIdx()-1] = 2.0;
        done = true;
        break;
        //case 98:
        //  atom->SetPartialCharge(3.0);
      default:
        break;
      }

      if (done)
        continue;
    
      if (type == 32) {
        int o_count = 0;
        bool sulfonamide = false;
        int s_count = 0;

        FOR_NBORS_OF_ATOM(nbr, &*atom) {
          FOR_NBORS_OF_ATOM(nbr2, &*nbr) {
            if (nbr2->IsOxygen() && (nbr2->GetValence() == 1))
              o_count++;
            if (nbr2->IsSulfur() && (nbr2->GetValence() == 1))
              s_count++;
            if (nbr2->IsNitrogen() && !IsAromatic(&*nbr2))
              sulfonamide = true;
          }
	
          if (nbr->IsCarbon())
            m_fCharges[atom->GetIdx()-1] = -0.5; // O2CM
	  
          if (nbr->IsNitrogen() && (o_count == 3))
            m_fCharges[atom->GetIdx()-1] = -1.0 / o_count; // O3N
	  
          if (nbr->IsSulfur() && !sulfonamide)
            if (((o_count + s_count) == 2) && (nbr->GetValence() == 3) && (nbr->BOSum() == 3))
              m_fCharges[atom->GetIdx()-1] = -0.5; // O2S
            else if ((o_count + s_count) == 3)
              m_fCharges[atom->GetIdx()-1] = -1.0 / 3.0; // O3S
            else if ((o_count + s_count) == 4)
              m_fCharges[atom->GetIdx()-1] = -0.5; // O4S
	  
          if (nbr->IsPhosphorus())
            if ((o_count + s_count) == 2)
              m_fCharges[atom->GetIdx()-1] = -0.5; // O2P
            else if ((o_count + s_count) == 3)
              m_fCharges[atom->GetIdx()-1] = -2.0 / 3.0; // O3P
            else if ((o_count + s_count) == 4)
              m_fCharges[atom->GetIdx()-1] = -0.25; // O4P
	  
          if (type == 77)
            m_fCharges[atom->GetIdx()-1] = -0.25; // O4CL
        }
      } else if (type == 61) {
        FOR_NBORS_OF_ATOM(nbr, &*atom)
          if (atom->GetBond(&*nbr)->IsTriple() && nbr->IsNitrogen())
            m_fCharges[atom->GetIdx()-1] = 1.0;
      } else if (type == 72) {
        int s_count = 0;

        FOR_NBORS_OF_ATOM(nbr, &*atom) {
          if (nbr->IsSulfur())
            s_count++;
	  
          if (nbr->IsPhosphorus() || nbr->IsSulfur()) {
            FOR_NBORS_OF_ATOM(nbr2, &*nbr)
              if ((nbr2->IsSulfur() || nbr2->IsOxygen()) && (nbr2->GetValence() == 1) && (atom->GetIdx() != nbr2->GetIdx()))
                m_fCharges[atom->GetIdx()-1] = -0.5;
          } else
            m_fCharges[atom->GetIdx()-1] = -1.0;

          if (nbr->IsCarbon())
            FOR_NBORS_OF_ATOM(nbr2, &*nbr)
              if (nbr2->IsSulfur() && (nbr2->GetValence() == 1) && (atom->GetIdx() != nbr2->GetIdx()))
                m_fCharges[atom->GetIdx()-1] = -0.5; // SSMO
	
          if (s_count >= 2)      
            m_fCharges[atom->GetIdx()-1] = -0.5; // SSMO
        }
      } else if (type == 76) {
       	vector<OBRing*> vr;
        vr = mol.GetSSSR();
        vector<OBRing*>::iterator ri;
        vector<int>::iterator rj;
        int n_count;

        for (ri = vr.begin();ri != vr.end();ri++) { // for each ring
          n_count = 0;

          if (IsAromatic(*ri) && (*ri)->IsMember(&*atom) && ((*ri)->Size() == 5)) {
            for(rj = (*ri)->_path.begin();rj != (*ri)->_path.end();rj++) // for each ring atom
              if (mol.GetAtom(*rj)->IsNitrogen())
                n_count++;
	    
            if (n_count > 1)
              m_fCharges[atom->GetIdx()-1] = -1.0 / n_count;
          }
        }
      } else if (type == 81) {
        m_fCharges[atom->GetIdx()-1] = 1.0;
        
        vector<OBRing*> vr;
        vr = mol.GetSSSR();
        vector<OBRing*>::iterator ri;
        vector<int>::iterator rj;
        for (ri = vr.begin();ri != vr.end();ri++) // for each ring
          if (IsAromatic(*ri) && (*ri)->IsMember(&*atom) && ((*ri)->Size() == 5)) {
            int n_count = 0;
            for(rj = (*ri)->_path.begin();rj != (*ri)->_path.end();rj++) // for each ring atom
              if (mol.GetAtom(*rj)->IsNitrogen() && (mol.GetAtom(*rj)->GetValence() == 3))
                n_count++;

            m_fCharges[atom->GetIdx()-1] = 1.0 / n_count; // NIM+
	    
            FOR_NBORS_OF_ATOM(nbr, &*atom)
              FOR_NBORS_OF_ATOM(nbr2, &*nbr)
              if (GetCachedType(&*nbr2) == 56)
                m_fCharges[atom->GetIdx()-1] = 1.0 / 3.0;
  	  
            FOR_NBORS_OF_ATOM(nbr, &*atom)
              FOR_NBORS_OF_ATOM(nbr2, &*nbr)
              if (GetCachedType(&*nbr2) == 55)
                m_fCharges[atom->GetIdx()-1] = 1.0 / (1.0 + n_count);
          }
      } 
    
    }

    return true;
  }

  bool MMFF94Common::SetPartialCharges(/*const*/ OBMol &mol)
  {
    m_pCharges.resize(mol.NumAtoms(), 0.0);
    double M, Wab, factor, q0a, q0b, Pa, Pb;

    FOR_ATOMS_OF_MOL (atom, mol) {
      int type = GetCachedType(&*atom);

      switch (type) {
      case 32:
      case 35:
      case 72:
        factor = 0.5;
        break;
      case 62:
      case 76:
        factor = 0.25;
        break;
      default:
        factor = 0.0;
        break;
      }
      
      M = GetCrd(type);
      q0a = m_fCharges.at(atom->GetIdx()-1);

      // charge sharing
      if (!factor)
        FOR_NBORS_OF_ATOM (nbr, &*atom)
          if (m_fCharges.at(nbr->GetIdx()-1) < 0.0)
            q0a += m_fCharges.at(nbr->GetIdx()-1) / (2.0 * nbr->GetValence());
      
      // needed for SEYWUO, positive charge sharing?
      if (type == 62)
        FOR_NBORS_OF_ATOM (nbr, &*atom)
          if (m_fCharges.at(nbr->GetIdx()-1) > 0.0)
            q0a -= m_fCharges.at(nbr->GetIdx()-1) / 2.0;
     
      q0b = 0.0;
      Wab = 0.0;
      Pa = Pb = 0.0;
      FOR_NBORS_OF_ATOM (nbr, &*atom) {
        int nbr_type = GetCachedType(&*nbr);

        q0b += m_fCharges.at(nbr->GetIdx()-1);
    
        int bondType = GetBondType(&*atom, &*nbr);
        bool bci_found = false;
        std::vector<OBParameterDB::Query> query;
        query.push_back( OBParameterDB::Query(0, OBVariant(bondType)) );
        query.push_back( OBParameterDB::Query(1, OBVariant(type), true) );
        query.push_back( OBParameterDB::Query(2, OBVariant(nbr_type), true) );
        bool swapped;
        std::vector<OBVariant> row = m_database->FindRow(MMFF94SimpleParameterDB::ChargeParameters, query, &swapped);

        if (row.size()) {
          if (swapped)
            Wab += row.at(3).AsDouble();
          else
            Wab -= row.at(3).AsDouble();
          bci_found = true;
        }

        if (!bci_found) {
          std::vector<OBParameterDB::Query> query_a;
          query.push_back( OBParameterDB::Query(0, OBVariant(type)) );
          std::vector<OBVariant> row_a = m_database->FindRow(MMFF94SimpleParameterDB::PartialBondChargeIncrements, query);
          if (row_a.size())
            Pa = row_a.at(1).AsDouble();

          std::vector<OBParameterDB::Query> query_b;
          query.push_back( OBParameterDB::Query(0, OBVariant(nbr_type)) );
          std::vector<OBVariant> row_b = m_database->FindRow(MMFF94SimpleParameterDB::PartialBondChargeIncrements, query);
          if (row_b.size())
            Pb = row_b.at(1).AsDouble();

          Wab += Pa - Pb;
        }
      }
      
      if (factor)
        m_pCharges[atom->GetIdx()-1] = (1.0 - M * factor) * q0a + factor * q0b + Wab;
      else 
        m_pCharges[atom->GetIdx()-1] = q0a + Wab;

    }
    
    return true;
  }

  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Calculate bond type, angle type, stretch-bend type, torsion type
  //
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
 
  //
  // MMFF part V - page 620
  //
  // BTij is 1 when:
  // a) single bond between atoms i and j, both i and j are not aromatic and both types have sbmb set in mmffprop.par, or
  // b) bewtween two aromatic atoms, but the bond is not aromatic (e.g. connecting bond in biphenyl)
  //
  int MMFF94Common::GetBondType(OBAtom* a, OBAtom* b)
  {
    OBMol *mol = a->GetParent();
    
    if (!mol->GetBond(a,b)->IsSingle())
      return 0;
    
    if (!IsAromatic(mol->GetBond(a,b)))
      if (HasAromSet(GetCachedType(a)) && HasAromSet(GetCachedType(b)))
        return 1;
      
    if (HasSbmbSet(GetCachedType(a)) && HasSbmbSet(GetCachedType(b)))
      return 1;
    
    return 0;
  }
  
  int MMFF94Common::GetAngleType(OBAtom* a, OBAtom* b, OBAtom *c)
  {
    int sumbondtypes;

    sumbondtypes = GetBondType(a,b) + GetBondType(b, c);

    if (a->IsInRingSize(3) && b->IsInRingSize(3) && c->IsInRingSize(3) && IsInSameRing(a, c))
      switch (sumbondtypes) {
      case 0:
        return 3; 
      case 1:
        return 5; 
      case 2:
        return 6; 
      }
    
    if (a->IsInRingSize(4) && b->IsInRingSize(4) && c->IsInRingSize(4) && IsInSameRing(a, c))
      switch (sumbondtypes) {
      case 0:
        return 4; 
      case 1:
        return 7; 
      case 2:
        return 8; 
      }
    
    return sumbondtypes;
  }
  
  int MMFF94Common::GetStrBndType(OBAtom* a, OBAtom* b, OBAtom *c)
  {
    int btab, btbc, atabc;
    bool inverse;

    btab = GetBondType(a, b);
    btbc = GetBondType(b, c);
    atabc = GetAngleType(a, b, c);

    if (GetCachedType(a) <= GetCachedType(c))
      inverse = false;
    else
      inverse = true;

    switch (atabc) {
    case 0:
      return 0;

    case 1:
      if (btab)
        if (!inverse)
          return 1;
        else
          return 2;
      if (btbc)
        if (!inverse)
          return 2;
        else
          return 1;

    case 2:
      return 3;

    case 3:
      return 5;

    case 4:
      return 4;

    case 5:
      if (btab)
        if (!inverse)
          return 6;
        else
          return 7;
      if (btbc)
        if (!inverse)
          return 7;
        else
          return 6;
      
    case 6:
      return 8;
      
    case 7:
      if (btab)
        if (!inverse)
          return 9;
        else
          return 10;
      if (btbc)
        if (!inverse)
          return 10;
        else
          return 9;
      
    case 8:
      return 11;
    }
    
    return 0;
  }
  
  //
  // MMFF part IV - page 609
  //
  // TTijkl = 1 when BTjk = 1
  // TTijkl = 2 when BTjk = 0 but BTij and/or BTkl = 1
  // TTijkl = 4 when i, j, k and l are all members of the same four-membered ring
  // TTijkl = 5 when i, j, k and l are members of a five-membered ring and at least one is a sp3-hybridized carbon (MMFF atom type 1)
  //
  int MMFF94Common::GetTorsionType(OBAtom* a, OBAtom* b, OBAtom *c, OBAtom *d)
  {
    OBMol *mol = a->GetParent();
    int btab, btbc, btcd;

    btab = GetBondType(a, b);
    btbc = GetBondType(b, c);
    btcd = GetBondType(c, d);
    
    if (btbc == 1)
      return 1;
    
    if (a->IsInRingSize(4) && b->IsInRingSize(4) && c->IsInRingSize(4) && d->IsInRingSize(4))
      if (IsInSameRing(a,b) && IsInSameRing(b,c) && IsInSameRing(c,d))
        return 4;
   
    if (mol->GetBond(b,c)->IsSingle()) {
      if (btab || btcd)
        return 2;
      /*
        unsigned int order1 = GetCXT(0, atoi(d->GetType()), atoi(c->GetType()), atoi(b->GetType()), atoi(a->GetType())); 
        unsigned int order2 = GetCXT(0, atoi(a->GetType()), atoi(b->GetType()), atoi(c->GetType()), atoi(d->GetType()));
    
        cout << "GetTorsionType(" << a->GetType() << ", " << b->GetType() << ", " << c->GetType() << ", " << d->GetType() << ")" << endl;
        cout << "    order1 = " << order1 << endl;
        cout << "    order2 = " << order2 << endl;
        cout << "    btab = " << btab << endl;
        cout << "    btbc = " << btbc << endl;
        cout << "    btcd = " << btcd << endl;
      */
    }
    
    if (a->IsInRingSize(5) && b->IsInRingSize(5) && c->IsInRingSize(5) && d->IsInRingSize(5)) {
      vector<OBRing*> vr;
      vr = mol->GetSSSR();
    
      if( !((GetCachedType(a) == 1) || (GetCachedType(b) == 1) || (GetCachedType(c) == 1) || (GetCachedType(d) == 1)) )
        return 0;

      vector<OBRing*>::iterator ri;
      vector<int>::iterator rj;
      for (ri = vr.begin();ri != vr.end();++ri) { // for each ring
        if (IsAromatic(*ri))
          continue;
	
        if ((*ri)->Size() != 5)
          continue;
        
        if (!(*ri)->IsMember(a) || !(*ri)->IsMember(b) || !(*ri)->IsMember(c) || !(*ri)->IsMember(d))
          continue;
	
        return 5;
      }
    }
 

    return 0;
  }

  // CXB = MC * (I * MA + J) + BTij
  unsigned int MMFF94Common::GetCXB(int type, int a, int b)
  {
    unsigned int cxb;
    cxb = 2 * (a * 136 + b) + type;
    return cxb;
  }
  
  // CXA = MC * (J * MA^2 + I * MA + K) + ATijk
  unsigned int MMFF94Common::GetCXA(int type, int a, int b, int c)
  {
    unsigned int cxa;
    cxa = 9 * (b * 18496 + a * 136 + c) + type;
    return cxa;
  }
  
  // CXS = MC * (J * MA^2 + I * MA + K) + STijk
  unsigned int MMFF94Common::GetCXS(int type, int a, int b, int c)
  {
    unsigned int cxs;
    cxs = 12 * (b * 18496 + a * 136 + c) + type;
    return cxs;
  }
  
  // CXO = J * MA^3 + I * MA^2 + K * MA + L
  unsigned int MMFF94Common::GetCXO(int a, int b, int c, int d)
  {
    unsigned int cxo;
    cxo = b * 2515456 + a * 18496 + c * 136 + d;
    return cxo;
  }
  
  // CXT = MC * (J * MA^3 + K * MA^2 + I * MA + L) + TTijkl
  unsigned int MMFF94Common::GetCXT(int type, int a, int b, int c, int d)
  {
    unsigned int cxt;
    cxt = 6 * (b * 2515456 + c * 18496 + a * 136 + d) + type;
    return cxt;
  }
  
  // CXQ = MC * (I * MA + J) + BTij
  unsigned int MMFF94Common::GetCXQ(int type, int a, int b)
  {
    unsigned int cxq;
    cxq = 2 * (a * 136 + b) + type;
    return cxq;
  }

  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  //
  //  Various tables & misc. functions
  //
  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
 
  // MMFF part V - TABLE I
  bool MMFF94Common::HasLinSet(int atomtype)
  {
    std::vector<OBParameterDB::Query> query;
    query.push_back( OBParameterDB::Query(0, OBVariant(atomtype)) );
    const std::vector<OBVariant> &row = m_database->FindRow(MMFF94SimpleParameterDB::AtomProperties, query);
    if (row.size())
      return row.at(7).AsBool();

    return false;
  }
 
  // MMFF part V - TABLE I
  bool MMFF94Common::HasPilpSet(int atomtype)
  {
    std::vector<OBParameterDB::Query> query;
    query.push_back( OBParameterDB::Query(0, OBVariant(atomtype)) );
    const std::vector<OBVariant> &row = m_database->FindRow(MMFF94SimpleParameterDB::AtomProperties, query);
    if (row.size())
      return row.at(4).AsBool();

    return false;
  }
  
  // MMFF part V - TABLE I
  bool MMFF94Common::HasAromSet(int atomtype)
  {
    std::vector<OBParameterDB::Query> query;
    query.push_back( OBParameterDB::Query(0, OBVariant(atomtype)) );
    const std::vector<OBVariant> &row = m_database->FindRow(MMFF94SimpleParameterDB::AtomProperties, query);
    if (row.size())
      return row.at(6).AsBool();

    return false;
  }
 
  // MMFF part V - TABLE I
  bool MMFF94Common::HasSbmbSet(int atomtype)
  {
    std::vector<OBParameterDB::Query> query;
    query.push_back( OBParameterDB::Query(0, OBVariant(atomtype)) );
    const std::vector<OBVariant> &row = m_database->FindRow(MMFF94SimpleParameterDB::AtomProperties, query);
    if (row.size())
      return row.at(8).AsBool();

    return false;
  }

  // MMFF part V - TABLE I
  int MMFF94Common::GetCrd(int atomtype)
  {
    std::vector<OBParameterDB::Query> query;
    query.push_back( OBParameterDB::Query(0, OBVariant(atomtype)) );
    const std::vector<OBVariant> &row = m_database->FindRow(MMFF94SimpleParameterDB::AtomProperties, query);
    if (row.size())
      return row.at(2).AsInt();

    return 0;
  }

  // MMFF part V - TABLE I
  int MMFF94Common::GetVal(int atomtype)
  {
    std::vector<OBParameterDB::Query> query;
    query.push_back( OBParameterDB::Query(0, OBVariant(atomtype)) );
    const std::vector<OBVariant> &row = m_database->FindRow(MMFF94SimpleParameterDB::AtomProperties, query);
    if (row.size())
      return row.at(3).AsInt();

    return 0;
  }

  // MMFF part V - TABLE I
  int MMFF94Common::GetMltb(int atomtype)
  {
    std::vector<OBParameterDB::Query> query;
    query.push_back( OBParameterDB::Query(0, OBVariant(atomtype)) );
    const std::vector<OBVariant> &row = m_database->FindRow(MMFF94SimpleParameterDB::AtomProperties, query);
    if (row.size())
      return row.at(5).AsInt();

    return 0;
  }
  
  // MMFF part I - TABLE IV
  int MMFF94Common::EqLvl2(int type)
  {
    std::vector<OBParameterDB::Query> query;
    query.push_back( OBParameterDB::Query(0, OBVariant(type)) );
    const std::vector<OBVariant> &row = m_database->FindRow(MMFF94SimpleParameterDB::AtomTypeLevels, query);
    if (row.size())
      return row.at(1).AsInt();

    return type; 
  }
  
  // MMFF part I - TABLE IV
  int MMFF94Common::EqLvl3(int type)
  {
    std::vector<OBParameterDB::Query> query;
    query.push_back( OBParameterDB::Query(0, OBVariant(type)) );
    const std::vector<OBVariant> &row = m_database->FindRow(MMFF94SimpleParameterDB::AtomTypeLevels, query);
    if (row.size())
      return row.at(2).AsInt();

    return type; 
  }
  
  // MMFF part I - TABLE IV
  int MMFF94Common::EqLvl4(int type)
  {
    std::vector<OBParameterDB::Query> query;
    query.push_back( OBParameterDB::Query(0, OBVariant(type)) );
    const std::vector<OBVariant> &row = m_database->FindRow(MMFF94SimpleParameterDB::AtomTypeLevels, query);
    if (row.size())
      return row.at(3).AsInt();

    return type; 
  }

  // MMFF part I - TABLE IV
  int MMFF94Common::EqLvl5(int type)
  {
    std::vector<OBParameterDB::Query> query;
    query.push_back( OBParameterDB::Query(0, OBVariant(type)) );
    const std::vector<OBVariant> &row = m_database->FindRow(MMFF94SimpleParameterDB::AtomTypeLevels, query);
    if (row.size())
      return row.at(4).AsInt();

    return type; 
  }
  
  // MMFF part V - TABLE VI
  double MMFF94Common::GetZParam(OBAtom* atom)
  {
    if (atom->IsHydrogen())
      return 1.395;
    if (atom->IsCarbon())
      return 2.494;
    if (atom->IsNitrogen())
      return 2.711;
    if (atom->IsOxygen())
      return 3.045;
    if (atom->GetAtomicNum() == 9) // F
      return 2.847;
    if (atom->GetAtomicNum() == 14) // Si
      return 2.350;
    if (atom->IsPhosphorus())
      return 2.350;
    if (atom->IsSulfur())
      return 2.980;
    if (atom->GetAtomicNum() == 17) // Cl
      return 2.909;
    if (atom->GetAtomicNum() == 35) // Br
      return 3.017;
    if (atom->GetAtomicNum() == 53) // I
      return 3.086;

    return 0.0;
  }
 
  // MMFF part V - TABLE VI
  double MMFF94Common::GetCParam(OBAtom* atom)
  {
    if (atom->GetAtomicNum() == 5) // B
      return 0.704;
    if (atom->IsCarbon())
      return 1.016;
    if (atom->IsNitrogen())
      return 1.113;
    if (atom->IsOxygen())
      return 1.337;
    if (atom->GetAtomicNum() == 14) // Si
      return 0.811;
    if (atom->IsPhosphorus())
      return 1.068;
    if (atom->IsSulfur())
      return 1.249;
    if (atom->GetAtomicNum() == 17) // Cl
      return 1.078;
    if (atom->GetAtomicNum() == 33) // As
      return 0.825;

    return 0.0;
  }
 
  // MMFF part V - TABLE X
  double MMFF94Common::GetUParam(OBAtom* atom)
  {
    if (atom->IsCarbon())
      return 2.0;
    if (atom->IsNitrogen())
      return 2.0;
    if (atom->IsOxygen())
      return 2.0;
    if (atom->GetAtomicNum() == 14) // Si
      return 1.25;
    if (atom->IsPhosphorus())
      return 1.25;
    if (atom->IsSulfur())
      return 1.25;
    
    return 0.0;
  }
  
  // MMFF part V - TABLE X
  double MMFF94Common::GetVParam(OBAtom* atom)
  {
    if (atom->IsCarbon())
      return 2.12;
    if (atom->IsNitrogen())
      return 1.5;
    if (atom->IsOxygen())
      return 0.2;
    if (atom->GetAtomicNum() == 14) // Si
      return 1.22;
    if (atom->IsPhosphorus())
      return 2.4;
    if (atom->IsSulfur())
      return 0.49;
    
    return 0.0;
  }
    
  // R Blom and A Haaland, J. Mol. Struct., 128, 21-27 (1985)
  double MMFF94Common::GetCovalentRadius(OBAtom* a) {

    switch (a->GetAtomicNum()) {
    case 1:
      return 0.33; // corrected value from MMFF part V
    case 5:
      return 0.81;
    case 6:
      return 0.77; // corrected value from MMFF part V
    case 7:
      return 0.73;
    case 8:
      return 0.72;
    case 9:
      return 0.74;
    case 13:
      return 1.22;
    case 14:
      return 1.15;
    case 15:
      return 1.09;
    case 16:
      return 1.03;
    case 17:
      return 1.01;
    case 31:
      return 1.19;
    case 32:
      return 1.20;
    case 33:
      return 1.20;
    case 34:
      return 1.16;
    case 35:
      return 1.15;
    case 44:
      return 1.46;
    case 50:
      return 1.40;
    case 51:
      return 1.41;
    case 52:
      return 1.35;
    case 53:
      return 1.33;
    case 81:
      return 1.51;
    case 82:
      return 1.53;
    case 83:
      return 1.55;
    default:
      return etab.GetCovalentRad(a->GetAtomicNum());
    }
  }
  
  double MMFF94Common::GetBondLength(OBAtom* a, OBAtom* b)
  {
    double rab;

    int bondType = GetBondType(a, b);
    int type_a = GetCachedType(a);
    int type_b = GetCachedType(b);

    std::vector<OBParameterDB::Query> query;
    query.push_back( OBParameterDB::Query(0, OBVariant(bondType)) );
    query.push_back( OBParameterDB::Query(1, OBVariant(type_a)) );
    query.push_back( OBParameterDB::Query(2, OBVariant(type_b)) );
    std::vector<OBVariant> row = m_database->FindRow(MMFF94SimpleParameterDB::BondParameters, query);
    if (!row.size()) {
      // try inverse order
      OBParameterDB::Query swap(query[1]);
      query[1] = query[2];
      query[2] = swap;
      row = m_database->FindRow(MMFF94SimpleParameterDB::BondParameters, query);
    }

    if (row.size())
      rab = row.at(4).AsDouble();
    else
      rab = GetRuleBondLength(a, b); 
  
    return rab;
  }
 
  // MMFF part V - page 625
  double MMFF94Common::GetRuleBondLength(OBAtom* a, OBAtom* b)
  {
    double r0ab, r0a, r0b, c, Xa, Xb;
    int Ha, Hb, BOab;
    r0a = GetCovalentRadius(a);
    r0b = GetCovalentRadius(b);
    Xa = etab.GetAllredRochowElectroNeg(a->GetAtomicNum());
    Xb = etab.GetAllredRochowElectroNeg(b->GetAtomicNum());
    
   
    if (a->IsHydrogen())
      r0a = 0.33;
    if (b->IsHydrogen())
      r0b = 0.33;
    
    if (a->IsHydrogen() || b->IsHydrogen())
      c = 0.050;
    else
      c = 0.085;

    if (GetMltb(GetCachedType(a) == 3))
      Ha = 1;
    else if ((GetMltb(GetCachedType(a)) == 1) || (GetMltb(GetCachedType(a)) == 2))
      Ha = 2;
    else
      Ha = 3;

    if (GetMltb(GetCachedType(b) == 3))
      Hb = 1;
    else if ((GetMltb(GetCachedType(b)) == 1) || (GetMltb(GetCachedType(b)) == 2))
      Hb = 2;
    else
      Hb = 3;

    BOab = a->GetBond(b)->GetBondOrder();
    if ((GetMltb(GetCachedType(a)) == 1) && (GetMltb(GetCachedType(b)) == 1))
      BOab = 4;
    if ((GetMltb(GetCachedType(a)) == 1) && (GetMltb(GetCachedType(b)) == 2))
      BOab = 5;
    if ((GetMltb(GetCachedType(a)) == 2) && (GetMltb(GetCachedType(b)) == 1))
      BOab = 5;
    if (IsAromatic(a->GetBond(b)))
      if (!HasPilpSet(GetCachedType(a)) && !HasPilpSet(GetCachedType(b)))
        BOab = 4;
      else
        BOab = 5;
     
    switch (BOab) {
    case 5:
      r0a -= 0.04;
      r0b -= 0.04;
      break;
    case 4:
      r0a -= 0.075;
      r0b -= 0.075;
      break;
    case 3:
      r0a -= 0.17;
      r0b -= 0.17;
      break;
    case 2:
      r0a -= 0.10;
      r0b -= 0.10;
      break;
    case 1:
      if (Ha == 1)
        r0a -= 0.08;
      if (Ha == 2)
        r0a -= 0.03;
      if (Hb == 1)
        r0b -= 0.08;
      if (Hb == 2)
        r0b -= 0.03;
    }
    
    /*
      cout << "Ha=" << Ha << "  Hb=" << Hb << "  BOab=" << BOab << endl;
      cout << "r0a=" << r0a << "  Xa=" << Xa << endl;
      cout << "r0b=" << r0b << "  Xb=" << Xb << endl;
      cout << "r0a + r0b=" << r0a +r0b << endl;
      cout << "c=" << c << "  |Xa-Xb|=" << fabs(Xa-Xb) << "  |Xa-Xb|^1.4=" << pow(fabs(Xa-Xb), 1.4) << endl;
    */
    r0ab = r0a + r0b - c * pow(fabs(Xa - Xb), 1.4) - 0.008; 

    return r0ab;
  }

  bool MMFF94Common::IsInSameRing(OBAtom* a, OBAtom* b)
  {
    OBMol *mol = a->GetParent();
    bool a_in, b_in;
    vector<OBRing*> vr;
    vr = mol->GetSSSR();
    
    vector<OBRing*>::iterator i;
    vector<int>::iterator j;
    
    for (i = vr.begin();i != vr.end();i++) {
      a_in = false;
      b_in = false;
      for(j = (*i)->_path.begin();j != (*i)->_path.end();j++) {
        if ((unsigned)(*j) == a->GetIdx())
          a_in = true;
        if ((unsigned)(*j) == b->GetIdx())
          b_in = true;
      }
      
      if (a_in && b_in)
        return true;
    }
    
    return false;
  }

  bool MMFF94Common::IsAromatic(OBRing *ring) const
  {
    vector<int>::iterator i;
    for (i = ring->_path.begin(); i != ring->_path.end(); ++i)
      if (!m_aromAtoms.at(*i-1))
        return false;

    return true;
  }

  int MMFF94Common::GetElementRow(OBAtom *atom)
  {
    int row;
    
    row = 0;

    if (atom->GetAtomicNum() > 2)
      row++;
    if (atom->GetAtomicNum() > 10)
      row++;
    if (atom->GetAtomicNum() > 18)
      row++;
    if (atom->GetAtomicNum() > 36)
      row++;
    if (atom->GetAtomicNum() > 54)
      row++;
    if (atom->GetAtomicNum() > 86)
      row++;
    
    return row;
  }

}
} // end namespace OpenBabel

//! \file forcefieldmmff94.cpp
//! \brief MMFF94 force field
