/*********************************************************************
MMFF94ElectroTermOpenCL - MMFF94 force field bond stratching term

Copyright (C) 2006-2008,2009 by Tim Vandermeersch
 
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

#include "vdw_opencl.h"
#include "parameter.h"
#include "common.h"

#include <openbabel/mol.h>

#include <OBLogFile>
#include <OBVectorMath>

using namespace std;

namespace OpenBabel {
namespace OBFFs {
 
  struct Parameter {
    enum { Type, Alpha_i, N_i, A_i, G_i, DA };
  };
 
  MMFF94VDWTermOpenCL::MMFF94VDWTermOpenCL(OBFunction *function, MMFF94Common *common) : OBFunctionTerm(function), m_common(common)
  {
    m_value = 999999.99;
  }
   
  void MMFF94VDWTermOpenCL::Compute(OBFunction::Computation computation)
  {
    //cout << "MMFF94VDWTermOpenCL::Compute" << endl;
    m_value = 0.0;
      
    // setup input (posistions + charges) to copy to device
    cl_float hostData[4*m_numAtoms];
    for (int i = 0; i < m_numAtoms; ++i) {
      const Eigen::Vector3d &atomPos = m_function->GetPositions().at(i);
      m_calcs[i].x = atomPos.x();
      m_calcs[i].y = atomPos.y();
      m_calcs[i].z = atomPos.z();
    }

    int p = 64;
    int globalWorkSize = m_numAtoms + p - m_numAtoms % p;
 
    try {
      // write positions (xyz) and charges (w) to device
      m_queue.enqueueWriteBuffer(m_devPos, CL_TRUE, 0, 8 * m_numAtoms * sizeof(cl_float), (void*)&m_calcs[0]);

      // execute kernel
      cl::KernelFunctor func = m_kernel.bind(m_queue, cl::NDRange(m_numAtoms), cl::NullRange);
      cl::LocalSpaceArg sharedPos = cl::__local(8 * p * sizeof(cl_float));
      func(m_devPos, m_devGrad, m_numAtoms, sharedPos);

      // read back gradients (xyz) and energies (w) from device
      m_queue.enqueueReadBuffer(m_devGrad, CL_TRUE, 0, 4 * m_numAtoms * sizeof(cl_float), (void*)hostData);

    } catch (cl::Error err) {
      std::cout << "  ERROR: " << err.what() << "(" << err.err() << ")" << std::endl;
    }

    // sum the energies
    double devEnergy = 0.0;
    for (int i = 0; i < m_numAtoms; i++) {
      int offset = 4*i;
      m_function->GetGradients()[i].x() += hostData[offset];
      m_function->GetGradients()[i].y() += hostData[offset+1];
      m_function->GetGradients()[i].z() += hostData[offset+2];
//      cout << i  << " = pos(" << hostData[offset] << ", " << hostData[offset+1] << ", " << hostData[offset+2] << ")  q = " << hostData[offset+3] << std::endl;
      devEnergy += hostData[offset+3];
    }
    //cout << "devEnergy = " << devEnergy << endl;

    // compute the self energy (i.e. i == j, bonded atoms, 1-3 (angle terminal atoms) and the 1-4 atoms contribution times 0.25
    double selfEnergy;
   // if (computation == OBFunction::Gradients)
   //   selfEnergy = ComputeSelfGradients();
   // else
      selfEnergy = ComputeSelfEnergy();

    m_value = devEnergy - selfEnergy;

    //cout << "E_VDW_dev = " << devEnergy << endl;
    //cout << "E_VDW_self = " << selfEnergy << endl;
  //  cout << "E_VDW = " << m_value << endl;
  }

  void MMFF94VDWTermOpenCL::InitSelfPairs(OBMol &mol)
  {
    unsigned int numAtoms = mol.NumAtoms();

    for (unsigned int i = 0; i < numAtoms; ++i) {
      for (unsigned int j = 0; j < numAtoms; ++j) {
        /*
        if (i == j) {
          m_selfPairs.push_back(std::pair<unsigned int, unsigned int>(i, j));
          continue;
        }
        */

        OBAtom *a = mol.GetAtom(i+1);
        OBAtom *b = mol.GetAtom(j+1);

        if (a->IsConnected(b))
          m_selfPairs.push_back(std::pair<unsigned int, unsigned int>(i, j));
        else if (a->IsOneThree(b))
          m_selfPairs.push_back(std::pair<unsigned int, unsigned int>(i, j));
      }
    }
  }

  double vdwEnergy(const MMFF94VDWTermOpenCL::Calculation &ai, const MMFF94VDWTermOpenCL::Calculation &aj)
  {
    Eigen::Vector3d r(ai.x - aj.x, ai.y - aj.y, ai.z - aj.z);
    double rab = r.norm();

    double R_AB, R_AB7, epsilon;
    if (ai.DA == 1) { // hydrogen bond donor
      R_AB = 0.5 * (ai.R_AA + aj.R_AA);
      double R_AB2 = R_AB * R_AB;
      double R_AB4 = R_AB2 * R_AB2;
      double R_AB6 = R_AB4 * R_AB2;

      if (aj.DA == 2) { // hydrogen bond acceptor
        epsilon = 0.5 * (181.16 * ai.alpha_G * aj.alpha_G) / (ai.sqrt_alpha_N + aj.sqrt_alpha_N) * (1.0 / R_AB6);
        // R_AB is scaled to 0.8 for D-A interactions. The value used in the calculation of epsilon is not scaled. 
        R_AB *= 0.8;
      } else
        epsilon = (181.16 * ai.alpha_G * aj.alpha_G) / (ai.sqrt_alpha_N + aj.sqrt_alpha_N) * (1.0 / R_AB6);

      R_AB2 = R_AB * R_AB;
      R_AB4 = R_AB2 * R_AB2;
      R_AB6 = R_AB4 * R_AB2;
      R_AB7 = R_AB6 * R_AB;
    } else if (aj.DA == 1) { // hydrogen bond donor
      R_AB = 0.5 * (ai.R_AA + aj.R_AA);
      double R_AB2 = R_AB * R_AB;
      double R_AB4 = R_AB2 * R_AB2;
      double R_AB6 = R_AB4 * R_AB2;

      if (ai.DA == 2) { // hydrogen bond acceptor
        epsilon = 0.5 * (181.16 * ai.alpha_G * aj.alpha_G) / (ai.sqrt_alpha_N + aj.sqrt_alpha_N) * (1.0 / R_AB6);
        // R_AB is scaled to 0.8 for D-A interactions. The value used in the calculation of epsilon is not scaled. 
        R_AB *= 0.8;
      } else
        epsilon = (181.16 * ai.alpha_G * aj.alpha_G) / (ai.sqrt_alpha_N + aj.sqrt_alpha_N) * (1.0 / R_AB6);

      R_AB2 = R_AB * R_AB;
      R_AB4 = R_AB2 * R_AB2;
      R_AB6 = R_AB4 * R_AB2;
      R_AB7 = R_AB6 * R_AB;
    } else {
      double g_AB = (ai.R_AA - aj.R_AA) / (ai.R_AA + aj.R_AA);
      double g_AB2 = g_AB * g_AB;
      R_AB =  0.5 * (ai.R_AA + aj.R_AA) * (1.0 + 0.2 * (1.0 - exp(-12.0 * g_AB2)));
      double R_AB2 = R_AB * R_AB;
      double R_AB4 = R_AB2 * R_AB2;
      double R_AB6 = R_AB4 * R_AB2;
      R_AB7 = R_AB6 * R_AB;
      epsilon = (181.16 * ai.alpha_G * aj.alpha_G) / (ai.sqrt_alpha_N + aj.sqrt_alpha_N) * (1.0 / R_AB6);
    }
    
    //cout << "R_AB = " << R_AB << endl;

    const double rab7 = rab*rab*rab*rab*rab*rab*rab;
    double erep = (1.07 * R_AB) / (rab + 0.07 * R_AB); //***
    double erep7 = erep*erep*erep*erep*erep*erep*erep;
    double eattr = (((1.12 * R_AB7) / (rab7 + 0.12 * R_AB7)) - 2.0);

    //cout << "epsilon = " << epsilon << endl;
    //cout << "erep7 = " << erep7 << endl;
    //cout << "eattr = " << eattr << endl;


    double e = 0.5 * epsilon * erep7 * eattr;

    return e;
  }

  double MMFF94VDWTermOpenCL::ComputeTotalEnergy()
  {
    double E = 0.0;
    for (int i = 0; i < m_numAtoms; ++i) {
      for (int j = 0; j < m_numAtoms; ++j) {
        E += vdwEnergy(m_calcs.at(i), m_calcs.at(j));
      }
    }

    return E;  
  }

  double MMFF94VDWTermOpenCL::ComputeSelfEnergy()
  {
    double E = 0.0;
    for (unsigned int k = 0; k < m_selfPairs.size(); ++k) {
      unsigned int i = m_selfPairs.at(k).first;
      unsigned int j = m_selfPairs.at(k).second;
      E += vdwEnergy(m_calcs.at(i), m_calcs.at(j));
    }

    //cout << "E_self = " << E << endl;
    return E;  
  }

  double MMFF94VDWTermOpenCL::ComputeSelfGradients()
  {
    double E = 0.0;
    /*
    for (unsigned int k = 0; k < m_selfPairs.size(); ++k) {
      unsigned int i = m_selfPairs.at(k).first;
      unsigned int j = m_selfPairs.at(k).second;
      const Eigen::Vector3d &ai = m_function->GetPositions().at(i);
      const Eigen::Vector3d &aj = m_function->GetPositions().at(j);     
      
      Eigen::Vector3d r = ai - aj;
      double dist = r.norm() + 0.05;
      double dist2 = dist * dist;

      double QiQj = charges.at(i) * charges.at(j);

      double e = QiQj / dist;
      double dE = 332.0716 * QiQj / (dist2 * dist);

      m_function->GetGradients()[i] += r * dE;
      m_function->GetGradients()[j] -= r * dE;

      E += e;
    }
*/
    //cout << "E_self = " << E << endl;
    return E;  
  }



  bool MMFF94VDWTermOpenCL::Setup(/*const*/ OBMol &mol)
  {
    OBLogFile *logFile = m_function->GetLogFile();
    OBParameterDB *database = m_function->GetParameterDB();

    if (logFile->IsLow())
      logFile->Write("SETTING UP VAN DER WAALS CALCULATIONS...\n");
    std::stringstream ss;

    m_numAtoms = mol.NumAtoms();
    for (unsigned int i = 0; i < m_numAtoms; ++i) {
      OBAtom *atom = mol.GetAtom(i+1);
    
      std::vector<OBParameterDB::Query> query;
      query.push_back( OBParameterDB::Query(Parameter::Type, OBVariant(m_common->GetCachedType(atom))) );
      std::vector<OBVariant> row = database->FindRow(MMFF94SimpleParameterDB::VanDerWaalsParameters, query);
 
      if (row.empty()) {
        if (logFile->IsLow()) {
          std::stringstream ss;
          ss << "   COULD NOT FIND VAN DER WAALS PARAMETERS FOR " << m_common->GetCachedType(atom) << " (MMFF94 Atom Type)..." << endl;
          logFile->Write(ss.str());
        }

        return false;
      }

      Calculation calc;

      double alpha = row.at(Parameter::Alpha_i).AsDouble();
      double Ai = row.at(Parameter::A_i).AsDouble();
      double Ni = row.at(Parameter::N_i).AsDouble();
      double Gi = row.at(Parameter::G_i).AsDouble();

      calc.R_AA = Ai * pow(alpha, 0.25);
      calc.sqrt_alpha_N = sqrt(alpha / Ni);
      calc.alpha_G = alpha * Gi;
      calc.DA = row.at(Parameter::DA).AsInt();

      m_calcs.push_back(calc);
    }

    InitSelfPairs(mol);

    try {
      m_context = cl::Context(CL_DEVICE_TYPE_GPU, 0, NULL, NULL);
      cl::vector<cl::Device> devices = m_context.getInfo<CL_CONTEXT_DEVICES>();
      ss << "  # OpenCL devices: " << devices.size() << std::endl;

      if (devices.empty())
        return false;

      // Print some info about the device
      ss << "  device 1:" << std::endl;
      cl_uint ret_uint;
      size_t ret_size;
      size_t ret_psize[3];
      devices[0].getInfo(CL_DEVICE_MAX_COMPUTE_UNITS, &ret_uint);
      ss << "    CL_DEVICE_MAX_COMPUTE_UNITS = " << ret_uint << std::endl;
      devices[0].getInfo(CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, &ret_uint);
      ss << "    CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS = " << ret_uint << std::endl;
      devices[0].getInfo(CL_DEVICE_MAX_WORK_ITEM_SIZES, &ret_psize);
      ss << "    CL_DEVICE_MAX_WORK_ITEM_SIZESS = (" << ret_psize[0] << ", " << ret_psize[1] << ", " << ret_psize[2] << ")" << std::endl;
      devices[0].getInfo(CL_DEVICE_MAX_WORK_GROUP_SIZE, &ret_size);
      ss << "    CL_DEVICE_MAX_WORK_GROUP_SIZE = " << ret_size << std::endl;

      // Open the kernel source file
      std::ifstream ifs;
      ifs.open("mmff94vdwkernel.cl");
      if (!ifs) {
        std::stringstream msg;
        msg << "Cannot open mmff94vdwkernel.cl..." << std::endl;
        logFile->Write(msg.str());
        return false;
      }

      // Read the file into a string
      std::string srcCode;
      std::string line;
      while (std::getline(ifs, line)) {
        srcCode += line + "\n";
      }
      ifs.close();

      // Create the OpenCL program 
      cl::Program::Sources source(1, std::make_pair(srcCode.c_str(), srcCode.size()));
      cl::Program program = cl::Program(m_context, source);
      program.build(devices);

      // Create the kernel
      m_kernel = cl::Kernel(program, "vdwKernel");
      m_queue = cl::CommandQueue(m_context, devices[0], 0);

      m_devPos = cl::Buffer(m_context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 8 * m_numAtoms * sizeof(cl_float), (void*)&m_calcs[0]);
      m_devGrad = cl::Buffer(m_context, CL_MEM_READ_WRITE, 4 * m_numAtoms * sizeof(cl_float));

    } catch (cl::Error err) {
      ss << "  ERROR: " << err.what() << "(" << err.err() << ")" << std::endl;
    }


    if (logFile->IsLow())
      logFile->Write(ss.str());

    return true;
  }       

}
} // end namespace OpenBabel

//! \file forcefieldmmff94.cpp
//! \brief MMFF94 force field
