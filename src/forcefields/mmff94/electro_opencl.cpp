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

#include "electro_opencl.h"
#include "parameter.h"
#include "common.h"

#include <openbabel/mol.h>

#include <OBLogFile>
#include <OBVectorMath>

using namespace std;

namespace OpenBabel {
namespace OBFFs {
 
  struct Parameter {
    enum { Class, TypeA, TypeB, TypeC, Ka, Theta0 };
  };
 
  MMFF94ElectroTermOpenCL::MMFF94ElectroTermOpenCL(OBFunction *function, MMFF94Common *common) : OBFunctionTerm(function), m_common(common)
  {
    m_value = 999999.99;
  }
   
  void MMFF94ElectroTermOpenCL::Compute(OBFunction::Computation computation)
  {
    cout << "MMFF94ElectroTermOpenCL::Compute" << endl;
    m_value = 0.0;

    /*
    try {
    cl::KernelFunctor func = m_kernel.bind(m_queue, cl::NDRange(1, 1), cl::NDRange(1, 1));

    func().wait();
    } catch (cl::Error err) {
      std::cout << "  ERROR: " << err.what() << "(" << err.err() << ")" << std::endl;
    }
    */



 /*
    unsigned int idxA, idxB;
    double rab, e;
    Eigen::Vector3d Fa, Fb;
    for (int i = 0; i < (int)m_calcs.size(); ++i) {
      idxA = m_calcs[i].idx1;
      idxB = m_calcs[i].idx2;
  
      if (computation == OBFunction::Gradients) {
        rab = VectorDistanceDerivative(m_function->GetPositions()[idxA], m_function->GetPositions()[idxB], Fa, Fb);
        rab += 0.05; // ??
        const double rab2 = rab * rab;
        const double dE = - (m_calcs[i].qq / rab2);
        Fa *= dE;
        Fb *= dE;
        m_function->GetGradients()[idxA] += Fa;
        m_function->GetGradients()[idxB] += Fb;
      } else {
        const Eigen::Vector3d ab = m_function->GetPositions()[idxA] - m_function->GetPositions()[idxB];
        rab = ab.norm();
        rab += 0.05; // ??
      }

      e = m_calcs[i].qq / rab;
      m_value += e;
    }
    */ 
    cout << "E_ele = " << m_value << endl;
 
  }

  void MMFF94ElectroTermOpenCL::InitSelfPairs(OBMol &mol)
  {
    unsigned int numAtoms = mol.NumAtoms();

    for (unsigned int i = 0; i < numAtoms; ++i) {
      for (unsigned int j = 0; j < numAtoms; ++j) {
        if (i == j) {
          m_selfPairs.push_back(std::pair<unsigned int, unsigned int>(i, j));
          continue;
        }

        OBAtom *a = mol.GetAtom(i+1);
        OBAtom *b = mol.GetAtom(j+1);

        if (a->IsConnected(b))
          m_selfPairs.push_back(std::pair<unsigned int, unsigned int>(i, j));
        else if (a->IsOneThree(b))
          m_selfPairs.push_back(std::pair<unsigned int, unsigned int>(i, j));
        else if (a->IsOneFour(b) )
          m_oneFourPairs.push_back(std::pair<unsigned int, unsigned int>(i, j));
      }
    }
  }

  double energy(const Eigen::Vector3d &ai, const Eigen::Vector3d &aj, double qi, double qj)
  {
    Eigen::Vector3d r = ai - aj;
    double dist = r.norm() + 0.05;
    double e = 0.5 * 332.0716 * qi * qj / dist;
  }

  double MMFF94ElectroTermOpenCL::ComputeTotalEnergy()
  {
    const std::vector<double> &charges = m_common->m_pCharges;
    int numAtoms = charges.size();

    double E = 0.0;
    for (int i = 0; i < numAtoms; ++i) {
      double e1 = 0.0;

      for (int j = 0; j < numAtoms; ++j) {
        const Eigen::Vector3d &ai = m_function->GetPositions().at(i);
        const Eigen::Vector3d &aj = m_function->GetPositions().at(j);
        e1 += energy(ai, aj, charges.at(i), charges.at(j));
      }

      cout << "e = " << e1 << endl;
      E += e1;
    }

    cout << "E_total = " << E << endl;
    return E;  
  }

  double MMFF94ElectroTermOpenCL::ComputeSelfEnergy(OBMol &mol)
  {
    const std::vector<double> &charges = m_common->m_pCharges;
    int numAtoms = charges.size();

    double E = 0.0;
    for (unsigned int k = 0; k < m_selfPairs.size(); ++k) {
      unsigned int i = m_selfPairs.at(k).first;
      unsigned int j = m_selfPairs.at(k).second;
      const Eigen::Vector3d &ai = m_function->GetPositions().at(i);
      const Eigen::Vector3d &aj = m_function->GetPositions().at(j);
          
      E += energy(ai, aj, charges.at(i), charges.at(j));
    }

    for (unsigned int k = 0; k < m_oneFourPairs.size(); ++k) {
      unsigned int i = m_oneFourPairs.at(k).first;
      unsigned int j = m_oneFourPairs.at(k).second;
      const Eigen::Vector3d &ai = m_function->GetPositions().at(i);
      const Eigen::Vector3d &aj = m_function->GetPositions().at(j);
          
      E += 0.25 * energy(ai, aj, charges.at(i), charges.at(j));
    }



    /*
    for (int i = 0; i < numAtoms; ++i) {
      double e1 = 0.0;

      for (int j = 0; j < numAtoms; ++j) {
        bool self = false, oneFour = false;
        if (i == j)
          self = true;

        OBAtom *a = mol.GetAtom(i+1);
        OBAtom *b = mol.GetAtom(j+1);

        if (a->IsConnected(b))
          self = true;
        else if (a->IsOneThree(b))
          self = true;
        else if (a->IsOneFour(b) ) {
          self = true;
          oneFour = true;
        }

        if (!self)
          continue;

        const Eigen::Vector3d &ai = m_function->GetPositions().at(i);
        const Eigen::Vector3d &aj = m_function->GetPositions().at(j);

        if (oneFour)
          e1 += 0.25 * energy(ai, aj, charges.at(i), charges.at(j)); // 1-4 are scaled by factor 0.75 --> self = 0.25
        else
          e1 += energy(ai, aj, charges.at(i), charges.at(j));
      }

      cout << "e = " << e1 << endl;
      E += e1;
    }
    */

    cout << "E_self = " << E << endl;
    return E;  
  }


  bool MMFF94ElectroTermOpenCL::Setup(/*const*/ OBMol &mol)
  {
    OBAtom *a, *b, *c;

    OBLogFile *logFile = m_function->GetLogFile();
    OBParameterDB *database = m_function->GetParameterDB();

    if (logFile->IsLow())
      logFile->Write("SETTING UP ELECTROSTATIC CALCULATIONS...\n");
    std::stringstream ss;

    InitSelfPairs(mol);

    try {
      m_context = cl::Context(CL_DEVICE_TYPE_GPU, 0, NULL, NULL);
      cl::vector<cl::Device> devices = m_context.getInfo<CL_CONTEXT_DEVICES>();
      ss << "  # OpenCL devices: " << devices.size() << std::endl;

      if (devices.empty())
        return false;

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

      // Read the kernel source file
      std::ifstream ifs;
      ifs.open("kernel.cl");
      if (!ifs) {
        std::stringstream msg;
        msg << "Cannot open kernel.cl..." << std::endl;
        logFile->Write(msg.str());
        return false;
      }

      std::string srcCode;
      std::string line;
      while (std::getline(ifs, line)) {
        srcCode += line + "\n";
      }
      
      ifs.close();

      cl::Program::Sources source(1, std::make_pair(srcCode.c_str(), srcCode.size()));
      cl::Program program = cl::Program(m_context, source);
      program.build(devices);

      m_kernel = cl::Kernel(program, "electrostaticKernel");
      m_queue = cl::CommandQueue(m_context, devices[0], 0);


      const std::vector<double> &charges = m_common->m_pCharges;
      int numAtoms = charges.size();

      int p = 4;
      int globalWorkSize = numAtoms + p - numAtoms % p;
      ss << "globalWorkSize = " << globalWorkSize << std::endl;

      cout << "---- TOTAL (CPU) ----" << endl;
      ComputeTotalEnergy();
      cout << "---- SELF (CPU) ----" << endl;
      ComputeSelfEnergy(mol);
      cout << "---------------------" << endl;
      
      cl_float positions[4*numAtoms];
      for (int i = 0; i < numAtoms; ++i) {
        int offset = 4 * i;
        const Eigen::Vector3d &atomPos = m_function->GetPositions().at(i);
        positions[offset  ] = atomPos.x();
        positions[offset+1] = atomPos.y();
        positions[offset+2] = atomPos.y();
        positions[offset+3] = charges.at(i);
      }

      cl::Buffer dev_pos(m_context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 4 * numAtoms * sizeof(cl_float), (void*)positions);
      cl::Buffer dev_grad(m_context, CL_MEM_READ_WRITE, 4 * numAtoms * sizeof(cl_float));

      cl::KernelFunctor func = m_kernel.bind(m_queue, cl::NDRange(numAtoms), cl::NullRange);
      cl::LocalSpaceArg dev_sharedPos = cl::__local(p * sizeof(cl_float));
      func(dev_pos, dev_grad, numAtoms, dev_sharedPos);

      m_queue.enqueueReadBuffer(dev_grad, CL_TRUE, 0, 4 * numAtoms * sizeof(cl_float), (void*)positions);

      for (int i = 0; i < numAtoms; i++) {
        int offset = 4*i;
        ss << i  << " = pos(" << positions[offset] << ", " << positions[offset+1]
           << ", " << positions[offset+2] << ")  q = " << positions[offset+3] << std::endl;
      }

    } catch (cl::Error err) {
      ss << "  ERROR: " << err.what() << "(" << err.err() << ")" << std::endl;
    }


    if (logFile->IsLow())
      logFile->Write(ss.str());



    // 
    // Electrostatic Calculations
    //
    /*
    if (logFile->IsLow())
      logFile->Write("SETTING UP ELECTROSTATIC CALCULATIONS...\n");
 
    Calculation elecalc;
    m_calcs.clear();
    
    FOR_PAIRS_OF_MOL(p, mol) {
      a = mol.GetAtom((*p)[0]);
      b = mol.GetAtom((*p)[1]);
      
      elecalc.qq = 332.0716 * m_common->GetPartialCharge(a->GetIdx()-1) * m_common->GetPartialCharge(b->GetIdx()-1);
      
      if (elecalc.qq) {
        elecalc.idx1 = a->GetIdx() - 1;
        elecalc.idx2 = b->GetIdx() - 1;
        
        // 1-4 scaling
        if (a->IsOneFour(b))
          elecalc.qq *= 0.75;
	  
        m_calcs.push_back(elecalc);
      }
    }
    */

    

    return true;


  }       

}
} // end namespace OpenBabel

//! \file forcefieldmmff94.cpp
//! \brief MMFF94 force field
