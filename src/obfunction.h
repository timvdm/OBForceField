/**********************************************************************
obfunction.h - Base class for force fields and scroring functions which
               depend on atom positions.
 
Copyright (C) 2008,2009 by Tim Vandermeersch
Copyright (C) 2009 by Frank Peters (some additions)
 
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

#ifndef OPENBABEL_OBFUNCTION_H
#define OPENBABEL_OBFUNCTION_H

#include <vector>
#include <Eigen/Core>

namespace OpenBabel {

  class OBMol;

namespace OBFFs {

  class OBLogFile;
  class OBFunctionTerm;
  class OBParameterDB;
  class OBFFType;
  class OBChargeMethod;

  /** @class OBFunction
   *  @brief Base class for functions (e.g. force fields, ...) of 3D variables (e.g. atom coordinates, ...).
   */
  class OBFunction 
  {
    public:
      enum Computation {
        Value,
        Gradients,
      };
      
      /**
       * Constructor
       */
      OBFunction();
      /**
       * Destructor
       */
      virtual ~OBFunction();
      /**
       * Get the name for this function ("MMFF94", "UFF", "My function 3", ...).
       */
      virtual std::string GetName() const = 0;
      /**
       * Setup this function.
       *
       * Atom positions and gradients will be handled by the OBFunction class. Subclasses can overload
       * this method to do their own setup but should always call OBFunction::Setup() to make sure the
       * positions are copied and gradients are set to zero.
       */
      virtual bool Setup(/*const*/ OBMol &mol);
      /**
       * Perform the specified OBFunction::Computation. 
       */
      virtual void Compute(Computation computation = Value) = 0;
      /**
       * Implemented by subclasses to return the current value (i.e. OBFunctionImpl).
       * Call Compute() before GetValue().
       */
      virtual double GetValue() const = 0;
      /**
       * Get the used unit (e.g. "kJ/mol", "kcal/mol") for this function.
       */
      virtual std::string GetUnit() const = 0;
      /**
       * Get the number of particles
       */
      unsigned int NumParticles() const { return m_positions.size(); }
      /** 
       * Get the atom positions.
       */
      std::vector<Eigen::Vector3d>&  GetPositions() { return m_positions; }
      const std::vector<Eigen::Vector3d>&  GetPositions() const { return m_positions; }
      /** 
       * Copy atom positions to molecule
       */
      bool CopyPositionsToMol(OBMol& mol) const;
      /** 
       * Get the atom gradients.
       */
      std::vector<Eigen::Vector3d>&  GetGradients() { return m_gradients; } 
      const std::vector<Eigen::Vector3d>&  GetGradients() const { return m_gradients; } 
      /**
       * @return True if this function has analytical gradients. 
       */
      virtual bool HasAnalyticalGradients() const { return false; }
      /**
       * Get the OBLogFile for this function.
       */
      OBLogFile* GetLogFile() { return m_logfile; }
      /**
       * Get the OBParameterDB for this function.
       */
      OBParameterDB* GetParameterDB() { return m_parameterDB; }
      /**
       * Set the OBParamterDB for this function.
       */
      void SetParameterDB(OBParameterDB *database);
      /**
       * Get the OBFFType for this function.
       */
      OBFFType* GetOBFFType() { return m_obffType; }
      /**
       * Set the OBFFType for this function.
       */
      void SetOBFFType(OBFFType *obfftype);
      /**
       * Get the OBChargeMethod for this function.
       */
      OBChargeMethod* GetOBChargeMethod() { return m_obChargeMethod; }
      /**
       * Set the OBChargeMethod for this function.
       */
      void SetOBChargeMethod(OBChargeMethod *obChargeMethod);
      /**
       * Add a term to this function.
       */
      void AddTerm(OBFunctionTerm *term);
      /**
       * Remove all terms from this function. The term instances are deleted.
       */ 
      void RemoveAllTerms(bool deleteTerms = true);
      /**
       * Get all terms (i.e. pointers to OBFunctionTerm objects) for this function.
       */
      const std::vector<OBFunctionTerm*>& GetTerms() const;

      std::string GetOptions() const;
      void SetOptions(const std::string &options);
    
      /*! Calculate the potential energy function derivative numerically with 
       *  repect to the coordinates of atom with index a (this vector is the gradient)
       *
       * \param a  provides coordinates
       * \param terms OBFF_ENERGY, OBFF_EBOND, OBFF_EANGLE, OBFF_ESTRBND, OBFF_ETORSION, 
       * OBFF_EOOP, OBFF_EVDW, OBFF_ELECTROSTATIC
       * \return the negative gradient of atom a
       */
      Eigen::Vector3d NumericalDerivative(unsigned int index);
      //! OB 3.0
      Eigen::Vector3d NumericalSecondDerivative(unsigned int index);
 
    protected:
      struct Option {
        Option(int _line, const std::string &_name, const std::string &_value) 
            : line(_line), name(_name), value(_value)
        {
        }
        int line; //!< the line number of the option
        std::string name; //!< the name of the option
        std::string value; //!< the value for the option
        std::string error; //!< subclasses can use this to report warnings/errors
      };
      virtual void ProcessOptions(std::vector<Option> &options) = 0;
      virtual std::string GetDefaultOptions() const = 0;

      OBLogFile *m_logfile;
      OBParameterDB *m_parameterDB;
      OBFFType *m_obffType;
      OBChargeMethod *m_obChargeMethod;
      std::string m_options;
      std::vector<OBFunctionTerm*> m_terms;
      std::vector<Eigen::Vector3d> m_positions;
      std::vector<Eigen::Vector3d> m_gradients;
  };

  class OBFunctionFactory
  {
    public:
      OBFunctionFactory();
      virtual std::string GetName() const = 0;
      virtual OBFunction* NewInstance() = 0;
      
      static std::vector<OBFunctionFactory*> GetFactories();
      static OBFunctionFactory* GetFactory(const std::string &name);
    protected:
      static std::vector<OBFunctionFactory*>& GetFactoriesRef();
  };

} // OBFFs
} // OpenBabel

#endif
