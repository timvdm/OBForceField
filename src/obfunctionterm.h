#ifndef OBFFS_FUNCTIONTERM_H
#define OBFFS_FUNCTIONTERM_H

#include <vector>

#include <OBVariant>
#include <OBFunction>

namespace OpenBabel {
namespace OBFFs {

  /**
   * 
   * @todo: explain rows & colums
   *
   *
   *
   */
  class OBFunctionTerm
  {
    public:
      /**
       * Constructor.
       */
      OBFunctionTerm(OBFunction *function);
      virtual ~OBFunctionTerm();
      /*
       * Get the name for this function term (e.g. "MMFF94 bond stretching", 
       * "UFF angle bending", ...).
       */
      virtual std::string GetName() const = 0;
      /**
       * Setup this term
       */
      virtual bool Setup(/*const*/ OBMol &mol) = 0;
      /**
       * Compute the value or gradients for this term.
       */
      virtual void Compute(OBFunction::Computation computation = OBFunction::Value) = 0;
      /**
       * Implemented by subclasses to return the current value (i.e. OBFunctionImpl). 
       * Call Compute() before GetValue().
       */
      virtual double GetValue() const = 0;
 
      /**
       * Get the the parameter data base for this term.
       */
      //virtual OBParameterDB* GetParameterDB() const = 0;

    protected:
      OBFunction *m_function;
  };

} // OBFFs
} // OpenBabel

#endif
