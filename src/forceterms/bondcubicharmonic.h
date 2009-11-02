#include <OBFunctionTerm>

namespace OpenBabel {
  namespace OBFFs {

    /**
     * The BondCubicHarmonicTerm is a bond stretching term using a quadratic function of the form
     *
     * \f[ E_{bond} = prefactor \times kb_{ij} \times \Delta r_{ij}^2 \times (1 + cs \times \Delta r_{ij} + cs2 \times \Delta r_{ij}^2) \f]
     * 
     * where \f$ \Delta r_{ij} = |r_{ij} - r_0|\f$. The \f$prefactor\f$, \f$cs\f$ and \f$cs2\f$ can be specified in the constructor.
     * The force constant \f$kb_{ij}\f$ and bond length \f$r_0\f$ are acquired from the database using the @p tableName, @p forceConstantColumn 
     * and @p bondLengthColumn specified in the constructor. 
     */
    class BondCubicHarmonicTerm : public OBFunctionTerm
    {
    public:
      struct Index
      {
	unsigned int iA, iB;
      };
      struct Parameter
      {
	double K, r0;
      };
      /**
       * Constructor specifying all the options for the term.
       *
       * @param function The function to which this term belongs.
       * @param prefactor The prefactor used in the function.
       * @param cs The cubic stretch constant used in the function.
       * @param cs The squared cubic stretch constant used in the function.
       * @param tableName The database table name from which to get the parameters.
       * @param forceConstantColumn The column in the database table contaning the bond stretching force constant (i.e. \f$kb_{ij}\f$).
       * @param bondLengthColumn The column in the database table containing the reference bond length (i.e. \f$r_0\f$).
       */
      BondCubicHarmonicTerm(OBFunction *function, double prefactor, double cs, double cs2, const std::string &tableName, int forceContantColumn, int bondLengthColumn);
      /**
       * Destructor.
       */
      ~BondCubicHarmonicTerm();
      /**
       * Get the name for this function term.
       */
      std::string GetName() const 
      { 
        return "Bond Cubic Harmonic";
      }
      /**
       * Setup this term.
       */
      bool Setup();
      /**
       * Perform the specified computation.
       */
      void Compute(OBFunction::Computation computation = OBFunction::Value);
      /**
       * Get the last computed value. Make sure to call Compute() before using GetValue().
       */
      double GetValue() const 
      { 
        return m_value; 
      }
    private:
      const std::string m_tableName; //!< The database table name (e.g. "Bond Parameters")
      const int m_forceConstantColumn; //!< The database table column containing \f$kb_{ij}\f$
      const int m_bondLengthColumn; //!< The database table column containing \f$r_0\f$
      unsigned int m_numBonds; //!< The number of bonds
      Parameter *m_calcs; //!< The parameters for the bond calculations
      Index *m_i; //!< The atom indexes for each bond
      double m_prefactor; //!< The \f$prefactor\f$ (see formula)
      double m_cs; //!< The \f$cs\f$ (see formula)
      double m_cs2; //!< The \f$cs2\f$ (see formula)
      double m_value; //!< The last computed value
    };

  } // OBFFs
} // OpenBabel
