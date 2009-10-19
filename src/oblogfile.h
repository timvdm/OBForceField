/**********************************************************************
obloggable.h - Class for handling force field logging.
 
Copyright (C) 2008,2009 by Tim Vandermeersch
 
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

#ifndef OB_LOGFILE_H
#define OB_LOGFILE_H

#include <openbabel/babelconfig.h>
#include <string>

#undef OBAPI
#define OBAPI

namespace OpenBabel {
namespace OBFFs {

  // inline if statements for logging.
#define IF_LOGLVL_LOW    if(GetLogLevel() >= LogLevel::Low)
#define IF_LOGLVL_MEDIUM if(GetLogLevel() >= LogLevel::Medium)
#define IF_LOGLVL_HIGH   if(GetLogLevel() >= LogLevel::High)

  namespace LogLevel
  {
    enum Level
    {
      None = 0,  //!< no output
      Low,       //!< SteepestDescent progress... (no output from Energy())
      Medium,    //!< individual energy terms
      High       //!< individual calculations and parameters
    };
  };

 
  class OBAPI OBLogFile
  {     
    public:
      enum LogLevel
      {
        None = 0,  //!< no output
        Low,       //!< SteepestDescent progress... (no output from Energy())
        Medium,    //!< individual energy terms
        High       //!< individual calculations and parameters
      };
 
      OBLogFile();
      /** 
       * @brief Destructor.
       */
      virtual ~OBLogFile(); 

      /////////////////////////////////////////////////////////////////////////
      // Logging                                                             //
      /////////////////////////////////////////////////////////////////////////

      //! \name Methods for logging
      //@{
      /**
       * Set the std::ostream to write log messages too.
       */
      bool SetOutputStream(std::ostream *pos);
      /** 
       * @brief Set the log level (OBFF_LOGLVL_NONE, OBFF_LOGLVL_LOW, OBFF_LOGLVL_MEDIUM, OBFF_LOGLVL_HIGH).
       * Inline if statements for logging are available: 
       * @code
       * #define IF_OBFF_LOGLVL_LOW    if(GetLogLevel() >= OBFF_LOGLVL_LOW)
       * #define IF_OBFF_LOGLVL_MEDIUM if(GetLogLevel() >= OBFF_LOGLVL_MEDIUM)
       * #define IF_OBFF_LOGLVL_HIGH   if(GetLogLevel() >= OBFF_LOGLVL_HIGH)
       * @endcode
       *
       * example:
       * @code
       * SetLogLevel(OBFF_LOGLVL_MEDIUM);
       * IF_OBFF_LOGLVL_HIGH {
       *   OBFFLog("this text will NOT be logged...\n");
       * }
       *
       * IF_OBFF_LOGLVL_LOW {
       *   OBFFLog"this text will be logged...\n");
       * }
       *
       * IF_OBFF_LOGLVL_MEDIUM {
       *   OBFFLog("this text will also be logged...\n");
       * }
       * @endcode
       */
      bool SetLogLevel(LogLevel level) { m_loglvl = level; }
      /** 
       * @return The log level.
       */ 
      int GetLogLevel() const { return m_loglvl; }
      /** 
       * @brief Print msg to the logfile regardless of current loglevel.
       * 
       * @param msg The message to print.
       */
      void Write(const std::string &msg);
      /** 
       * @brief Print msg to the logfile if the current loglevel is .
       * 
       * @param lvl The minimum required loglevel.
       * @param msg The message to print.
       */
      void Write(LogLevel lvl, const std::string &msg);
      std::ostream& operator()() { return *m_logos; }
      //@}
      
      bool IsNone() const { return (m_loglvl >= None) ? true : false; }
      bool IsLow() const { return (m_loglvl >= Low) ? true : false; }
      bool IsMedium() const { return (m_loglvl >= Medium) ? true : false; }
      bool IsHigh() const { return (m_loglvl >= High) ? true : false; }

    protected:
      // some data in class for fast inline Get/Set functions
      LogLevel      m_loglvl; //!< Log level for output
      std::ostream* m_logos; //!< Output for logfile 
  }; // class OBForceField

}// namespace OpenBabel
}

#endif   // OB_FORCEFIELD_H

//! @file obloggable.h
//! @brief OBLoggable class
