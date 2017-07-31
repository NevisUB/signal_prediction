/**
 * \file sp_base.h
 *
 * \ingroup SP
 * 
 * \brief Class definition file of sp_base
 *
 * @author Kazu - Nevis 2015
 */

/** \addtogroup SP

    @{*/

#ifndef __SP_BASE_H__
#define __SP_BASE_H__

#include <vector>
#include "sp_logger.h"

namespace sp {
    
  class sp_base {
    
  public:
    
    sp_base(const std::string logger_name="sp_base")
      : _logger(nullptr)
    { _logger = &(::sp::logger::get(logger_name)); }
    
    sp_base(const sp_base &original) : _logger(original._logger) {}
    
    virtual ~sp_base(){};

    inline const sp::logger& logger() const
    { return *_logger; }
    
    void set_verbosity(::sp::msg::Level_t level)
    { _logger->set(level); }

    const std::string& name() const
    { return logger().name(); }
    
  private:
    
    sp::logger *_logger;   //! don't save
    
  };
}
#endif

/** @} */ // end of doxygen group
