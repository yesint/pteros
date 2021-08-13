/*
 * This file is a part of
 *
 * ============================================
 * ###   Pteros molecular modeling library  ###
 * ============================================
 *
 * https://github.com/yesint/pteros
 *
 * (C) 2009-2021, Semen Yesylevskyy
 *
 * All works, which use Pteros, should cite the following papers:
 *
 *  1.  Semen O. Yesylevskyy, "Pteros 2.0: Evolution of the fast parallel
 *      molecular analysis library for C++ and python",
 *      Journal of Computational Chemistry, 2015, 36(19), 1480–1488.
 *      doi: 10.1002/jcc.23943.
 *
 *  2.  Semen O. Yesylevskyy, "Pteros: Fast and easy to use open-source C++
 *      library for molecular analysis",
 *      Journal of Computational Chemistry, 2012, 33(19), 1632–1636.
 *      doi: 10.1002/jcc.22989.
 *
 * This is free software distributed under Artistic License:
 * http://www.opensource.org/licenses/artistic-license-2.0.php
 *
*/


#pragma once

#include "spdlog/spdlog.h"
#include "fmt/ostream.h"
#include "spdlog/sinks/stdout_sinks.h"
#include <iostream>

namespace pteros {

decltype(spdlog::details::registry::instance()) get_registry();

class Log
{
public:
  static Log& instance()
  {
    static Log *instance = new Log();    
    return *instance;
  }

  std::shared_ptr<spdlog::logger> logger;
  std::shared_ptr<spdlog::sinks::stdout_sink_st> console_sink;
  std::string generic_pattern;
  spdlog::level::level_enum default_level;

private:
  Log() {
      generic_pattern = "[%n]\t(%l)\t%v";
      default_level = spdlog::level::info;
      try {
          console_sink = std::make_shared<spdlog::sinks::stdout_sink_st>();
          logger = std::make_shared<spdlog::logger>("pteros", console_sink);
          logger->set_pattern(generic_pattern);
          logger->set_level(default_level);
          spdlog::register_logger(logger);
      } catch (const spdlog::spdlog_ex& ex) {
          std::cout << "Log initialization failed: " << ex.what() << std::endl;
      }
  }
};

std::shared_ptr<spdlog::logger> create_logger(const std::string& name);

void set_log_level(const std::string& lev);

}

#define LOG() (Log::instance().logger)



