/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2017, Semen Yesylevskyy
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of Artistic License:
 *
 * Please note, that Artistic License is slightly more restrictive
 * then GPL license in terms of distributing the modified versions
 * of this software (they should be approved first).
 * Read http://www.opensource.org/licenses/artistic-license-2.0.php
 * for details. Such license fits scientific software better then
 * GPL because it prevents the distribution of bugged derivatives.
 *
*/

#ifndef LOGGING_H
#define LOGGING_H

#include "spdlog/spdlog.h"
#include <iostream>

namespace pteros {

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

private:
  Log() {
      generic_pattern = "[%n]\t(%l)\t%v";
      try {
          console_sink = std::make_shared<spdlog::sinks::stdout_sink_st>();
          logger = std::make_shared<spdlog::logger>("pteros", console_sink);
          logger->set_pattern(generic_pattern);
      } catch (const spdlog::spdlog_ex& ex) {
          std::cout << "Log initialization failed: " << ex.what() << std::endl;
      }
  }
};

}

#define LOG() (Log::instance().logger)

#endif
