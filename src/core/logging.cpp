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

#include "pteros/core/logging.h"
#include "pteros/core/pteros_error.h"

std::shared_ptr<spdlog::logger> pteros::create_logger(const std::string &name)
{
    auto log = std::make_shared<spdlog::logger>(name, Log::instance().console_sink);
    log->set_pattern(Log::instance().generic_pattern);    
    log->set_level(Log::instance().default_level);
    get_registry().register_logger(log);
    return log;
}


decltype(spdlog::details::registry::instance()) pteros::get_registry(){return spdlog::details::registry::instance();}

void pteros::set_log_level(const std::string &lev)
{
    if(lev=="trace") Log::instance().default_level = spdlog::level::trace;
    else if(lev=="debug") Log::instance().default_level = spdlog::level::debug;
    else if(lev=="info") Log::instance().default_level = spdlog::level::info;
    else if(lev=="warn") Log::instance().default_level = spdlog::level::warn;
    else if(lev=="err") Log::instance().default_level = spdlog::level::err;
    else if(lev=="critical") Log::instance().default_level = spdlog::level::critical;
    else if(lev=="off") Log::instance().default_level = spdlog::level::off;
    else throw Pteros_error("Invalid log level '{}'",lev);
    spdlog::set_level(Log::instance().default_level);
}
