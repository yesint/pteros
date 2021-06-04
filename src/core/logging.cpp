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
    else throw PterosError("Invalid log level '{}'",lev);
    spdlog::set_level(Log::instance().default_level);
}




