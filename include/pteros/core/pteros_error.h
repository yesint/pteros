/*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   Pteros   ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2009-2013, Semen Yesylevskyy
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

#ifndef PERROR_H
#define PERROR_H

#include "spdlog/fmt/fmt.h"

namespace pteros {

/// Represents an error in the Pteros code. Used for all Pteros-related exceptions.
class Pteros_error {
public:

    Pteros_error(const Pteros_error& p){
        text = p.text;
    }

    Pteros_error(): text("") { }

    /// Constructs an exception object with text message
    Pteros_error(std::string s): text(s) { }

    template <typename... Args>
    Pteros_error(const char *format, const Args & ... args) {
        text = fmt::format(format, args...);
    }


    /// Return error message as string
    std::string what() const {
        return text;
    }

private:
    std::string text;
};

}
#endif
