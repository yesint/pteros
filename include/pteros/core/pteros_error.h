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

#include <sstream>
#include <iostream>
#include <exception>

namespace pteros {

using namespace std;

/// Represents an error in the Pteros code. Used for all Pteros-related exceptions.
class Pteros_error {
public:

    Pteros_error(const Pteros_error& p){
        text.str(p.text.str());
    }

    Pteros_error(){
        text.str("");
    }

    /// Constructs an exception object with text message
    Pteros_error(string s){
        text.str(s);
    }

    /** Operator << allows constructing error strings on the fly
         like this:
            \code throw Pteros_error("Wrong number ") << 5
                    << " should be between "
                    << 7 << " and " << 10;
            \endcode
        */
    template<class T>
    Pteros_error& operator<<(T data){
        text << data; //Collect data
        return *this;
    };

    /// Print error message
    void print() const {
        cout << endl << "PTEROS terminated due to the following error:"
             << endl << text.str() << endl;
    }

    /// Return error message as string
    std::string what() const {
        return text.str();
    }

private:
    std::stringstream text;
};

}
#endif
