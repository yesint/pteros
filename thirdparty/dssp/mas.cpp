// Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
//   Distributed under the Boost Software License, Version 1.0.
//       (See accompanying file LICENSE_1_0.txt or copy at    
//             http://www.boost.org/LICENSE_1_0.txt)      
// 
// mas.cpp - simple attempt to write a multiple sequence alignment application

#include "mas.h"
#include <algorithm>

using namespace std;

// 22 real letters and 1 dummy
const char kResidues[] = "ABCDEFGHIKLMNPQRSTVWYZX";
const uint8 kResidueNrTable[] = {
//	A   B   C   D   E   F   G   H   I       K   L   M   N       P   Q   R   S   T  U=X  V   W   X   Y   Z
	0,  1,  2,  3,  4,  5,  6,  7,  8, 23,  9, 10, 11, 12, 23, 13, 14, 15, 16, 17, 22, 18, 19, 22, 20, 21
};

sequence encode(const string& s)
{
	sequence result(s.length(), 0);
	for (unsigned int i = 0; i < s.length(); ++i)
		result[i] = is_gap(s[i]) ? '-' : ResidueNr(s[i]);
	return result;
}

string decode(const sequence& s)
{
	string result(s.length(), 0);
	for (unsigned int i = 0; i < s.length(); ++i)
		result[i] = s[i] >= 23 ? '.' : kResidues[s[i]];
	return result;
}
