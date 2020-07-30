/**
 *  @file RxnSpecies.h
 *
 */

// Copyright 2001  California Institute of Technology


#ifndef CKR_RXNSPECIES_H
#define CKR_RXNSPECIES_H

#include <string>
#include <vector>
//#include "Cantera.h"

using namespace std;

namespace ckr {

    typedef vector_int        group_t;
    typedef vector<group_t>   grouplist_t; 


/**
 * A class for species in a reaction.
 *
 */
class RxnSpecies {
public:
    string         name;        //!< The name of the object.
    double         number;      //!< The number of units (molecules, etc.).
    grouplist_t    groups;
};

}

#endif
