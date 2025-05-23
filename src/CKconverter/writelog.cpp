/**
 *  @file converters/writelog.cpp
 *
 */

// Copyright 2001  California Institute of Technology


// turn off warnings about truncating long names under Windows
#ifdef _WIN32
#pragma warning(disable:4786)
#endif

#include <stdio.h>
#include "writelog.h"

namespace ckr {


    /// format a string for the log file message printed when starting a new task
    string newTask(string msg) {
        string s = "\n";
        s += msg + "...\n";
        return s;
    }


    /// print the falloff parameters for a pressure-dependent reaction
    bool writeFalloff(int type, const vector_fp& c, ostream& log) {

        int nonStdFlag = 0;
        log.precision(6);
        log.width(0);
        log.flags(ios::uppercase);

        //    bool ok = true;
        switch (type) {

        case Lindemann:
            log << "   Lindemann falloff function" << endl;
            return true;

        case Troe:

            log << "   Troe falloff function: " << endl;
            if (c.size() == 3) {
                log << "      alpha, T***, T* = (" << c[0] << ", " << c[1] 
                    << ", " << c[2] << ")" << endl;
            }
            else if (c.size() == 4) {
                log << "      alpha, T***, T*, T** = (" << c[0] << ", " << c[1] 
                    << ", " << c[2] << ", " << c[3] << ")" << endl;
            }
            else {
				for (size_t n = 0; n < c.size(); n++) {
					log << c[n] << ", "; log << endl;
				}
                log << "###### ERROR #####   incorrect number of parameters" << endl;
                return false;
            }

            if(c[0] <= 0.0 || c[0] >= 1.0) {
              nonStdFlag += 1;
            }
            if(c[1] <= 0.0) {
              nonStdFlag += 2;
            }
            if(c[2] <= 0.0) {
              nonStdFlag += 4;
            }
            if(c.size() == 4) {
              if(c[3] <= 0.0) {
                nonStdFlag += 8;
              }
            }
            if(nonStdFlag > 0) {
              log << "# INFO: Reaction may not have the intended centering function." << endl;
              log << "#       Troe centering fucntion Fc(T) being used:" << endl;
              log << "#           Fc(T) = (" << 1.0-c[0] << ")exp(-T/(" << c[1]<< "))" << endl;
              log << "#                   +(" << c[0] << ")exp(-T/(" << c[2] <<"))" << endl;
              if(c.size() == 4) {
                log << "#                 +exp(" << -c[3] << "/T)" << endl;
              }  
              log << "#" << endl;
              log << "#     In some mechanisms, non-positive parameters were" << endl
                  << "#     used as flags for alternative forms for the" << endl
                  << "#     centering function Fc.  Please verify that the above" << endl
                  << "#     definition is consistent with your expectations." << endl;
            }
            return true;

        case SRI:
            log << "   SRI falloff function: " << endl;
            if (c.size() == 3) {
                log << "      a, b, c = (" << c[0] << ", " << c[1] 
                    << ", " << c[2] << ")" << endl;
            }
            else if (c.size() == 5) {
                log << "      a, b, c, d, e = (" << c[0] << ", " << c[1] 
                    << ", " << c[2] << ", " << c[3] << ", " << c[4] 
                    << ")" << endl;
            }
            else {
				for (size_t n = 0; n < c.size(); n++) {
					log << c[n] << ", "; log << endl;
				}
                log << "##### ERROR #####  incorrect number of parameters" << endl;
                return false;
            }
            return true;

        default:
            log << "unknown falloff type: " << type << endl;
            return false;
        }
    }


    /// print the rate coefficient parameters
    bool writeRateCoeff(const RateCoeff& k, ostream& log) {

        log.precision(10);
        log.width(0);
        log.flags(ios::uppercase);
        int n;

        bool ok = true;
        int nb;

        switch (k.type) {

        case Arrhenius:
	  //printf("Arrhenius type?\n");
            log <<" A, n, E = (" << k.A << ", " << k.n 
                << ", " << k.E << ")" << endl;
            break;
    
        case LandauTeller:
            log << "A, n, E, B, C = (" << k.A << ", " << k.n
                << ", " << k.E << ", " << k.B << ", " << k.C 
                << ") *** Landau-Teller ***" << endl;
            break;
    
        case Jan:
            log <<" A, n, E = (" << k.A << ", " << k.n 
                << ", " << k.E << ") *** JAN *** " << endl;
            nb = k.b.size();
            for (n = 0; n < nb; n++) {
                log << "   b" << n+1 << ": " << k.b[n] << endl;
            }
            if (nb != 9) log 
                      << "warning: number of b coefficients should be 9." 
                      << endl;
            break;

        case Fit1:
            log <<" A, n, E = (" << k.A << ", " << k.n 
                << ", " << k.E << ") *** FIT1 *** " << endl;
            nb = k.b.size();
            for (n = 0; n < nb; n++) {
                log << "   b" << n+1 << ": " << k.b[n] << endl;
            }
            if (nb != 9) log 
                      << "warning: number of b coefficients should be 4." 
                      << endl;
            break;
        case PLogInterpolation:
          log << "   Default rate coeff: A, n, E = (" << k.A << ", " << k.n 
	      << ", " << k.E << ") only used if PLOG is ignored" << endl;
          nb = k.pres_pts_plog.size();
          for(n=0; n<nb; ++n) {
            log << "   P = " << k.pres_pts_plog[n] 
                << " atm rate coeff: A, n, E = (" << k.A_plog[n] << ", " 
                << k.n_plog[n] << ", " << k.E_plog[n] << ")" << endl;
          }
	  break;

        default:
            log << "unknown rate coefficient type: " << k.type << endl;
            ok = false;
        }
        return ok;
    }

    /**
     * Write onto an output stream the chemical equation for a reaction.
     */
    void printReactionEquation(ostream& f, const Reaction& r) {
        //        r.write(f);
        f << reactionEquation(r);
    }


    /**
     * Write to a string the chemical equation for a reaction.
     */
    string reactionEquation(const Reaction& r) {
        string s = "";
        int nr = static_cast<int>(r.reactants.size());
        int np = static_cast<int>(r.products.size());
        int k;
        double m;
        char buf[20];

        for (k = 0; k < nr; k++) {
            m =  r.reactants[k].number;
            if (m != 1.0) {
                sprintf(buf,"%g",m);
                s += string(buf);
                s += " ";
            }
            s +=  r.reactants[k].name;
            if (k < nr - 1) s += " + ";
        }

        if (r.isFalloffRxn) s += " (+ " + r.thirdBody + ")";
        else if (r.isThreeBodyRxn) s += " + " + r.thirdBody;
        if (r.isReversible) s += " <=> ";
        else s += " => ";

        for (k = 0; k < np; k++) {
            m =  r.products[k].number;
            if (m != 1.0) {
                sprintf(buf,"%g",m);
                s += string(buf);
                s += " ";
            }
            s += r.products[k].name;
            if (k < np - 1) s += " + ";
        }
        if (r.isFalloffRxn) s += " (+ " + r.thirdBody + ")";
        else if (r.isThreeBodyRxn) s += " + " + r.thirdBody;
        return s;
    }



    /** 
     * 
     * Write a summary of the properties of one species to the log file.
     *  @param log     log file output stream
     *  @param spec    instance of Species class
     */

    void writeSpeciesData(ostream& log, const Species& spec) {
        
        if (!spec.id.empty()) 
            log << endl << "   id/date: " << spec.id << endl;
        else 
            log << " ... " << endl;
        
        log << "   phase: " 
            << spec.phase << endl 
            << "   composition: (";
        
        for (size_t ie = 0; ie < spec.elements.size(); ie++) {
            if (!spec.elements[ie].name.empty()) {
                log.flags(ios::fixed);
                log.precision(0);
                if (ie > 0) log << ", ";
                log << spec.elements[ie].number << " " 
                    << spec.elements[ie].name;
            }
        } 
        log << ")";

        log.flags(ios::showpoint | ios::fixed);
        log.precision(2);
        log << endl << "   Tlow, Tmid, Thigh: (" << spec.tlow << ", " << 
            spec.tmid << ", " << spec.thigh << ")" << endl << endl;
        log << "   coefficients (low, high):" << endl;
        log.flags(ios::scientific | ios::uppercase | ios::internal);
        log.precision(8);
        for (int j = 0; j < 7; j++) { 
            log << "   a" << j + 1;
            log.setf(ios::showpos);
            log << "  \t" << spec.lowCoeffs[j] 
                << "  \t" << spec.highCoeffs[j] << endl;
            log.unsetf(ios::showpos);
        }
        log << endl;
    }

    /**
     * Write to the log file the real order rate of progress.
     */
    void writeRealOrder(const Reaction& rxn, ostream& log) {
        
        log.precision(15);
        log.width(0);
        log.flags(ios::uppercase);

        const int nreactants = rxn.reactants.size();
        const int nproducts  = rxn.products.size();
        const int nford      = rxn.fwdOrder.size();
        const int nrord      = rxn.revOrder.size();

        std::map<string, double>::iterator it;
        log << "    forward reaction concentration powers:" << endl;

        // report all the specified real-order coefficients
        for(it=rxn.fwdOrder.begin(); it != rxn.fwdOrder.end(); it++) {
            log << "        " << "[" << it->first << "]**"
                << it->second << "    (FORD)" << endl;
        }

        // check the reactant list to see if any species are not found in
        // the fwdOrder map
        for(int j=0; j<nreactants; ++j) {

            it = rxn.fwdOrder.find(rxn.reactants[j].name);
            if(it == rxn.fwdOrder.end()) {
                // reactant species not found in FORD list
                log << "        [" << rxn.reactants[j].name << "]**"
                    << rxn.reactants[j].number << "    (STOIC)" << endl;

            }
        }

        if(rxn.isReversible) {
            log << "    reverse reaction concentration powers:" << endl;

            // report all the specified real-order coefficients
            for(it=rxn.revOrder.begin(); it != rxn.revOrder.end(); it++) {
                log << "        " << "[" << it->first << "]**"
                    << it->second << "    (RORD)" << endl;
            }

            // check the reactant list to see if any species are not found in
            // the revOrder map
            for(int j=0; j<nproducts; ++j) {

                it = rxn.revOrder.find(rxn.products[j].name);
                if(it == rxn.revOrder.end()) {
                    // product species not found in RORD list
                    log << "        [" << rxn.products[j].name << "]**"
                        << rxn.products[j].number << "    (STOIC)" << endl;

                }
            }

            // TODO: Add equilibrium check/warning 
        }
    }


}


