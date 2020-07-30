/**
 *  @file Reaction.h
 *
 */

// Copyright 2001  California Institute of Technology


#ifndef CKR_REACTION_H
#define CKR_REACTION_H

#include <string>
#include <vector>
#include <map>
using namespace std;

#include "ckr_defs.h"
#include "ckr_utils.h"
#include "RxnSpecies.h"


namespace ckr {


    /** @name Reaction Types
     */

    //@{
    const int Elementary = 2000;  ///< elementary, pressure-independent reaction
    const int ThreeBody = 2001;   ///< three-body reaction
    const int Falloff = 2002;     ///< falloff reaction
    const int ChemAct = 2003;     ///< chemical activation reaction

#ifdef ENABLE_PLASMA_REACTIONS
    const int CrossSection = 2004;   ///< momentum transfer cross section
    const int CollisionFreq = 2005;  ///< momentum transfer collision frequency
#endif

    const int PLogInterpolation = 10000;
    //@}

    /**
     *  Reaction rate coefficient class. For simplicity, the RateCoeff
     *  class has public members representing the parameters of all rate
     *  coefficient types supported by Chemkin. The 'type' attribute
     *  specifies the rate coefficient type.
     */
    class RateCoeff {
    public:
        
        /**
         *  Constructor. Construct a default rate coefficient object. The
         *  type is set to Arrhenius, and all parameter values are set to 
         *  zero.
         */
    RateCoeff() : A(-1.0e+300), n(0.0), E(0.0), // MJM mod
            B(0.0), C(0.0), type(Arrhenius) {}

        /**
         *  Copy constructor
         */
        RateCoeff(const RateCoeff& k) : A(k.A), n(k.n), E(k.E), 
	B(k.B), C(k.C), type(k.type), b(k.b), A_plog(k.A_plog),
        n_plog(k.n_plog), E_plog(k.E_plog), pres_pts_plog(k.pres_pts_plog) {}

        /// Destructor. Does nothing.
        ~RateCoeff() {}

        RateCoeff& operator=(const RateCoeff& k) {
            if (this == &k) return *this;
            A = k.A;
            n = k.n;
            E = k.E;
            B = k.B;
            C = k.C;
            type = k.type;
            b = k.b;
            A_plog = k.A_plog;
            n_plog = k.n_plog;
            E_plog = k.E_plog;
            pres_pts_plog = k.pres_pts_plog;
            return *this;
        }

        // Modified Arrhenius parameters, common to all types.
        double A;          ///< pre-exponential factor (all types)
        double n;          ///< temperature exponent   (all types)
        double E;          ///< activation energy      (all types)

        // Landau-Teller parameters
        double B;          ///< Landau-Teller B parameter
        double C;          ///< Landau-Teller C parameter

        int type;          ///< rate coefficient type

        // Coefficients for the special forms allowed by Chemkin-III
        vector_fp b;  ///< coefficients for JAN or FIT1 form

        // coefficients for the pressure logarithmic interpolation
        vector_fp A_plog;
        vector_fp n_plog;
        vector_fp E_plog;
        vector_fp pres_pts_plog; 
    };


    //////////////////////////////////////////////////////////////////////////
    
    /**
     *   Specifies the units for all reactions.  
     */
    
    class ReactionUnits {
    public:
        int ActEnergy;    ///< Activation energy unit flag 
        int Quantity;     ///< Moles or molecules unit flag
    };
    

    //////////////////////////////////////////////////////////////////////////


    /// A class for reactions.

    // Note: if you add data items to this class, be sure to update
    // the copy constructor and the assignment operator !

    class Reaction {
    public:
        
        /// a list of auxiliary data values
        typedef vector_fp auxdata;
        
        /// Construct an empty Reaction object
        Reaction() : type(Elementary), 
                     isFalloffRxn(false), 
                     isChemActRxn(false),
                     isThreeBodyRxn(false), 
                     isReversible(false),
                     isDuplicate(false),
	             isPLogInterpolation(false), 
                     duplicate(0), 
                     thirdBody("<none>"), 
                     number(0),
                     falloffType(Lindemann) {} 

        /// Copy constructor
        Reaction(const Reaction& r) : type(r.type), 
                                    isFalloffRxn(r.isFalloffRxn), 
                                    isChemActRxn(r.isChemActRxn),
                                    isThreeBodyRxn(r.isThreeBodyRxn), 
                                    isReversible(r.isReversible),
                                    isDuplicate(r.isDuplicate),
	                            isPLogInterpolation(r.isPLogInterpolation), 
                                    duplicate(r.duplicate),
                                    thirdBody(r.thirdBody), 
                                    number(r.number), 
                                    reactants(r.reactants), 
                                    fwdOrder(r.fwdOrder),
                                    products(r.products),
                                    e3b(r.e3b),
                                    kf(r.kf), 
                                    kf_aux(r.kf_aux), 
                                    krev(r.krev), 
                                    falloffType(r.falloffType), 
                                    falloffParameters(r.falloffParameters),
                                    otherAuxData(r.otherAuxData),
                                    lines(r.lines), comment(r.comment) {}

        /// Destructor
        virtual ~Reaction() {}
        
        Reaction& operator=(const Reaction& r);

        int type;                 ///< Reaction type.
        
        bool isFalloffRxn;        ///< True if reaction is a falloff reaction.
        bool isChemActRxn;        ///< True if reaction is a chemical activation reaction.
        bool isThreeBodyRxn;      ///< True if reaction is a three-body reaction.
        bool isReversible;        ///< True if reaction is reversible.
        bool isDuplicate;         ///< True if reaction is declared to be a duplicate;
        bool isPLogInterpolation; ///< True if reaction is declared to be a
                                  ///< pressure interpolation reaction
        /**
         * reaction number this one is a duplicate to (declared or not). 
         * If the reaction is not a duplicate, the value is zero.
         */
        int duplicate;
        
        /**
         *  For pressure-dependent reactions (including three-body ones)
         *  this string contains either "M" if all species may act as
         *  third body collision partners, or a species name if only 
         *  one species does.
         */
        string thirdBody;

        /// Reaction number.
        int number;
        
        /** 
         * list of species that participate as reactants,
         * and their stoichiometric coefficients
         */
        vector<RxnSpecies> reactants;
        
        mutable map<string, double> fwdOrder;

        /** 
         * list of species that participate as products,
         * and their stoichiometric coefficients
         */
        vector<RxnSpecies> products;
        
        
        /**
         * map from species names to enhanced third-body collision efficiencies
         */
        mutable map<string, double> e3b;
        
        /**
         *  Forward rate coefficient. For falloff reactions, this is the 
         *  high-pressure rate coefficient, and for chemical activation
         *  reactions it is the low-pressure one.
         */
        RateCoeff kf;
        
        /**
         *  For pressure-dependent reactions, the rate coefficient for the
         *  opposite pressure limit as kf (
         */
        RateCoeff kf_aux;
        
        /// Reverse rate coefficient. Empty unless REV auxiliary data given.
        RateCoeff krev;


        int falloffType;        
        vector_fp falloffParameters;
        
        /**
         * auxiliary data not handled elsewhere.
         */
        mutable map<string, auxdata> otherAuxData;

        /**
         * input file lines
         */ 
        vector<string> lines;

        /**
         * comments
         */
        vector<string> comment;

        // methods

        double stoichCoefficient(const string& s) const;        
        bool operator==(const Reaction& r) const;
        void write(ostream& s) const;

    };

    /// a list of Reaction objects
    typedef vector<Reaction>     reactionList;

    Reaction forwardReaction(const Reaction& rxn);
    Reaction reverseReaction(const Reaction& rxn);
}




#endif





