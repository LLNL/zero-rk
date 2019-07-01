
//spify Example:
//
//  See README and appExampleIFP.cpp for
//  details.
//
//  The structure of this file should not require
//  modification.  Only the name of the class.

#ifndef zerork_reactorIFP_H
#define zerork_reactorIFP_H

#include "spify/spify.h"

class zerork_reactorIFP : public spify::parser
{
    public:
        zerork_reactorIFP(const std::string& fileName);
        ~zerork_reactorIFP(){};
    private:
        zerork_reactorIFP(); //disallow default constructor
    public:
        const bool& cuda() const;
        const int& cudaDevice() const;
        const std::string& mechFile() const;
        const std::string& thermFile() const;
        const std::string& mechLogFile() const;
        const std::string& outFile() const;
        const int& nReactors() const;
        const int& nReactorsMax() const;
        const int& nReactorsMin() const;
        const int& nMatrixReactors() const;
        const int& nMatrixThreads() const;
        const std::map<std::string,double>& fuelComp() const;
        const std::map<std::string,double>& oxidizerComp() const;
        const std::string& reactorType() const;
        const double& deltaTign() const;
        const double& reactorTime() const;
        const double& printDt() const;
        const double& maxDtInternal() const;
        const int& maxSteps() const;
        const double& relTol() const;
        const double& absTol() const;
        const double& refTemp() const;
        const double& TMin() const;
        const double& TMax() const;
        const double& pMin() const;
        const double& pMax() const;
        const double& phiMin() const;
        const double& phiMax() const;
        const double& egrMin() const;
        const double& egrMax() const;
        const double& precThresh() const;
        const int& krylovDim() const;
        const bool& oneStepMode() const;
        const bool& printAllSteps() const;
        const bool& doFakeUpdate() const;
        const bool& doILU() const;
        const int& precThreshType() const;
        const double& partialPivotThresh() const;
        const int& permutationType() const;
        const double& permThresh() const;
        const bool& strictSamePattern() const;
        const double& maxGammaOrderChange() const;
        const bool& cvStabLimDet() const;
        const int& maxNumNonLinIters() const;
        const int& maxOrd() const;
        const double& cvEpsLin() const;
        const double& cvNlConvCoeff() const;
        const bool& constantPressure() const;
        const bool& reInitOnPrint() const;
        const bool& dumpStateOnPrint() const;
        const std::vector<std::string>& stateFiles() const;
        const std::string& linearSolver() const;
        const bool& cusolverRf() const;
        const int& CUSOLVER_FACTOR_ALG() const;
        const int& CUSOLVER_SOLVE_ALG() const;
        const bool& cusolverSp() const;
        const double& logiA() const;
        const double& logiK() const;
        const double& logiQ() const;
        const double& logiM() const;
        const double& logiB() const;
        const double& logirNu() const;
        const double& min_scaled_dt() const;

    private:
        bool m_cuda;
        int m_cudaDevice;
        std::string m_mechFile;
        std::string m_thermFile;
        std::string m_mechLogFile;
        std::string m_outFile;
        int m_nReactors;
        int m_nReactorsMax;
        int m_nReactorsMin;
        int m_nMatrixReactors;
        int m_nMatrixThreads;
        std::map<std::string,double> m_fuelComp;
        std::map<std::string,double> m_oxidizerComp;
        std::string m_reactorType;
        double m_deltaTign;
        double m_reactorTime;
        double m_printDt;
        double m_maxDtInternal;
        int m_maxSteps;
        double m_relTol;
        double m_absTol;
        double m_refTemp;
        double m_TMin;
        double m_TMax;
        double m_pMin;
        double m_pMax;
        double m_phiMin;
        double m_phiMax;
        double m_egrMin;
        double m_egrMax;
        double m_precThresh;
        int m_krylovDim;
        bool m_oneStepMode;
        bool m_printAllSteps;
        bool m_doFakeUpdate;
        bool m_doILU;
        int m_precThreshType;
        double m_partialPivotThresh;
        int m_permutationType;
        double m_permThresh;
        bool m_strictSamePattern;
        double m_maxGammaOrderChange;
        bool m_cvStabLimDet;
        int m_maxNumNonLinIters;
        int m_maxOrd;
        double m_cvEpsLin;
        double m_cvNlConvCoeff;
        bool m_constantPressure;
        bool m_reInitOnPrint;
        bool m_dumpStateOnPrint;
        std::vector<std::string> m_stateFiles;
        std::string m_linearSolver;
        bool m_cusolverRf;
        int m_CUSOLVER_FACTOR_ALG;
        int m_CUSOLVER_SOLVE_ALG;
        bool m_cusolverSp;
        double m_logiA;
        double m_logiK;
        double m_logiQ;
        double m_logiM;
        double m_logiB;
        double m_logirNu;
        double m_min_scaled_dt;

};

#endif
