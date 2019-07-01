
#include "zerork_reactorIFP.h"

//See README for details on implementing new parser.

zerork_reactorIFP::zerork_reactorIFP(const std::string& fileName)
: spify::parser(fileName)
{
    addParameter
    (
        (new spify::scalar_param<bool>("cuda"))
        ->defaultValue(1)
        ->shortDesc("Enable CUDA?")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<int>("cudaDevice"))
        ->defaultValue(-1)
        ->shortDesc("CUDA Device to use for Simulation")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<std::string>("mechFile"))
        ->shortDesc("Chemkin Format Mechansim File")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<std::string>("thermFile"))
        ->shortDesc("Chemkin Format Thermodynamics File")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<std::string>("mechLogFile"))
        ->defaultValue("/dev/null")
        ->shortDesc("Mechanism Parser Log File")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<std::string>("outFile"))
        ->shortDesc("Ignition Delay Output File")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<int>("nReactors"))
        ->boundMin(1)
        ->shortDesc("Number of reactors to solve")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<int>("nReactorsMax"))
        ->defaultValue(8192)
        ->boundMax(8192)
        ->boundMin(1)
        ->shortDesc("Maximum number of reactors to solve")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<int>("nReactorsMin"))
        ->defaultValue(256)
        ->boundMin(1)
        ->shortDesc("Minimum number of reactors to solve")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<int>("nMatrixReactors"))
        ->defaultValue(8192000)
        ->boundMin(1)
        ->shortDesc("Number of reactors to solve")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<int>("nMatrixThreads"))
        ->defaultValue(1)
        ->boundMax(32)
        ->boundMin(1)
        ->shortDesc("Number of threads to use in matrix factor")  
    ); 
    
    addParameter
    (
        (new spify::map_param<std::string,double>("fuelComp"))
        ->shortDesc("Fuel composition map")  
    ); 

    addParameter
    (
        (new spify::map_param<std::string,double>("oxidizerComp"))
        ->shortDesc("Oxidizer composition map")  
    ); 

    std::vector<std::string> reactorType_disv;
    reactorType_disv.push_back("CV");
    reactorType_disv.push_back("CP");
    addParameter
    (
        (new spify::scalar_param<std::string>("reactorType"))
        ->defaultValue("CV")
        ->discreteValues(reactorType_disv)  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("deltaTign"))
        ->defaultValue(400.0)
        ->boundMax(100000.0)
        ->boundMin(50.0)  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("reactorTime"))
        ->shortDesc("Total integration time")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("printDt"))
        ->shortDesc("Print Interval")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("maxDtInternal"))
        ->defaultValue(0.05)
        ->shortDesc("Maximum Internal Integrator Step Size")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<int>("maxSteps"))
        ->defaultValue(1000000)
        ->shortDesc("Maximum Number of Integrator Steps")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("relTol"))
        ->defaultValue(1e-08)
        ->shortDesc("Relative Integrator Tolerance")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("absTol"))
        ->defaultValue(1e-20)
        ->shortDesc("Absolute Integrator Tolerance")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("refTemp"))
        ->defaultValue(1000.0)
        ->shortDesc("Ignition Delat Metric: Maximum species concentration")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("TMin"))
        ->shortDesc("Array of initial temperatures to sweep over")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("TMax"))
        ->shortDesc("Array of initial temperatures to sweep over")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("pMin"))
        ->shortDesc("Array of initial temperatures to sweep over")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("pMax"))
        ->shortDesc("Array of initial temperatures to sweep over")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("phiMin"))
        ->shortDesc("Array of initial temperatures to sweep over")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("phiMax"))
        ->shortDesc("Array of initial temperatures to sweep over")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("egrMin"))
        ->shortDesc("Array of initial temperatures to sweep over")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("egrMax"))
        ->shortDesc("Array of initial temperatures to sweep over")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("precThresh"))
        ->shortDesc("Preconditioner threshold values")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<int>("krylovDim"))
        ->defaultValue(5)
        ->shortDesc("Integrator max krylov dimension")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<bool>("oneStepMode"))
        ->defaultValue(0)
        ->shortDesc("Switch: integrate in one step mode.")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<bool>("printAllSteps"))
        ->defaultValue(0)
        ->shortDesc("Switch: print every integrator step.")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<bool>("doFakeUpdate"))
        ->defaultValue(0)
        ->shortDesc("Switch: fakeupdate.")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<bool>("doILU"))
        ->defaultValue(0)
        ->shortDesc("Switch: ILU.")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<int>("precThreshType"))
        ->defaultValue(1)
        ->shortDesc("Switch: Preconditioner Threshold Type.")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("partialPivotThresh"))
        ->defaultValue(0.0)
        ->shortDesc("partial pivot threshold.")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<int>("permutationType"))
        ->defaultValue(1)
        ->shortDesc("perm type.")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("permThresh"))
        ->defaultValue(0.3)
        ->shortDesc("perm threshold.")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<bool>("strictSamePattern"))
        ->defaultValue(0)
        ->shortDesc("Switch")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("maxGammaOrderChange"))
        ->defaultValue(3.0)
        ->shortDesc("??")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<bool>("cvStabLimDet"))
        ->defaultValue(0)
        ->shortDesc("??")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<int>("maxNumNonLinIters"))
        ->defaultValue(3)
        ->shortDesc("??")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<int>("maxOrd"))
        ->defaultValue(5)
        ->shortDesc("??")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("cvEpsLin"))
        ->defaultValue(0.1)
        ->shortDesc("??")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("cvNlConvCoeff"))
        ->defaultValue(0.05)
        ->shortDesc("??")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<bool>("constantPressure"))
        ->defaultValue(0)
        ->shortDesc("Constant Pressure instead of constant volume")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<bool>("reInitOnPrint"))
        ->defaultValue(0)
        ->shortDesc("Re-initialize CVODE at every print interval")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<bool>("dumpStateOnPrint"))
        ->defaultValue(0)
        ->shortDesc("Dump system state at every print interval")  
    ); 
    
    std::vector<std::string> stateFiles_defv;
    addParameter
    (
        (new spify::vector_param<std::string>("stateFiles"))
        ->shortDesc("Files to read initial states from.")
        ->defaultValue(stateFiles_defv)  
    ); 
    
    std::vector<std::string> linearSolver_disv;
    linearSolver_disv.push_back("IterativeSparse");
    linearSolver_disv.push_back("DirectDense");
    linearSolver_disv.push_back("DirectDenseDVD");
    linearSolver_disv.push_back("DirectSparse");
    linearSolver_disv.push_back("DirectSparseDVD");
    linearSolver_disv.push_back("IterativeSparseCusolverRf");
    linearSolver_disv.push_back("IterativeSparseCusolverRfBlocked");
    linearSolver_disv.push_back("IterativeSparseCusolverRfBatched");
    linearSolver_disv.push_back("IterativeSparseCusolverSpBatched");
    addParameter
    (
        (new spify::scalar_param<std::string>("linearSolver"))
        ->defaultValue("IterativeSparse")
        ->shortDesc("CVODE Nonlinear solver")
        ->discreteValues(linearSolver_disv)  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<bool>("cusolverRf"))
        ->defaultValue(1)
        ->shortDesc("Enable cusolverRf Sparse Solver")  
    ); 
    
    std::vector<int> CUSOLVER_FACTOR_ALG_disv;
    CUSOLVER_FACTOR_ALG_disv.push_back(0);
    CUSOLVER_FACTOR_ALG_disv.push_back(1);
    CUSOLVER_FACTOR_ALG_disv.push_back(2);
    CUSOLVER_FACTOR_ALG_disv.push_back(3);
    addParameter
    (
        (new spify::scalar_param<int>("CUSOLVER_FACTOR_ALG"))
        ->defaultValue(2)
        ->shortDesc("cusolverRf Factorization algorithm")
        ->discreteValues(CUSOLVER_FACTOR_ALG_disv)  
    ); 
    
    std::vector<int> CUSOLVER_SOLVE_ALG_disv;
    CUSOLVER_SOLVE_ALG_disv.push_back(0);
    CUSOLVER_SOLVE_ALG_disv.push_back(1);
    CUSOLVER_SOLVE_ALG_disv.push_back(2);
    addParameter
    (
        (new spify::scalar_param<int>("CUSOLVER_SOLVE_ALG"))
        ->defaultValue(0)
        ->shortDesc("cusolverRf Back Substition algorithm")
        ->discreteValues(CUSOLVER_SOLVE_ALG_disv)  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<bool>("cusolverSp"))
        ->defaultValue(1)
        ->shortDesc("Enable cusolverSp Sparse Solver")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("logiA"))
        ->defaultValue(0.015)
        ->shortDesc("logistical function A factor")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("logiK"))
        ->defaultValue(2.35)
        ->shortDesc("logistical function K factor")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("logiQ"))
        ->defaultValue(3.2)
        ->shortDesc("logistical function Q factor")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("logiM"))
        ->defaultValue(-2.2)
        ->shortDesc("logistical function M factor")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("logiB"))
        ->defaultValue(1.0)
        ->shortDesc("logistical function B factor")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("logirNu"))
        ->defaultValue(1.0)
        ->shortDesc("logistical function rNu factor")  
    ); 
    
    addParameter
    (
        (new spify::scalar_param<double>("min_scaled_dt"))
        ->defaultValue(1e-05)
        ->shortDesc("minimum scaled time step for grouping reactors")  
    ); 
    
    validateParameters();

    //Assign values to member data.
    get("cuda",m_cuda);
    get("cudaDevice",m_cudaDevice);
    get("mechFile",m_mechFile);
    get("thermFile",m_thermFile);
    get("mechLogFile",m_mechLogFile);
    get("outFile",m_outFile);
    get("nReactors",m_nReactors);
    get("nReactorsMax",m_nReactorsMax);
    get("nReactorsMin",m_nReactorsMin);
    get("nMatrixReactors",m_nMatrixReactors);
    get("nMatrixThreads",m_nMatrixThreads);
    get("fuelComp",m_fuelComp);
    get("oxidizerComp",m_oxidizerComp);
    get("reactorType",m_reactorType);
    get("deltaTign",m_deltaTign);
    get("reactorTime",m_reactorTime);
    get("printDt",m_printDt);
    get("maxDtInternal",m_maxDtInternal);
    get("maxSteps",m_maxSteps);
    get("relTol",m_relTol);
    get("absTol",m_absTol);
    get("refTemp",m_refTemp);
    get("TMin",m_TMin);
    get("TMax",m_TMax);
    get("pMin",m_pMin);
    get("pMax",m_pMax);
    get("phiMin",m_phiMin);
    get("phiMax",m_phiMax);
    get("egrMin",m_egrMin);
    get("egrMax",m_egrMax);
    get("precThresh",m_precThresh);
    get("krylovDim",m_krylovDim);
    get("oneStepMode",m_oneStepMode);
    get("printAllSteps",m_printAllSteps);
    get("doFakeUpdate",m_doFakeUpdate);
    get("doILU",m_doILU);
    get("precThreshType",m_precThreshType);
    get("partialPivotThresh",m_partialPivotThresh);
    get("permutationType",m_permutationType);
    get("permThresh",m_permThresh);
    get("strictSamePattern",m_strictSamePattern);
    get("maxGammaOrderChange",m_maxGammaOrderChange);
    get("cvStabLimDet",m_cvStabLimDet);
    get("maxNumNonLinIters",m_maxNumNonLinIters);
    get("maxOrd",m_maxOrd);
    get("cvEpsLin",m_cvEpsLin);
    get("cvNlConvCoeff",m_cvNlConvCoeff);
    get("constantPressure",m_constantPressure);
    get("reInitOnPrint",m_reInitOnPrint);
    get("dumpStateOnPrint",m_dumpStateOnPrint);
    get("stateFiles",m_stateFiles);
    get("linearSolver",m_linearSolver);
    get("cusolverRf",m_cusolverRf);
    get("CUSOLVER_FACTOR_ALG",m_CUSOLVER_FACTOR_ALG);
    get("CUSOLVER_SOLVE_ALG",m_CUSOLVER_SOLVE_ALG);
    get("cusolverSp",m_cusolverSp);
    get("logiA",m_logiA);
    get("logiK",m_logiK);
    get("logiQ",m_logiQ);
    get("logiM",m_logiM);
    get("logiB",m_logiB);
    get("logirNu",m_logirNu);
    get("min_scaled_dt",m_min_scaled_dt);

}

const bool& zerork_reactorIFP::cuda() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"cuda\" before validation."));
    }
    return m_cuda;
}

const int& zerork_reactorIFP::cudaDevice() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"cudaDevice\" before validation."));
    }
    return m_cudaDevice;
}

const std::string& zerork_reactorIFP::mechFile() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"mechFile\" before validation."));
    }
    return m_mechFile;
}

const std::string& zerork_reactorIFP::thermFile() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"thermFile\" before validation."));
    }
    return m_thermFile;
}

const std::string& zerork_reactorIFP::mechLogFile() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"mechLogFile\" before validation."));
    }
    return m_mechLogFile;
}

const std::string& zerork_reactorIFP::outFile() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"outFile\" before validation."));
    }
    return m_outFile;
}

const int& zerork_reactorIFP::nReactors() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"nReactors\" before validation."));
    }
    return m_nReactors;
}

const int& zerork_reactorIFP::nReactorsMax() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"nReactorsMax\" before validation."));
    }
    return m_nReactorsMax;
}

const int& zerork_reactorIFP::nReactorsMin() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"nReactorsMin\" before validation."));
    }
    return m_nReactorsMin;
}

const int& zerork_reactorIFP::nMatrixReactors() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"nMatrixReactors\" before validation."));
    }
    return m_nMatrixReactors;
}

const int& zerork_reactorIFP::nMatrixThreads() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"nMatrixThreads\" before validation."));
    }
    return m_nMatrixThreads;
}

const std::map<std::string,double>& zerork_reactorIFP::fuelComp() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"fuelComp\" before validation."));
    }
    return m_fuelComp;
}

const std::map<std::string,double>& zerork_reactorIFP::oxidizerComp() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"oxidizerComp\" before validation."));
    }
    return m_oxidizerComp;
}

const std::string& zerork_reactorIFP::reactorType() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"reactorType\" before validation."));
    }
    return m_reactorType;
}

const double& zerork_reactorIFP::deltaTign() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"deltaTign\" before validation."));
    }
    return m_deltaTign;
}

const double& zerork_reactorIFP::reactorTime() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"reactorTime\" before validation."));
    }
    return m_reactorTime;
}

const double& zerork_reactorIFP::printDt() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"printDt\" before validation."));
    }
    return m_printDt;
}

const double& zerork_reactorIFP::maxDtInternal() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"maxDtInternal\" before validation."));
    }
    return m_maxDtInternal;
}

const int& zerork_reactorIFP::maxSteps() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"maxSteps\" before validation."));
    }
    return m_maxSteps;
}

const double& zerork_reactorIFP::relTol() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"relTol\" before validation."));
    }
    return m_relTol;
}

const double& zerork_reactorIFP::absTol() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"absTol\" before validation."));
    }
    return m_absTol;
}

const double& zerork_reactorIFP::refTemp() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"refTemp\" before validation."));
    }
    return m_refTemp;
}

const double& zerork_reactorIFP::TMin() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"TMin\" before validation."));
    }
    return m_TMin;
}

const double& zerork_reactorIFP::TMax() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"TMax\" before validation."));
    }
    return m_TMax;
}

const double& zerork_reactorIFP::pMin() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"pMin\" before validation."));
    }
    return m_pMin;
}

const double& zerork_reactorIFP::pMax() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"pMax\" before validation."));
    }
    return m_pMax;
}

const double& zerork_reactorIFP::phiMin() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"phiMin\" before validation."));
    }
    return m_phiMin;
}

const double& zerork_reactorIFP::phiMax() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"phiMax\" before validation."));
    }
    return m_phiMax;
}

const double& zerork_reactorIFP::egrMin() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"egrMin\" before validation."));
    }
    return m_egrMin;
}

const double& zerork_reactorIFP::egrMax() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"egrMax\" before validation."));
    }
    return m_egrMax;
}

const double& zerork_reactorIFP::precThresh() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"precThresh\" before validation."));
    }
    return m_precThresh;
}

const int& zerork_reactorIFP::krylovDim() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"krylovDim\" before validation."));
    }
    return m_krylovDim;
}

const bool& zerork_reactorIFP::oneStepMode() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"oneStepMode\" before validation."));
    }
    return m_oneStepMode;
}

const bool& zerork_reactorIFP::printAllSteps() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"printAllSteps\" before validation."));
    }
    return m_printAllSteps;
}

const bool& zerork_reactorIFP::doFakeUpdate() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"doFakeUpdate\" before validation."));
    }
    return m_doFakeUpdate;
}

const bool& zerork_reactorIFP::doILU() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"doILU\" before validation."));
    }
    return m_doILU;
}

const int& zerork_reactorIFP::precThreshType() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"precThreshType\" before validation."));
    }
    return m_precThreshType;
}

const double& zerork_reactorIFP::partialPivotThresh() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"partialPivotThresh\" before validation."));
    }
    return m_partialPivotThresh;
}

const int& zerork_reactorIFP::permutationType() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"permutationType\" before validation."));
    }
    return m_permutationType;
}

const double& zerork_reactorIFP::permThresh() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"permThresh\" before validation."));
    }
    return m_permThresh;
}

const bool& zerork_reactorIFP::strictSamePattern() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"strictSamePattern\" before validation."));
    }
    return m_strictSamePattern;
}

const double& zerork_reactorIFP::maxGammaOrderChange() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"maxGammaOrderChange\" before validation."));
    }
    return m_maxGammaOrderChange;
}

const bool& zerork_reactorIFP::cvStabLimDet() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"cvStabLimDet\" before validation."));
    }
    return m_cvStabLimDet;
}

const int& zerork_reactorIFP::maxNumNonLinIters() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"maxNumNonLinIters\" before validation."));
    }
    return m_maxNumNonLinIters;
}

const int& zerork_reactorIFP::maxOrd() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"maxOrd\" before validation."));
    }
    return m_maxOrd;
}

const double& zerork_reactorIFP::cvEpsLin() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"cvEpsLin\" before validation."));
    }
    return m_cvEpsLin;
}

const double& zerork_reactorIFP::cvNlConvCoeff() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"cvNlConvCoeff\" before validation."));
    }
    return m_cvNlConvCoeff;
}

const bool& zerork_reactorIFP::constantPressure() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"constantPressure\" before validation."));
    }
    return m_constantPressure;
}

const bool& zerork_reactorIFP::reInitOnPrint() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"reInitOnPrint\" before validation."));
    }
    return m_reInitOnPrint;
}

const bool& zerork_reactorIFP::dumpStateOnPrint() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"dumpStateOnPrint\" before validation."));
    }
    return m_dumpStateOnPrint;
}

const std::vector<std::string>& zerork_reactorIFP::stateFiles() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"stateFiles\" before validation."));
    }
    return m_stateFiles;
}

const std::string& zerork_reactorIFP::linearSolver() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"linearSolver\" before validation."));
    }
    return m_linearSolver;
}

const bool& zerork_reactorIFP::cusolverRf() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"cusolverRf\" before validation."));
    }
    return m_cusolverRf;
}

const int& zerork_reactorIFP::CUSOLVER_FACTOR_ALG() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"CUSOLVER_FACTOR_ALG\" before validation."));
    }
    return m_CUSOLVER_FACTOR_ALG;
}

const int& zerork_reactorIFP::CUSOLVER_SOLVE_ALG() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"CUSOLVER_SOLVE_ALG\" before validation."));
    }
    return m_CUSOLVER_SOLVE_ALG;
}

const bool& zerork_reactorIFP::cusolverSp() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"cusolverSp\" before validation."));
    }
    return m_cusolverSp;
}

const double& zerork_reactorIFP::logiA() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"logiA\" before validation."));
    }
    return m_logiA;
}

const double& zerork_reactorIFP::logiK() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"logiK\" before validation."));
    }
    return m_logiK;
}

const double& zerork_reactorIFP::logiQ() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"logiQ\" before validation."));
    }
    return m_logiQ;
}

const double& zerork_reactorIFP::logiM() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"logiM\" before validation."));
    }
    return m_logiM;
}

const double& zerork_reactorIFP::logiB() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"logiB\" before validation."));
    }
    return m_logiB;
}

const double& zerork_reactorIFP::logirNu() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"logirNu\" before validation."));
    }
    return m_logirNu;
}

const double& zerork_reactorIFP::min_scaled_dt() const
{
    if(!m_parametersValidated)
    {
       throw(std::logic_error("zerork_reactorIFP: attempted to access parameter "
                             "\"min_scaled_dt\" before validation."));
    }
    return m_min_scaled_dt;
}




#ifdef ZERORK_REACTORIFP_MAKEMASTER
#include <unistd.h>
#include <string.h>
#include <stdio.h>

int main(int argc, char **argv)
{
    FILE * tmpFile;
    char tmpFileName[32];

    memset(tmpFileName,0,sizeof(tmpFileName));
    strncpy(tmpFileName,"ymlfXXXXXX",10);
    mkstemp(tmpFileName);
    tmpFile = fopen(tmpFileName,"w");
    fprintf(tmpFile,"printMasterFileTo: zerork_reactorIFPMaster.yml\n");
    fclose(tmpFile);

    try //This will fail.
    {
        zerork_reactorIFP inputData(tmpFileName);
    }
    catch (...)
    {
        remove(tmpFileName);
    }
    return(0);
}
#endif

