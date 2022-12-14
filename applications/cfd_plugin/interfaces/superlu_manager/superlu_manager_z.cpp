#include "superlu_manager_z.h"


superlu_manager_z::superlu_manager_z() :
  m_last_factor_nnz(-1),
  m_factored(false)
{
  zCreate_CompCol_Matrix(&m_M, 1, 1, 1, NULL, NULL, NULL, SLU_NC,SLU_Z,SLU_GE);
  zCreate_Dense_Matrix(&m_B,1,1,NULL,1, SLU_DN,SLU_Z,SLU_GE);
  zCreate_Dense_Matrix(&m_X,1,1,NULL,1, SLU_DN,SLU_Z,SLU_GE);
  set_default_options(&m_options);
  m_options.ColPerm = MY_PERMC;
  m_options.RowPerm = LargeDiag_MC64;
  m_options.Equil   = YES;

  StatInit(&m_stats);
}


superlu_manager_z::~superlu_manager_z()
{
  StatFree(&m_stats);
  Destroy_SuperMatrix_Store(&m_B);
  Destroy_SuperMatrix_Store(&m_X);
  if(m_factored) {
    Destroy_SuperNode_Matrix(&m_L);
    Destroy_CompCol_Matrix(&m_U);
  }
  Destroy_SuperMatrix_Store(&m_M);
}

void superlu_manager_z::setup_memory(int n, int nnz, const int* indexes, const int* sums, const doublecomplex* values) {
  m_M.nrow = n;
  m_M.ncol = n;
  ((NCformat *)m_M.Store)->nnz = nnz;
  ((NCformat *)m_M.Store)->nzval = const_cast<doublecomplex*>(values);
  ((NCformat *)m_M.Store)->rowind = const_cast<int*>(indexes);
  ((NCformat *)m_M.Store)->colptr = const_cast<int*>(sums);

  m_colPermutation.assign(n,0);
  m_rowPermutation.assign(n,0);
  m_colElimTree.assign(n,0);
  m_R.assign(n,0.0);
  m_C.assign(n,0.0);

  m_B.nrow=n;
  m_X.nrow=n;
  m_B.ncol=1;
  m_X.ncol=1;
  ((DNformat *)m_B.Store)->lda = n;
  ((DNformat *)m_X.Store)->lda = n;
}


int superlu_manager_z::factor(const std::vector<int>& indexes,
                              const std::vector<int>& sums,
                              const std::vector<doublecomplex>& values,
                              matrix_t type)
{
  int n = sums.size()-1;
  int nnz = indexes.size();

  if(values.size() != nnz) {
    return 1;
  }
  return this->factor(n, nnz, &(indexes[0]), &(sums[0]), &(values[0]), type);
}

int superlu_manager_z::factor(int n,
                              int nnz,
                              const int* indexes,
                              const int* sums,
                              const doublecomplex* values,
                              matrix_t type)
{
  if(nnz != m_last_factor_nnz || n != m_last_factor_n) {
    setup_memory(n,nnz, indexes, sums, values);
  }

  if(m_factored) {
    Destroy_SuperNode_Matrix(&m_L);
    Destroy_CompCol_Matrix(&m_U);
  }

  m_B.ncol=0;
  m_options.Fact=DOFACT;
  get_perm_c(MMD_AT_PLUS_A, &m_M, &m_colPermutation[0]);

  int lwork=0;      // SuperLU allocates its own memeory
  void *work=NULL;  // ''
  double rpg,rcond; // reciprocal pivot growth, recip. condition number
  double *ferr = NULL; //Not factoring so un-used
  double *berr = NULL;
  int flag;
  zgssvx(&m_options,&m_M, &m_colPermutation[0],&m_rowPermutation[0],
         &m_colElimTree[0],m_equed,&m_R[0],&m_C[0],&m_L, &m_U, work, lwork,
         &m_B, &m_X, &rpg, &rcond, ferr, berr,
#if SUPERLU_MAJOR_VERSION == 5
         &m_Glu,
#endif
         &m_mem_usage, &m_stats, &flag);

  m_factored=true;
  m_last_factor_n = n;
  m_last_factor_nnz = nnz;

  return flag;
}


int superlu_manager_z::refactor(const std::vector<doublecomplex>& values) {
  int nnz = values.size();
  return this->refactor(nnz, &(values[0]));
}

int superlu_manager_z::refactor(int nnz, const doublecomplex* values)
{
  if(nnz != m_last_factor_nnz) {
    return 1;
  }
  ((NCformat *)m_M.Store)->nzval = const_cast<doublecomplex*>(values);
  Destroy_SuperNode_Matrix(&m_L);
  Destroy_CompCol_Matrix(&m_U);

  m_B.ncol=0;
  m_options.Fact=SamePattern;

  int lwork=0;      // SuperLU allocates its own memeory
  void *work=NULL;  // ''
  double rpg,rcond; // reciprocal pivot growth, recip. condition number
  double *ferr = NULL; //Not factoring so un-used
  double *berr = NULL;
  int flag;
  zgssvx(&m_options,&m_M, &m_colPermutation[0],&m_rowPermutation[0],
         &m_colElimTree[0],m_equed,&m_R[0],&m_C[0],&m_L, &m_U, work, lwork,
         &m_B, &m_X, &rpg, &rcond, ferr, berr,
#if SUPERLU_MAJOR_VERSION == 5
         &m_Glu,
#endif
         &m_mem_usage, &m_stats, &flag);

  if(m_options.DiagPivotThresh > 0.0) {
    for(int j = 0; j < m_last_factor_n; ++j) {
      if(m_rowPermutation[j] == -1 ) {
        flag = 1;
      }
    }
  }

  return flag;
}


int superlu_manager_z::solve(const std::vector<doublecomplex> rhs, std::vector<doublecomplex>* soln) {
  if(rhs.size() != soln->size()) {
    return 1;
  }
  int n = rhs.size();
  return this->solve(n, &(rhs[0]), &((*soln)[0]));
}

int superlu_manager_z::solve(int n, const doublecomplex* rhs, doublecomplex* soln) {
  if(n != m_last_factor_n) {
    return 1;
  }
  int flag;
  int lwork=0;
  void *work=NULL;
  double rpg,rcond;
  double ferr[1],berr[1]; // length is the number of RHS

  m_options.Fact=FACTORED;
  m_B.ncol=1;

  ((DNformat *)m_B.Store)->nzval = const_cast<doublecomplex*>(rhs);
  ((DNformat *)m_X.Store)->nzval = soln;

  //backsolve with expert driver function
  zgssvx(&m_options,&m_M,&m_colPermutation[0],&m_rowPermutation[0],
         &m_colElimTree[0], m_equed, &m_R[0], &m_C[0], &m_L, &m_U,
	 work, lwork, &m_B, &m_X, &rpg, &rcond, ferr, berr,
#if SUPERLU_MAJOR_VERSION == 5
         &m_Glu,
#endif
         &m_mem_usage, &m_stats, &flag);

  return flag;
}
