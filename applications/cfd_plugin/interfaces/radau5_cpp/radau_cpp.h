#ifndef RADAU_H
#define RADAU_H

#include <vector>

#include "sundials/sundials_nvector.h"

namespace radau_cpp {

typedef int (*radau_deriv_fcn)(double x, N_Vector y, N_Vector dy,
                               void* user_data);
typedef int (*radau_jac_fcn)(double x, N_Vector y, N_Vector fy,
                             void* user_data);
typedef int (*radau_jac_decomp_fcn)(double hji, void* user_data);
typedef int (*radau_jac_solve_fcn)(double x, N_Vector y, N_Vector fy,
                                   N_Vector r, N_Vector z, void* user_data);
typedef int (*radau_jac_complex_decomp_fcn)(int k, double alpha, double beta,
                                            void* user_data);
typedef int (*radau_jac_complex_solve_fcn)(int k, N_Vector ax, N_Vector bx,
                                           void* user_data);
typedef int (*radau_output_fcn)(int nsteps, double x, double h, N_Vector y,
                                N_Vector dy, void* user_data);

class radau {
 public:
  radau(int, N_Vector);
  virtual ~radau();

  void set_nmax(int);
  void set_uround(double);

  void set_hinit(double);
  void set_hmax(double);
  void set_thet(double);
  void set_step_size_params(double, double);
  void set_order_params(double, double, double, double);
  void set_safe(double);
  void set_quot_factors(double, double);
  void set_nsmin(int);
  void set_nsmax(int);
  void set_nsus(int);
  void set_nit(int);
  void set_startn(bool);
  void set_pred(bool);
  void set_tolerances(N_Vector, N_Vector);
  void set_tolerances(double, double);
  void set_deriv_fcn(radau_deriv_fcn);
  void set_jac_fcn(radau_jac_fcn);
  void set_jac_decomp_fcn(radau_jac_decomp_fcn);
  void set_jac_solve_fcn(radau_jac_solve_fcn);
  void set_jac_complex_decomp_fcn(radau_jac_complex_decomp_fcn);
  void set_jac_complex_solve_fcn(radau_jac_complex_solve_fcn);
  void set_output_fcn(radau_output_fcn, void*);
  void set_user_data(void*);
  int solve(double* x_in, double xend, N_Vector y);
  void get_integrator_stats(int* nfcn, int* njac, int* nstep, int* naccpt,
                            int* nrejct, int* ndec, int* nsol);
  void contra(const double x, N_Vector yout);

 private:
  void coercv(int);
  int slvrar(double fac1, N_Vector z1, N_Vector f1);
  int slvrai(int k, double alpha, double beta, N_Vector z2, N_Vector z3,
             N_Vector f2, N_Vector f3, N_Vector cc);
  int slvrad(double fac1, double alpha, double beta, N_Vector z1, N_Vector z2,
             N_Vector z3, N_Vector f1, N_Vector f2, N_Vector f3, N_Vector cc);
  double estrav(N_Vector y);
  double estrad(N_Vector y);

  int n;
  int nfcn;
  int njac;
  int nstep;
  int naccept;
  int nreject;
  int ndec;
  int nsol;
  // -------- NUMBER MAXIMAL AND MINIMAL OF STAGES  NS
  int ns;
  int nsmin;
  int nsmax;
  int nsus;
  int nit;
  // -------- NMAX , THE MAXIMAL NUMBER OF STEPS -----
  int nmax;
  // -------- STARTN  SWITCH FOR STARTING VALUES OF NEWTON ITERATIONS
  bool startn;
  // -------- PRED   STEP SIZE CONTROL
  bool pred;
  // -------- UROUND   SMALLEST NUMBER SATISFYING 1.D0+UROUND>1.D0  double
  // uround;
  double uround;
  // -------- MAXIMAL STEP SIZE
  double h;
  double hacc;
  double hmax;
  // ------ THET     DECIDES WHETHER THE JACOBIAN SHOULD BE RECOMPUTED;
  double thet;

  // --- QUOT1 AND QUOT2: IF QUOT1 < HNEW/HOLD < QUOT2, STEP SIZE = CONST.
  double quot1;
  double quot2;

  // -------  FACL,FACR     PARAMETERS FOR STEP SIZE SELECTION
  double facl;
  double facr;

  // ------- SAFETY FACTOR FOR STEP SIZE PREDICTION
  double safe;
  // -------- PARAMETERS FOR ORDER SELECTION STRATEGY
  double vitu;
  double vitd;
  double hhou;
  double hhod;

  double u1;
  double c[8];
  double dd[7];
  double alph[3];
  double beta[3];
  double alphn[3];
  double betan[3];

  constexpr static double t311 = 0.9123239487089294279155e-01;
  constexpr static double t312 = -0.1412552950209542084280e+00;
  constexpr static double t313 = -0.3002919410514742449186e-01;
  constexpr static double t321 = 0.2417179327071070189575e+00;
  constexpr static double t322 = 0.2041293522937999319960e+00;
  constexpr static double t323 = 0.3829421127572619377954e+00;
  constexpr static double t331 = 0.9660481826150929361906e+00;
  constexpr static double ti311 = 0.4325579890063155351024e+01;
  constexpr static double ti312 = 0.3391992518158098695428e+00;
  constexpr static double ti313 = 0.5417705399358748711865e+00;
  constexpr static double ti321 = -0.4178718591551904727346e+01;
  constexpr static double ti322 = -0.3276828207610623870825e+00;
  constexpr static double ti323 = 0.4766235545005504519601e+00;
  constexpr static double ti331 = -0.5028726349457868759512e+00;
  constexpr static double ti332 = 0.2571926949855605429187e+01;
  constexpr static double ti333 = -0.5960392048282249249688e+00;

  constexpr static double t511 = -0.1251758622050104589014e-01;
  constexpr static double t512 = -0.1024204781790882707009e-01;
  constexpr static double t513 = 0.4767387729029572386318e-01;
  constexpr static double t514 = -0.1147851525522951470794e-01;
  constexpr static double t515 = -0.1401985889287541028108e-01;
  constexpr static double t521 = -0.1491670151895382429004e-02;
  constexpr static double t522 = 0.5017286451737105816299e-01;
  constexpr static double t523 = -0.9433181918161143698066e-01;
  constexpr static double t524 = -0.7668830749180162885157e-02;
  constexpr static double t525 = 0.2470857842651852681253e-01;
  constexpr static double t531 = -0.7298187638808714862266e-01;
  constexpr static double t532 = -0.2305395340434179467214e+00;
  constexpr static double t533 = 0.1027030453801258997922e+00;
  constexpr static double t534 = 0.1939846399882895091122e-01;
  constexpr static double t535 = 0.8180035370375117083639e-01;
  constexpr static double t541 = -0.3800914400035681041264e+00;
  constexpr static double t542 = 0.3778939022488612495439e+00;
  constexpr static double t543 = 0.4667441303324943592896e+00;
  constexpr static double t544 = 0.4076011712801990666217e+00;
  constexpr static double t545 = 0.1996824278868025259365e+00;
  constexpr static double t551 = -0.9219789736812104884883e+00;
  constexpr static double ti511 = -0.3004156772154440162771e+02;
  constexpr static double ti512 = -0.1386510785627141316518e+02;
  constexpr static double ti513 = -0.3480002774795185561828e+01;
  constexpr static double ti514 = 0.1032008797825263422771e+01;
  constexpr static double ti515 = -0.8043030450739899174753e+00;
  constexpr static double ti521 = 0.5344186437834911598895e+01;
  constexpr static double ti522 = 0.4593615567759161004454e+01;
  constexpr static double ti523 = -0.3036360323459424298646e+01;
  constexpr static double ti524 = 0.1050660190231458863860e+01;
  constexpr static double ti525 = -0.2727786118642962705386e+00;
  constexpr static double ti531 = 0.3748059807439804860051e+01;
  constexpr static double ti532 = -0.3984965736343884667252e+01;
  constexpr static double ti533 = -0.1044415641608018792942e+01;
  constexpr static double ti534 = 0.1184098568137948487231e+01;
  constexpr static double ti535 = -0.4499177701567803688988e+00;
  constexpr static double ti541 = -0.3304188021351900000806e+02;
  constexpr static double ti542 = -0.1737695347906356701945e+02;
  constexpr static double ti543 = -0.1721290632540055611515e+00;
  constexpr static double ti544 = -0.9916977798254264258817e-01;
  constexpr static double ti545 = 0.5312281158383066671849e+00;
  constexpr static double ti551 = -0.8611443979875291977700e+01;
  constexpr static double ti552 = 0.9699991409528808231336e+01;
  constexpr static double ti553 = 0.1914728639696874284851e+01;
  constexpr static double ti554 = 0.2418692006084940026427e+01;
  constexpr static double ti555 = -0.1047463487935337418694e+01;

  constexpr static double t711 = -0.2153754627310526422828e-02;
  constexpr static double t712 = 0.2156755135132077338691e-01;
  constexpr static double t713 = 0.8783567925144144407326e-02;
  constexpr static double t714 = -0.4055161452331023898198e-02;
  constexpr static double t715 = 0.4427232753268285479678e-02;
  constexpr static double t716 = -0.1238646187952874056377e-02;
  constexpr static double t717 = -0.2760617480543852499548e-02;
  constexpr static double t721 = 0.1600025077880428526831e-02;
  constexpr static double t722 = -0.3813164813441154669442e-01;
  constexpr static double t723 = -0.2152556059400687552385e-01;
  constexpr static double t724 = 0.8415568276559589237177e-02;
  constexpr static double t725 = -0.4031949570224549492304e-02;
  constexpr static double t726 = -0.6666635339396338181761e-04;
  constexpr static double t727 = 0.3185474825166209848748e-02;
  constexpr static double t731 = -0.4059107301947683091650e-02;
  constexpr static double t732 = 0.5739650893938171539757e-01;
  constexpr static double t733 = 0.5885052920842679105612e-01;
  constexpr static double t734 = -0.8560431061603432060177e-02;
  constexpr static double t735 = -0.6923212665023908924141e-02;
  constexpr static double t736 = -0.2352180982943338340535e-02;
  constexpr static double t737 = 0.4169077725297562691409e-03;
  constexpr static double t741 = -0.1575048807937684420346e-01;
  constexpr static double t742 = -0.3821469359696835048464e-01;
  constexpr static double t743 = -0.1657368112729438512412e+00;
  constexpr static double t744 = -0.3737124230238445741907e-01;
  constexpr static double t745 = 0.8239007298507719404499e-02;
  constexpr static double t746 = 0.3115071152346175252726e-02;
  constexpr static double t747 = 0.2511660491343882192836e-01;
  constexpr static double t751 = -0.1129776610242208076086e+00;
  constexpr static double t752 = -0.2491742124652636863308e+00;
  constexpr static double t753 = 0.2735633057986623212132e+00;
  constexpr static double t754 = 0.5366761379181770094279e-02;
  constexpr static double t755 = 0.1932111161012620144312e+00;
  constexpr static double t756 = 0.1017177324817151468081e+00;
  constexpr static double t757 = 0.9504502035604622821039e-01;
  constexpr static double t761 = -0.4583810431839315010281e+00;
  constexpr static double t762 = 0.5315846490836284292051e+00;
  constexpr static double t763 = 0.4863228366175728940567e+00;
  constexpr static double t764 = 0.5265742264584492629141e+00;
  constexpr static double t765 = 0.2755343949896258141929e+00;
  constexpr static double t766 = 0.5217519452747652852946e+00;
  constexpr static double t767 = 0.1280719446355438944141e+00;
  constexpr static double t771 = -0.8813915783538183763135e+00;
  constexpr static double ti711 = -0.2581319263199822292761e+03;
  constexpr static double ti712 = -0.1890737630813985089520e+03;
  constexpr static double ti713 = -0.4908731481793013119445e+02;
  constexpr static double ti714 = -0.4110647469661428418112e+01;
  constexpr static double ti715 = -0.4053447889315563304175e+01;
  constexpr static double ti716 = 0.3112755366607346076554e+01;
  constexpr static double ti717 = -0.1646774913558444650169e+01;
  constexpr static double ti721 = -0.3007390169451292131731e+01;
  constexpr static double ti722 = -0.1101586607876577132911e+02;
  constexpr static double ti723 = 0.1487799456131656281486e+01;
  constexpr static double ti724 = 0.2130388159559282459432e+01;
  constexpr static double ti725 = -0.1816141086817565624822e+01;
  constexpr static double ti726 = 0.1134325587895161100083e+01;
  constexpr static double ti727 = -0.4146990459433035319930e+00;
  constexpr static double ti731 = -0.8441963188321084681757e+01;
  constexpr static double ti732 = -0.6505252740575150028169e+00;
  constexpr static double ti733 = 0.6940670730369876478804e+01;
  constexpr static double ti734 = -0.3205047525597898431565e+01;
  constexpr static double ti735 = 0.1071280943546478589783e+01;
  constexpr static double ti736 = -0.3548507491216221879730e+00;
  constexpr static double ti737 = 0.9198549132786554154409e-01;
  constexpr static double ti741 = 0.7467833223502269977153e+02;
  constexpr static double ti742 = 0.8740858897990081640204e+02;
  constexpr static double ti743 = 0.4024158737379997877014e+01;
  constexpr static double ti744 = -0.3714806315158364186639e+01;
  constexpr static double ti745 = -0.3430093985982317350741e+01;
  constexpr static double ti746 = 0.2696604809765312378853e+01;
  constexpr static double ti747 = -0.9386927436075461933568e+00;
  constexpr static double ti751 = 0.5835652885190657724237e+02;
  constexpr static double ti752 = -0.1006877395780018096325e+02;
  constexpr static double ti753 = -0.3036638884256667120811e+02;
  constexpr static double ti754 = -0.1020020865184865985027e+01;
  constexpr static double ti755 = -0.1124175003784249621267e+00;
  constexpr static double ti756 = 0.1890640831000377622800e+01;
  constexpr static double ti757 = -0.9716486393831482282172e+00;
  constexpr static double ti761 = -0.2991862480282520966786e+03;
  constexpr static double ti762 = -0.2430407453687447911819e+03;
  constexpr static double ti763 = -0.4877710407803786921219e+02;
  constexpr static double ti764 = -0.2038671905741934405280e+01;
  constexpr static double ti765 = 0.1673560239861084944268e+01;
  constexpr static double ti766 = -0.1087374032057106164456e+01;
  constexpr static double ti767 = 0.9019382492960993738427e+00;
  constexpr static double ti771 = -0.9307650289743530591157e+02;
  constexpr static double ti772 = 0.2388163105628114427703e+02;
  constexpr static double ti773 = 0.3927888073081384382710e+02;
  constexpr static double ti774 = 0.1438891568549108006988e+02;
  constexpr static double ti775 = -0.3510438399399361221087e+01;
  constexpr static double ti776 = 0.4863284885566180701215e+01;
  constexpr static double ti777 = -0.2246482729591239916400e+01;

  // double rtol;
  // double atol;
  N_Vector rtol;
  N_Vector atol;

  bool caljac;
  bool first;
  bool last;
  bool reject;
  double err;
  double erracc;
  double x;

  void check_tolerances(void);
  void create_work_vectors(N_Vector y);
  void destroy_work_vectors();

  // Work arraeys
  N_Vector y;
  N_Vector y0;
  N_Vector dy;
  N_Vector scal;
  N_Vector tmp1;
  N_Vector tmp2;
  N_Vector tmp3;
  std::vector<N_Vector> ff;
  std::vector<N_Vector> zz;
  std::vector<N_Vector> cont;
  std::vector<N_Vector> tmp;

  radau_deriv_fcn deriv_fcn;
  radau_jac_fcn jac_fcn;
  radau_jac_decomp_fcn jac_decomp_fcn;
  radau_jac_solve_fcn jac_solve_fcn;
  radau_jac_complex_decomp_fcn jac_complex_decomp_fcn;
  radau_jac_complex_solve_fcn jac_complex_solve_fcn;
  radau_output_fcn output_fcn;

  void* user_data;
  void* output_fcn_data;
};

}  // end namespace radau_cpp

#endif
