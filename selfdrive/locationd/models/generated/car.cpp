#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with sympy 1.9                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_6029616048175551447) {
   out_6029616048175551447[0] = delta_x[0] + nom_x[0];
   out_6029616048175551447[1] = delta_x[1] + nom_x[1];
   out_6029616048175551447[2] = delta_x[2] + nom_x[2];
   out_6029616048175551447[3] = delta_x[3] + nom_x[3];
   out_6029616048175551447[4] = delta_x[4] + nom_x[4];
   out_6029616048175551447[5] = delta_x[5] + nom_x[5];
   out_6029616048175551447[6] = delta_x[6] + nom_x[6];
   out_6029616048175551447[7] = delta_x[7] + nom_x[7];
   out_6029616048175551447[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_114715398326382947) {
   out_114715398326382947[0] = -nom_x[0] + true_x[0];
   out_114715398326382947[1] = -nom_x[1] + true_x[1];
   out_114715398326382947[2] = -nom_x[2] + true_x[2];
   out_114715398326382947[3] = -nom_x[3] + true_x[3];
   out_114715398326382947[4] = -nom_x[4] + true_x[4];
   out_114715398326382947[5] = -nom_x[5] + true_x[5];
   out_114715398326382947[6] = -nom_x[6] + true_x[6];
   out_114715398326382947[7] = -nom_x[7] + true_x[7];
   out_114715398326382947[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_2186784832712934011) {
   out_2186784832712934011[0] = 1.0;
   out_2186784832712934011[1] = 0;
   out_2186784832712934011[2] = 0;
   out_2186784832712934011[3] = 0;
   out_2186784832712934011[4] = 0;
   out_2186784832712934011[5] = 0;
   out_2186784832712934011[6] = 0;
   out_2186784832712934011[7] = 0;
   out_2186784832712934011[8] = 0;
   out_2186784832712934011[9] = 0;
   out_2186784832712934011[10] = 1.0;
   out_2186784832712934011[11] = 0;
   out_2186784832712934011[12] = 0;
   out_2186784832712934011[13] = 0;
   out_2186784832712934011[14] = 0;
   out_2186784832712934011[15] = 0;
   out_2186784832712934011[16] = 0;
   out_2186784832712934011[17] = 0;
   out_2186784832712934011[18] = 0;
   out_2186784832712934011[19] = 0;
   out_2186784832712934011[20] = 1.0;
   out_2186784832712934011[21] = 0;
   out_2186784832712934011[22] = 0;
   out_2186784832712934011[23] = 0;
   out_2186784832712934011[24] = 0;
   out_2186784832712934011[25] = 0;
   out_2186784832712934011[26] = 0;
   out_2186784832712934011[27] = 0;
   out_2186784832712934011[28] = 0;
   out_2186784832712934011[29] = 0;
   out_2186784832712934011[30] = 1.0;
   out_2186784832712934011[31] = 0;
   out_2186784832712934011[32] = 0;
   out_2186784832712934011[33] = 0;
   out_2186784832712934011[34] = 0;
   out_2186784832712934011[35] = 0;
   out_2186784832712934011[36] = 0;
   out_2186784832712934011[37] = 0;
   out_2186784832712934011[38] = 0;
   out_2186784832712934011[39] = 0;
   out_2186784832712934011[40] = 1.0;
   out_2186784832712934011[41] = 0;
   out_2186784832712934011[42] = 0;
   out_2186784832712934011[43] = 0;
   out_2186784832712934011[44] = 0;
   out_2186784832712934011[45] = 0;
   out_2186784832712934011[46] = 0;
   out_2186784832712934011[47] = 0;
   out_2186784832712934011[48] = 0;
   out_2186784832712934011[49] = 0;
   out_2186784832712934011[50] = 1.0;
   out_2186784832712934011[51] = 0;
   out_2186784832712934011[52] = 0;
   out_2186784832712934011[53] = 0;
   out_2186784832712934011[54] = 0;
   out_2186784832712934011[55] = 0;
   out_2186784832712934011[56] = 0;
   out_2186784832712934011[57] = 0;
   out_2186784832712934011[58] = 0;
   out_2186784832712934011[59] = 0;
   out_2186784832712934011[60] = 1.0;
   out_2186784832712934011[61] = 0;
   out_2186784832712934011[62] = 0;
   out_2186784832712934011[63] = 0;
   out_2186784832712934011[64] = 0;
   out_2186784832712934011[65] = 0;
   out_2186784832712934011[66] = 0;
   out_2186784832712934011[67] = 0;
   out_2186784832712934011[68] = 0;
   out_2186784832712934011[69] = 0;
   out_2186784832712934011[70] = 1.0;
   out_2186784832712934011[71] = 0;
   out_2186784832712934011[72] = 0;
   out_2186784832712934011[73] = 0;
   out_2186784832712934011[74] = 0;
   out_2186784832712934011[75] = 0;
   out_2186784832712934011[76] = 0;
   out_2186784832712934011[77] = 0;
   out_2186784832712934011[78] = 0;
   out_2186784832712934011[79] = 0;
   out_2186784832712934011[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_5798323304716649818) {
   out_5798323304716649818[0] = state[0];
   out_5798323304716649818[1] = state[1];
   out_5798323304716649818[2] = state[2];
   out_5798323304716649818[3] = state[3];
   out_5798323304716649818[4] = state[4];
   out_5798323304716649818[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_5798323304716649818[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_5798323304716649818[7] = state[7];
   out_5798323304716649818[8] = state[8];
}
void F_fun(double *state, double dt, double *out_2117145641404181425) {
   out_2117145641404181425[0] = 1;
   out_2117145641404181425[1] = 0;
   out_2117145641404181425[2] = 0;
   out_2117145641404181425[3] = 0;
   out_2117145641404181425[4] = 0;
   out_2117145641404181425[5] = 0;
   out_2117145641404181425[6] = 0;
   out_2117145641404181425[7] = 0;
   out_2117145641404181425[8] = 0;
   out_2117145641404181425[9] = 0;
   out_2117145641404181425[10] = 1;
   out_2117145641404181425[11] = 0;
   out_2117145641404181425[12] = 0;
   out_2117145641404181425[13] = 0;
   out_2117145641404181425[14] = 0;
   out_2117145641404181425[15] = 0;
   out_2117145641404181425[16] = 0;
   out_2117145641404181425[17] = 0;
   out_2117145641404181425[18] = 0;
   out_2117145641404181425[19] = 0;
   out_2117145641404181425[20] = 1;
   out_2117145641404181425[21] = 0;
   out_2117145641404181425[22] = 0;
   out_2117145641404181425[23] = 0;
   out_2117145641404181425[24] = 0;
   out_2117145641404181425[25] = 0;
   out_2117145641404181425[26] = 0;
   out_2117145641404181425[27] = 0;
   out_2117145641404181425[28] = 0;
   out_2117145641404181425[29] = 0;
   out_2117145641404181425[30] = 1;
   out_2117145641404181425[31] = 0;
   out_2117145641404181425[32] = 0;
   out_2117145641404181425[33] = 0;
   out_2117145641404181425[34] = 0;
   out_2117145641404181425[35] = 0;
   out_2117145641404181425[36] = 0;
   out_2117145641404181425[37] = 0;
   out_2117145641404181425[38] = 0;
   out_2117145641404181425[39] = 0;
   out_2117145641404181425[40] = 1;
   out_2117145641404181425[41] = 0;
   out_2117145641404181425[42] = 0;
   out_2117145641404181425[43] = 0;
   out_2117145641404181425[44] = 0;
   out_2117145641404181425[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_2117145641404181425[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_2117145641404181425[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2117145641404181425[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2117145641404181425[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_2117145641404181425[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_2117145641404181425[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_2117145641404181425[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_2117145641404181425[53] = -9.8000000000000007*dt;
   out_2117145641404181425[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_2117145641404181425[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_2117145641404181425[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2117145641404181425[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2117145641404181425[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_2117145641404181425[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_2117145641404181425[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_2117145641404181425[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2117145641404181425[62] = 0;
   out_2117145641404181425[63] = 0;
   out_2117145641404181425[64] = 0;
   out_2117145641404181425[65] = 0;
   out_2117145641404181425[66] = 0;
   out_2117145641404181425[67] = 0;
   out_2117145641404181425[68] = 0;
   out_2117145641404181425[69] = 0;
   out_2117145641404181425[70] = 1;
   out_2117145641404181425[71] = 0;
   out_2117145641404181425[72] = 0;
   out_2117145641404181425[73] = 0;
   out_2117145641404181425[74] = 0;
   out_2117145641404181425[75] = 0;
   out_2117145641404181425[76] = 0;
   out_2117145641404181425[77] = 0;
   out_2117145641404181425[78] = 0;
   out_2117145641404181425[79] = 0;
   out_2117145641404181425[80] = 1;
}
void h_25(double *state, double *unused, double *out_7879684886680657494) {
   out_7879684886680657494[0] = state[6];
}
void H_25(double *state, double *unused, double *out_5321245738136129279) {
   out_5321245738136129279[0] = 0;
   out_5321245738136129279[1] = 0;
   out_5321245738136129279[2] = 0;
   out_5321245738136129279[3] = 0;
   out_5321245738136129279[4] = 0;
   out_5321245738136129279[5] = 0;
   out_5321245738136129279[6] = 1;
   out_5321245738136129279[7] = 0;
   out_5321245738136129279[8] = 0;
}
void h_24(double *state, double *unused, double *out_7049418576302359049) {
   out_7049418576302359049[0] = state[4];
   out_7049418576302359049[1] = state[5];
}
void H_24(double *state, double *unused, double *out_3148596139130629713) {
   out_3148596139130629713[0] = 0;
   out_3148596139130629713[1] = 0;
   out_3148596139130629713[2] = 0;
   out_3148596139130629713[3] = 0;
   out_3148596139130629713[4] = 1;
   out_3148596139130629713[5] = 0;
   out_3148596139130629713[6] = 0;
   out_3148596139130629713[7] = 0;
   out_3148596139130629713[8] = 0;
   out_3148596139130629713[9] = 0;
   out_3148596139130629713[10] = 0;
   out_3148596139130629713[11] = 0;
   out_3148596139130629713[12] = 0;
   out_3148596139130629713[13] = 0;
   out_3148596139130629713[14] = 1;
   out_3148596139130629713[15] = 0;
   out_3148596139130629713[16] = 0;
   out_3148596139130629713[17] = 0;
}
void h_30(double *state, double *unused, double *out_4360242752360011888) {
   out_4360242752360011888[0] = state[4];
}
void H_30(double *state, double *unused, double *out_6208807994081805582) {
   out_6208807994081805582[0] = 0;
   out_6208807994081805582[1] = 0;
   out_6208807994081805582[2] = 0;
   out_6208807994081805582[3] = 0;
   out_6208807994081805582[4] = 1;
   out_6208807994081805582[5] = 0;
   out_6208807994081805582[6] = 0;
   out_6208807994081805582[7] = 0;
   out_6208807994081805582[8] = 0;
}
void h_26(double *state, double *unused, double *out_3532692919509055069) {
   out_3532692919509055069[0] = state[7];
}
void H_26(double *state, double *unused, double *out_1579742419262073055) {
   out_1579742419262073055[0] = 0;
   out_1579742419262073055[1] = 0;
   out_1579742419262073055[2] = 0;
   out_1579742419262073055[3] = 0;
   out_1579742419262073055[4] = 0;
   out_1579742419262073055[5] = 0;
   out_1579742419262073055[6] = 0;
   out_1579742419262073055[7] = 1;
   out_1579742419262073055[8] = 0;
}
void h_27(double *state, double *unused, double *out_1480643369266808348) {
   out_1480643369266808348[0] = state[3];
}
void H_27(double *state, double *unused, double *out_8383571305882230493) {
   out_8383571305882230493[0] = 0;
   out_8383571305882230493[1] = 0;
   out_8383571305882230493[2] = 0;
   out_8383571305882230493[3] = 1;
   out_8383571305882230493[4] = 0;
   out_8383571305882230493[5] = 0;
   out_8383571305882230493[6] = 0;
   out_8383571305882230493[7] = 0;
   out_8383571305882230493[8] = 0;
}
void h_29(double *state, double *unused, double *out_1323456700915679203) {
   out_1323456700915679203[0] = state[1];
}
void H_29(double *state, double *unused, double *out_8349810040957770090) {
   out_8349810040957770090[0] = 0;
   out_8349810040957770090[1] = 1;
   out_8349810040957770090[2] = 0;
   out_8349810040957770090[3] = 0;
   out_8349810040957770090[4] = 0;
   out_8349810040957770090[5] = 0;
   out_8349810040957770090[6] = 0;
   out_8349810040957770090[7] = 0;
   out_8349810040957770090[8] = 0;
}
void h_28(double *state, double *unused, double *out_810739376155975909) {
   out_810739376155975909[0] = state[0];
}
void H_28(double *state, double *unused, double *out_3267411023888239516) {
   out_3267411023888239516[0] = 1;
   out_3267411023888239516[1] = 0;
   out_3267411023888239516[2] = 0;
   out_3267411023888239516[3] = 0;
   out_3267411023888239516[4] = 0;
   out_3267411023888239516[5] = 0;
   out_3267411023888239516[6] = 0;
   out_3267411023888239516[7] = 0;
   out_3267411023888239516[8] = 0;
}
void h_31(double *state, double *unused, double *out_2070398505765154130) {
   out_2070398505765154130[0] = state[8];
}
void H_31(double *state, double *unused, double *out_5351891700013089707) {
   out_5351891700013089707[0] = 0;
   out_5351891700013089707[1] = 0;
   out_5351891700013089707[2] = 0;
   out_5351891700013089707[3] = 0;
   out_5351891700013089707[4] = 0;
   out_5351891700013089707[5] = 0;
   out_5351891700013089707[6] = 0;
   out_5351891700013089707[7] = 0;
   out_5351891700013089707[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_6029616048175551447) {
  err_fun(nom_x, delta_x, out_6029616048175551447);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_114715398326382947) {
  inv_err_fun(nom_x, true_x, out_114715398326382947);
}
void car_H_mod_fun(double *state, double *out_2186784832712934011) {
  H_mod_fun(state, out_2186784832712934011);
}
void car_f_fun(double *state, double dt, double *out_5798323304716649818) {
  f_fun(state,  dt, out_5798323304716649818);
}
void car_F_fun(double *state, double dt, double *out_2117145641404181425) {
  F_fun(state,  dt, out_2117145641404181425);
}
void car_h_25(double *state, double *unused, double *out_7879684886680657494) {
  h_25(state, unused, out_7879684886680657494);
}
void car_H_25(double *state, double *unused, double *out_5321245738136129279) {
  H_25(state, unused, out_5321245738136129279);
}
void car_h_24(double *state, double *unused, double *out_7049418576302359049) {
  h_24(state, unused, out_7049418576302359049);
}
void car_H_24(double *state, double *unused, double *out_3148596139130629713) {
  H_24(state, unused, out_3148596139130629713);
}
void car_h_30(double *state, double *unused, double *out_4360242752360011888) {
  h_30(state, unused, out_4360242752360011888);
}
void car_H_30(double *state, double *unused, double *out_6208807994081805582) {
  H_30(state, unused, out_6208807994081805582);
}
void car_h_26(double *state, double *unused, double *out_3532692919509055069) {
  h_26(state, unused, out_3532692919509055069);
}
void car_H_26(double *state, double *unused, double *out_1579742419262073055) {
  H_26(state, unused, out_1579742419262073055);
}
void car_h_27(double *state, double *unused, double *out_1480643369266808348) {
  h_27(state, unused, out_1480643369266808348);
}
void car_H_27(double *state, double *unused, double *out_8383571305882230493) {
  H_27(state, unused, out_8383571305882230493);
}
void car_h_29(double *state, double *unused, double *out_1323456700915679203) {
  h_29(state, unused, out_1323456700915679203);
}
void car_H_29(double *state, double *unused, double *out_8349810040957770090) {
  H_29(state, unused, out_8349810040957770090);
}
void car_h_28(double *state, double *unused, double *out_810739376155975909) {
  h_28(state, unused, out_810739376155975909);
}
void car_H_28(double *state, double *unused, double *out_3267411023888239516) {
  H_28(state, unused, out_3267411023888239516);
}
void car_h_31(double *state, double *unused, double *out_2070398505765154130) {
  h_31(state, unused, out_2070398505765154130);
}
void car_H_31(double *state, double *unused, double *out_5351891700013089707) {
  H_31(state, unused, out_5351891700013089707);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
