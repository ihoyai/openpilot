#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_4230691914079283868);
void live_err_fun(double *nom_x, double *delta_x, double *out_74140843859662143);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_5960310541044414783);
void live_H_mod_fun(double *state, double *out_388164681519952464);
void live_f_fun(double *state, double dt, double *out_5945865858291546245);
void live_F_fun(double *state, double dt, double *out_4275040546717650741);
void live_h_4(double *state, double *unused, double *out_2049851605477456516);
void live_H_4(double *state, double *unused, double *out_8819864892898239864);
void live_h_9(double *state, double *unused, double *out_7309656425162044801);
void live_H_9(double *state, double *unused, double *out_2339660245546864282);
void live_h_10(double *state, double *unused, double *out_3857645970319454097);
void live_H_10(double *state, double *unused, double *out_834033692795469833);
void live_h_12(double *state, double *unused, double *out_3865856022971558194);
void live_H_12(double *state, double *unused, double *out_4607422772779349957);
void live_h_31(double *state, double *unused, double *out_2476635822166799026);
void live_H_31(double *state, double *unused, double *out_6260217123438704376);
void live_h_32(double *state, double *unused, double *out_7890212855320444047);
void live_H_32(double *state, double *unused, double *out_4319378786105695575);
void live_h_13(double *state, double *unused, double *out_8684233441508878315);
void live_H_13(double *state, double *unused, double *out_5175995149130776497);
void live_h_14(double *state, double *unused, double *out_7309656425162044801);
void live_H_14(double *state, double *unused, double *out_2339660245546864282);
void live_h_33(double *state, double *unused, double *out_3859204275898927551);
void live_H_33(double *state, double *unused, double *out_3109660118799846772);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}