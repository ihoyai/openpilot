#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_6029616048175551447);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_114715398326382947);
void car_H_mod_fun(double *state, double *out_2186784832712934011);
void car_f_fun(double *state, double dt, double *out_5798323304716649818);
void car_F_fun(double *state, double dt, double *out_2117145641404181425);
void car_h_25(double *state, double *unused, double *out_7879684886680657494);
void car_H_25(double *state, double *unused, double *out_5321245738136129279);
void car_h_24(double *state, double *unused, double *out_7049418576302359049);
void car_H_24(double *state, double *unused, double *out_3148596139130629713);
void car_h_30(double *state, double *unused, double *out_4360242752360011888);
void car_H_30(double *state, double *unused, double *out_6208807994081805582);
void car_h_26(double *state, double *unused, double *out_3532692919509055069);
void car_H_26(double *state, double *unused, double *out_1579742419262073055);
void car_h_27(double *state, double *unused, double *out_1480643369266808348);
void car_H_27(double *state, double *unused, double *out_8383571305882230493);
void car_h_29(double *state, double *unused, double *out_1323456700915679203);
void car_H_29(double *state, double *unused, double *out_8349810040957770090);
void car_h_28(double *state, double *unused, double *out_810739376155975909);
void car_H_28(double *state, double *unused, double *out_3267411023888239516);
void car_h_31(double *state, double *unused, double *out_2070398505765154130);
void car_H_31(double *state, double *unused, double *out_5351891700013089707);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}