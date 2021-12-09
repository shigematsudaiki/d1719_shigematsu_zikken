/*
 * MPAをArduinoからコントロールするにあたって
 * 2DOF legged modelをクラス化する試み
 * ver.1.0 2019.02.04
 * 中西 大輔
 * 
 */

#include <math.h>
#include <arduino.h>
#include "MPA.h"
#include "Fourier_Series_Class.h"

#ifndef ONE_DOF_LEGGED_MODL_H_INCLUDE
#define ONE_DOF_LEGGED_MODL_H_INCLUDE

//#define PI 3.141592653589793   

class Robot_Class{

  public:
    double M = 2.0;       // 代表の重さ
    double g = 9.80665;
    double Lu = 0.500;
    double Lo = 0.500;
    double L1 = 0.180;
    double L2 = 0.100;
    double L3 = 0.080;
    double L4 = 0.110;
    double r1 = 0.065;    //プーリ半径
    double r2 = 0.040;    //プーリ半径


    // 初期姿勢，ワイヤ長関係 
    double R0 = 800/1000;
    double theta0 = 0*pi/180;
    double beta0;
            

  public:
    Robot_Class();

  
    double calculate_Lh(double R, double gamma1, double lmh, double theta *MPA_h);
    double calculate_Lr(double R, double gamma2, double lmr, double theta *MPA_r);
    double calculate_lmh(double R);
    double calculate_lmr(double R, double theta, double eat1);
    double calculate_beta(double R, double dzeta1, double dzeta2);
    double calculate_eta1(double R, double theta, double dzeta1, double dzeta2);
    double calculate_eta2(double R, double theta, double dzeta1, double dzeta2);
    double calculate_dzeta1(double R);
    double calculate_dzeta2(double R);
    double calculate_P(double A, double phi);
    double calculate_vm_r(double theta, double dthetadt, double R, double dRdt, double GrR, double GrT);
    double calculate_vm_h(double theta, double dthetadt, double R, double dRdt, double GhR, double GhT);
    double calculate_Vm_h(double theta, double dthetadt, double R, double dRdt);
    double calculate_Vm_r(double theta, double dthetadt, double R, double dRdt);
    double calculate_la(double R, double theta);
    double calculate_lb(double R, double theta, double eta1);
    double calculate_gamma1(double R, double theta, double gamma3, double eta2, double lmh, double la);
    double calculate_gamma2(double R, double theta, double lmr, double lb);
    double calculate_gamma3(double R, double theta, double lmh, double la);
    double calculate_gamma4(double R, double theta);
    double calculate_gamma5(double R, double theta);
    double calculate_Grx(double R, double theta, double eta1, double eta2, double beta, double gamma1, double gamma2, double gamma3, double gamma4, double gamma5);
    double calculate_Ghx(double R, double theta, double eta1, double eta2, double beta, double gamma1, double gamma2, double gamma3, double gamma4, double gamma5);
    double calculate_Gry(double R, double theta, double eta1, double eta2, double beta, double gamma1, double gamma2, double gamma3, double gamma4, double gamma5);
    double calculate_Ghy(double R, double theta, double eta1, double eta2, double beta, double gamma1, double gamma2, double gamma3, double gamma4, double gamma5);
    double calculate_GrR(double R, double theta, double Grx, double Ghx, double Gry, double Ghy);
    double calculate_GhR(double R, double theta, double Grx, double Ghx, double Gry, double Ghy);
    double calculate_GrT(double R, double theta, double Grx, double Ghx, double Gry, double Ghy);
    double calculate_GhT(double R, double theta, double Grx, double Ghx, double Gry, double Ghy);

};











#endif
