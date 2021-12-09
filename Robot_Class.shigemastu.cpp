#include "Robot_Class.h"
#include "MPA.h"
#include "Fourier_Series_Class.h"


// コンストラクタ
Robot_Class::Robot_Class() {
  beta0 = calculate_beta(R0);
  eta10 = calculate_eta(R0, theta0);
  eta20 = calculate_eta(R0, theta0);
  lMh0 = calculate_lmh(R0);
  lMr0 = calculate_lmr(R0, theta0);
  gamma10 = calculate_gamma(R0, theta0);
  gamma20 = calculate_gamma(R0, theta0);
  gamma30 = calculate_gamma(R0, theta0);
  gamma40 = calculate_gamma(R0, theta0);
  gamma50 = calculate_gamma(R0, theta0);
};

double Robot_Class::calculate_beta(double R, double dzeta1, double dzeta2); {
  return beta = pi-(dzeta1+dzeta2);
}
double Robot_Class::calculate_Lh(double R, double gamma1, double lmh, double theta *MPA_h) {
  double Lh0 = MPA_h->get_init_length();
  return Lh = r1*(gamma10-gamma1) + lmh - lMh0 + Lh0;
}

double Robot_Class::calculate_Lr(double R, double gamma2, double lmr, double theta *MPA_r, double eta2) {
  double Lr0 = MPA_r->get_init_length();
  return Lr = r2*(eta20-eta2) + r2*(gamma20-gamma2) + (lmr - lMr0) + Lr0;
}

double Robot_Class::calculate_lmh(double R) {
  double la = sqrt((Lu*Lo^2+Lu*L3^2-L3*Lo^2-L3*Lu^2+L3*R^2)/Lu);
  return lmh = sqrt(la^2-r1^2);
} 

double Robot_Class::double calculate_lmr(double R, double theta, double eat1) {
  double lb = sqrt(L2^2+Lo^2 - 2*L2*Lo*cos(eta1)); 
  return lmr = sqrt(lb^2-r2^2);
}

double Robot_Class::calculate_eta1(double R, double theta, double dzeta1, double dzeta2, double beta) {
  if theta < 0 && dzeta1 < abs(theta)
    return eta1 = pi/2 + abs(theta+dzeta1);
  else
    return eta1 = pi/2-(dzeta1+theta);
}

double Robot_Class::calculate_eta2(double R, double theta, double dzeta1, double dzeta2, double beta) {
  return eta2 = beta-eta1;
}

double Robot_Class::calculate_dzeta1(double R) {
  return dzeta1  = acos((R^2+Lo^2-Lu^2)/(2*R*Lo));
}

double Robot_Class::calculate_dzeta2(double R) {
  return dzeta2  = acos((R^2+Lu^2-Lo^2)/(2*R*Lu));
}

double Robot_Class::calculate_P(double A, double phi) {
  double a1 = A(1);
  double a2 = A(2);
  return P = a1 + a2*(1-sin(phi));
}

double Robot_Class::calculate_vm_r(double theta, double dthetadt, double R, double dRdt, double GrR, double GrT) {
  return vm_r = GrR*dRdt + GrT*R*dthetadt;
}
double Robot_Class::calculate_vm_h(double theta, double dthetadt, double R, double dRdt, double GhR, double GhT) {
  return vm_h = GhR*dRdt + GhT*R*dthetadt;
}

double Robot_Class::calculate_Vm_h(double theta, double dthetadt, double R, double dRdt) {
  return Vm_h = -1*(L3*R*dRdt/(Lu*sqrt(Lo^2 + L3^2 - (L3*Lo^2/Lu) - L3*Lu - r1^2 + (L3*R^2)/Lu)));
}

double Robot_Class::calculate_Vm_r(double theta, double dthetadt, double R, double dRdt) {
  double Lr_R = -1*R*L2*Lo*sin(theta-0.5*acos((Lo^2+Lu^2-R^2)/(2*Lo*Lu)))/(2*Lo*Lu*sqrt(L2^2+Lo^2-r2^2-2*L2*Lo*cos(theta-0.5*acos((Lo^2+Lu^2-R^2)/(2*Lo*Lu))))*sqrt(1-((Lo^2+Lu^2-R^2)/(2*Lo*Lu))^2));
  double Lr_theta = L2*Lo*sin(theta-0.5*acos((Lo^2+Lu^2-R^2)/(2*Lo*Lu)))/sqrt(L2^2+Lo^2-r2^2-2*L2*Lo*cos(theta-0.5*acos((Lo^2+Lu^2-R^2)/(2*Lo*Lu))));
  double Lr_t = Lr_R*dRdt + Lr_theta*dthetadt;
  return Vm_r = -1*Lr_t;
} 

double Robot_Class::calculate_la(double R, double theta) {
  return la = sqrt((Lu*Lo^2+Lu*L3^2-L3*Lo^2-L3*Lu^2+L3*R^2)/Lu);
}

double Robot_Class::calculate_lb(double R, double theta, double eta1) {
  return lb = sqrt(L2^2+Lo^2 - 2*L2*Lo*cos(eta1));
}

double Robot_Class::calculate_gamma1(double R, double theta, double gamma3, double eta2, double lmh, double la) {
  return gamma1 = pi - (eta2+gamma3);
}

double Robot_Class::calculate_gamma2(double R, double theta, double lmr, double lb) {
  return gamma2 = pi - (atan(r2/lmr) + acos((L2^2+lb^2-Lo^2)/(2*L2*lb)));
}

double Robot_Class::calculate_gamma3(double R, double theta, double lmh, double la) {
  return gamma3 = atan(r1/lmh) + acos((L3^2+la^2-Lo^2)/(2*L3*la));
}

double Robot_Class::calculate_gamma4(double R, double theta) {
  return gamma4 = asin(r2/L4);
}

double Robot_Class::calculate_gamma5(double R, double theta) {
  return gamma5 = asin(r1/L1);
}
double Robot_Class::calculate_Grx(double R, double theta, double eta1, double eta2, double beta, double gamma1, double gamma2, double gamma3, double gamma4, double gamma5) {
  return Grx = (-cos(eta1)*cos(eta2)*sin(gamma4-eta2) + cos(eta1)*cos(eta2)*(sin(gamma4-eta2)+sin(gamma2)) +sin(eta1)*cos(eta2)*cos(gamma4-eta2)-sin(eta1)*cos(eta2)*(cos(gamma4-eta2)+cos(gamma2)) +(L4*cos(eta1)*sin(gamma4))/Lo)/(sin(eta1+eta2));
}

double Robot_Class::calculate_Ghx(double R, double theta, double eta1, double eta2, double beta, double gamma1, double gamma2, double gamma3, double gamma4, double gamma5) {
  return Ghx = ( cos(eta1)*cos(eta2)*sin(gamma3+eta2)+sin(eta1)*cos(eta2)*cos(gamma3+eta2)-(L3*cos(eta1)*sin(gamma3))/Lo)/(sin(eta1+eta2));
}

double Robot_Class::calculate_Gry(double R, double theta, double eta1, double eta2, double beta, double gamma1, double gamma2, double gamma3, double gamma4, double gamma5) {
  return Gry = (-cos(eta1)*sin(eta2)*sin(gamma4-eta2)+cos(eta1)*sin(eta2)*(sin(gamma4-eta2)+sin(gamma2))+sin(eta1)*sin(eta2)*cos(gamma4-eta2)-sin(eta1)*sin(eta2)*(cos(gamma4-eta2)+cos(gamma2))-(L4*sin(eta1)*sin(gamma4))/Lo)/(sin(eta1+eta2));
} 

double Robot_Class::calculate_Ghy(double R, double theta, double eta1, double eta2, double beta, double gamma1, double gamma2, double gamma3, double gamma4, double gamma5) {
  return Ghy = (cos(eta1)*sin(eta2)*sin(gamma3+eta2)+sin(eta1)*sin(eta2)*cos(gamma3+eta2)+(L3*sin(eta1)*sin(gamma3))/Lo)/(sin(eta1+eta2));

double Robot_Class::calculate_GrR(double R, double theta, double Grx, double Ghx, double Gry, double Ghy) {
  return GrR = Grx*sin(theta) - Gry*cos(theta);
}

double Robot_Class::calculate_GhR(double R, double theta, double Grx, double Ghx, double Gry, double Ghy) {
  return GhR = Ghx*sin(theta) - Ghy*cos(theta);
}

double Robot_Class::calculate_GrT(double R, double theta, double Grx, double Ghx, double Gry, double Ghy) {
  return GrT = Grx*cos(theta) + Gry*sin(theta);
}

double Robot_Class::calculate_GhT(double R, double theta, double Grx, double Ghx, double Gry, double Ghy) {
  return GhT = Ghx*cos(theta) + Ghy*sin(theta);
}



