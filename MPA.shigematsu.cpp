#include <arduino.h>
#include "MPA.h"

//コンストラクタ

MPA_Class::MPA_Class(double init_length, double init_diameter, double rubber_thickness, int mesh_type)
{
  // MPAの定数パラメータを初期化していく
  L0 = init_length;
  D0 = init_diameter;
  t0 = rubber_thickness;
  double b_m, L_m, D_m, n_m;

  if(mesh_type==1){
    // type1:現行MPAスリーブ
    b_m  = 54.6/1000.0;                    //実測繊維長さ
    L_m = 51/1000.0;                      //実測メッシュ長さ
    D_m = 10/1000.0;                      //実測メッシュ直径
  }else if(mesh_type==2){
    // type2:大成製作所メッシュ
    b_m = 36.55/1000.0;                   //実測繊維長さ
    L_m = 34.5/1000.0;                    //実測メッシュ長さ
    D_m = 12/1000.0;                      //実測メッシュ直径
  }

  n_m = pow((pow(b_m,2)-pow(L_m,2)),0.5)/(D_m*PI);   //推定巻数
  theta = acos(pow((pow(b_m,2)-pow((D0*n_m*PI),2)),0.5)/b_m);     //初期メッシュ角
  n = n_m*L0/pow((pow(b_m,2)-pow((D0*n_m*PI),2)),0.5);       //巻数
  b = L0/cos(theta);                          //繊維長さ
  Vr  = PI*(D0*t0 - pow(t0,2))*L0;         //ゴム管の体積

};

double MPA_Class::get_initial_length(){
  return L0;
};

double MPA_Class::calc_tension(double L, double V, double P){
  check_MPA_length(L, 0);
  check_MPA_pressure(P, 0);
  double dVdL, dWdL, fm;
  dVdL = calc_dVdL(L);
  dWdL = calc_dWdL(L);
  fm = -P*dVdL + Vr*dWdL + V*dfmdv;
  if(fm > 0){
    return fm;
  }else{
    while(1){
      //Serial.println("Actuator tension is Negative! : in calc_tension");
      return 0;
      //delay(1000);
    };  
  };
};


double MPA_Class::calc_dVdL(double L){
  check_MPA_length(L, 0);
  //dV/dL:内側MPA長さ変化に対するMPA体積変化率
  return (pow(b,2)-3*pow(L,2))/(4*PI*pow(n,2));
};

double MPA_Class::calc_dWdL(double L){
  // MPA長さが正のとき，それぞれの値を計算．負のときエラーを出す．
  check_MPA_length(L, 0);
  double TermC1, TermC2;
  //TermC1:内側MPAに関して，dWfdLのC1がかかった項  
  TermC1 = 2*L/pow(L0,2) - 2*L/pow((D0*n*PI),2) - 2*pow((L0*D0*n*PI),2)*(pow(b,2)-2*pow(L,2))/(pow(L,3)*pow((pow(b,2)-pow(L,2)),2));  
  //TermC2:内側MPAに関して，dWfdLのC2がかかった項
  TermC2 = -2*pow(L0,2)/pow(L,3) + 2*pow((D0*n*PI),2)*L/(pow((pow(b,2)-pow(L,2)),2)) + (2*L*pow(b,2)-4*pow(L,3))/pow((L0*D0*n*PI),2);
  //dWdL
  return C1*TermC1 + C2*TermC2;
};

double MPA_Class::calc_TermC1(double L){
  double TermC1;
  TermC1 = 2*L/pow(L0,2) - 2*L/pow((D0*n*PI),2) - 2*pow((L0*D0*n*PI),2)*(pow(b,2)-2*pow(L,2))/(pow(L,3)*pow((pow(b,2)-pow(L,2)),2));
  //dWdL
  return TermC1;
};

double MPA_Class::calc_TermC2(double L){
  double TermC2;
  //TermC2:内側MPAに関して，dWfdLのC2がかかった項
  TermC2 = -2*pow(L0,2)/pow(L,3) + 2*pow((D0*n*PI),2)*L/(pow((pow(b,2)-pow(L,2)),2)) + (2*L*pow(b,2)-4*pow(L,3))/pow((L0*D0*n*PI),2);
  //dWdL
  return TermC2;
}

double MPA_Class::get_Vr(){
  return Vr;
};

// MPA長から算出される，張力が0になるために最低限必要な印加圧力を計算する
double  MPA_Class::calc_P_limit(double L){
  // 長さが正の時，値を計算する．またL < L0の場合P_limitは負になるので，0を代入．
  check_MPA_length(L, 0);
  double dVdL, dWdL, P_limit;
  dVdL = calc_dVdL(L);
  dWdL = calc_dWdL(L);
  
  if(Vr*dWdL/dVdL>0){
    return (Vr*dWdL/dVdL>0);
  }else{
    return 0;
  };
};

void MPA_Class::check_MPA_length(double L, int no){
  if(L<0){
    while(1){
      Serial.print("Actuator Length is Negative ! : in ");
      Serial.println(no);
      delay(1000);
    };    
  }else if(L>b){
    while(1){
      Serial.print("Actuator Length exceeds Fiber Length !: in ");
      Serial.println(no);
      delay(1000);
    };  
  };
};


void MPA_Class::check_MPA_pressure(double P, int no){
  if(P<0){
    while(1){
      Serial.print("Actuator pressure is Negative! : in ");
      Serial.println(no);
      delay(1000);
    };    
  };
};

double MPA_Class::get_initial_rength(){
  return L0;
}

double MPA_Class::get_dfmdv(){
  return L0;
}
