/*
 * MPAをArduinoからコントロールするにあたって
 * 長さ，収縮速度，印加圧力から発生張力を計算する関数群を
 * 標準化したMPAライブラリ，そのヘッダファイル．
 * 
 * ver.1.0 2019.02.03
 * 中西 大輔
 * 
 */

#ifndef MPA_H_INCLUDE
#define MPA_H_INCLUDE


//#define pi 3.14159265359



class MPA_Class{

  private:
    double L0; // 初期長さ
    double D0; // 直直径
    double t0; // ゴム厚
    double b;  // 繊維長さ
    double n;  // 巻き数
    double theta;  // 編み角
    double Vr; // ゴム部体積

    double R_air = 8.3143;             //気体定数[J/(k*mol)]
    double T = 28+271.15;              //空気温度（絶対温度）
    double P0 = 101325;                //大気圧[Pa]
    double C1 = 192000.0;  //旧MPAの値を暫定的に使用
    double C2 = 1300.0;  //旧MPAの値を暫定的に使用
    double dfmdv = -35.0;    // 減衰係数


  public:
    // コンストラクタ．MPA_Class(init_length：初期長さ, init_diameter：初期直径, rubber_thickness：ゴム厚, mesh_type：メッシュタイプ)
    MPA_Class(double init_length, double init_diameter, double rubber_thickness, int mesh_type);

    double get_initial_length();
    
    // 張力を算出する．calc_tension(obj：プライベート変数, L：長さ, V：収縮速度, P：印加圧力)
    double calc_tension(double L, double V, double P);

    double calc_dWdL(double L);

    double calc_dVdL(double L);

    double calc_TermC1(double L);

    double calc_TermC2(double L);

    double get_Vr();
    double get_initial_rength();
    double get_dfmdv();
    
    double calc_P_limit(double L);

    void check_MPA_length(double L, int no);
    void check_MPA_pressure(double P, int no);


};



#endif

