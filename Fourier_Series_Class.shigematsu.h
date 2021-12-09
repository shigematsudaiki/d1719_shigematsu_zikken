/*
 * MPAをArduinoからコントロールするにあたって
 * 印加圧力とか，目標姿勢とか，フーリエ級数展開で表現する
 * いろんな波を表現するためのクラス
 * ver.1.0 2019.02.04
 * 中西 大輔
 * 
 */

#include <math.h>
#include <arduino.h>


#ifndef Fourier_Series_Class_INCLUDE
#define Fourier_Series_Class_INCLUDE

#define DATA_NUM 3

class Fourier_Series_Class{

  public:
    double center;
    double sin_amplitude[DATA_NUM];
    double sin_frequency[DATA_NUM];
    double cos_amplitude[DATA_NUM];
    double cos_frequency[DATA_NUM];

  public:
    Fourier_Series_Class(double center_ref, double *sin_amplitude_ref, double *sin_frequency_ref, double *cos_amplitude_ref, double *cos_frequency_ref);
    double calc_waveform_r(double t);
    double calc_waveform_drdt(double t);
    double calc_waveform_drddt(double t);
    
};

#endif


//        function [r, drdt, drddt] = calc_waveform(obj, t)
//            r = obj.center;
//            drdt = 0;
//            drddt = 0;
//            
//            for i=1:length(obj.sin_amplitude)
//                r = r + obj.sin_amplitude(i)*sin(obj.sin_frequency(i)*t);
//                drdt = drdt + obj.sin_amplitude(i)*obj.sin_frequency(i)*cos(obj.sin_frequency(i)*t);
//                drddt = drddt - obj.sin_amplitude(i)*obj.sin_frequency(i)^2*sin(obj.sin_frequency(i)*t);
//            end
//            for i=1:length(obj.cos_amplitude)
//                r = r + obj.cos_amplitude(i)*cos(obj.cos_frequency(i)*t);
//                drdt = drdt - obj.cos_amplitude(i)*obj.cos_frequency(i)*sin(obj.cos_frequency(i)*t);
//                drddt = drddt - obj.cos_amplitude(i)*obj.cos_frequency(i)^2*cos(obj.cos_frequency(i)*t);
//            end
//            
//            
//        end
//    end
//end
