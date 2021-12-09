 

#include "Fourier_Series_Class.h"


Fourier_Series_Class::Fourier_Series_Class(double center_ref, double *sin_amplitude_ref, double *sin_frequency_ref, double *cos_amplitude_ref, double *cos_frequency_ref){
  center = center_ref;
  for(int i=0;i<DATA_NUM;i++){
    sin_amplitude[i] = sin_amplitude_ref[i];
    sin_frequency[i] = sin_frequency_ref[i];
    cos_amplitude[i] = cos_amplitude_ref[i];
    cos_frequency[i] = cos_frequency_ref[i];
  };
}

double Fourier_Series_Class::calc_waveform_r(double t){
  double r = center;
  for(int i=0;i<DATA_NUM;i++){
    r = r + sin_amplitude[i]*sin(sin_frequency[i]*t) + cos_amplitude[i]*cos(cos_frequency[i]*t);
  };
  return r;
}

double Fourier_Series_Class::calc_waveform_drdt(double t){
  double drdt = 0;
  for(int i=0;i<DATA_NUM;i++){
    drdt = drdt + sin_amplitude[i]*sin_frequency[i]*cos(sin_frequency[i]*t) - cos_amplitude[i]*cos_frequency[i]*sin(cos_frequency[i]*t);
  };
  return drdt;
}

double Fourier_Series_Class::calc_waveform_drddt(double t){
  double drddt = 0;
  for(int i=0;i<DATA_NUM;i++){
    drddt = drddt - sin_amplitude[i]*pow(sin_frequency[i],2)*sin(sin_frequency[i]*t) - cos_amplitude[i]*pow(cos_frequency[i],2)*cos(cos_frequency[i]*t);
  };
  return drddt;
}

