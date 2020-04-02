#include <iostream>
#include <string>
#include <math.h>

double CalculateMean(double values[], int length) {
  double sum = 0;
  for(int i=0; i<length; i++)
    sum+=values[i];
  return (sum / length);
}

double CalculateVariance(double values[], int length) {
  double mean = CalculateMean(values, length);
  double temp = 0;
  for(int i=0; i<length; i++){
    temp += (values[i]-mean)*(values[i]-mean);
  }
  return temp / length;
}

double CalculateSampleVariance(double values[], int length) {
  double mean = CalculateMean(values, length);
  double temp = 0;
  for(int i=0; i<length; i++)
    temp += (values[i]-mean)*(values[i]-mean);
  return temp / (length - 1);
}

double GetStandardDeviation(double values[], int length) {
  return sqrt(CalculateVariance(values, length));
}

double GetSampleStandardDeviation(double values[], int length) {
  return sqrt(CalculateSampleVariance(values, length));
}

int main() {
  double sample[] = {
    15.17, 16.94, 14.94, 14.99, 13.77, 13.75,
    12.67, 12.14, 12.59, 12.48, 14.81, 14.29,
    12.74, 12.52, 11.65, 12.24, 11.42, 12.25,
    12.72, 11.64, 11.09, 11.22, 11.50, 11.36,
    11.84, 12.18, 11.04, 10.90, 11.80, 11.84,
  };
  int sample_length = sizeof(sample) / sizeof(sample[0]); // usual trick

  double mean =  CalculateMean(sample, sample_length);
  double variance = CalculateVariance(sample, sample_length);
  double samplevariance =  CalculateSampleVariance(sample, sample_length);
  double sampledevi =  GetStandardDeviation(sample, sample_length);
  double devi =  GetSampleStandardDeviation(sample, sample_length);

  std::cout << "Total Numbers\t\t\t: " << sample_length << "\n"
            << "Mean\t\t\t\t: " << mean << "\n"
            << "Population Variance\t\t: " << variance << "\n"
            << "Sample variance\t\t\t: " << samplevariance << "\n"
            << "Population Standard Deviation\t: " << devi << "\n"
            << "Sample Standard Deviation\t: " << sampledevi << "\n";

  return 0;
}

