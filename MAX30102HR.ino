#include <Wire.h>
#include "MAX30105.h"

#include "heartRate.h"

#include "arduinoFFT.h"
#include "KickSort.h"

MAX30105 particleSensor;

const byte RATE_SIZE = 4; //Increase this for more averaging. 4 is good.
byte rates[RATE_SIZE]; //Array of heart rates
byte rateSpot = 0;
long lastBeat = 0; //Time at which the last beat occurred

float beatsPerMinute;
int beatAvg;


#include <FIR.h>
#define FILTERTAPS 35
FIR fir;

float data[51] = { 0.2, -0.01757378,  0.25261   ,  0.50501472,  0.28169386,
                   0.73591371,  0.67214444,  0.63916107,  1.04340099,  0.75089683,
                   0.97104406,  1.1070092 ,  0.79944834,  1.15681282,  0.95424177,
                   0.83291634,  1.1026499 ,  0.68138434,  0.80743739,  0.79725188,
                   0.39288974,  0.65097057,  0.32531128,  0.14508086,  0.32099195,
                   -0.17069916, -0.07172355, -0.14896219, -0.55850617, -0.30387694,
                   -0.64596943, -0.77404044, -0.57958674, -1.0231335 , -0.8365037 ,
                   -0.86670651, -1.16875508, -0.81455928, -1.07313216, -1.05899594,
                   -0.76797161, -1.09233578, -0.76336743, -0.70354745, -0.86712553,
                   -0.40092413, -0.57401086, -0.4319664 , -0.07473948, -0.32007654,
                   0.0936594
                 };

#define movingAvg_HR_DC_step 100
float movingAvgArray[movingAvg_HR_DC_step] = {0};
int index_AvgArray = 0;
float  movingAvg_forDC = 100000;
float movingAverage_smoother = 0;
float finalSignal = 0;
long previousIrValue= 0;

#define finalHR_AvgSteps 3  //min = 3    if finalHR_AvgSteps < 3 : line of "int result_HR = ..." should be changed

float array4Avg_HR[finalHR_AvgSteps] = {0};
int counter_avg = 0;
int numberOfNoSkins = 0;

arduinoFFT FFT = arduinoFFT(); /* Create FFT object */
/*
  These values can be changed in order to evaluate the functions
*/

#define SAMPLES 512             //Must be a power of 2
#define SAMPLING_FREQUENCY 500

unsigned int sampling_period_us;
unsigned long microseconds;

double vReal[SAMPLES];
double vImag[SAMPLES];


void setup()
{
  Serial1.begin(38400);
  Serial1.println("Initializing...");

  // Initialize sensor
  if (!particleSensor.begin(Wire, I2C_SPEED_FAST)) //Use default I2C port, 400kHz speed
  {
    Serial1.println("MAX30105 was not found. Please check wiring/power. ");
    while (1);
  }
  Serial1.println("Place your index finger on the sensor with steady pressure.");


  particleSensor.setup(); //Configure sensor with default settings
  particleSensor.setPulseAmplitudeRed(0x0A); //Turn Red LED to low to indicate sensor is running
  particleSensor.setPulseAmplitudeGreen(0); //Turn off Green LED

  float coef[FILTERTAPS] = {0.0011,   -0.0048,    0.0093  , -0.0088  , -0.0013 ,   0.0141  , -0.0136  , -0.0053 ,   0.0216  , -0.0079 ,  -0.0287 ,   0.0401,
    0.0074  , -0.0705  ,  0.0502   , 0.1007  , -0.2996,    0.3919  , -0.2996 ,   0.1007   , 0.0502 ,  -0.0705 ,   0.0074 ,   0.0401,
   -0.0287 ,  -0.0079 ,   0.0216  , -0.0053   ,-0.0136  ,  0.0141   ,-0.0013 ,  -0.0088 ,   0.0093  , -0.0048,    0.0011};
  
  fir.setCoefficients(coef);

    //declare gain coefficient to scale the output back to normal
  float gain = 1; // set to 1 and input unity to see what this needs to be
  fir.setGain(gain);
  
  sampling_period_us = round(1000000 * (1.0 / (SAMPLING_FREQUENCY / 10)));

}

void loop()
{
    for(int i=0; i<SAMPLES; i++)
    {
        microseconds = micros();    //Overflows after around 70 minutes!
        long irValue = (particleSensor.getIR());
        
        if (irValue<50000){ //if the sensor is not on the skin
          finalSignal = 0;
          movingAverage_smoother = 0;
          movingAvg_forDC = 0;
          numberOfNoSkins ++;
          }
        else if (previousIrValue < 50000 && irValue > 50000){
          movingAverage_smoother = irValue;
          movingAvg_forDC = irValue;
        }
        else{
          long lowPassed = lowPassFIRFilter(irValue);  //Sparkfun Max30102 library, lowPassFIRFilter function
          movingAverage_smoother = movingAverage_smoother * 0.6 + 0.4*lowPassed;
          movingAvg_forDC = 0.95 * movingAvg_forDC + 0.05 * movingAverage_smoother;
        }
        
             
        if (abs(movingAverage_smoother -  movingAvg_forDC) < 2000){
          finalSignal = movingAverage_smoother -  movingAvg_forDC;
          }
        else{
         // finalSignal = movingAvg_forDC; finalSignal = finalSiganl! no change
          }
       
        vReal[i] = finalSignal;
        vImag[i] = 0;
        previousIrValue = irValue;
        
        while(micros() < (microseconds + sampling_period_us)){
        }
    }
    if (numberOfNoSkins > 100){
      Serial1.println("No skin");
      }
    else{
      FFT.Windowing(vReal, SAMPLES, FFT_WIN_TYP_HAMMING, FFT_FORWARD);
      FFT.Compute(vReal, vImag, SAMPLES, FFT_FORWARD);
      FFT.ComplexToMagnitude(vReal, vImag, SAMPLES);
      vReal[0] = 0;
      vReal[1] = 0;vReal[2] = 0;vReal[3] = 0;vReal[4] = 0;
      vReal[5] = 0;
      vReal[6] = 0;vReal[7] = 0;
      vReal[8] = 0;
      double peak = FFT.MajorPeak(vReal, SAMPLES, SAMPLING_FREQUENCY); //if other parts of the code takes noticable time, the SAMPLING_FREQUENCY will be changed and it should be considered.
      array4Avg_HR[counter_avg] = peak / 10;  //because at first we said freq = 500 but it is 50
      counter_avg ++;
    }
    numberOfNoSkins = 0;


    if (counter_avg == finalHR_AvgSteps){
    
        KickSort<float>::bubbleSort(array4Avg_HR, finalHR_AvgSteps);
        Serial1.println(array4Avg_HR[0]);
        Serial1.println(array4Avg_HR[1]);
        Serial1.println(array4Avg_HR[2]);
        int result_HR = (array4Avg_HR[2]*0.6 + array4Avg_HR[1]*0.4);
        counter_avg = 0;

        Serial1.println(result_HR);
    } 
     
    
}
