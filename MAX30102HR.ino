/*
  Optical Heart Rate Detection (PBA Algorithm) using the MAX30105 Breakout
  By: Nathan Seidle @ SparkFun Electronics
  Date: October 2nd, 2016
  https://github.com/sparkfun/MAX30105_Breakout

  This is a demo to show the reading of heart rate or beats per minute (BPM) using
  a Penpheral Beat Amplitude (PBA) algorithm.

  It is best to attach the sensor to your finger using a rubber band or other tightening
  device. Humans are generally bad at applying constant pressure to a thing. When you
  press your finger against the sensor it varies enough to cause the blood in your
  finger to flow differently which causes the sensor readings to go wonky.

  Hardware Connections (Breakoutboard to Arduino):
  -5V = 5V (3.3V is allowed)
  -GND = GND
  -SDA = A4 (or SDA)
  -SCL = A5 (or SCL)
  -INT = Not connected

  The MAX30105 Breakout can handle 5V or 3.3V I2C logic. We recommend powering the board with 5V
  but it will also run at 3.3V.
*/

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

// Make two instances of the FIR filter. The first is a 13 element float
// and the second is a 10 point float. The 13 element filter is a 2 Hz low-pass
// for a 10 SPS signal and the 10 element is a moving average over one second.

/*
FIR<float, 13> fir_lp;
FIR<float, 50> fir_avg;
FIR<float, 27> fir_hp;
*/
// Some data that is a 0.2 Hz sine wave with a 0.3 Hz sine wave noise on it.
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
/*
  float coef_lp[13] = { 660, 470, -1980, -3830, 504, 10027, 15214,
                        10027, 504, -3830, -1980, 470, 660
                      };

  // For a moving average we use all ones as coefficients.
  float coef_avg[50] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};

  float coef_hp[27] = {-0.0040 ,   0.0097 ,  -0.0074 ,  -0.0094  ,  0.0265  , -0.0173   ,-0.0183,    0.0356,    0.0059 ,  -0.0713 ,   0.0580 ,   0.0920,
   -0.2981 ,  0.3952   ,-0.2981   , 0.0920  ,  0.0580 ,  -0.0713  ,  0.0059,    0.0356 ,  -0.0183,   -0.0173  ,  0.0265   ,-0.0094, -0.0074    ,0.0097  , -0.0040};
      
  // Set the coefficients
  fir_lp.setFilterCoeffs(coef_lp);
  fir_avg.setFilterCoeffs(coef_avg);
    fir_hp.setFilterCoeffs(coef_hp);

  // Set the gain
  Serial1.print("Low Pass Filter Gain: ");
  Serial1.println(fir_lp.getGain());
  Serial1.print("Moving Average Filter Gain: ");
  Serial1.println(fir_avg.getGain());
    fir_hp.getGain();
*/
 /* float coef[FILTERTAPS] ={-0.0040 ,   0.0097 ,  -0.0074 ,  -0.0094  ,  0.0265  , -0.0173   ,-0.0183,    0.0356,    0.0059 ,  -0.0713 ,   0.0580 ,   0.0920,
   -0.2981 ,  0.3952   ,-0.2981   , 0.0920  ,  0.0580 ,  -0.0713  ,  0.0059,    0.0356 ,  -0.0183,   -0.0173  ,  0.0265   ,-0.0094, -0.0074    ,0.0097  , -0.0040};
  */
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
  //long irValue = particleSensor.getIR();
  /*
    //Serial1.print("FIR:");
    float FIR_result = fir_lp.processReading(irValue);

    /////////my moving average//////////
    movingAvgArray[index_AvgArray] = irValue;
    index_AvgArray++;
    for (int i = 0; i<movingAvg_HR_DC_step; i++){
    sum = sum + movingAvgArray[i];
    }
    movingAvg_forDC = sum / movingAvg_HR_DC_step;
    if (index_AvgArray >= 100){
    index_AvgArray = 0;
    }
  */
  
    //movingAvg_forDC = 0.99 * movingAvg_forDC + 0.01 * FIR_result;
    //Serial1.println(irValue );
    
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
      double peak = FFT.MajorPeak(vReal, SAMPLES, SAMPLING_FREQUENCY);
      array4Avg_HR[counter_avg] = peak / 10;  //because at first we said freq = 500 but it is 50
      counter_avg ++;
    }
    numberOfNoSkins = 0;


        /*
    for(int i=0; i<(SAMPLES/2); i++)
    {

        //Serial.print((i * 1.0 * SAMPLING_FREQUENCY) / SAMPLES, 1);
        //Serial.print(" ");
        Serial1.println(vReal[i], 1);    //View only this line in serial plotter to visualize the bins
    }
     Serial1.print("Peak: ");
    Serial1.println(peak);     //Print out what frequency is the most dominant.
    */
    
    if (counter_avg == finalHR_AvgSteps){
    
        KickSort<float>::bubbleSort(array4Avg_HR, finalHR_AvgSteps);
        Serial1.println(array4Avg_HR[0]);
        Serial1.println(array4Avg_HR[1]);
        Serial1.println(array4Avg_HR[2]);
        int result_HR = (array4Avg_HR[0]*0.6 + array4Avg_HR[1]*0.4);// + array4Avg_HR[2]*0.2) * 60;
        counter_avg = 0;

        Serial1.println(result_HR);
    } 
     
    
}
