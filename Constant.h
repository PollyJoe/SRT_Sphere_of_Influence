#ifndef CONSTANT
#define CONSTANT
#include<math.h>

const int TreeNum = 5;
const int DebrisNum = 123;

//extern double debris_data[123][8];   //全局碎片信息
extern double tsp_data[130][2];         //全局碎片信息

const double Day2Second = 86400.0;      //天转s
const double mu =1.9891*6.67259*1E19   ;//太阳引力常数
const double req = 6378137.0;			//地球半径(m)
const double J2 = 1.08262668e-3;		//J2摄动
const double MJD_Init= 23467.0;         //初始时刻

//PI
const double pi = 3.1415926535897932384626433832795;	//圆周率
const double DPI = 3.1415926535897932384626433832795;
const double D2PI = 6.283185307179586476925286766559;
const double D2R = 0.017453292519943295769236907684886;
const double R2D = 57.295779513082320876798154814105;


const double EPSILON = 1.0e-14;
const double ee = 6.69437999014E-3;//地球椭圆形子午圈的偏心率平方
const double eep = 6.73949674228E-3;//地球椭圆形子午圈的偏心率平方


#endif

