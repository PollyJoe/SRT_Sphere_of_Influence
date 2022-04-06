#include <iostream>
#include <cmath>

#include "Constant.h"
#include "Ephemeris.h"
#include "Lambert.hpp"
#include "PSO_DE_NLOPT.hpp"
#include "OrbitFun.h"
using namespace std;
//注意：兰伯特求解器的输入为两次位置矢量，输出为轨道的半长轴、离心率和两个速度


//规定变轨次数
const int TransNum = 3;

//设置星历
string dir_Earth("/Users/polly/Desktop/SRT/Ephemeris/Earth_1.txt");
string dir_Mars("/Users/polly/Desktop/SRT/Ephemeris/Mars_1.txt");

//题目规定数据
const int five_years = 5 * 365 * 24 * 3600;
const double M_Mars = 6.4171E23;
ephemeris Earth(dir_Earth, 24 * 3600);                                                              //读取地球星历得到地球的位置和速度
ephemeris Mars(dir_Mars, 24 * 3600);                                                                //读取火星星历得到火星的位置和速度
const double Au = 149597870700;                                                                     //日地距离（m）
const double Esph = 9.25E8;                                                                         //地球影响球半径（m）
const double Msph = pow(M_Mars / 1.989E30, 2 / 5) * Au;                                             //火星影响球半径（m）
const double emu = 3.98600441800e+14;                                                               //地球引力常数GM

//定义所需变量
double** r;                                                                                         //通过二级索引定义位矢
double*** v;                                                                                        //通过三级索引定义每次变轨前后的速度矢量
double total_t;                                                                                     //定义总时间
double* delta_t;                                                                                    //定义每次变轨所用的时间
double* a;                                                                                          //使用一级索引定义每次变轨行进轨道的半长轴
double* e;                                                                                          //使用一级索引定义每次变轨行进轨道的离心率

          
/****************************************************常用工具函数***************************************************************/
/************************************包括求速度、初始化变量、以及用来进行优化的函数*********************************************/

//设置一个初始化所有变量的函数
void init()
{
    int i, j;
    r = new double* [TransNum + 1];
    v = new double** [TransNum + 1];
    delta_t = new double[TransNum];
    a = new double[TransNum];
    e = new double[TransNum];
    for (i = 0; i < TransNum + 1; i++)
    {
        r[i] = new double[3];
        v[i] = new double* [2];
        for (j = 0; j < 2; j++)
        {
            v[i][j] = new double[3];
        }
    }
}
//设置一个求速度变量的函数
double delta_velocity(double*** v)
{
    double result = 0;
    int turn, i;
    for (turn = 0; turn < TransNum + 1; turn++)
    {
        double temp = 0;
        for (i = 0; i < 3; i++)
        {
            temp += pow(v[turn][0][i] - v[turn][1][i], 2);
        }
        result += sqrt(temp);
    }
    return result;
}

/*********************************************************求解模块******************************************************************/
/************************************************主要用于还原实际变轨情景***********************************************************/


double Multi_Transfer(const std::vector<double>& X, std::vector<double>& grad, void* f_data)  //对多变轨模块的优化
{
    
    double* para = (double*)f_data; //类型强制转换
   
    int j;
    double sum = 0;
    for (j = 0; j < TransNum + 1; j++)                                                  //构建一个和，便于使用边界约束
    {
        sum += X[j];
    }
    if (sum > 1)                                                                        //利用惩罚函数限制在五年之内
    {
        return 1E20 * pow(sum, 2);
    }
    init();
    int i;

    Earth.get_position(X[0] * five_years);                                               //得到出发时地球的位置和速度，前者用于固定一个出发点，后者用于计算出发时的速度变量
    Earth.get_velocity(X[0] * five_years);
    //求出变轨总时间
    total_t = 0;
    for (i = 0; i < TransNum + 1; i++)                                                   //优化变量：时间；X[0]*五年是出发时的时刻（因为不是t=0时立刻出发)
    {
        total_t += X[i] * five_years;
    }

    Mars.get_position(total_t);                                                          //得到到达时火星的位置和速度，前者用于固定一个到达点，后者用于计算到达时的速度变量
    Mars.get_velocity(total_t);

    r[0][0] = Earth.r[0];                                                                //初始化出发点
    r[0][1] = Earth.r[1];
    r[0][2] = Earth.r[2];
    v[0][0][0] = Earth.v[0];                                                             //初始化出发加速前速度
    v[0][0][1] = Earth.v[1];
    v[0][0][2] = Earth.v[2];

    r[TransNum][0] = Mars.r[0];                                                          //得到终点和终点的末速度
    r[TransNum][1] = Mars.r[1];
    r[TransNum][2] = Mars.r[2];
    v[TransNum][1][0] = Mars.v[0];                                                       //得到终点的末速度
    v[TransNum][1][1] = Mars.v[1];
    v[TransNum][1][2] = Mars.v[2];
 
    if (TransNum > 1)
    {
        for (i = 1; i < TransNum; i++)                                                   //优化变量：变轨位置
        {
            r[i][0] = cos(pi * X[TransNum - 2 + 3 * i]) * 1.6 * Au;
            r[i][1] = cos(pi * X[TransNum - 2 + 3 * i + 1]) * 1.6 * Au;
            r[i][2] = cos(pi * X[TransNum - 2 + 3 * i + 2]) * 1.6 * Au;
        }
    }

    for (i = 0; i < TransNum; i++)                                                       //接入lambert求解器,遍历求解
    {
        delta_t[i] = X[i + 1] * five_years;
        double* v1 = v[i][1];
        double* v2 = v[i + 1][0];
        const double* R1 = r[i];
        const double* R2 = r[i + 1];
        int flag = 1;
        int way = 0; int N = 0; int branch = 0; int Maxiter = 60; double tol = 1.0e-11;
        const double unith[3] = { 0,0,1 };
        lambert(&v1[0], &v2[0], a[i], e[i], &R1[0], &R2[0], delta_t[i], &unith[0], flag, mu, way, N, branch, Maxiter, tol);
    }

    double obj = delta_velocity(v);

    return obj;
}



/**************************************************************优化模块******************************************************************/
/***********************该优化器可以随时调用优化各个函数，容器内种子数量现阶段足够，不必反复使用测试算例即可优化*************************/
void optimization(double (*ObjFun)(const std::vector<double>& X, std::vector<double>& grad, void* f_data),int vnum,int psonum)
{
    std::vector<double> X = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
    double f = 1000.0;
    PSO(ObjFun, nullptr, X, f, vnum, 1000, psonum, 1000);
}





/************************************************************影响球模块****************************************************/
/*
输入：
输出：
*/
double Earth_Sphere(const std::vector<double>& X, std::vector<double>& grad, void* f_data)//地心坐标系
{
    int e_trans = 3;//球内总变轨次数
    double **re;//每次变轨的坐标
    double ***ve;//每次变轨前后的速度，第二角标0表示变轨前，1表示变轨后
                                                                                            
    double total_te=0; //变轨总时间
    double *delta_te;//每次变轨的时间
    double *ae;  //每次轨道的半长轴和偏心率
    double *ee;
    //变量初始化
    int i, j;
    re = new double* [e_trans + 1];
    ve = new double** [e_trans + 1];
    delta_te = new double[e_trans];
    ae = new double[e_trans];
    ee = new double[e_trans];
    for (i = 0; i < e_trans + 1; i++)
    {
        re[i] = new double[3];
        ve[i] = new double* [2];
        for (j = 0; j < 2; j++)
        {
            ve[i][j] = new double[3];
        }
    }
    //设置停泊轨道半径
    double R0 = 1E8;
    //设置停泊轨道初位置
    re[0][0] = R0 * cos(X[0]) * cos(X[1]);
    re[0][1] = R0 * cos(X[0]) * sin(X[1]);
    re[0][2] = R0 * sin(X[0]);
   //设置变轨时间
    for (int k = 0; k < e_trans + 1; k++)
    {
        delta_te[k] = X[k + 2] * 3600;
    }
    //求解总时间
    for(int i=0;i<e_trans;i++)
    {
        total_te+=delta_te[i];
    }
    
    //求解此时地球所处的位置
    const double rv0[6]={r[0][0],r[0][1],r[0][2],v[0][0][0],v[0][0][1],v[0][0][2]};
    double rv1[6];//实际用来求交点坐标的rv参数
    int flag;
    rv02rvf(flag,rv1,rv0,total_te,mu);
    
    //求解轨迹与影响球的交点并且从日心系转换为地心系
    
    
    
    double obj;
    return obj;
}















int main()
{
    optimization(Multi_Transfer,4*TransNum-2,10000);
    return 0;
}
