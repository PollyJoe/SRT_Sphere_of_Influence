//
//  Ephemeris.cpp
//  Lambert's problem (trans version)
//
//  Created by 周懿 on 2022/1/10.
//

#include "Ephemeris.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

ephemeris::ephemeris(string dir,double scale){
    
    r = new double[3];
    v = new double[3];
    
    int round = get_total(*this,dir);
    this->scale = scale;

    r_list = new double*[round];
    for(int i = 0;i<round;i++){
        r_list[i] = new double[3];
    }

    v_list = new double*[round];
    for(int i = 0;i<round;i++){
        v_list[i] = new double[3];
    }

    get_v_r(*this,dir);

}
    

void ephemeris::get_velocity(double time){
    int t = int(time/(scale));
    v[0] = v_list[t][0];
    v[1] = v_list[t][1];
    v[2] = v_list[t][2];
}



void ephemeris::get_position(double time){
    int t = int(time/(scale));
    r[0] = r_list[t][0];
    r[1] = r_list[t][1];
    r[2] = r_list[t][2];
}


int get_total(ephemeris a, string dir){
    fstream fin(dir); //打开文件
    string readline;
    string begin = "$$SOE";
    string end = "$$EOE";
    int flag = 0;
    int round = 0;
    while (getline(fin, readline)) //逐行读取，直到结束
    {
        if(readline==begin){
            flag = 1;
            continue;
        }
        if(readline==end){
            flag = 0;
            break;
        }
        if(flag)
            round++;
    }
    return round;
}


void get_v_r(ephemeris a,string dir) {
    fstream fin(dir); //打开文件
    FILE *fp=fopen("/Users/zhouyi/Desktop/SRT/Other files/Earth.txt","r");
    if(fp=NULL)
    {
        cout<<"Error!"<<endl;
    }
    else
    {
        cout<<"Normal"<<endl;
    }
    string readline;
    string begin = "$$SOE";
    string end = "$$EOE";
    int flag = 0;
    int round = -1;

    round = -1;
    while (getline(fin, readline)) //逐行读取，直到结束
    {
        if(readline==begin){
            flag = 1;
            continue;
        }
        if(readline==end){
            flag = 0;
            break;
        }
        if(flag){
            round++;
            istringstream input(readline);
            string substr;
            char temp_char;
            input>>substr>>substr>>substr>>substr>>a.r_list[round][0]>>temp_char>>a.r_list[round][1]>>temp_char>>a.r_list[round][2]>>temp_char>>a.v_list[round][0]>>temp_char>>a.v_list[round][1]>>temp_char>>a.v_list[round][2];

            a.r_list[round][0] = 1000* a.r_list[round][0];
            a.r_list[round][1] = 1000* a.r_list[round][1];
            a.r_list[round][2] = 1000* a.r_list[round][2];
            a.v_list[round][0] = 1000* a.v_list[round][0];
            a.v_list[round][1] = 1000* a.v_list[round][1];
            a.v_list[round][2] = 1000* a.v_list[round][2];

        }
    }
}
