//
//  ephemeris.h
//  Lambert's problem (trans version)
//
//  Created by 周懿 on 2022/1/10.
//

#ifndef Ephemeris_hpp
#define Ephemeris_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

class ephemeris{
public:
    ephemeris(string dir,double scale);
    void get_velocity(double time);
    void get_position(double time);
    friend void get_v_r(ephemeris a, string dir);
    friend int get_total(ephemeris a, string dir);
public:
    double *r ;
    double *v ;
    double **r_list;
    double **v_list;
    double scale;
};

void get_v_r(ephemeris a,string dir);

int get_total(ephemeris a, string dir);

#endif
