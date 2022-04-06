//
//  random_threads.hpp
//  Lambert's Problem (Self-Version)
//
//  Created by 周懿 on 2021/11/22.
//

#ifndef RANDOM_THREADS_H
#define RANDOM_THREADS_H

#include<random>

//在[min,max)范围内生成随机实数
 double realRand(const double& min, const double& max);
//在[min,max]范围内生成随机整数
int intRand(const int& min, const int& max);

#endif

