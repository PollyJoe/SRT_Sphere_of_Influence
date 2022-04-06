# Sphere of Influence



## Part 1 Physics Background ：

Instead of just considering the gravity of the sun, consider the gravity of the Earth and Mars. Perturbation equation can be used to solve the range of planetary gravitational action within the acceptable range of error, which is approximately a sphere, namely the sphere of influence. The whole orbit is divided into three sections, the earth's gravitational sphere, the solar system cruise, Mars's gravitational sphere. Among them, the error of the gravitational sphere segment is smaller than that of the cruise segment, so it is considered to optimize the cruise segment first to get the rough trajectory, and then conduct trajectory correction when approaching the planet through the gravitational sphere segment.



## Part 2  Process of Solution

### Section 1  Trajectory optimization of cruise section

#### 1. Idealization Model : 
Treating the two planets as two single points, where the condition of Lambert's Problem is satisfied. Here we can use Lambert Solver. 

#### 2. Optimization Objective :
$$
\Delta V_{opt}=\sum_{i=1}^n||\Delta v_i||
$$

$||\Delta v_i||$  is the norm of velocity change vector for each orbital change. 

##### >>>	Solution of $||\Delta v_i||$ :

**Once $r_i$ and $r_{i+1}$ are fixed, we have $\textcolor{red}{Lambert\space Solver}$ to solve $v_i$ and $v_{i+1}$, then we can figure out the difference of the two velocity vectors, the norm of the vector difference is $||\Delta v_i||$.**

#### 3. Optimization Variables :

##### (1) $N_{trans}$

   := The number of  orbital change in cruise section $\textcolor{red}{(manually\space fixed, not\space in\space PSO)} $     

##### (2) { $r_i$ } 

:= Array of position vectors of each orbital change  (PSO) 

##### (3) { $\Delta t_i$ } 

:= Array of time of each orbital change (PSO)

#### 4. Constraint :

$t_{total}$ := The total time of the whole launch mission $ \textcolor{red}{( including\space both\space the\space waiting\space time\space and\space launching\space time)}$

We have the requirment :
$$
t_{total}\leq 5\space years.
$$

### Section 2  Trajectory optimization of in-Sphere flight section

#### 1. Earth section

##### (1)  Optimization Objective :

$$
\Delta V_{sph}=\sum_{i=1}^{m}||\Delta v_{sph_i}||
$$

$||\Delta v_{sph_i}||$ is the norm of velocity change for each orbital change in the sphere.

Solution is almost the same as that in section 1. 

##### (2) Optimization Variables :

###### a. $n_{Etrans}$

  := The number of orbital change in Earth-Sphere section $\textcolor{red}{(manually\space fixed, not\space in\space PSO)} $     

###### b. { $r_{e_i}$ }

 := Array of position vectors of each orbital change in Earth-Sphere ( Earth-centered coodination system ,PSO )

###### c. { $\Delta t_{e_i}$ } 

:= Array of time of each orbital change (PSO)

##### (3) Constraint :

$$
||r_{e_i}||\leq R_{earth} 
$$

##### (4) Formula :

###### a. Give the transfer time array  { $\Delta t_{sph_i}$ } :

Since in a normal launch mission, it costs only several minutes to escape from the gravity field of Earth, and given that the seed $X[i]$ is in the range of [0,1], we can give the time array :
$$
\Delta t_{e_i}=X[i]*3600 \space \space (unit :s)
$$
This is required to be done first because we need the time array to figure out the intersection of the sphere surface and cruise trajectory. 

###### b. Give the position array { $r_{e_i}$ } :

**(1) We use double-dimensional array to mark those positions :**

```c++
double re[e_trans][3];//Define an array to save the position vectors.
```



**(2) We need to choose a parking orbit in order to give initial velocity and choose position for departure :**
$$
R_{park}=1*10^8(m)
$$
$\textcolor{yellow}{p.s.\space The \space radius \space R_{park} \space is \space randomly \space chosen, \space it \space depends \space on \space the \space launch \space mission. }$

**(3) Then we can use spherical coordinate frame to give the position and velocity for departure :**

```c++
#include <cmath>
//Use spherical coordinate system
re[0][0]=R0*cos(X[0])*cos(X[1]);
re[0][1]=R0*cos(X[0])*sin(X[1]);
re[0][2]=R0*sin(X[0]);
ve[0][0][0]=;
ve[0][0][1]=;
ve[0][0][2]=;
```



**(4) Then is the most difficult part ---- solving the intersection of sphere surface and trajectory, because this is going to be formed as the position of:**

```c++
void rv02rvf(int& flag, double* rv1, const double* rv0, double dt, double mu);
//input initial r amd v (rv0), output final r and v (rv1), using two-body model
```

**We use $\textcolor{red}{rv02rvf\space function}$ to solve the intersection via tracing back from the first cruise-orbital-change position vector. The solving process is listed below : **

Step 1 : Give an array of time { $t_{e_i}$ } .

```c++
delta_te[i]=X[m]*3600;
```

$\textcolor{yellow}{p.s.\space Since \space the \space time \space to \space escape\space  from\space the \space gravity\space field \space of\space Earth \space is \space usually \space several\space minutes.\space }$

$\textcolor{yellow}{Therefore\space we \space assume \space the \space mission\space costs \space less \space than\space 1\space hour.}$

Step 2 : Solve the Earth's position when the areo leave the sphere.

```c++
double rv0[6]={r[0][0],r[0][1],r[0][2],v[0][0][0],v[0][0][1],v[0][0][2]};
//This is the starting position and velocity of Earth we got from cruise section.
double rv1[6];//The Earth's position and velocity we need
int flag;
rv02rvf(flag,rv1,rv0,delta_te[i],mu);
```

Here we solved $r_{Earth}$.

Step 3 : Solve the intersection between the trajectory and the surface of sphere.

We use interpolation to figure out a numerical solution. 

1. Make 1000 attempts, use rv02rvf to solve the result ( *rv1 );

2. Figure out the distance between the result position and the position of Earth $||l_{rv-earth}||$; 

3. Compare the distance $||l_{rv-earth}||$ and the radius of sphere of influence  $r_{Esph}$  :

   a. If $||l_{rv-earth}||> r_{Esph} $, the areo is not in the sphere, go to the next attempt;

   b. If $||l_{rv-earth}||< r_{Esph} $, the areo is in the sphere, stop attempt, take the data ( assume is rv[k] ) and the one in the last trial ( rv[k-1] ) .

4. Take the equinoctical point of rv[k] and rv[k-1] as the velocity and position vector when the areo is at the surface of the sphere. We will define it as resph[3] and vesph[3] ( the statistics are taken from the equinotical vector ) . 

   $\textcolor{red}{Caution: Pay\space attention\space to\space the\space coodinate \space system \space change}$

   $\textcolor{yellow}{p.s.\space In \space (5) \space we \space need \space to \space transform \space them \space from \space Solar \space coordinate \space system \space into \space Earth-centered \space coordinate \space system.}$

**(5) Now we can give the position and velocity when the areo is at the surface via the data we got from (4) :**



