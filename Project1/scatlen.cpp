
#include <stdio.h>
#include <math.h>
#include<iostream>
#include <fstream>
#include <cmath>

using namespace std;
// Ruge-Kuttaのきざみ幅
double h = 0.0005;

//Gaussian width
double R0 = 0.015 * sqrt(2);

extern double k0;
double k0;

double  kcotdel(double rho, double drho, double k) {
    return k * (k + drho * cos(k * rho) / sin(k * rho)) / (k * cos(k * rho) / sin(k * rho) - drho);
}
// kは入射波のエネルギーの波数。kは0にもっていく。散乱長を求めるため。

double V(double x,double v)
{
    return v/(R0*R0) * exp(-x * x / (R0 * R0));
}
double ddphi(double x, double phi, double z, double k,double v)    // dz/dx = (V(x)-E)phi
{
    return((V(x,v) - k * k) * phi);
}




/* ---------------------------------------------------------------------
   ルンゲ・クッタ法
   --------------------------------------------------------------------- */

double rungeKuttaMethod(double x0, double phi0, double z0, double h, double k,double v) {
    int i, j, loop;
    double x, phi, dphi, p[4], q[4];
    loop = 200 / h;
    x = x0; phi = phi0, dphi = z0;  // 波動関数初期化
    
    //p is discretized phi
    //q is descretized dphi
    for (i = 1; i <= loop; i++) {
        // 4次のRunge-kutta法 の連立
        p[0] = h * dphi;
        q[0] = h * ddphi(x, phi, dphi, k,v);
        p[1] = h * (dphi + q[0] / 2);
        q[1] = h * ddphi(x + h / 2, phi + p[0] / 2, dphi + q[0] / 2, k,v);
        p[2] = h * (dphi + q[1] / 2);
        q[2] = h * ddphi(x + h / 2, phi + p[1] / 2, dphi + q[1] / 2, k,v);
        p[3] = h * (dphi + q[2]);
        q[3] = h * ddphi(x + h, phi + p[2], dphi + q[2], k,v);
        phi += p[0] / 6 + p[1] / 3 + p[2] / 3 + p[3] / 6;
        dphi += q[0] / 6 + q[1] / 3 + q[2] / 3 + q[3] / 6;
        x = x0 + i * h;
        
    }
    // 1/aを返す。
    return -kcotdel(loop * h, dphi / phi, k);
    x = x0; phi = phi0, dphi = z0;
}

int main(void)
{
    double k = 0.001;
    //std::ofstream writing_file;
    //std::string filename = "scatlendata0015.txt";
    //writing_file.open(filename, std::ios::out);

    double phi0, dphi0;
    //x0＝0で波動関数phi0=0が境界条件(一般的な3次元schrodinger eqの扱いにより。)
    phi0 = 0.0000;
    dphi0 = 0.1;
    /*
    for (double v = -2.0; v >= -3.2; v-=0.0002) {
        cout << v << endl;
        double inva = rungeKuttaMethod(0.0, phi0, dphi0, h, 0.01, v);
        writing_file << v << " " <<inva<< std::endl;
    }
    */
    //for (double inv_scat = 2.2; inv_scat < 10.01; inv_scat += 0.2) {
        double inv_scat = 7.6;
        double v = -1.0;
        double vlow=-9.4;
        double vhigh=-2.0;
        double err;
        while (1) {
             err= rungeKuttaMethod(0.0, phi0, dphi0, h, 0.01, v) - inv_scat;
             //cout << err;
            if (abs(err) < 0.001) break;
            else if (err<=0) {
                vhigh = v;
                v = (v + vlow) / 2;
            }
            else if (err>0) {
                vlow = v;
                v = (v + vhigh) / 2;
            }

        }
        //cout << inv_scat << endl;
        std::cout << v << ",";
    //}

    //writing_file.close();
    

    return 0;
}