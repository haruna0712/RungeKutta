

#include <stdio.h>
#include <math.h>
#include<iostream>
#include <fstream>

using namespace std;
// きざみ幅
double h = 0.0005;

double k = 0.0001;
//Gaussian width
double R0 = 0.025 * sqrt(2);

extern double V0;
double V0;

double  kcotdel(double rho, double drho) {
    return k * (k + drho * cos(k * rho) / sin(k * rho)) / (k * cos(k * rho) / sin(k * rho) - drho);
}
// kは入射波のエネルギーの波数。kは0にもっていく。散乱長を求めるため。

double dphi(double x, double phi, double z)  // dphi/dx = z
{
    return(z);
}
double V(double x)
{
    //return -1.0/(r0*sqrt(6.28))*exp(-x*x/(2*r0*r0));
    return V0 * exp(-x * x / (R0 * R0));
    //return 0.5 * x * x;
}
double dz(double x, double phi, double z)    // dz/dx = (V(x)-E)phi
{
    return((V(x) - k * k) * phi);
    //return((V(x) - 0.5) * phi);
}




/* ---------------------------------------------------------------------
   ルンゲ・クッタ法
   --------------------------------------------------------------------- */

double rungeKuttaMethod(double x0, double phi0, double z0, double h) {
    int i, j, n, loop, color;
    double x, phi, z, p[4], q[4], phi_old[2];
    loop = 200 / h;
    x = x0; phi = phi0, z = z0;  // 波動関数初期化


    //z is a one derivative of phi
    //p is discretized phi
    //q is descretized z
    for (i = 1; i <= loop; i++) {
        // 4次のRunge-kutta法 の連立
        p[0] = h * dphi(x, phi, z);
        q[0] = h * dz(x, phi, z);
        p[1] = h * dphi(x + h / 2, phi + p[0] / 2, z + q[0] / 2);
        q[1] = h * dz(x + h / 2, phi + p[0] / 2, z + q[0] / 2);
        p[2] = h * dphi(x + h / 2, phi + p[1] / 2, z + q[1] / 2);
        q[2] = h * dz(x + h / 2, phi + p[1] / 2, z + q[1] / 2);
        p[3] = h * dphi(x + h, phi + p[2], z + q[2]);
        q[3] = h * dz(x + h, phi + p[2], z + q[2]);
        phi += p[0] / 6 + p[1] / 3 + p[2] / 3 + p[3] / 6;
        z += q[0] / 6 + q[1] / 3 + q[2] / 3 + q[3] / 6;
        x = x0 + i * h;
        if (i % 1000 == 0) {
            //printf("i = %d\n", i);
            //printf("%2.15lf\n", phi);
            //printf("rho = %2.15lf\n", z/phi);
            //printf("kcotdel = %2.15lf\n", kcotdel(i*h,z/phi));
            //printf("%2.15lf\n", atan(z/k/phi)+k*i*h);
        };
    }

    return kcotdel(loop * h, z / phi);
    x = x0; phi = phi0, z = z0;

    /*
    for (i = 0; i < 50000; i++) {
        printf("phi = %2.15lf\n", phiarray[i]);
    }*/

}

int main(void)
{
    std::ofstream writing_file;
    std::string filename = "scatlendata0015.txt";
    writing_file.open(filename, std::ios::out);
    double phi0, z0;

    phi0 = 0.0000;
    z0 = 0.1;
    //65.5~sqrt(2.68/(0.025*0.025))=sqrt(2.68)/k0
    // //163.8~sqrt(2.68/(0.01*0.01))=sqrt(2.68)/0.01
    // 109.2~sqrt(2.68/(0.015*0.015))=sqrt(2.68)/0.015
    //DrawPotential(xmax);
    for (V0 = -20000; V0 <= 20000; V0 += 25) {
        //cout <<(k0 * R0) * (k0 * R0)/2 << " ";
        cout << V0*R0 * R0 << ",";
    }
    cout << endl;
    cout << endl;
    for (V0 = -20000; V0 <= 20000; V0 += 25) {
        //cout <<(k0 * R0) * (k0 * R0)/2 << " ";
        //cout << (k0 * R0) * (k0 * R0) / 2 << ",";
        double r = rungeKuttaMethod(0.0, phi0, z0, h);
        cout << 1.0/r<<",";
        writing_file << V0* R0 *R0 << "," << 1.0/r << std::endl;
    }
    //std::string writing_text = "

    writing_file.close();
    //cpgclos();

    return 0;
}