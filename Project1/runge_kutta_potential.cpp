/* =====================================================================
   ルンゲ・クッタ法で一次元調和振動子ポテンシャルの固有値を求める
   Last modified: Oct 20 2003
   ===================================================================== */

#include <stdio.h>
#include <math.h>


   // 物理設定
double E = 0.0;
double dE = -1.0;
int parity = 1;

// きざみ幅
double h = 0.0001;

int linestyle = 1;
int linewidth = 1;

/* ---------------------------------------------------------------------
   potential function
   --------------------------------------------------------------------- */
/*
double V(double x)  // V(x)=1/2*x^2
{
    return (pow(x, 2) / 2);
}
*/
double V(double x)
{
    double R0 = 0.0001;
    return -1/sqrt(2*3.141592) / R0  * exp(-x * x / (R0 * R0)/2);
}


/* ---------------------------------------------------------------------
   導関数の定義
   --------------------------------------------------------------------- */

double dphi(double x, double phi, double z)  // dphi/dx = z
{
    return(z);
}
double dz(double x, double phi, double z)    // dz/dx = (V(x)-E)phi
{
    return((V(x) - E) * phi);
}



/* ---------------------------------------------------------------------
   ルンゲ・クッタ法
   --------------------------------------------------------------------- */

void rungeKuttaMethod(double x0, double phi0, double z0, double h) {
    int i, n, loop, color;
    double x, phi, z, k[4], q[4], phi_old[2];
    loop = 200 / h; color = 9;

    for (n = 1; n <= 14; n++) {
        dE = -1 * dE / 10;
        x = x0; phi = phi0, z = z0;  // 波動関数初期化
        linewidth = 2;

        while (1) {
            for (i = 1; i <= loop; i++) {
                // 4次のRunge-kutta法 の連立
                k[0] = h * dphi(x, phi, z);
                q[0] = h * dz(x, phi, z)*2;
                k[1] = h * dphi(x + h / 2, phi + k[0] / 2, z + q[0] / 2);
                q[1] = h * dz(x + h / 2, phi + k[0] / 2, z + q[0] / 2)*2;
                k[2] = h * dphi(x + h / 2, phi + k[1] / 2, z + q[1] / 2);
                q[2] = h * dz(x + h / 2, phi + k[1] / 2, z + q[1] / 2)*2;
                k[3] = h * dphi(x + h, phi + k[2], z + q[2]);
                q[3] = h * dz(x + h, phi + k[2], z + q[2])*2;


                phi += k[0] / 6 + k[1] / 3 + k[2] / 3 + k[3] / 6;
                z += q[0] / 6 + q[1] / 3 + q[2] / 3 + q[3] / 6;
                x = x0 + i * h;

            }

            // 今回と前回のphiを保存
            phi_old[0] = phi_old[1]; phi_old[1] = phi;

            // 今回と前回のphiが異符号ならwhile(1)を抜ける
            // i=1ではphi_oldの値が無意味なので比較しない
            if ((i >= 2) && ((phi_old[0] > 0) && (phi < 0))) break;
            if ((i >= 2) && ((phi_old[0] < 0) && (phi > 0))) break;

            //printf("%1.10\n", E);  // モニタ
            E += dE;  // 固有値変更
            linestyle = 2;
            linewidth = 1;
            x = x0; phi = phi0, z = z0;
        }
        linestyle = 1;
        linewidth = 2;
    }

    // 波動関数の描画
    x = x0; phi = phi0, z = z0;
    for (i = 1; i <= loop; i++) {
        k[0] = h * dphi(x, phi, z);
        q[0] = h * dz(x, phi, z);
        k[1] = h * dphi(x + h / 2, phi + k[0] / 2, z + q[0] / 2);
        q[1] = h * dz(x + h / 2, phi + k[0] / 2, z + q[0] / 2);
        k[2] = h * dphi(x + h / 2, phi + k[1] / 2, z + q[1] / 2);
        q[2] = h * dz(x + h / 2, phi + k[1] / 2, z + q[1] / 2);
        k[3] = h * dphi(x + h, phi + k[2], z + q[2]);
        q[3] = h * dz(x + h, phi + k[2], z + q[2]);
        
        phi += k[0] / 6 + k[1] / 3 + k[2] / 3 + k[3] / 6;
        z += q[0] / 6 + q[1] / 3 + q[2] / 3 + q[3] / 6;
        x = x0 + i * h;
    }

    // 求まった固有値を表示
    printf("E = %1.10e\n", E);
}



/* ---------------------------------------------------------------------
   調和振動子ポテンシャルの描画(青線)  V(x)=1/2(x)
   --------------------------------------------------------------------- */




/* ---------------------------------------------------------------------
   main
   --------------------------------------------------------------------- */

int main(void)
{
    double phi0, z0;
    double xmax, bottom, top;

    // グラフ描画設定
    xmax = 8.0;
    bottom = 2.0;
    top = 2.0;

    //cpgopen( "/xserv" );

    if (parity == -1) {  // 奇関数
        phi0 = 0.0; z0 = 1.0;
    }
    else {              // 偶関数
        phi0 = 1.0; z0 = 0.0;
    }
    rungeKuttaMethod(0., phi0, z0, h);


    return 0;
}


