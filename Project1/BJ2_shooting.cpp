#include <stdio.h>
#include <math.h>
#define N 500
#define NARRAY 2*N+1
#define F(x,e) (-3*(pow(1/cosh((x)),2))-(e))

int main() {
    
    double eigenf[NARRAY], xc[NARRAY];
    double x, xstart, xend, dx;
    double e;
    double u, umax, p, unew, pnew;
    double dudx1, dudx2, dudx3, dudx4;
    double dpdx1, dpdx2, dpdx3, dpdx4;
    double uproduct, uprevious;
    double pproduct, pprevious;
    int i, j, igraph;

    xstart = -50.0;
    xend = 0;
    dx = (xend - xstart) / N;


    for (j = 0; j < 10000; j++) {
        e = -1.7 + 0.000001 * j;
        x = xstart;
        u = 1.0e-5;
        xc[0] = x;
        eigenf[0] = u;

        umax = 0;

        p = u * sqrt(F(xstart, e));
        for (i = 1; i <= N; i++) {
            dudx1 = p;
            dpdx1 = F(x, e) * u;
            dudx2 = p + dpdx1 * dx / 2;
            dpdx2 = F(x + dx / 2, e) * (u + dudx1 * dx / 2);
            dudx3 = p + dpdx2 * dx / 2;
            dpdx3 = F(x + dx / 2, e) * (u + dudx2 * dx / 2);
            dudx4 = p + dpdx3 * dx;
            dpdx4 = F(x + dx, e) * (u + dudx3 * dx);

            unew = u + (dudx1 + 2 * dudx2 + 2 * dudx3 + dudx4) * dx / 6;
            pnew = p + (dpdx1 + 2 * dpdx2 + 2 * dpdx3 + dpdx4) * dx / 6;

            pprevious = p;
            u = unew;
            p = pnew;
            x = x + dx;
            xc[i] = x;
            eigenf[i] = u;

            if (umax < fabs(u)) {
                umax = fabs(u);
            }
        }
        if (j > 0) {
            pproduct = pprevious * p;
            if (pproduct <= 0.0) {
                break;
            }
        }
    }
    printf("energy = %12.5e\n", e);
    

    return 0;
}