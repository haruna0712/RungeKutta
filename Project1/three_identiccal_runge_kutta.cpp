/* =====================================================================
   �����Q�E�N�b�^�@�ňꎟ�����a�U���q�|�e���V�����̌ŗL�l�����߂�
   Last modified: Oct 20 2003
   ===================================================================== */
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>

   // �����ݒ�
double E = 0.0;
double dE = -1.0;


   /* ---------------------------------------------------------------------
      main
      --------------------------------------------------------------------- */

int main(void)
{
    
    const int N = 20000;

    using namespace std;
    vector<double> potential;
    std::fstream data_file_read("C:/Users/harun/experiment/discrete.csv", ios_base::in | ios_base::binary);
    data_file_read.seekg(0);
    for (int i = 0; i < N; ++i) {
        double val;
        data_file_read.read(reinterpret_cast<char*>(&val), sizeof(val));
        potential.push_back(val);
    }
    for (int i = 0; i < potential.size(); i+=100) {
        cout << potential[i] << endl;
    }
    double phi0, z0;

    data_file_read.close();

    phi0 = 0.0; z0 = 0.0000001; double x0 = 0.0001;
    
    int i, n;
    double x, phi, z, k[4], q[4], phi_old[2];

    double step_size_RK = 0.001;
    //runge-kutta method
    for (n = 1; n <= 14; n++) {
        dE = -1 * dE / 10;
        x = x0; phi = phi0, z = z0;  // �g���֐�������

        while (1) {
            for (i = 0; i < N/2; i++) {
                // 4����Runge-kutta�@ �̘A��
                k[0] = step_size_RK * z;
                q[0] = step_size_RK * ((potential[i*2] - E) * phi) ;
                k[1] = step_size_RK * (z + q[0] / 2);
                q[1] = step_size_RK * ((potential[i * 2+1] - E) * (phi + k[0] / 2)) ;
                k[2] = step_size_RK * (z + q[1] / 2);
                q[2] = step_size_RK * ((potential[i * 2+1] - E) * (phi + k[1] / 2)) ;
                k[3] = step_size_RK * (z + q[2]);
                q[3] = step_size_RK * ((potential[i * 2+2] - E) * (phi + k[2] / 2)) ;


                phi += k[0] / 6 + k[1] / 3 + k[2] / 3 + k[3] / 6;
                z += q[0] / 6 + q[1] / 3 + q[2] / 3 + q[3] / 6;
                x = x0 + i * step_size_RK;

            }

            // ����ƑO���phi��ۑ�
            phi_old[0] = phi_old[1]; 
            phi_old[1] = phi;

            // ����ƑO���phi���ٕ����Ȃ�while(1)�𔲂���
            // i=1�ł�phi_old�̒l�����Ӗ��Ȃ̂Ŕ�r���Ȃ�
            if ((i >= 2) && ((phi_old[0] > 0) && (phi < 0))) break;
            if ((i >= 2) && ((phi_old[0] < 0) && (phi > 0))) break;

            cout << E <<endl;
            E += dE;  // �ŗL�l�ύX
            x = x0; phi = phi0, z = z0;
        }
    }
    
    // ���܂����ŗL�l��\��
    printf("E = %1.10e\n", E);
    for (int i = 0; i < 2 * N; ++i) {
        //cout << to[i] << '\t' << to[i];
    }

    



    return 0;
}


