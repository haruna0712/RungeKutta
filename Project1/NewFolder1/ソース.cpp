/* =====================================================================
     	   �����Q�E�N�b�^�@�ňꎟ�����a�U���q�|�e���V�����̌ŗL�l�����߂�
     	   Last modified: Oct 20 2003
     	   ===================================================================== */

	#include <stdio.h>
	#include <math.h>
	//#include "cpgplot.h"

	// �����ݒ�
	double E = -5.0;
	double dE = 1.0;
	int parity = -1;

	// �����ݕ�
	double h = 0.001;



	/* ---------------------------------------------------------------------
	   potential function
	   --------------------------------------------------------------------- */
    double r0 = 0.001;
	double V(double x)  // V(x)=1/2*x^2
	{
    	    //return -1.0/(r0*sqrt(6.28))*exp(-x*x/(2*r0*r0));
        return -3.0 / (pow(cosh(x), 2));
    	}



	/* ---------------------------------------------------------------------
	   ���֐��̒�`
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
	   �����Q�E�N�b�^�@
	   --------------------------------------------------------------------- */

	void rungeKuttaMethod(double x0, double phi0, double z0, double h) {
    	  int i, j, n, loop, color;
    	  double x, phi, z, k[4], q[4], phi_old[2];
    	  loop = 15 / h; color = 9;
    
        	  for (n = 1; n <= 14; n++) {
        	    dE = -1 * dE / 10;
        	    x = x0; phi = phi0, z = z0;  // �g���֐�������
        
            	    while (1) {
            	      for (i = 1; i <= loop; i++) {
                		// 4����Runge-kutta�@ �̘A��
                    	k[0] = h * dphi(x, phi, z);
                		q[0] = 2*h * dz(x, phi, z);
                		k[1] = h * dphi(x + h / 2, phi + k[0] / 2, z + q[0] / 2);
                		q[1] = 2*h * dz(x + h / 2, phi + k[0] / 2, z + q[0] / 2);
                		k[2] = h * dphi(x + h / 2, phi + k[1] / 2, z + q[1] / 2);
                		q[2] = 2*h * dz(x + h / 2, phi + k[1] / 2, z + q[1] / 2);
                		k[3] = h * dphi(x + h, phi + k[2], z + q[2]);
                		q[3] = 2*h * dz(x + h, phi + k[2], z + q[2]);
                		phi += k[0] / 6 + k[1] / 3 + k[2] / 3 + k[3] / 6;
                		z += q[0] / 6 + q[1] / 3 + q[2] / 3 + q[3] / 6;
                		x = x0 + i * h;
                
            }
            
                	      // ����ƑO���phi��ۑ�
                	      phi_old[0] = phi_old[1]; phi_old[1] = phi;
            
                	      // ����ƑO���phi���ٕ����Ȃ�while(1)�𔲂���
                	      // i=1�ł�phi_old�̒l�����Ӗ��Ȃ̂Ŕ�r���Ȃ�
                	      if ((i >= 2) && ((phi_old[0] > 0) && (phi < 0))) break;
            	      if ((i >= 2) && ((phi_old[0] < 0) && (phi > 0))) break;
            
                	      printf("%2.15lf\n", E);  // ���j�^
            	      E += dE;  // �ŗL�l�ύX
            	      x = x0; phi = phi0, z = z0;
            
        }
        
    }
    
        	  // ���܂����ŗL�l��\��
        	  printf("E = %2.15lf\n", E);
    
        	  // �g���֐��̕`��
        	  x = x0; phi = phi0, z = z0;
              /*
    	  for (i = 1; i <= loop; i++) {
        	    k[0] = h * dphi(x, phi, z);
        	    q[0] = h * dz(x, phi, z);
       	        k[1] = h * dphi(x + h / 2, phi + k[0] / 2, z + q[0] / 2);
        	    q[1] = h * dz(x + h / 2, phi + k[0] / 2, z + q[0] / 2);
        	    k[2] = h * dphi(x + h / 2, phi + k[1] / 2, z + q[1] / 2);
        	    q[2] = h * dz(x + h / 2, phi + k[1] / 2, z + q[1] / 2);
        	    k[3] = h * dphi(x + h, phi + k[2], z + q[2]);
        	    q[3] = h * dz(x + h, phi + k[2], z + q[2]);
        	    //cpgsci(color); cpgmove(x, phi); cpgdraw(x + h / 6, phi + k[0] / 6);
        	    //cpgsci(color); cpgmove(x + h / 6, phi + k[0] / 6); cpgdraw(x + h / 2, phi + k[1] / 2);
        	    //cpgsci(color); cpgmove(x + h / 2, phi + k[1] / 2); cpgdraw(x + h * 5 / 6, phi + k[2] * 5 / 6);
            //cpgsci(color); cpgmove(x + h * 5 / 6, phi + k[2] * 5 / 6);  cpgdraw(x + h, phi + k[3]);
            phi += k[0] / 6 + k[1] / 3 + k[2] / 3 + k[3] / 6;
            z += q[0] / 6 + q[1] / 3 + q[2] / 3 + q[3] / 6;
           x = x0 + i * h;
        
    }*/
    
}


	/* ---------------------------------------------------------------------
	   ���a�U���q�|�e���V�����̕`��(��)  V(x)=1/2(x)
	   --------------------------------------------------------------------- 
	void DrawPotential(double xmax)
	{
    	  double i;
    	  cpgsci(4); cpgslw(1);
    	  cpgmove(-xmax, pow(xmax, 2) / 2);
    	  for (i = -xmax; i < xmax; i += 0.1) {
        129	    cpgdraw(i, pow(i, 2) / 2);
        130
    }
    	}*/



/* ---------------------------------------------------------------------
   main
   --------------------------------------------------------------------- */

	int main(void)
	{
    	  double phi0, z0;
    	  double xmax, bottom, top;
    
        	  // �O���t�`��ݒ�
          /*
        	  xmax = 12.0;
      bottom = 2.0;
    	  top = 2.0;
    */
        	  //cpgopen("/xserv");
      //cpgpap(8.0, .6);
      //cpgenv(0, xmax, -bottom, top, 0, 0);
    
        	  if (parity == -1) {  // ��֐�
        	    phi0 = 0.0; z0 = 1.0;
        
    }
    	  else {              // ��֐�
        	    phi0 = 1.0; z0 = 0.0;
        
    }
    
          //DrawPotential(xmax);
    	  rungeKuttaMethod(0.0, phi0, z0, h);
    
        	  //cpgclos();
    
          return 0;
    }