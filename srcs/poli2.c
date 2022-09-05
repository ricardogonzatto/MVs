#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <mps/mps.h>
#include <complex.h>
#define o 1326
#define als_numero 200  // Numero de als gerados
#define l 50
#define lmax 50


// Função para calcular os coeficientes
void coefi_pol(double al_real[], double al_imag[], mpf_t coef_real[], mpf_t coef_imag[])
{
    mpz_t int_binom;
    mpf_t float_binom;
    mpf_t raiz_binom;
    mpf_t mp_al_real[o];
	mpf_t mp_al_imag[o];
    mpf_t coef_negativ;
    int index;                         // levar em conta que esse valor do index pode explodir para listas grandes

    for(int i = 0; i < o; i++)
    { //Convertendo os als double em float do gmp
        mpf_inits(mp_al_real[i], mp_al_imag[i], NULL);

        mpf_set_d(mp_al_real[i], al_real[i]);
        mpf_set_d(mp_al_imag[i], al_imag[i]);
    }

    mpz_inits(int_binom, NULL);
    mpf_inits(float_binom, raiz_binom, coef_negativ, NULL);

    for(int i = 0; i < l; i++)
    {
        mpf_inits(coef_real[l], coef_imag[l], NULL);

        mpz_bin_uiui(int_binom, (2*l), l);     // Calculando Binomio
        mpf_set_z(float_binom, int_binom);     // Convertendo inteiro em float para tirar a raiz
        mpf_sqrt(raiz_binom, float_binom);     // Tirando a raiz

        mpf_mul(coef_real[l], raiz_binom, mp_al_real[l]);  // Multiplicando o al pela raiz do binom
        mpf_mul(coef_imag[l], raiz_binom, mp_al_imag[l]);

        for(int m = 1; m <= l; m++)
        {
            index = ((m * ((2*lmax) + 1 - m))/ 2) + l;       //calculando a posição dos alms

            mpf_inits(coef_real[l-m], coef_imag[l-m], coef_real[l+m], coef_imag[l+m], NULL);

            mpz_bin_uiui(int_binom, (2*l), (m+l));        // Calculando o binômio
            mpf_set_z(float_binom, int_binom);            // Convertendo inteiro em float para tirar a raiz
            mpf_sqrt(raiz_binom, float_binom);            // Tirando a raiz

            mpf_mul(coef_real[l-m], raiz_binom, mp_al_real[index]);
            mpf_set_si(coef_negativ, pow(-1,m));
            mpf_mul(coef_real[l-m], coef_real[l-m], coef_negativ);

            mpf_mul(coef_imag[l-m], raiz_binom, mp_al_imag[index]);
            mpf_set_si(coef_negativ, pow(-1,m) * (-1));
            mpf_mul(coef_imag[l-m], coef_imag[l-m], coef_negativ);

            mpf_mul(coef_real[l+m], raiz_binom, mp_al_real[index]);
            mpf_mul(coef_imag[l+m], raiz_binom, mp_al_imag[index]);

        }
    }

    mpz_clears(int_binom, NULL);
	mpf_clears(float_binom, raiz_binom, NULL);
}

// Função para extrair as raízes
void raizes_pol(mpf_t coef_real[], mpf_t coef_imag[], double* raiz_real, double* raiz_imag)
{
    mpq_t rat_coef_real[(2*l)+1], rat_coef_imag[(2*l)+1];

    for(int i = 0; i < (2*l)+1; i++)   
    { //convertendo float em rational (pois até então o algoritmo só resolve com racional)
        mpq_inits(rat_coef_real[i], rat_coef_imag[i], NULL);
        mpq_set_f(rat_coef_real[i], coef_real[i]);   
        mpq_set_f(rat_coef_imag[i], coef_imag[i]);
    }

    mps_monomial_poly *p;
    mps_context *s;
    s = mps_context_new ();
    p = mps_monomial_poly_new (s, 2*l);
    mps_context_select_algorithm(s, MPS_ALGORITHM_SECULAR_GA);

    for (int i = 0; i < ((2*l)+1); i++){
        mps_monomial_poly_set_coefficient_q(s, p, i, rat_coef_real[i], rat_coef_imag[i]);  // o primeiro número é a ordem do coeficiente
                                                                                            // o segundo o coeficiente real                                                                        // o terceiro o coeficiente imaginário
    }

    mps_context_set_input_poly (s, MPS_POLYNOMIAL (p));
    cplx_t *results = cplx_valloc (2*l);
    mps_mpsolve (s);
    mps_context_get_roots_d (s, &results, NULL);

    for(int i = 0; i < (2*l); i++){
        raiz_real[i] = cplx_Re(results[i]);
        raiz_imag[i] = cplx_Im(results[i]);
    }
}

// Função para encontrar as coordenadas dos vetores multipolares 
void coord_pol(double raiz_real[], double raiz_imag[], double* theta, double* phi)
{
//    double theta[2*l], phi[2*l], R[2*l];
    double R[2*l];
    static double complex z[2*l];

    for(int i = 0; i < (2*l); i++)
    {
        z[i] = raiz_real[i] + (raiz_imag[i] * I);
        R[i] = cabs(z[i]);
        phi[i] = M_PI + carg(z[i]);
        theta[i] = 2 * atan(1/R[i]);
    }
}

void eta_phi_pol(double theta[], double phi[], double* eta, double* varphi)
{
    for(int i = 0; i < (2*l); i++)
    {
        eta[i] = (1-cos(theta[i]))/2;
        varphi[i] = phi[i]/(2*M_PI);
    }
}

// Função para encontrar as coordenadas dos vetores de Fréchet
void frechet_pol(double theta[], double phi[], double theta_esfera[], double phi_esfera[], double* frechet_vec_eta, double* frechet_vec_varphi)
{
    double frechet_vec_theta, frechet_vec_phi;
    double vec_test[2];

    for(int i = 0; i < 768; i++ ) //aqui coloquei até 768 pq usei o nside = 8 que divide a esfera em 768 pixeis
    {
        vec_test[0] = 0;

        for(int j = 0; j < (2*l); j++)
        {
            vec_test[0] += pow(acos(cos(theta[j]) * cos(theta_esfera[i])) 
            + ((sin(theta[j]) * sin(theta_esfera[i])) * cos(phi[j] - phi_esfera[i])), 2);
        }

        if (i == 0)
        {
           frechet_vec_theta = theta_esfera[i];
           frechet_vec_phi = phi_esfera[i];
           
           vec_test[1] = vec_test[0]; 
        }
        if (i != 0 && vec_test[1] > vec_test[0])
        {
           frechet_vec_theta = theta_esfera[i];
           frechet_vec_phi = phi_esfera[i];
          
           vec_test[1] = vec_test[0];
        }   
    } 

    *frechet_vec_eta = (1-cos(frechet_vec_theta))/2;
    *frechet_vec_varphi = frechet_vec_phi/(2*M_PI);
}

int main ()
{
    FILE *file;
    double al_real[o], al_imag[o];
    mpf_t coef_real[(2*l) + 1], coef_imag[(2*l) + 1];
    double raiz_real[2*l], raiz_imag[2*l];
    double theta[2*l], phi[2*l], eta[2*l], varphi[2*l];
    char filename[400];
    double theta_esfera[768], phi_esfera[768];
    double frechet_vec_eta, frechet_vec_varphi;

    // importanto as coordenadas para encontrarmos os vetores de Fréchet
    file = fopen("/home/ricardo/Documentos/EstudosMestrado/Programas/EstudoMest7/datas/ThetaPhi.dat", "r");
    for(int i = 0; i < 768; i++)
    {
	    fscanf(file, "%lf %lf\n", &theta_esfera[i], &phi_esfera[i]); 
    }
    fclose(file);

    for(int i = 0; i < als_numero; i++)
    {
        sprintf(filename,"/home/ricardo/Documentos/EstudosMestrado/Programas/EstudoMest7/datas/ALS/AL%d.dat", i);
		file = fopen(filename, "r");

	    for(int j = 0; j < o; j++)
        {
		    fscanf(file, "%lf %lf\n", &al_real[j], &al_imag[j]);
	    }
	    fclose(file); 

        coefi_pol(al_real, al_imag, coef_real, coef_imag);
        raizes_pol(coef_real, coef_imag, raiz_real, raiz_imag);
        coord_pol(raiz_real, raiz_imag, theta, phi);
        eta_phi_pol(theta, phi, eta, varphi);
        frechet_pol(theta, phi, theta_esfera, phi_esfera, &frechet_vec_eta, &frechet_vec_varphi);

        if (i == 0)
        {
            file = fopen("/home/ricardo/Documentos/EstudosMestrado/Programas/EstudoMest7/datas/frechet.dat", "w");
            fprintf(file, "%f %f\n", frechet_vec_eta, frechet_vec_varphi);
            fclose(file);

            file = fopen("/home/ricardo/Documentos/EstudosMestrado/Programas/EstudoMest7/datas/EtaePhi.dat", "w");
            for(int j = 0; j < (2*l); j++) 
            {
                fprintf(file, "%lf %lf\n", eta[j], varphi[j]);
            }
            fclose(file);
        }

        else{
            file = fopen("/home/ricardo/Documentos/EstudosMestrado/Programas/EstudoMest7/datas/frechet.dat", "a");
            fprintf(file, "%f %f\n", frechet_vec_eta, frechet_vec_varphi);
            fclose(file);

            file = fopen("/home/ricardo/Documentos/EstudosMestrado/Programas/EstudoMest7/datas/EtaePhi.dat", "a");
            for(int j = 0; j < (2*l); j++)
            {
                fprintf(file, "%lf %lf\n", eta[j], varphi[j]);
            }
            fclose(file);
        }
    }
}



//    for(int i = 0; i < (2*l); i++)
//    {
//        printf("%f %f\n", eta[i], varphi[i]);    
//        gmp_printf("%Ff %Ff\n", coef_real[i], coef_imag[i]);
//    }