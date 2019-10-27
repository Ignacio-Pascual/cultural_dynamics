#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "funciones_4_con_switch.h"


int rk4sys(int longs, double h, double *X, double *eX, int nsteps_d, 
double *Fit, int numtr, double Lambda, double Mu, double Beta, int K, 
double **a, int longa, int Cuantos, double **b, 
double diferencia, int **dis, int tlim, int col_salida_1, int col_salida_2, 
int col_llegada_1, int col_llegada_2, int condini, int comienzo, int final, 
char *encadenados_pob, char *encadenados_fit, char *encadenados_teq, int iniant, double alpha){

	int alo=2;//numero de pasos hasta equilibrio
	
	int iw;
	int k;
	double sumaS;
	
	FILE *archivopob, *archivofit;
	FILE *archivoteq;
	
	
	archivoteq = fopen(encadenados_teq, "w");
	archivopob = fopen(encadenados_pob, "w");
	archivofit = fopen(encadenados_fit, "w");	
	
	for (k = 0; k < nsteps_d; k++){
		//Runge-Kutta:
		alo=estabiliza(longs, h, X, eX, Fit, numtr, Lambda, Mu, Beta, K, 
		diferencia,dis,tlim, 
		col_salida_1, col_salida_2, col_llegada_1, col_llegada_2, 
		k, condini, comienzo, final, iniant, alpha);		
		fprintf(archivoteq, "%d\n", alo);
		fflush(archivoteq);
		alo=2;
		
		
		sumaS=0.0;
		for(iw=0;iw<longs;iw++){
			sumaS += eX[iw];
		}
		printf("sum(eX[i]) = %f\n", sumaS);
		
		for(iw=0;iw<longs;iw++){
			fprintf(archivopob, "%.8lf\t", eX[iw]/sumaS);
			fflush(archivopob);
			fprintf(archivofit, "%.8lf\t", Fit[iw]);
			fflush(archivofit);
		}
		
		fprintf(archivopob, "\n");
		fflush(archivopob);
		fprintf(archivofit, "\n");
		fflush(archivofit);
		
		//actualiza matriz de fitness:
		matriz_switch_suma(a,b,longa,numtr);
		
		
		
		//actualiza fitness sin cambiar agujeros:
		Cuantos=fitness_aguj_anteriores(Fit, a, K, longs, numtr);
		
		if(iniant%2==0){
			llenar_vector_igual(X, Fit, Cuantos, longs);
			expon_vector(eX, X, Fit, longs);
		}
		
		
		
	}//termina for(..;k<nsteps_d;..)
	
	printf("\n termina for con %d pasos\n",k);
	

	//Runge-Kutta:
	alo=estabiliza(longs, h, X, eX, Fit, numtr, Lambda, Mu, Beta, K, 
	diferencia,dis, tlim, 
	col_salida_1, col_salida_2, col_llegada_1, col_llegada_2, 
	k, condini, comienzo, final, iniant, alpha);
	fprintf(archivoteq, "%d\n", alo);
	fflush(archivoteq);


	sumaS=0.0;
	for(iw=0;iw<longs;iw++){
		sumaS += eX[iw];
	}
	for(iw=0;iw<longs;iw++){
		fprintf(archivopob, "%.8lf\t", eX[iw]/sumaS);
		fflush(archivopob);
		fprintf(archivofit, "%.8lf\t", Fit[iw]);
		fflush(archivofit);
	}
	fprintf(archivopob, "\n");
	fflush(archivopob);
	fprintf(archivofit, "\n");
	fflush(archivofit);	
	
	fclose(archivopob);
	fclose(archivofit);
	fclose(archivoteq);
	
	
	return(1);
	
}


