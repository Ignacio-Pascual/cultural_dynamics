#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "funciones_4_con_switch.h"


//diferencia es variable de convergencia
int estabiliza(int longs, double h, double *X, double *eX, double *Fit, 
int numtr, double Lambda, double Mu, double Beta, int K, 
double diferencia, int **dis, int tlim, 
int col_salida_1, int col_salida_2, int col_llegada_1, int col_llegada_2, 
int k_de_rk, int condini, int comienzo, int final, int iniant, double alpha){
	
	
	double suma=0.0;
	
	int klim=0;
	
	
	double maxima_dif=1.0;//mide maximo de las diferencias de la poblacion (convergencia)
	double *eX_ant;
	
	int RK_order=4;
	int i, iw, jw;
	int k;
	double *Y, **Frk;
	double *Fbeta;
	
	
	eX_ant = (double *) calloc((longs+1),sizeof(double));
	Frk = (double **) calloc((RK_order+1),sizeof(double *));
	Y = (double *) calloc((longs+1),sizeof(double));
	Fbeta = (double *) calloc((longs+1),sizeof(double));
	
	if (NULL == eX_ant || NULL == Frk || NULL == Y || NULL == Fbeta){
		printf("ERROR: no hay memoria \n");}
	for (i = 0; i < RK_order + 1; i++){
		Frk[i] = (double *) calloc((longs+1),sizeof(double));
		if (NULL == Frk[i]) {printf("ERROR: no hay memoria Frk[%d]! \n",i);}
	}
	if(NULL == Frk){printf("ERROR: no hay memoria Frk \n");}
    
	
	//VARIABLE PARA CRITERIO DE CONVERGENCIA:
	klim=(int) ceil(tlim/h);
	
	
	for(i=0;i<longs;i++){Fbeta[i]=pow(Fit[i], Beta);}
	
	
	k=0;//inicializo contador de pasos hasta convergencia
	while(maxima_dif > diferencia){
		
		//Paso 1 Runge-Kutta:
		suma=suma_pob(eX, longs);		
		
		if(iniant==3 || iniant==2){//Desde aqui si lleva homofilia
			
			for(i=0;i<longs;i++){
				Frk[1][i] = actualizacion_dham_alpha(i, Fbeta, X, longs, numtr, Lambda, Mu, Beta, suma, dis,alpha);
			}
			for (i = 0; i <= longs; i++){
				Y[i] = X[i] + 0.5 * h * Frk[1][i];
				eX[i]=exp(Y[i]);
			}
			
			//Paso 2 Runge-Kutta:
			suma=suma_pob(eX, longs);
			for(i=0;i<longs;i++){
				Frk[2][i] = actualizacion_dham_alpha(i, Fbeta, Y, longs, numtr, Lambda, Mu, Beta, suma, dis,alpha);
			}
			for (i = 0; i <= longs; i++){
				Y[i] = X[i] + 0.5 * h * Frk[2][i];
				eX[i]=exp(Y[i]);
			}
			
			//Paso 3 Runge-Kutta:
			suma=suma_pob(eX, longs);
			for(i=0;i<longs;i++){
				Frk[3][i] = actualizacion_dham_alpha(i, Fbeta, Y, longs, numtr, Lambda, Mu, Beta, suma, dis,alpha);
			}
			for (i = 0; i <= longs; i++){
				Y[i] = X[i] + h * Frk[3][i];
				eX[i]=exp(Y[i]);
			}
			
			//Paso 4 Runge-Kutta:
			suma=suma_pob(eX, longs);
			for(i=0;i<longs;i++){
				Frk[4][i] = actualizacion_dham_alpha(i, Fbeta, Y, longs, numtr, Lambda, Mu, Beta, suma, dis,alpha);
			}
			
			//Actualiza Runge-Kutta:
			for (i = 0; i <= longs; i++){
				X[i] = X[i] + (h/6) * (Frk[1][i] + 2 * Frk[2][i] + 2 * Frk[3][i] + Frk[4][i]);
				eX[i] = exp(X[i]);
			}
			
			
		}//Hasta aqui si lleva homofilia
		
		
		
		else{//Desde aqui si no lleva homofilia
			
			//Paso 1 Runge-Kutta:
			for(i=0;i<longs;i++){
				Frk[1][i] = actualizacion(i, Fbeta, X, longs, numtr, Lambda, Mu, Beta, suma, dis);
			}
			for (i = 0; i <= longs; i++){
				Y[i] = X[i] + 0.5 * h * Frk[1][i];
				eX[i]=exp(Y[i]);
			}
			
			//Paso 2 Runge-Kutta:
			suma=suma_pob(eX, longs);
			for(i=0;i<longs;i++){
				Frk[2][i] = actualizacion(i, Fbeta, Y, longs, numtr, Lambda, Mu, Beta, suma, dis);
			}
			for (i = 0; i <= longs; i++){
				Y[i] = X[i] + 0.5 * h * Frk[2][i];
				eX[i]=exp(Y[i]);
			}
			
			//Paso 3 Runge-Kutta:
			suma=suma_pob(eX, longs);
			for(i=0;i<longs;i++){
				Frk[3][i] = actualizacion(i, Fbeta, Y, longs, numtr, Lambda, Mu, Beta, suma, dis);
			}
			for (i = 0; i <= longs; i++){
				Y[i] = X[i] + h * Frk[3][i];
				eX[i]=exp(Y[i]);
			}
			
			//Paso 4 Runge-Kutta:
			suma=suma_pob(eX, longs);
			for(i=0;i<longs;i++){
				Frk[4][i] = actualizacion(i, Fbeta, Y, longs, numtr, Lambda, Mu, Beta, suma, dis);
			}
			
			//Actualiza Runge-Kutta:
			for (i = 0; i <= longs; i++){
				X[i] = X[i] + (h/6) * (Frk[1][i] + 2 * Frk[2][i] + 2 * Frk[3][i] + Frk[4][i]);
				eX[i] = exp(X[i]);
			}
			
			
		}//Hasta aqui si no lleva homofilia
		
		
		
		if(k%klim==0){
			suma=0.0;//suma exponencial poblacion
			for(iw=0;iw<longs;iw++){suma += eX[iw];}
			
			
			maxima_dif=0.0;//reinicializa maxima_dif
			
			for(jw=0;jw<longs;jw++){
				if(maxima_dif < fabs(eX_ant[jw]-eX[jw]) ){
					maxima_dif = fabs(eX_ant[jw]-eX[jw]);}
			}//maxima_dif (norma infinito) no normalizada.
			//maxima_dif normalizada (sumaS es constante en cada paso)
			maxima_dif = maxima_dif/suma;
			
			
			//actualiza eX_ant:
			for(jw=0;jw<longs;jw++){eX_ant[jw] = eX[jw];}
			
		}//TERMINA CRITERIO DE CONVERGENCIA
		
		k++;//contador de pasos
		
		
		
	}//termina while(maxim_dif > diferencia)
	
	printf("funcion estabiliza termina for con maxima_dif=%f en %d pasos\n", maxima_dif, k);
	
	
	for (i = 0; i < RK_order +1; i++){free(Frk[i]);}
	free(Frk); free(Y); free(eX_ant); free(Fbeta);
	
	
	
	return(k);
	
}

