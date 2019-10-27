#include <stdlib.h>
#include <time.h>  
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "funciones_4_con_switch.h"

/* Parametros */

static int N = 6;
static int K = 1;
static double beta=1.0; 
static int agujeros_deseados=24;
static double lambda = 0.999; //lambda = lambda/N;

static int ini0ant1=1;	//0: inipo con lambda+mu=1. 
			//1: antpo con lambda+mu=1
			//2: inipo con lambda+mu=1 y distancia hamming
			//3: antpo con lambda+mu=1 y distancia hamming
			//5: antpo con lambda+mu=1 y distancia hamming

//para reordenar columnas:
static int Col_salida_1=2;
static int Col_salida_2=3;
static int Col_llegada_1=0;
static int Col_llegada_2=2;

static double alpha=1.0; //exponente distancia hamming

//auxiliares:
static int condini=0;
static int comienzo=0;
static int final=0;

//para comprobar convergencia:
static double h=0.1;
static int dif_base=1;
static int dif_exponente=-4;
static int tlim=100;

static double eps=0.01; //variacion maxima de la matriz random por paso


int main()
{
	double mu = 0.01;
	
	double f_l=0.0;
	
	double diferencia=dif_base*pow(10,dif_exponente);
	
	int hola=2;
	
	int Cuantos=0;
	
	double nsteps_dos=0.0;
	
	int **dist;
	int *aa,*bb;
	
	double *eS;
	double *S, *F;
	double **A,**B;
	
	int i,j; 
	int iw,jw;
	int longS, longA;
		
	int cadlarga=160;
	char encadenados_pob[200], encadenados_fit[200], cadena[170], cadena2[25];
	char encadenados_teq[200];
	
	
	lambda = lambda/N;
	
	
	longS = (int)pow(2,N);
	longA = (int)pow(2,K+1);
	
	if (K == 0) {
		printf("ERROR: K=0 no se usa con este archivo \n");
		return(0);
	}
	
	
	if(mu != 0.01){
		printf("ERROR: mu se define como 1 - lambda \n"); return(0);
	}
	if(lambda > 1.0/N){
		printf("ERROR: lambda entre 0 y 1 (es %lf)\n", lambda); return(0);
	}
	if(lambda < 0.0){
		printf("ERROR: lambda entre 0 y 1 (es %lf)\n", lambda); return(0);
	}
	if(ini0ant1 != 0 && ini0ant1 != 1 && ini0ant1 != 2 && ini0ant1 != 3 && ini0ant1 != 5){
		printf("ERROR: escribe variable ini0ant1=0,1,2,3,5 para no sobreescribir \n");return(0);}
	
	mu=1-lambda*N;
	printf("mu es %f y (lambda/N)*N es %f \n", mu, lambda*N); 
	
	
	if(ini0ant1==3 || ini0ant1==2){printf("alpha -el exponente de homofilia- es %.1f\n", alpha);}
	
	printf("beta es %f \n", beta); printf("N es %d y K es %d \n", N, K);
	printf("agujeros_deseados es %d \n", agujeros_deseados);
	printf("longS es 2^N %d \t", longS); printf("longA es 2^(K+1)=%d \n \n", longA);
	
	
	A = (double **) calloc((longA+1),sizeof(double *));
	B = (double **) calloc((longA+1),sizeof(double *));
	if(NULL == A || NULL == B){
        printf("ERROR: couldn't allocate memory! \n"); return(0);}
        for(i = 0; i < longA; i++){
		A[i] = (double *) calloc(N+1,sizeof(double));
		B[i] = (double *) calloc(N+1,sizeof(double));
		if(NULL == A[i] || NULL == B[i]){
			printf("ERROR: couldn't allocate memory! \n"); return(0);}
	}
	if(NULL == A || NULL == B){printf("ERROR: couldn't allocate memory! \n"); return(0);}
	
	dist = (int **) calloc((longS+1),sizeof(int *));
	if(NULL == dist){printf("ERROR: couldn't allocate memory! \n"); return(0);}
	for(i = 0; i < longS; i++){
		dist[i] = (int *) calloc(longS+1,sizeof(int));
		if(NULL == dist[i]){printf("ERROR: couldn't allocate memory! \n"); return(0);}
	}
	if(NULL == dist){printf("ERROR: couldn't allocate memory! \n"); return(0);}
	
	S = (double *) calloc((longS+1),sizeof(double));	
	F = (double *) calloc((longS+1),sizeof(double));
	eS = (double *) calloc((longS+1),sizeof(double));
	if(NULL == dist || NULL == F || NULL == S || NULL == eS){
		printf("ERROR: couldn't allocate memory! \n"); return(0);}
	
	
	srand(10);
	srand48(10);
	
	
	aa = (int *) calloc((longS+1),sizeof(int));
	bb = (int *) calloc((longS+1),sizeof(int));
	for(i=0;i<longS;i++){
		for(j=0;j<longS;j++){
			dist[i][j]= distancia_hamming(i,j,aa, bb, N);}
	}
	free(aa);
	free(bb);
	
	//matriz de epistasis inicial:
	llenar_matriz_rand( A, longA, N);

	//intercambia columnas de matriz de epistasis inicial:
	intercambiar_columnas(A, longA, N, Col_salida_1, Col_salida_2);
	
	//matriz de epistasis final:
	llenar_matriz_rand( B, longA, N);

	//intercambia columnas de matriz de epistasis final:
	intercambiar_columnas(B, longA, N, Col_llegada_1, Col_llegada_2);
    
    
    
	//CALCULA EL PARAMETRO f_l:
	f_l=fitness_parametro_fl(F, A, B, K, longS, N, agujeros_deseados);

    
	for(i=0;i<longS;i++){F[i]=0.0;}

	printf("f_l es %.16f para hacer agujeros=%d\n",f_l,agujeros_deseados);


	//CREA UN NUMERO DE AGUJEROS SOLO SI LO SON EN AMBAS MATRICES:    
	Cuantos=fitness_aguj_ambas(F, A, B, K, longS, N, f_l+0.000000005);/*agujeros*/
	printf("Hay agujeros=%d\n",Cuantos);
	
	if(Cuantos!=agujeros_deseados){
		printf("ERROR AL CALCULAR PARAMETRO f_l PARA CALCULAR agujeros_deseados\n"); return(0);}
	
	
	nsteps_dos=matriz_switch_paso( A, B, eps, longA, N);
	printf("eps dado era %f, para el cual nsteps_dos=%f \n \n", eps, nsteps_dos);
	
	//nueva matriz B se usara para actualizar la matriz en cada paso:
	matriz_cambioresta_steps(A, B, nsteps_dos, longA, N);
	
	//vector de poblacion inicial:
	llenar_vector_igual(S, F, Cuantos, longS);
	expon_vector(eS, S, F, longS);
	
	
	
	
	
	
	memset(encadenados_pob, 0, sizeof encadenados_pob);	memset(encadenados_fit, 0, sizeof encadenados_fit);
	memset(encadenados_teq, 0, sizeof encadenados_teq);
	memset(cadena, 0, sizeof cadena); memset(cadena2, 0, sizeof cadena2);
	
	
	
	
	sprintf(cadena2, "N%d", N); strncat(cadena, cadena2, cadlarga);
	
	sprintf(cadena2, "_K%d", K); strncat(cadena, cadena2, cadlarga);
	
	sprintf(cadena2, "_hole%d", Cuantos); strncat(cadena, cadena2, cadlarga);
	
	sprintf(cadena2, "_beta%.2lf", beta); strncat(cadena, cadena2, cadlarga);
	
	sprintf(cadena2, "_col%d-%d_%d-%d", Col_salida_1, Col_salida_2, Col_llegada_1, Col_llegada_2);
	strncat(cadena,cadena2,cadlarga);
	
	if(ini0ant1==0){
		sprintf(cadena2, "_inipo");
		strncat(cadena,cadena2,cadlarga);
		sprintf(cadena2, "_lam%.3lf",lambda*N);
		strncat(cadena,cadena2,cadlarga);
	}
	if(ini0ant1==1){
		sprintf(cadena2, "_antpo"); strncat(cadena,cadena2,cadlarga);
		sprintf(cadena2, "_lam%.3lf",lambda*N);
		strncat(cadena,cadena2,cadlarga);
	}
	if(ini0ant1==2){
		sprintf(cadena2, "_inipo");	strncat(cadena,cadena2,cadlarga);
		sprintf(cadena2, "_lam%.3lf",lambda*N);
		strncat(cadena,cadena2,cadlarga);
		sprintf(cadena2, "_ham%.0lf", alpha); strncat(cadena,cadena2,cadlarga);
	}
	if(ini0ant1==3){
		sprintf(cadena2, "_antpo");
		strncat(cadena,cadena2,cadlarga);
		sprintf(cadena2, "_lam%.3lf",lambda*N);
		strncat(cadena,cadena2,cadlarga);
		sprintf(cadena2, "_ham%.0lf", alpha);
		strncat(cadena,cadena2,cadlarga);
	}
	if(ini0ant1==4){printf("ERROR inipo e histeresis no tiene sentido\n"); return(0);}
	if(ini0ant1==5){
		printf("Histeresis: camino de ida \n");
		sprintf(cadena2, "_antpo");
		strncat(cadena,cadena2,cadlarga);
		sprintf(cadena2, "_lam%.3lf",lambda*N);
		strncat(cadena,cadena2,cadlarga);
	}
	
	
	sprintf(cadena2, "_eps%.3lf", eps); strncat(cadena, cadena2, cadlarga);
	
	sprintf(encadenados_pob, "res/pob_");
	sprintf(encadenados_fit, "res/fit_"); sprintf(encadenados_teq, "teq/teq_");
	
	strncat(encadenados_pob, cadena, cadlarga);	strcat(encadenados_pob, ".dat");
	
	strncat(encadenados_fit, cadena, cadlarga);	strcat(encadenados_fit, ".dat");
	
	strncat(encadenados_teq, cadena, cadlarga);	strcat(encadenados_teq, ".dat");
	
	
	
	
	
	
	
	hola=rk4sys(longS, h, S, eS, nsteps_dos, F, N, lambda, mu, beta, K, 
	A, longA, Cuantos, B, diferencia, dist, tlim, 
	Col_salida_1, Col_salida_2, Col_llegada_1, Col_llegada_2, condini, comienzo, final, 
	encadenados_pob, encadenados_fit, encadenados_teq, ini0ant1, alpha);
	printf("termina principal \n");
	
	
	
	
	
	//esta parte es solo para ver histeresis:
	if(ini0ant1==5){
		printf("empieza histeresis\n");		
		
		srand(10);
		srand48(10);

		//matriz de epistasis inicial:		
		llenar_matriz_rand( A, longA, N);

		//intercambia columnas de matriz de epistasis inicial:
		intercambiar_columnas(A, longA, N, Col_salida_1, Col_salida_2);
		
		//matriz de epistasis final:	
		llenar_matriz_rand( B, longA, N);

		//intercambia columnas de matriz de epistasis final:
		intercambiar_columnas(B, longA, N, Col_llegada_1, Col_llegada_2);
		
		
		
		//CALCULA EL PARAMETRO f_l:		
		f_l=fitness_parametro_fl(F, A, B, K, longS, N, agujeros_deseados);
		
		for(i=0;i<longS;i++){F[i]=0.0;}
		printf("f_l es %.16f para hacer agujeros=%d\n",f_l,agujeros_deseados);
		
		//CREA UN NUMERO DE AGUJEROS SOLO SI LO SON EN AMBAS MATRICES 
		//en histeresis se intercambian las matrices
		Cuantos=fitness_aguj_ambas(F, B, A, K, longS, N, f_l+0.000000005);
		printf("Hay agujeros=%d\n",Cuantos);
		
		if(Cuantos!=agujeros_deseados){
			printf("ERROR AL CALCULAR PARAMETRO f_l PARA CALCULAR agujeros_deseados\n"); return(0);}
		
		
		nsteps_dos=matriz_switch_paso( A, B, eps, longA, N);
		printf("eps dado era %f, para el cual nsteps_dos=%f \n \n", eps, nsteps_dos);
		
		//nueva matriz B se usara para actualizar la matriz en cada paso:
		//en histeresis se intercambian las matrices
		matriz_cambioresta_steps(B, A, nsteps_dos, longA, N);
		
		
		
		memset(encadenados_pob, 0, sizeof encadenados_pob);
		memset(encadenados_fit, 0, sizeof encadenados_fit);
		memset(encadenados_teq, 0, sizeof encadenados_teq);
		memset(cadena, 0, sizeof cadena);
		memset(cadena2, 0, sizeof cadena2);
		
		sprintf(cadena2, "N%d", N); strncat(cadena, cadena2, cadlarga);
		
		sprintf(cadena2, "_K%d", K); strncat(cadena, cadena2, cadlarga);
		
		sprintf(cadena2, "_hole%d", Cuantos); strncat(cadena, cadena2, cadlarga);
		
		sprintf(cadena2, "_beta%.2lf", beta); strncat(cadena, cadena2, cadlarga);
		
		sprintf(cadena2, "_col%d-%d_%d-%d", Col_salida_1, Col_salida_2, Col_llegada_1, Col_llegada_2);
		strncat(cadena,cadena2,cadlarga);
		
		
		if(ini0ant1==5){
		sprintf(cadena2, "_antpo"); strncat(cadena,cadena2,cadlarga);
		sprintf(cadena2, "_lam%.3lf",lambda*N);

		strncat(cadena,cadena2,cadlarga);
		sprintf(cadena2, "_his"); strncat(cadena,cadena2,cadlarga);
		}
		
		sprintf(cadena2, "_eps%.3lf", eps); strncat(cadena, cadena2, cadlarga);
		
		sprintf(encadenados_pob, "res/pob_"); sprintf(encadenados_fit, "res/fit_");
		sprintf(encadenados_teq, "teq/teq_");
		
		strncat(encadenados_pob, cadena, cadlarga);
		strcat(encadenados_pob, ".dat");
		
		strncat(encadenados_fit, cadena, cadlarga);
		strcat(encadenados_fit, ".dat");
		
		strncat(encadenados_teq, cadena, cadlarga);
		strcat(encadenados_teq, ".dat");
		
		
		//en histeresis se intercambian las matrices
		hola=rk4sys(longS, h, S, eS, nsteps_dos, F, N, lambda, mu, beta, K, 
		B, longA, Cuantos, A, diferencia, dist, tlim, 
		Col_salida_1, Col_salida_2, Col_llegada_1, Col_llegada_2, condini, comienzo, final, 
		encadenados_pob, encadenados_fit, encadenados_teq, ini0ant1, alpha);
		
		printf("termina histeresis \n");
		
	}//termina if(histeresis)
	
	
	
	
	
	
	
	
	
	for(i = 0; i < longA; i++){free(A[i]);}
	free(A); free(S); free(eS); free(F);
	for(i = 0; i < longA; i++){free(B[i]);}
	free(B); free(dist);
	
	printf("termina hola=%d \n", hola);
	
	return(0);
}




















