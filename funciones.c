#include <stdlib.h>
#include <time.h> 
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "funciones_4_con_switch.h"


//llenar una matriz de numeros random
void llenar_matriz_rand(double **enc, int fila, int columna) {
	
	int i, j;
	
	for(i=0;i< fila;i++){
		for(j=0;j<columna;j++){
			enc[i][j] = drand48(); //con numeros entre 0 y 1
			}
		}
	
	
} 


double matriz_switch_paso(double **enc, double **enc_dos, double eps, int fila, int columna){
	
	int i,j;
	double maxim=fabs(enc[0][0]-enc_dos[0][0]);
	
	double nsteps_switch;
	
	for(i=0;i<fila;i++){
		for(j=0;j<columna;j++){
			if(maxim<fabs(enc[i][j]-enc_dos[i][j])){maxim=fabs(enc[i][j]-enc_dos[i][j]);}
		}
	}
	printf("maxima diferencia en valor absoluto entre las dos matrices es maxim=%f \n", maxim);
	
	nsteps_switch =ceil(maxim/eps);
	
	
	return(nsteps_switch);
	
}

void matriz_cambioresta_steps(double **enc, double **enc_dos, double numsteps, int fila, int columna){
	
	int i,j;
	
	for(i=0;i<fila;i++){
		for(j=0;j<columna;j++){
			enc_dos[i][j]=((enc_dos[i][j]-enc[i][j])/numsteps);
		}
	}
	
}

void matriz_switch_suma(double **enc, double **enc_dos, int fila, int columna){
	
	int i,j;
	
	for(i=0;i<fila;i++){
		for(j=0;j<columna;j++){
			enc[i][j]=enc[i][j]+enc_dos[i][j];
		}
	}
	
	
}

//llenar un vector de numeros igualmente distribuidos
void llenar_vector_igual(double *s, double *f, int cuantos, int longvec){
	int i;
	double anchura=0.0;
	double sumaini=0.0;
	
	anchura=1/(double)(longvec-cuantos);
	
	for(i=0;i<longvec;i++){
		if(f[i]==0.0) {s[i]=-1/0.0;}
		else{s[i]=anchura;}
	}
	
	printf("llenando vector inicial s[i] = %.4lf \t", anchura);
	for(i=0;i<longvec;i++){
		sumaini+=exp(s[i]);
	}
	printf("suma(exp(s[i])) inicial = %f \n", sumaini);
	
}

	

int pasa_puntero_entero(int *p, int numtr){
	//p es el vector a convertir. numtr es el numero de traits: N
	int i;
	int g;
	g=0;
	for(i=0;i<numtr;i++){
		g = g + p[i]*(pow(2,(numtr-i-1)));
	}
	
	return(g);
}


void pasa_entero_puntero(int *pu, int j, int numtr){
	// j es el numero int asignado a un vector, numtr es el numero de 
	// traits: N para los vectores completos (pero no para los auxiliares)
	// guardar la conversion binaria en pu. OJO: definirlo fuera
	
	int aux, i, pot; //auxiliares
	
	
	pot=2;
	for(i=0;i<numtr;i++){
		aux=(j%pot);
		pu[numtr-1-i]=aux;
		j=(j-aux)/2; //actualizar j para el siguiente paso
		
	}
	
}


int distancia_hamming(int m,int l,int *p, int *q, int numtr){
	// m y l son los enteros asignados a dos vectores
	// *p y *q seran los vectores binarios de longitud numtr=N correspondientes
	
	int distancia, i;
	distancia=0;
	
	
	for(i=0;i<numtr;i++){
		p[i]=0; q[i]=0;
	}
	
	pasa_entero_puntero(p, m, numtr);
	pasa_entero_puntero(q, l, numtr);
	for(i=0;i<numtr;i++){
		if(p[i]!=q[i]) distancia++;
	}
	
	return(distancia);
	
}


double actualizacion(int l, double *fbeta, double *s, int longs, 
int numtr, double Lambda, double Mu, double Beta, double suma,int **dis){
	// l es el paso de actualizacion en el que estoy (vector que estoy actualizando)
	// su integer correspondiente lo guardo en vectorori. los demas con 
	// los que comparo los guardo en vectorcomp
	
	
	int jw;
	int i; //contador de traits
	int j; //contador de vertices
	
	int liaux=0;//liaux es el estado l (que se llama j en main), pero con el trait i cambiado

	int *vectorori, *vectorcomp, *p, *q;
	int *vectororili;
	
	
	double flux=0.0;
	double influxint=0.0, outfluxint=0.0;
	double influxmut=0.0, outfluxmut=0.0;
	double finf=0.0,fout=0.0, fdiv=0.0;
	//fdiv es la fitness por el que se divide en ambos casos (cambia en cada trait)
	
	if(fbeta[l]==0.0){return(0.0);}//(se anula influx porque f[l]=0.0 y outflux porqe x[l]=0.0 -> consistente)
	
	vectorori= (int *) calloc((numtr+1),sizeof(int));
	//vectorori(=l) es el estado que vamos a actualizar con esta funcion
	vectorcomp= (int *) calloc((numtr+1),sizeof(int));
	//vectorcomp son los vectores con los que interacciona (o no). 
	//los que suman o restan poblacion al interaccionar con vectorori
	
	vectororili= (int *) calloc((numtr+1),sizeof(int));//asociado a liaux
	
	//vectores auxiliares para calcular distancias de hamming
	p= (int *) calloc((numtr+1),sizeof(int));
	q= (int *) calloc((numtr+1),sizeof(int));
	
	pasa_entero_puntero(vectorori, l, numtr);
	
	//en el for: 
	//i corresponde a los traits de vectorori, el estado que estoy actualizando
	//j corresponde a vectorcomp: los estados con los que se interacciona
	for(i=0;i<numtr;i++){
		
		//vectororili: vectorori (asociado a l) con trait i cambiado
		for(jw=0;jw<numtr;jw++){vectororili[jw]=0;}
		for(jw=0;jw<numtr;jw++){vectororili[jw]=vectorori[jw];}
		vectororili[i]=1-vectororili[i];
		
		liaux=pasa_puntero_entero(vectororili,numtr);//liaux es el entero asociado a vectororili
		
		
		fdiv=fbeta[liaux] + fbeta[l];
		finf=fbeta[l]/fdiv;          
		fout= fbeta[liaux]/fdiv;     
		
		
		for(j=0;j<longs;j++){
			
			pasa_entero_puntero(vectorcomp, j, numtr);
			
			//Dentro de for(j...) vectorori, vectorcomp y demas siguen evaluados en el trait i
			//(j es el nodo)

			//if para LAS CONTRIBUCIONES INFLUX por INTERACCION:
			if(vectorori[i] == vectorcomp[i]){
				influxint += dis[liaux][j]*finf*(exp(s[liaux]+s[j]-s[l]));
			}
			
			
			
			//else para LAS CONTRIBUCIONES OUTFLUX por INTERACCION:
			else{
				outfluxint +=  dis[l][j]*fout*(exp(s[j]));
			}
			
			
		} //termina el for(j=0;j<longs;j++)
		
		
		//mutacion: van fuera de for(j=0;...)   (se suman solo en i)
		influxmut += finf*(exp(s[liaux]-s[l]));

		outfluxmut += fout;

		
		
		
	} //termina el for (i=0;i<numrt;i++)
	
	
	influxint=(Lambda*influxint)/suma;
	outfluxint=(Lambda*outfluxint)/suma;
	
	influxmut=influxmut*Mu;
	outfluxmut=outfluxmut*Mu;
	
	flux=influxint+influxmut-outfluxint-outfluxmut;

	
	
	free(vectorori); free(vectorcomp); free(p); free(q); free(vectororili);
	
	return(flux);
}

void expon_vector(double *es, double *s, double *f, int longvec){
	int i;

	
	for(i=0;i<longvec;i++){
		es[i]=exp(s[i]);		
	}
	
	
}

double suma_pob(double *eX, int longs){
	int i;
	double Suma=0.0;
	
	for(i=0;i<longs;i++){
		Suma+=eX[i];
	}
	
	return(Suma);
}

void llenar_vector_unos(double *s, double *f, int longvec){
	int i;
	
	for(i=0;i<longvec;i++){
		if(f[i]==0.0) {s[i]=-1/0.0;}
		else{s[i]=1.0;}
	}

}



int fitness_aguj_ambas(double *f, double **a, double **b, int k, 
int longs, int numtr, double f_l){
	//longs es longS=2^N, k es el "numero de interaccion/epistasis",
	//**a es la matriz de epistasis de salida y **b de llegada, en *f guardo
	//el resultado. f_l es parametro agujeros
	
	int cuantos=0;
	int i, j, jw, t; 
	
	int kfila;
	int *pmuyaux, *paux;
	double faux=0.0, f_Baux=0.0;
	double *f_B;
	
	f_B = (double *) calloc((longs+1),sizeof(double));
	if(f_B == NULL){printf("ERROR: couldn't allocate memory! \n");}
	
	paux = (int *) calloc((numtr+1),sizeof(int));
	pmuyaux = (int *) calloc((k+1+1),sizeof(int));
	if(paux == NULL || pmuyaux == NULL){
		printf("ERROR: couldn't allocate memory! \n");}
	
	for(i=0;i<longs;i++){f[i]=0.0; f_B[i]=0.0;}
	
	
	//FITNESS DE LA MATRIZ A y B 	
	for(i=0;i<longs;i++){
		pasa_entero_puntero(paux, i, numtr); //guarda en paux el vector asociado a i
		
		for(j=0;j<numtr;j++){
			
			for(t=0;t<k+1;t++){
				pmuyaux[t]=paux[(j+t)%numtr];
			}
			
			kfila = pasa_puntero_entero( pmuyaux , k+1);
			//interacciona 1 trait con otros k traits=k+1 
			faux=faux+a[kfila][j];
			f_Baux=f_Baux+b[kfila][j];
		}
		
		
		for(jw=0;jw<numtr;jw++){paux[jw]=0;}
		for(jw=0;jw<k+1;jw++){pmuyaux[jw]=0;}

		f[i]=faux;
		f_B[i]=f_Baux;
		faux=0.0;
		f_Baux=0.0;
		
		
	}//acaba for(i=0;i<longs;i++) que calcula fitness provisional (suma de las columnas de la matriz)
	
	free(pmuyaux);
	free(paux);
	
	
	//NORMALIZA LA FITNESS ENTRE 0 Y 1:
	for(i=0;i<longs;i++){
		f[i]=f[i]/numtr;
		f_B[i]=f_B[i]/numtr;
	}
	
	for(i=0;i<longs;i++){
		if(f[i]<f_l || f_B[i]<f_l){
			f[i]=0.0;
			f_B[i]=0.0;
			cuantos++;
		}
	}
	
	free(f_B);
	
	
	return(cuantos);
}


int fitness_aguj_anteriores(double *f, double **a, int k, int longs, int numtr){
	// longs es longS=2^N, donde N=numtr=numero de traits, 
	//k es el "numero de interaccion/epistasis", **a es la matriz de epistasis,
	//*f es la fitness (y donde guardo valores nuevos)
	
	int cuantos=0;
	int i, j, jw, t;
	
	int kfila;
	int *pmuyaux, *paux;
	double faux=0.0;
	
	paux = (int *) calloc((numtr+1),sizeof(int));
	pmuyaux = (int *) calloc((k+1+1),sizeof(int));
	
	//NO INICIALIZAR A CERO LAS FITNESS, PERDERIA LOS AGUJEROS QUE TIENEN QUE COINCIDIR EN AMBAS MATRICES
	
	for(i=0;i<longs;i++){
		
		if(f[i]!=0.0){
			pasa_entero_puntero(paux, i, numtr); //guarda en paux el vector asociado a i
			
			for(j=0;j<numtr;j++){
				
				for(t=0;t<k+1;t++){
					pmuyaux[t]=paux[(j+t)%numtr];
				}
				
				kfila = pasa_puntero_entero( pmuyaux , k+1);
				//interacciona 1 trait con otros k traits=k+1
				faux=faux+a[kfila][j];
			}
			
			for(jw=0;jw<numtr;jw++){
				paux[jw]=0;
			}
			for(jw=0;jw<k+1;jw++){
				pmuyaux[jw]=0;
			}
			
			
			f[i]=faux;
			faux=0.0;
			
		}//acaba el calculo de las fitness tales que f[i]!=0.0
		
		//else{no hace nada, continua siendo f[i]=0.0;}
		
		
	}//acaba for(i=0;i<longs;i++) que calcula fitness provisional (suma de las columnas de la matriz)
	
	free(pmuyaux);
	free(paux);
	
	//NORMALIZA LA FITNESS ENTRE 0 Y 1:
	for(i=0;i<longs;i++){
		f[i]=f[i]/numtr;
	}
	
	
	//comprueba que hay el numero correcto de agujeros
	for(j=0;j<longs;j++){
			if(f[j]==0.0){
				cuantos++;
			}
		}
	
	return(cuantos);
}


void intercambiar_columnas(double **enc, int n_filas, int n_columnas, int col_1, int col_2){

	int i;
	double *ayuda;
	
	ayuda = (double *) calloc(n_filas+1,sizeof(double));
	if (ayuda == NULL) {
		printf("ERROR: couldn't allocate memory! \n");
	}
	
	for(i=0;i<n_filas;i++){
		ayuda[i]=enc[i][col_1];
		enc[i][col_1]=enc[i][col_2];
		enc[i][col_2]=ayuda[i];
	}
	
	free(ayuda);
}



//intercambia las filas que suman mayor y menor fitness con la fila cero o ultima, segun corresponda (variable orden)
void colocar_filas(double **enc, int fila, int columna, int orden){
	
	int i,j;
	double *sumarcol,*ayuda;
	double maxim=0.0,minim=1.0;
	int i_max=-1, i_min=-1;
	int orden_baja=-1, orden_alta=-1;
	
	sumarcol = (double *) calloc(fila+1,sizeof(double));
	if (NULL == sumarcol) {
        printf("ERROR: couldn't allocate memory! \n");
    }
	for(i=0;i<fila;i++){
		sumarcol[i]=0.0;
	}
	
	ayuda = (double *) calloc(columna+1,sizeof(double));
	if (NULL == ayuda) {
        printf("ERROR: couldn't allocate memory! \n");
    }
	
	
	for(i=0;i< fila;i++){
		for(j=0;j<columna;j++){
			sumarcol[i]+=enc[i][j];
		}
	}
	

	minim=columna*2;//cada columna es un valor entre 1 y 0
	for(i=0;i<fila;i++){
		if(maxim<sumarcol[i]){
			maxim=sumarcol[i];
			i_max=i;
		}
		if(minim>sumarcol[i]){
			minim=sumarcol[i];
			i_min=i;
		}
	}
	//en este for encontramos la fila con mayor fitness y la de menor fitness
	
	//recolocamos las filas segun convenga
	if(orden==0){
		orden_alta=0;
		orden_baja=fila-1;
	}
	if(orden==1){
		orden_baja=0;
		orden_alta=fila-1;
	}
	
	printf("i_max es %d con maxim=%f\n",i_max,maxim);
	printf("i_min es %d con minim=%f\n",i_min,minim);
	printf("orden_alta es %d y orden_baja es %d\n",orden_alta,orden_baja);
	
	for(j=0;j<columna;j++){
		ayuda[j]=enc[orden_baja][j];
		enc[orden_baja][j]=enc[i_min][j];
		enc[i_min][j]=ayuda[j];
		
		ayuda[j]=enc[orden_alta][j];
		enc[orden_alta][j]=enc[i_max][j];
		enc[i_max][j]=ayuda[j];
	}
	
	
	free(sumarcol);
	free(ayuda);
	
}

//intercambia la fila que suma mayor fitness con la fila cero o ultima, segun corresponda (variable orden)
void colocar_filas_solo_una(double **enc, int fila, int columna, int orden){
	
	int i,j;
	double *sumarcol,*ayuda;
	double maxim=0.0,minim=1.0;
	int i_max=-1, i_min=-1;
	int orden_baja=-1, orden_alta=-1;
	
	sumarcol = (double *) calloc(fila+1,sizeof(double));
	if (NULL == sumarcol) {
        printf("ERROR: couldn't allocate memory! \n");
    }
	for(i=0;i<fila;i++){
		sumarcol[i]=0.0;
	}
	
	ayuda = (double *) calloc(columna+1,sizeof(double));
	if (NULL == ayuda) {
        printf("ERROR: couldn't allocate memory! \n");
    }
	
	
	for(i=0;i< fila;i++){
		for(j=0;j<columna;j++){
			sumarcol[i]+=enc[i][j];
		}
	}
	
	minim=columna*2;//cada columna es un valor entre 1 y 0 (podria poner =sumarcol[0])
	for(i=0;i<fila;i++){
		if(maxim<sumarcol[i]){
			maxim=sumarcol[i];
			i_max=i;
		}
		if(minim>sumarcol[i]){
			minim=sumarcol[i];
			i_min=i;
		}
	}
	//en este for encontramos la fila con mayor fitness y la de menor fitness
	
	//recolocamos las filas segun convenienga
	if(orden==0){
		orden_alta=0;
	}
	if(orden==1){
		orden_alta=fila-1;
	}
	
	printf("i_max es %d con maxim=%f\n",i_max,maxim);
	printf("i_min es %d con minim=%f\n",i_min,minim);
	printf("orden_alta es %d\n",orden_alta);
	
	for(j=0;j<columna;j++){
		ayuda[j]=enc[orden_alta][j];
		enc[orden_alta][j]=enc[i_max][j];
		enc[i_max][j]=ayuda[j];
	}
	
	
	free(sumarcol);
	free(ayuda);
	
}

//intercambia la fila que elijamos con la fila cero o ultima (segun variable orden)
void colocar_elige_una_filas(double **enc, int fila, int columna, int orden, int quefila){
	
	int i,j;
	double *ayuda;
	int orden_alta=-1;
	
	ayuda = (double *) calloc(columna+1,sizeof(double));
	if (NULL == ayuda) {
        printf("ERROR: couldn't allocate memory! \n");
    }
	
	//recolocamos las filas segun convenienga
	if(orden==0){
		orden_alta=0;
	}
	if(orden==1){
		orden_alta=fila-1;
	}
	
	printf("orden_alta es %d\n",orden_alta);
	
	for(j=0;j<columna;j++){
		ayuda[j]=enc[orden_alta][j];
		enc[orden_alta][j]=enc[quefila][j];
		enc[quefila][j]=ayuda[j];
	}
	
	free(ayuda);
	
}

//intercambia dos filas que elijamos
void colocar_elige_dos_filas(double **enc, int fila, int columna, int quefila_1, int quefila_2){
	int j;
	double *ayuda;
	
	ayuda = (double *) calloc(columna+1,sizeof(double));
	if (NULL == ayuda) {
        printf("ERROR: couldn't allocate memory! \n");
    }
	
	
	for(j=0;j<columna;j++){
		ayuda[j]=enc[quefila_1][j];
		enc[quefila_1][j]=enc[quefila_2][j];
		enc[quefila_2][j]=ayuda[j];
	}
	
	free(ayuda);
	
}





double fitness_parametro_fl(double *f, double **a, double **b, int k, 
int longs, int numtr, int agujeros){
	//longs es longS=2^N, k es el "numero de interaccion/epistasis",
	//**a es la matriz de epistasis de salida y **b de llegada, en *f guardo
	//el resultado. f_l es parametro agujeros
	
	double necesitas=0.0;
	int flag_para=0;
	
	double minimo_prov=0.0;
	int nodo=0, nodo_B=0, nose=0, verdaderos=0, verdaderos2=0;
	
	int i, j, jw, t; 
	
	int kfila;
	int *pmuyaux, *paux;
	double faux=0.0, f_Baux=0.0;
	double *f_B;//para comparar con las fitness de la matriz de llegada para hacer agujeros
	
	
	f_B = (double *) calloc((longs+1),sizeof(double));
	paux = (int *) calloc((numtr+1),sizeof(int));
	pmuyaux = (int *) calloc((k+1+1),sizeof(int));
	if (f_B == NULL || paux == NULL || pmuyaux == NULL) {printf("ERROR: couldn't allocate memory! \n");}
		
	
	for(i=0;i<longs;i++){
		f[i]=0.0; f_B[i]=0.0;}
	
	//FITNESS DE LA MATRIZ A y B:
	for(i=0;i<longs;i++){
		pasa_entero_puntero(paux, i, numtr); //guarda en paux el vector asociado a i
		
		for(j=0;j<numtr;j++){
			
			for(t=0;t<k+1;t++){
				pmuyaux[t]=paux[(j+t)%numtr];
			}
			
			
			kfila = pasa_puntero_entero( pmuyaux , k+1);
			//interacciona 1 trait con otros k traits=k+1 
			faux=faux+a[kfila][j];
			f_Baux=f_Baux+b[kfila][j];
		}
		
		for(jw=0;jw<numtr;jw++){paux[jw]=0;}
		for(jw=0;jw<k+1;jw++){pmuyaux[jw]=0;}
		
		f[i]=faux;
		f_B[i]=f_Baux;
		faux=0.0;
		f_Baux=0.0;
		
		
	}//acaba for(i=0;i<longs;i++) que calcula fitness provisional (suma de las columnas de la matriz)
	
	free(pmuyaux);
	free(paux);
	
	
	//NORMALIZA LA FITNESS ENTRE 0 Y 1:
	for(i=0;i<longs;i++){
		f[i]=f[i]/numtr;
		f_B[i]=f_B[i]/numtr;
	}
	

	//CALCULA AGUJEROS DE UNO EN UNO:	
	verdaderos=0;
	while(verdaderos<agujeros){
		nose=0; minimo_prov=0.0;
		while(minimo_prov==0.0){
			minimo_prov=f[nose];  nodo=nose;
			nose++;
			if(nose==longs){minimo_prov=1.0; flag_para=1; break;}
		}
		if(flag_para==1){break;}
		
		
		//SE CALCULA CUAL ES EL SIGUIENTE AGUJERO:
		for(j=0;j<longs;j++){
			if(minimo_prov>f[j] && f[j]!=0.0){
				nodo=j; nodo_B=-1;
				minimo_prov=f[j];
			}
			if(minimo_prov>f_B[j] && f_B[j]!=0.0){
				nodo_B=j; nodo=-1;
				minimo_prov=f_B[j];
			}
		}//termina for(j<longs)
		
		if(nodo_B>=0 && nodo<0){//if(f_B[nodo_B]<f[nodo])
			necesitas=f_B[nodo_B];
			f_B[nodo_B]=0.0;
			f[nodo_B]=0.0;
		}
		else if(nodo>=0 && nodo_B<0){//if(f_B[nodo_B]>f[nodo])
			necesitas=f[nodo];
			f[nodo]=0.0;
			f_B[nodo]=0.0;
		}
		
		
		//COMPROBAMOS QUE SE HAN CREADO LOS AGUJEROS NECESARIOS:
		verdaderos=0;
		for(i=0;i<longs;i++){
			if(f[i]==0.0){verdaderos++;}
		}
		
	}//termina while(verdaderos<agujeros)
	
	verdaderos=0;verdaderos2=0;
	for(i=0;i<longs;i++){
			if(f[i]==0.0){verdaderos++;}
			if(f_B[i]==0.0){verdaderos2++;}
		}
	if(verdaderos!=verdaderos2){printf("\n\nERROR\nERROR, AGUJEROS EN fitness de A=%d y AGUJEROS EN fitness de B=%d\nERROR\n\n",
		verdaderos,verdaderos2);  flag_para=3;
		necesitas=1;
		}
	
	if(flag_para==1){printf("\n\nERROR\nERROR, HAY QUE HACER los %d AGUJEROS EN fitness de B\nERROR\n\n",agujeros);
		necesitas=1;}
	if(flag_para==2){printf("\n\nERROR\nERROR, HAY QUE HACER los %d AGUJEROS EN fitness de A\nERROR\n\n",agujeros);
		necesitas=1;}
			
	free(f_B);
	
	
	return(necesitas);
}




int maximos_locales(int longs,int nsteps_d,double *Fit,int numtr,int K,
double **a,int longa,double Eps,int Cuantos, double **b, int **dis,
int col_salida_1, int col_salida_2, int col_llegada_1, int col_llegada_2, int semilla){
	
	int i,j,iw,jw;
	int no;
	int cadlarga=160;
	char encadenados_MAX[200], cadena[170], cadena2[25];
	FILE *archivoMAX;
	
	memset(encadenados_MAX, 0, sizeof encadenados_MAX);
	memset(cadena, 0, sizeof cadena);
	memset(cadena2, 0, sizeof cadena2);
	
	
	sprintf(cadena2, "N%d", numtr);
	strncat(cadena, cadena2, cadlarga);
	
	sprintf(cadena2, "_K%d", K);
	strncat(cadena, cadena2, cadlarga);
	
	sprintf(cadena2, "_hole%d", Cuantos);
	strncat(cadena, cadena2, cadlarga);
	if(K==0){
		sprintf(cadena2, "_sem%.3d", semilla);
		strncat(cadena,cadena2,cadlarga);
	}
	else{
		sprintf(cadena2, "_col%d-%d_%d-%d", col_salida_1, col_salida_2, col_llegada_1, col_llegada_2);
		strncat(cadena,cadena2,cadlarga);
	}
	
	sprintf(cadena2, "_eps%.3lf", Eps);
	strncat(cadena, cadena2, cadlarga);
	
	sprintf(encadenados_MAX, "maxlocal/max_");
	
	strncat(encadenados_MAX, cadena, cadlarga);
	strcat(encadenados_MAX, ".dat");
	
	
	archivoMAX = fopen(encadenados_MAX, "w");
	
	
	for (jw = 0; jw < nsteps_d; jw++){
	//SE CREA LA LINEA QUE TIENE EL VALOR DE LA FITNESS CORRESPONDIENTE SOLO SI ES UN MAXIMO LOCAL
		for(i=0;i<longs;i++){
			no=0;
			for(j=0;j<longs;j++){
			//dis no tiene en cuenta los agujeros, pero mirar todos los casos en//
			//que dis[i][j]==1 no cambia quienes son los maximos locales 	  //
				if(dis[i][j]==1){
					if(Fit[i]<Fit[j]){
						no=1;
					}
				}
			}//termina for(j=0)
			if(no==0){
				//Fit[i] es maximo local
				fprintf(archivoMAX, "%.8lf\t", Fit[i]);
				fflush(archivoMAX);
			}
			else{
				//Fit[i] NO es maximo local
				fprintf(archivoMAX, "%.8lf\t", 0.0);
				fflush(archivoMAX);
			}
		}//termina for(i=0)
		fprintf(archivoMAX, "\n");
		fflush(archivoMAX);
		
		//actualiza matriz A de epistasis:
		matriz_switch_suma(a,b,longa,numtr);
		
		//nueva fitness sin hacer agujeros nuevos:
		Cuantos=fitness_aguj_anteriores(Fit, a, K, longs, numtr);
		
	}//termina for(..;jw<nsteps_d;..)
	
	
//se calculan los maximos locales del ultimo paso:
	for(i=0;i<longs;i++){
		no=0;
		for(j=0;j<longs;j++){
			if(dis[i][j]==1){
				if(Fit[i]<Fit[j]){
					no=1;
				}
			}
		}//termina for(j=0)
		if(no==0){
			//Fit[i] es maximo local
			fprintf(archivoMAX, "%.8lf\t", Fit[i]);
			fflush(archivoMAX);
		}
		else{
			//Fit[i] NO es maximo local
			fprintf(archivoMAX, "%.8lf\t", 0.0);
			fflush(archivoMAX);
		}
	}//termina for(i=0)
	fprintf(archivoMAX, "\n");
	fflush(archivoMAX);
	
	
	fclose(archivoMAX);
	
	
	
	
	
	return(0);
}



void juzga_guardar_maximo(int *S_A_max, int fila_actual, double *S_A, int longs){
	
	double prov=0.0;
	int nodo;
	int i;
	
	for(i=0;i<longs;i++){
		if(S_A[i]>prov){
			prov=S_A[i];
			nodo=i;
		}
	}
	
	S_A_max[fila_actual]=nodo;
	
}

void juzga_comparar_maximos(int *S_A_max, double *S_B, int dis_salto, 
int filas_S_A_max, int numtr, int col_B_i, int col_B_j, 
int K, int Cuantos, double Beta){
	//dis_salto=N-2, por ejemplo, para pedir la distancia a la que damos por valido el salto
	//filas_S_A_max=(int)(numtr*(numtr-1))/2
	//col_B_i, col_B_j son las columnas que se intercambian en B.
	
	double prov=0.0;
	int nodo, dist_entr_pobla, *p, *q;
	int i,j, longs, no_flag;
	int cadlarga=160;
	FILE *archivoju;
	char cadenadef[200], cadena[200], cadena2[25];
	
	longs = (int)pow(2,numtr);
	
	for(i=0;i<longs;i++){
		if(S_B[i]>prov){
			prov=S_B[i];
			nodo=i;
		}
	}
	
	
	memset(cadena, 0, sizeof cadena);
	memset(cadena2, 0, sizeof cadena2);
	memset(cadenadef, 0, sizeof cadenadef);
	sprintf(cadena2, "N%d", numtr);
	strncat(cadena, cadena2, cadlarga);
	sprintf(cadena2, "_K%d", K);
	strncat(cadena, cadena2, cadlarga);
	sprintf(cadena2, "_hole%d", Cuantos);
	strncat(cadena, cadena2, cadlarga);
	sprintf(cadena2, "_beta%.2lf", Beta);
	strncat(cadena, cadena2, cadlarga);
	sprintf(cadena2, "_colB%d-%d", col_B_i, col_B_j);
	strncat(cadena,cadena2,cadlarga);
	sprintf(cadenadef, "juzga/");
	strncat(cadenadef, cadena, cadlarga);

	strcat(cadenadef, ".dat");
	
	//vectores auxiliares para calcular distancias de hamming
	p= (int *) calloc((numtr+1),sizeof(int));
	if (p == NULL) {printf("ERROR: couldn't allocate memory! \n");}
	q= (int *) calloc((numtr+1),sizeof(int));
	if (q == NULL) {printf("ERROR: couldn't allocate memory! \n");}
	
	no_flag=0;
	for(i=0;i<filas_S_A_max;i++){
		
		dist_entr_pobla=distancia_hamming(nodo,S_A_max[i],p,q,numtr);
		
		if(dist_entr_pobla>=dis_salto){
			no_flag++;
			printf("colB%d-%d\t", col_B_i, col_B_j);// entra %d veces//, no_flag+1
		}
		if(no_flag==1){
			archivoju=fopen(cadenadef, "w");
			fprintf(archivoju, "%s %d %s %d\n", "ColB", col_B_i, "-", col_B_j );
			fflush(archivoju);
		}
		if(dist_entr_pobla>=dis_salto){
			//graba esta combinacion de columnas de A y B a traves de la fila de la matriz:
			printf("ha entrado con fila %d de S_A_max",i);
			fprintf(archivoju, "%s %d %s %d\n", "ColA correspondiente a fila:", i,"con un salto de distancia", dist_entr_pobla);
			fflush(archivoju);
		}
		
		
	}
	
	if(archivoju!=NULL){fclose(archivoju);}
		
	free(p);
	free(q);
	
	
}


void matriz_switch_suma_peso(double **enc, double **enc_dos, int fila, 
int columna, int peso){
	
	int i,j;
	
	for(i=0;i<fila;i++){
		for(j=0;j<columna;j++){
			enc[i][j]=enc[i][j]+peso*enc_dos[i][j];
		}
	}
	
	
}



double actualizacion_dham_alpha(int l, double *fbeta, double *s, int longs, 
int numtr, double Lambda, double Mu, double Beta, double suma,int **dis, double alpha){
	// l es el paso de actualizacion en el que estoy (vector que estoy actualizando)
	// su integer correspondiente lo guardo en vectorori. los demas con 
	// los que comparo los guardo en vectorcomp
	
	
	int jw; 
	int i; //contador de traits
	int j; //contador de vertices
	
	int liaux=0;//liaux es el estado l (que se llama j en el main), con el trait i cambiado
	
	int *vectorori, *vectorcomp, *p, *q;
	int *vectororili;
	
	
	double flux=0.0;
	double influxint=0.0, outfluxint=0.0;
	double influxmut=0.0, outfluxmut=0.0;
	double finf=0.0,fout=0.0, fdiv=0.0;
	//fdiv es el fitness por el que se divide en ambos casos (cambia en cada trait)
	
	double factor_normaliza=pow(1+alpha,1+alpha)/pow(alpha,alpha);//normaliza factor homofilia
	
	if(fbeta[l]==0.0){return(0.0);}//tener en cuenta que 
	//se anula influx porqe f[l]=0.0 y outflux porqe x[l]=0.0 (-> es consistente)
	
	vectorori= (int *) calloc((numtr+1),sizeof(int));
	//vectorori(=l) es el estado que vamos a actualizar con esta funcion
	vectorcomp= (int *) calloc((numtr+1),sizeof(int));
	//vectorcomp son los vectores con los que interacciona (o no). los que 
	//suman o restan poblacion al interaccionar con vectorori
	
	vectororili= (int *) calloc((numtr+1),sizeof(int));//asociado a liaux
	
	//vectores auxiliares para calcular distancias de hamming
	p= (int *) calloc((numtr+1),sizeof(int));
	q= (int *) calloc((numtr+1),sizeof(int));
	
	pasa_entero_puntero(vectorori, l, numtr);
	
	//en el for, i corresponde a los traits de vectorori, el estado que estoy actualizando
	//j corresponde a vectorcomp: los estados con los que se interacciona
	for(i=0;i<numtr;i++){
		
		//vectororili: vectorori (asociado a l) pero con trait i cambiado
		for(jw=0;jw<numtr;jw++){vectororili[jw]=0;}
		for(jw=0;jw<numtr;jw++){vectororili[jw]=vectorori[jw];}
		vectororili[i]=1-vectororili[i];
		
		liaux=pasa_puntero_entero(vectororili,numtr);//liaux es el entero asociado a vectororili
		
		
		fdiv=fbeta[liaux] + fbeta[l];
		finf=fbeta[l]/fdiv;          
		fout= fbeta[liaux]/fdiv;     
		
		
		for(j=0;j<longs;j++){
			
			pasa_entero_puntero(vectorcomp, j, numtr);
			
			//Dentro de for(j...) vectorori, vectorcomp y demas siguen evaluados en i, que es el trait. 
			//(no equivocarse con j que es el nodo)
			//este if PARA LAS CONTRIBUCIONES INFLUX por INTERACCION
			if(vectorori[i] == vectorcomp[i]){
				
				influxint += factor_normaliza*pow(1-(double)(dis[liaux][j])/(double)numtr , alpha)*dis[liaux][j]*finf*(exp(s[liaux]+s[j]-s[l]));//alpha ham
				
			}
			
			
			
			//el else para cambios en CONTRIBUCIONES OUTFLUX por INTERACCION
			else{
				
				outfluxint += factor_normaliza*pow(1-(double)(dis[l][j])/(double)numtr , alpha)*dis[l][j]*fout*(exp(s[j]));//alpha ham
			}
			
			
		} //termina el for(j=0;j<longs;j++)
		
		
		//estos (por mutacion) van fuera de for(j=0;...) porque se suman solo en i:
		influxmut += finf*(exp(s[liaux]-s[l]));
		outfluxmut += fout;
		
		
	} //termina el for (i=0;i<numrt;i++)
	
	
	influxint=(Lambda*influxint)/suma;
	outfluxint=(Lambda*outfluxint)/suma;
	
	influxmut=influxmut*Mu;
	outfluxmut=outfluxmut*Mu;
	
	flux=influxint+influxmut-outfluxint-outfluxmut;
	
	
	free(vectorori); free(vectorcomp); free(p); free(q); free(vectororili);
	
	return(flux);
}


