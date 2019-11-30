# cultural_dynamics

Using this code you can compute the stationary states of the dynamics described in 
'Epistasis between cultural traits causes paradigm shifts in cultural evolution'.
To do this, you may want to set different values of the parameters defined in the main file: main.c


N is an integer that corresponds to the length of the vector of the NK model. 
In the cultural model N is the number of attributes that can be changed. They are defined in binary form.


K measures the degree of dependence between cultural traits, which in the cultural model is called epistasis. 
If k = 0 each trait contributes additively and independently to the fitness. 
If k > 0 changing a trait affects the contribution of other k traits to the fitness.


beta measures individuals’ intolerance to inconsistencies: 
the larger beta the more prone they are to adopt traits that increase internal consistency.
A value of beta close to zero implies greater reluctance to change a cultural attribute, even if the final state has greater fitness.


agujeros_deseados allows to introduce epistasis in an alternative way: 
those nodes whose fitnees is below a certain threshold are removed. 


lambda tunes the rate of changes in a cultural state that occur through pairwise meetings.
It is defined so that lambda + mu = 1, where mu corresponds to innovation.


The following variables modify the initial (resp. ending) fitness landscapes. 
This is done by swapping two columns of the inital (resp. ending) epistasis matrix. 
Col_salida_1 and Col_salida_2 are the columns that are swapped in the initial epistasis matrix.
Col_llegada_1 and Col_llegada_2 are the columns that are swapped in the ending epistasis matrix.
These parameters range between 0 and N-1.


alpha measures the influence of homophily. The larger alpha, the more influence homophily has in the cultural dynamics.


ini0ant1 allows to obtain different results:
- ini0ant1=1 allows to obtain the successive stationary states between two fitness landscapes (varying slowly from 
the initial landscape A to the ending landscape B).
- ini0ant1=3 must be selected if you want to set alpha different to 1. Again, it allows to obtain the successive stationary states 
between two fitness landscapes (varying slowly from the initial landscape A to the ending landscape B).
- ini0ant1=5 allows to obtain both, the successive stationary states between two fitness landscapes (varying slowly from 
the initial landscape A to the ending landscape B) and the successive stationary states when the variation between 
the two fitness landscapes is done reversely (varying slowly from the ending landscape B to the initial landscape A). 
Select this option to check for hysteresis.


The results are saved in the folder 'res/' (it must be created in the same path where the code is run).
Time to convergence is saved in the folder 'teq/' (it must be created in the same path where the code is run).




