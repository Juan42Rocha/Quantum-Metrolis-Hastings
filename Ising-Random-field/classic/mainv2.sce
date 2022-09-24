// Este codigo varia el pareameto B, dado una T fija
stacksize('max')
exec('/home/baltazar/Documentos/Maestria/Tesis/1step/scilab/finlibv.sci',-1)

grand('setsd',getdate('s')); ///rand('seed',getdate('s'));

n=10;                                      //Tamanio de la cadena
Temp=0.1;                                 //Temp. de evol
B=1.5;                                     //parametro de desorden, puede variar
J=1.0;                                     //parametro de intereaccion

mts=1000;
iteraciones= mts*n;                          //numero de iteraciones
naprom=5000;//10*10^2;                           //numero de realizaciones a prom

printf("# n=%i, MontStep=%i, Temp=%.2f, reali=%i \n",n,mts,Temp,naprom)
printf('#iteraciones  energiaProm  DominiosProm  MagnetizacionProm \n')

iter=[1:iteraciones];
//Blist=[5.0, 5.0/2, 5.0/3, 5.0/4, 5.0/5,  5.0/7, 5.0/10]         // v1
//Blist=[ 2.0, 1.0/1, 1.0/2, 1.0/4, 1.0/8, 1.0/16, 1.0/32, 1/64 ];  // v2
//Blist=linspace(1/64,0.3,8);                                       // v3
//Blist=ones(5,1).*(1.0/64.0);
//Blist=[0:0.05:1];           //v4
//Blist=[1:2:9]*0.0125;  // 5vlas   [0:0.0125:0.25]           // v5
//Blist=[11:2:31]*0.0125  //v6
//Blist=[0:9]*0.5;         // v7
//Blist=[0:10]*0.25;           // v8
//Blist=[1:0.25:5];           // v9
//Blist=[0:50]*0.1;           // v10

Blist = [0.001  0.01  0.05  0.1  0.5  1] //  v100

//pesos=grand(1,n,'nor',0,1);          //pesos gaussianos
//savematfile("pesosN10.mat","pesos");

loadmatfile('pesosN10.mat');
p=pesos;
pesos=0;
pesos=p(1:n);

for Bcont=1:length(Blist)                         //ciclo para variar B
    B=Blist(Bcont);
    printf("# B=%.4f\n",B)

    energia=zeros(1,iteraciones);
    dominios=zeros(1,iteraciones);
    mag=zeros(1,iteraciones);


    for k=1:naprom                                //ciclos para realizaciones 

    	//S0=arrayrand(n,[-1,1]);                 //arreglo inicial
    	S0=[1 -1 1 -1 1 -1 1 -1 1 -1 ];
   	// pesos=grand(1,n,'nor',0,1);            //pesos gaussiano
    	H0=Energia(S0,J,B*pesos);                 //Energia de S0

    	for i=1:iteraciones                       // ciclo de iteraciones

            S=variaS(S0);                         //Crea nuevo arreglo
        
	        H=Energia(S,J,B*pesos);               //Energia de S
            difE=H-H0;

            energia(i)=energia(i)+H0;                    //Guarda ops energia
	        dominios(i)=dominios(i)+ndominios(S0);       //Guarda ops dominios
	        mag(i)=mag(i)+abs(mean(S0));                 //Guarda ops magnetizacion

  
	        propcamb=grand(1,1,'def');
		    if(difE<0 | propcamb<(1/(1+exp(difE/Temp))))//Si el nuevo estado tiene menor energia nos lo quedamos
                S0=S;H0=H;
            end
    	
	
	    end // fin del ciclo de interaciones
	
	

       	//  MS0f(Bcont,:)=S0;
    end // fin del ciclo de naprom numero de arreglos a promediar

    //MS0f(Bcont,:)=S0;

    menergia=energia/naprom;
    mdominios=dominios/naprom;
    mmag=mag/naprom;
    //vmag=variance(mag,1);


    printf('%f  %f  %f  %f  \n',iter',menergia',mdominios',mmag')
    printf('\n\n\n')

end  // fin del ciclo que cambi el valor de B

printf("#Este archivo deberia tener los valores de B\n")
printf("# Bvals")
printf("   B= %f ",Blist')
mfprintf(0,"%g \n", Bcont/length(Blist))

exit
