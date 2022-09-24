// Este codigo varia la temperatura dado un  B fijo 
clear 
stacksize('max')
exec('/home/baltazar/Documentos/Maestria/Tesis/1step/scilab/finlib.sci',-1)
grand('setsd',getdate('s'));

n=20;                                    //Tamanio de la cadena
Temp=3.2;                                  //Temp. de evol
B=0.3;                                     //parametro de desorden, puede variar
J=1.0;                                     //parametro de intereaccion

iteraciones=500*n;                          //numero de iteraciones
naprom=4*5*10^2;                           //numero de realizaciones a prom

printf("# n=%i, MontStep=%i, B=%.5f, reali=%i \n",n,50,B,iteraciones)
printf('#iteraciones  energiaProm  DominiosProm  MagnetizacionProm  r\n')

iter=[1:iteraciones];
//Tlist=[0.5:0.5:9.0]      //v1
//Tlist=[1.0:2.0:40]      //v2
Tlist=[0.01:0.01:0.2];   // v3


for Tcont=1:length(Tlist)
Temp=Tlist(Tcont);
printf("# Temp=%.4f\n",Temp)

energia=zeros(1,iteraciones);
dominios=zeros(1,iteraciones);
mag=zeros(1,iteraciones);


for k=1:naprom                             //ciclos para realizaciones 

    S0=arrayrand(n,[-1,1]);                //arreglo inicial
    pesos=B*grand(1,n,'nor',0,1);          //pesos gaussianos
    H0=Energia(S0,J,pesos);                //Energia de S0

    for i=1:iteraciones                    

        S=variaS(S0);                      //Crea nuevo arreglo
        
        H=Energia(S,J,pesos);              //Energia de S
        difE=H0-H;

        energia(i)=energia(i)+H0;                   //Guarda ops energia
  		dominios(i)=dominios(i)+ndominios(S0);       //Guarda ops dominios
		mag(i)=mag(i)+abs(mean(S0));             //Guarda ops magnetizacion

  
        propcamb=grand(1,1,'def');
	
        if(difE>0 | propcamb<(1/(1+exp(difE/Temp))) //Si el nuevo estado tiene menor energia nos lo quedamos
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

//savematfile("arreglsofinalesYpesosN100MS500R2000patodos.mat","MS0f","pesos","Blist");

//plot2d(menergia)
//ylabel('Energia-promedio')
//xlabel('Iteraciones')
exit
