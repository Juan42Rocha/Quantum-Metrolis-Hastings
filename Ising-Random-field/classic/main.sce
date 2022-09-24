clear 
stacksize('max')
exec('/home/baltazar/Documentos/Maestria/Tesis/1step/scilab/finlibv.sci',-1)
grand('setsd',getdate('s'))  //rand('seed',getdate('s'));

n=10;                                     //Tamanio de la cadena
Temp=0.0;                                  //Temp. de evol
B=0*1.0/10.0;                                //parametro de desorden, puede variar
J=1.0;                                     //parametro de intereaccion

iteraciones=300*n;                          //numero de iteraciones
naprom=4*(5*10^2);                         //numero de realizaciones a prom

printf("# n=%i, MontStep=%i, Temp=%.2f, itera=%i \n",n,400,Temp,iteraciones)
printf('#iteraciones  energiaProm  DominiosProm  MagnetizacionProm  \n')

iter=[1:iteraciones];

energia=zeros(1,iteraciones);
dominios=zeros(1,iteraciones);
mag=zeros(1,iteraciones);

for k=1:naprom                             //ciclos para realizaciones 

    S0=[1 1 1 1 1 -1 -1 -1 -1 -1 ];//arrayrand(n,[-1,1]);                //arreglo inicial
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
	
        if(difE>0 | propcamb<0.5)//exp(difE/Temp))//Si el nuevo estado tiene menor energia nos lo quedamos
            S0=S;H0=H;
        end
  
    end
    
   
    //MS0f(k,:)=S0;
end

menergia=energia/naprom;
mdominios=dominios/naprom;
mmag=mag/naprom;
//vmag=variance(mag,1);


printf('%f  %f  %f  %f \n',iter',menergia',mdominios',mmag')





//plot2d(menergia)
//ylabel('Energia-promedio')
//xlabel('Iteraciones')
exit


//f=gcf();
//f.color_map=jetcolormap(64);
//grayplot([900:1000],[1:100],MS0f(900:1000,:) )
//colorbar(0,1)
