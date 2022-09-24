// Copyright (C) 2017 - Corporation - Gustavo Montes
//
// About your license if you have any
//
// Date of creation: 10/02/2017
//
// ==================================
// Load quantum information libraries
// ==================================
//exec('/home/juan/Documentos/Maestria/Tesis/2step/quinf.sci',-1) 
exec('/home/baltazar/Documentos/Maestria/Tesis/2step/LCFuct.sci',-2)
exec('/home/baltazar/Documentos/Maestria/Tesis/2step/Map3qb.sci',-2)
//exec('/home/juan/Documentos/Maestria/Tesis/2step/LCFuct.sci',-2)
//exec('/home/juan/Documentos/Maestria/Tesis/2step/Map3qb.sci',-2)


// ===  Use memory ===
stacksize('max');

grand('setsd',getdate('s'));

// ==========================
//   Definition of parameters
// ==========================
Nqb   = 10;             // Numerode qbits
dim   = 2^Nqb;       // Dimension del vec
nstep = 300;         // numero de pasos mc

T= 0.0;              // Factor de escala
B= 0.000;              // Temp
//dir   = './p_0/Rho'; // directory to save the density matrix
Ncase='Scomp_10_'
// ==================================================================
// Mode selection:
// set 'ModeSelecat' to one of the following options
//    0   : initial congiguration with zero magnetization choosen
//          randomly 
//    1   : initial configuration defined in 'binaryState'
//
//    2   : maximal mixture  id / dim
//
ModeSelect  = 1;
//binaryState = '0101010101';
binaryState = '0000011111';

// ===============================
//   Definition of some  variables
// ===============================
maxNumberStates = 1;

//pesos=grand(1,Nqb,'nor',0,1);   // pesos del campo aleatorio 
loadmatfile('pesosN10.mat');
pesos=pesos.*[1 1 1 1 1 1 1 1 1 1];
CampB=B*pesos;
// =========================
//
so1 = swapOperation(Nqb, 1);
so2 = so1';
//so2 = swapOperation(Nqb, 2);

basisDomain = getDomainOfBasis(Nqb);
HH=2*basisDomain-Nqb-B*getEnergyBField(Nqb,pesos);

// ==================================================================
// Modes configuration
// ==================================================================

if (ModeSelect==0) then

   infoState    = getInitialStates(Nqb, 1);
   initialState = infoState.matrix;
   filename     = './data/'+'N'+string(Nqb)+'_B'+string(B)+'_T'+string(T)+'_Mts'+string(nstep)+'_'+infoState.binary+'_'+string(Nqb)+'.ou';
elseif (ModeSelect==1) then
   initialState = bin2projector(binaryState, Nqb);
   filename     = './data/'+Ncase+'N'+string(Nqb)+'_B'+string(B)+'_T'+string(T)+'_Mts'+string(nstep)+'_'+binaryState+'_'+string(Nqb)+'.ou';
 //   filename='./data/prueba.dat';
else
   maxMixture   = eye(dim, dim) / dim;
   initialState = maxMixture; 
   filename     = './data/IS_maxMixture_'+string(Nqb)+'.ou';
end 

// ==================================================================
//                                           Here start the evolution
// ==================================================================
for j=1:maxNumberStates  

    // ==============================================================
    // Choose initial pure states with same magnetization 

    //
             Rho = initialState;

        trRho(1) = getTrace(Rho);
       domain(1) = getDomain(Rho, basisDomain);
       Purity(1) = getPurity(Rho);
    Coherence(1) = getCoherence(Rho);

       Energy(1) = diag(Rho)'*HH; 
    state(1,:)   = diag(Rho)' ;
	//saveDensityMatrix(Rho, dir, 1);

    for n = 2:nstep // loop for time svolution
       tmp = 0;
       for idx = 1:Nqb  // loop for convex operation
          // idx    = selectRandomqb(Nqb);
	   bnk = CampB(idx);
	
          if (idx==1 ) then
             tmp1     = so2*Rho*so1;
             lstKraus = MapNqb(Nqb, idx+1,bnk,T);
             tmp1     = applyMap(lstKraus, tmp1);
             tmp1     = so1*tmp1*so2;
             //Rho      = so1*Rho*so2;
			elseif (idx==Nqb) then
			   tmp2     = so1*Rho*so2;
			   lstKraus = MapNqb(Nqb, idx-1,bnk,T);
			   tmp2     = applyMap(lstKraus, tmp2);
			   tmp2     = so2*tmp2*so1;
			   //Rho      = so2*Rho*so1;
			else  
			   lstKraus = MapNqb(Nqb, idx,bnk,T);
			   tmp      = tmp + applyMap(lstKraus, Rho);
			end
		 //pause
		 end // end convex operation
		 tmps = tmp + tmp1 + tmp2;
		 Rho = tmps / Nqb;
		 
		 // compute observables
		 trRho(n)     = getTrace(Rho);
		 domain(n)    = getDomain(Rho, basisDomain);
		 Purity(n)    = getPurity(Rho);
		 Coherence(n) = getCoherence(Rho);
		 
		 Energy(n)    = diag(Rho)'*HH; 
		 state(n,:)   = diag(Rho)';

		 //saveDensityMatrix(Rho, dir, n);
		//pause
	  end // end time evolution
	   
	   // ===== save Purity =====
	   unt = file('open', filename, 'unknown')
	   write(unt, [trRho, domain/Nqb, Purity, Coherence, Energy], '(9f12.6)')
	   write(unt,' ')
	   write(unt,' ')
	   write(unt,' ')
	   write(unt,'# Rho')
	   write(unt, state' , '('+string(nstep)+'f10.6)')
	   write(unt,' ')
	   write(unt,' ')
	   write(unt,' ')
	   write(unt,'# Pesos ')
	   write(unt, pesos' , '('+string(Nqb)+'f10.6)')
	   file('close', unt)

  end // end loop for realization

  disp("===================================================")
  disp(" DONE !!!")
  quit();




