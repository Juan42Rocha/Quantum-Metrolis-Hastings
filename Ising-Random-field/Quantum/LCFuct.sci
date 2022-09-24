// =====================================================
// Functions for N-qubit Linear Chain         
// =====================================================

// ==================================
// Here we definy global variables
// ==================================

// Quantum Swap gate
   Swapa = [1 0 0 0;
           0 0 1 0;
	   0 1 0 0;
	   0 0 0 1];

// ==================================
// n-qubit basis
// __________________________________
// Nqb: number of qubits
// i  : i-th element of the basis
//    : 0 -> 00...0, 1-> 00...1,...,
//    : n-1-> 11...1
// ==================================
function [vec] = nqbBasis(Nqb, i)
   dimH   = 2^Nqb;
   m      = i+1;
   vec    = zeros(dimH,1); 
   vec(m) = 1.0;
endfunction 

// ==================================
// Projector
// ==================================
function [Pn] = getProjector(Nqb, i)
   ket  = nqbBasis(Nqb,i);
   bra  = ket';
   Pn   = ket*bra;
endfunction

// ==============================================
// Generate single qubit Hadamard gate
// ==============================================
 function [pm]= H()
   isq2= 1.0/sqrt(2.0);
   pm=[isq2, isq2; isq2, -isq2];
 endfunction

// ==============================================
// Generate single qubit NOT gate
// ==============================================
function [mtx] = Sx()
  mtx = [0 1; 1 0];
endfunction




// ==================================================================
// MapNqb:
// This fuction is the extension of Map3qb, to be applied to a linear
// N qubit chain
// ==================================================================
function [lst] = MapNqb(Nqb,indx,bnk,T)
   L  = eye(2^(indx-2), 2^(indx-2));
   R  = eye(2^(Nqb-1-indx), 2^(Nqb-1-indx));
   //
   lst0 = Map3qb(bnk,T);
   lst  = list();
   for element = lst0
     new = kron(L, kron(element, R));
     lst($+1) = new;
   end
endfunction

// ==================================================================
// applyMap: 
// transform a density matrix rho ==> rho'=Sum_j Kj rho K_j^dgr
// INPUT DATA:
// KrausOperators : a 'list(K1,...,Km)' containing the Kraus elements
// rho            : initial density matrix
// OUTPUT DATA:
// newRho         : new density matrix   
// ==================================================================
function [newRho] = applyMap(KrausOperators, rho)
   newRho = 0.0;
   for element = KrausOperators
      //tmp    = element * rho;
      newRho = newRho  + element*rho*(element');
   end
endfunction

// ==============================================
// getBinaryPureState: 
// ______________________________________________ 
// rst   : random number bewteen: m <= rst <= n 
// bst   : binary representation of "rst" at the
//         Nqb basis. 
// qst   : projector: qst = |bst><bst|
// ==============================================
function [bst, qst] = getRandomBaseState(Nqb)
   m   = 0;
   n   = 2^(Nqb)-1;
   rst = floor(m + rand(1,1)*(n-m+1));
   bst = dec2bin(rst,Nqb);
   qst = getProjector(Nqb,rst);
endfunction 

// ==============================================
// getMagneticMomentum
// ______________________________________________ 
// Return the total magnetic momentum of a chain
// of N-qubits.
// ==============================================
function [gmM] = getMagneticMomentum(binState)
   gmM = 0;
   for j = 1:length(binState)
      if part(binState,j)=='0' then
         gmM = gmM + 1;
      else 
         gmM = gmM - 1;
      end 
   end 
endfunction 


// =====================================================
// getDomainNumber
// -----------------------------------------------------
// Return the number of domains of the binary represent-
// ation of one state of the basis.
// =====================================================
function [dnum] = getDomainNumber(binState)
   dnum = 0;
   dim  = length(binState);
   for j = 1:dim-1 
      if (part(binState, j)<>part(binState, j+1)) then
         dnum = dnum + 1;
      end 
   end 
   if (part(binState, dim)<>part(binState, 1)) then
      dnum = dnum +1;
   end
   //if dnum==0 then
   //   dnum = 1;
   //end
endfunction

// =====================================================
// getDomainOfBasis
// -----------------------------------------------------
// Returns an array with the number of domains for each
// state of the basis {|00...0>, |0...01>,...,|11...1>}
// -----------------------------------------------------
function [dom] = getDomainOfBasis(Nqb)
   dimH = 2^Nqb;
   for j=0:dimH-1 
      str = dec2bin(j, Nqb);
      dom(j+1) = getDomainNumber(str);
      clear str;
   end
endfunction 

// ====================================================================
// getInitialStates: 
// This function returns a list of initial quantum pure states from the
// basis B:={|00...0>,|00...01>,...,|11...1>} 
// input data:
// Nqb          : # of qubits in the chain.
// numberStates : # of states in the list.
// output data:
// lqs  : variable of kind struct with the following fields
//      : lqs := {binary, magnetization, quantumState}.
//      : quantumstate ==> vector representation of the state |psi>.
//      : binary       ==> binary representation of the state |psi>.
//      : magnetization==> "0" or "1", if Nqb is "even" or "odd". 
// --------------------------------------------------------------------
function [lqs] = getInitialStates(Nqb, numberStates)
   // ==============================================
   // choose magnetization 0 if Nqb is even and 1 if
   // Nqb iss odd.
   // |0> ==> spin up
   // |1> ==> spin down
   // ----------------------------------------------
   if (modulo(Nqb, 2)==0) then 
      numZeros = Nqb/2;      
      mu       = 0;   // for magnetization 0
   else
      numZeros = (Nqb + 1)/2;
      mu       = 1;   // for magnetization 1
   end
   // ==============================================
   vc = zeros(1, Nqb);
   vc(numZeros+1:Nqb) = 1;
   // ==============================================
   // Make different permutations of the vector vc
   // ----------------------------------------------
   mvc = grand(numberStates, 'prm', vc);

   // ============================================================================
   // Save initial states in a variable of kind list called "initialState"
   // ____________________________________________________________________________
   dimH = 2^(Nqb);
   for i = 1:numberStates 
      // ======================================
      // Binary representation of initial state
      // 
      bn  = ' ';
      for j = 1:Nqb
         bn = bn + string(mvc(i,j));
      end 
      
      // ===============================================
      // vector representation for the initial states
      // 
      j_idx     = bin2dec(bn);
      //rows_cols = [j_idx+1, 1]; // sparse version 
      //non_zeros = [1];
      //mn        = [dimH, 1];
      //   qSt    = sparse(rows_cols, non_zeros, mn);
      //2
      //qSt       = getBasisNqubits(j_idx+1, Nqb);
      qSt = getProjector(Nqb, j_idx);
      // ===================================================================
      // Save initial states in a variable of kind struct with the following
      // fields:
      // ___________________________________________________________________
      //
      lqs(i) = struct('binary', bn, 'magnetization', mu,'matrix', qSt);
   end

endfunction 

// ==================================================================
// bin2projector:
//
function [prj] = bin2projector(bn, Nqb)
   dcimal = bin2dec(bn);
      prj = getProjector(Nqb, dcimal);
endfunction

//function = printInitialState(listStates)
//   disp("===================================================")
//   print(%io(2),"initial states with equal magnetic momentum")
//   disp("===================================================")
//  
//   printf("#   binary      magnetization \n")
//   for j = 1:indk
//       printf("%2d  " + listStates.binary(j) + "     %3d \n",..
//       j, listStates.magnetization(j))
//   end 
//endfunction 


// ==============================================
// projectorToBinary
// ______________________________________________ 
//
function [binary] = projectorToBinary(projector, Nqb)
   sp = sparse(projector);
   [indx, mn, v] = spget(sp);
   binary = dec2bin(indx(1)-1,Nqb);
endfunction 


// ==============================================
// swapOperation
// ______________________________________________ 
// 
// ==============================================
function [so] = swapOperation(Nqb, idx)
   so = 1;
   select idx 
   case 1 then
      for j=1:Nqb-1
          L  = eye(2^(Nqb-j-1),2^(Nqb-j-1));    
          R  = eye(2^(j-1),2^(j-1));
          so = so * kron(L,kron(Swapa,R));
      end 
   case 2
      for j=1:Nqb-1
          L  = eye(2^(j-1),2^(j-1));
          R  = eye(2^(Nqb-j-1),2^(Nqb-j-1));    
          so = so * kron(L,kron(Swapa,R));
      end 
   else 
      print(%io(2),"Choose a valid case: idx = 1 or 2")
   end 
endfunction  

// ============================================================================
// Swap operations.
// Apply successively swap operations to permute the first element 
// of the chain with the last one.
// input data  : 
//   Nqb       : # of qubits in the chain
// output data : 
//   prm       : sparse matrix that change the first and the last element of
//             : the chain 
// ____________________________________________________________________________
//
function [prm] = permuteChain(Nqb)
   dimH     = 2^Nqb;
   prm      = speye(dimH, dimH);
   
   for j = 1:(Nqb-1)
      spid = speye(2^(j-1), 2^(j-1));
      spzr = spzeros(2^(j-1), 2^(j-1));
   
       B11 = [spid,spzr;spzr,spzr];
       B12 = [spzr,spzr;spid,spzr];
       B21 = [spzr,spid;spzr,spzr];
       B22 = [spzr,spzr;spzr,spid];
   
      spswid = [B11,B12;B21,B22];
   
      sp_tmp = spzeros(dimH, dimH);
      for k = 1:2^(Nqb-j-1)
          m = (2^(j+1))*k-(2^(j+1)-1);
          n = (2^(j+1))*k;
          //
          sp_tmp(m:n,m:n) = spswid;
      end 
      prm = prm * sp_tmp;
      prm = full(prm);
   end 
endfunction




// ==============================================
// This function generate an uniformly distribu-
// ted integer in the range [m,n] (including the 
// endpoints).
// To avoid the endpoints of the chain we 
// choose m=2 and n=Nqb-1. 
// Nqb    : integer, # of qubits in the chain.
// nrqn   : integer, random qubit selected.
// ==============================================
function [nrqb] = selectRandomqb(Nqb)
   m = 1;
   n = Nqb;
   nrqb = floor(m + rand(1,1)*(n-m+1));
endfunction

// ==============================================
// getTrace:
// tr[Rho]
// ==============================================
function [tr] = getTrace(Rho)
   tr = trace(Rho);
endfunction 

// ==============================================
// getPurity:
// tr[Rho^2]
// ==============================================
function [P] = getPurity(Rho)
   P = trace(Rho*Rho);
endfunction 

// ==============================================
// getCoherence:
// 2 * ( Sum_{j/=k} |Rho_{jk}| )  
// ==============================================
function [cohe] = getCoherence(Rho)
   tmp  = abs(Rho);
   cohe = max(0, (sum(tmp)-1.0));
endfunction 

// ==============================================
// getDomains: 
// Compute <dom> = Sum_{j} #(j)* rho_{jj}
// #(j) : number of domains of j-th basis confi-
// uration  (0 <= j <= 2^N - 1)
// ==============================================
function [dm] = getDomain(Rho, BasisDomain)  
   dRho = diag(Rho);
   dm   = dRho' * BasisDomain;
endfunction

// ==============================================
// saveDensityMatrix:
// save the matrix elements of Rho
// ==============================================
function [] = saveDensityMatrix(Rho, dir, label)
   filename = dir+string(label)+'.ou';
   unt = file('open', filename, 'unknown')
   dim = size(Rho, 1);
   for i = 1:dim
      write(unt, Rho(i, :), '('+string(dim)+'e16.6)');
   end
   file('close', unt)
endfunction

// ==============================================
// getDiagonalElements:
// This function returns all diagonal elments of
// Rho different from zero.
// ==============================================
function [indx, dElem] = getDiagonalElements(Rho)
   diagonal         = diag(Rho);
   [diagonal,vindx] = gsort(diagonal); // sort in decreasing order  
   m = 1; // local index
   while %t,
      if diagonal(m)==0 | m>length(diagonal) then
         break
      end
      dElem(m) = diagonal(m);
      indx(m)  = vindx(m);
      m = m + 1;
   end 
endfunction 



// ==============================================
// Funciones agregadas para H.B
// ==============================================


// ==============================================
// getDeltaEnergy(s,B,nk,idx):
// This function returns the change fo energy if
// we changes the spine idx
// ==============================================
function dE=getDeltaEnergy(s,b,nk,idx,Nqb)
	
	s  = strtod(part(s,idx));
	
	if (idx==1 ) then
		sm = strtod(part(s,Nqb));
		sp = strtod(part(s,idx,2));
	    
	elseif (idx==Nqb) then
		sm = strtod(part(s,Nqb-1));
		sp = strtod(part(s,1));
 	else
		sm = strtod(part(s,idx-1));
		sp = strtod(part(s,idx+1));
 	end 

	dE = 2*s*(sm+sp)+2*b*nk*s;

endfunction


// ==============================================
// This function help us to construct the hamiltonian
// ==============================================
function [pHam] = PHamConst(Nqb, indx)
	L  = eye(2^(indx-1), 2^(indx-1));
   	R  = eye(2^(Nqb-1-indx), 2^(Nqb-1-indx));
   	//
   	ss=kron(Sx(),Sx());
   	pHam = kron(L, kron(ss, R));
   
   	if(indx==Nqb) then
   		C  = eye(2^(indx-2), 2^(indx-2));
		pHam= kron(Sx(),kron(C,Sx()));
   	end

endfunction


// ==============================================
// This function construct the hamiltonian
// ==============================================
function [pHam] = Hamil(Nqb,b,Pes)
    pHam=0;
	for i=1:Nqb
		L  = eye(2^(i-1), 2^(i-1));
   		R  = eye(2^(Nqb-i), 2^(Nqb-i));
   		//
		pHam=pHam-PHamConst(Nqb,i)+b*Pes(i)*kron(L,kron(Sx(),R));

	end


endfunction

// Function for energy of magnetic field interaction
function [Esb]=getEnergyBField(Nqb,Pesos);
    dimH = 2^Nqb;
	Esb=zeros(dimH,1);
	for j=0:dimH-1 
        str = dec2bin(j, Nqb);		
		for k=1:Nqb
			s=part(str,k);
			if(s=='0') then
				Esb(j+1)=Esb(j+1)+Pesos(k);
			else
				Esb(j+1)=Esb(j+1)-Pesos(k);
			end
		end
      clear str;
    end

endfunction
