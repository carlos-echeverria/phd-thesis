function [ A, b, dyi, dxi ] = block_Schwarz_2D1L_get_linsys( N, xi, yi, C, I, f, EQparams, scaled, disc )
%BLOCK_SCHWARZ_2D1L_GET_LINSYS. 
%   Based on the equation parameters of a convection-diffusion problem:
%
%  -eps*Delta(u) + alpha.Nabla(u) + beta*u = f in Omega,  u = g on Gamma,
%
%   and using the mesh parameters of the discretized domain, this function 
%   generates the corresponding linear system obtained by  using either 
%   upwind or central difference operators for approximating the 
%   first and second order derivatives of the equation. This code uses a
%   for loop for constructing the matrix "row by row".
%
%   function call: 
%
%  [ A, b, dyi,dxi ] = block_Schwarz_2D1L_get_linsys( N, xi, yi, C, I, f, g, EQparams, scaled, disc )
%
%   input:
%          N: number of global intervals in each direction of the mesh.
%          xi: row vector with the x-coordinates of the nodes in the mesh.
%          yi: row vector with the y-coordinates of the nodes in the mesh.
%           C: matirx with the physical coordinates of each node.
%           I: matrix with the numbering of each node.
%           f: discrete right hand side of PDE  with boundary conditions.
%    EQparams: diffusion coefficient, components of convection, etc.
%      scaled: scaling parameter - 0:no, 1:yes.
%        disc: discretization parameter - 1:upwind, 2:central.
%
%   output:
%           A: discrete convection-diffusion matrix.
%           b: discrete right hand side. 
%         dyi: vector with distances to the neighbours in x-direction.
%         dxi: vector with distances to the neighbours in y-direction.
%
% Written by Carlos Echeverria on August 8, 2018.
% Based on a code written by Petr Tichy before 2016.

  epsi = EQparams.epsi;  
    wx = EQparams.wx;   
    wy = EQparams.wy;  
  beta = EQparams.beta; 

n = N-1;  % global number of nodes in each direction
K = n^2;  % total number of nodes (unknowns)  
      
    for i=1:N,
        dxi(i) = xi(i+1)-xi(i);
        dyi(i) = yi(i+1)-yi(i);
    end;

    
    Au  = spalloc(K,K,5*K);
   eps2 = 2 * epsi;
   
 for k=1:K,
        
 %  determine the physical coordinates of node k
        i = C(k,1);
        j = C(k,2);
   

 %  determine distances to neighbours
        hs = dyi(i);
        hw = dxi(j);
        hn = dyi(i+1);
        ho = dxi(j+1);
        
%         hs = 1; % <--- uncomment for discrete laplacian
%         hw = 1;
%         hn = 1;
%         ho = 1;

 %  determine x- and y-size of currents 5-point stencil
        dx = hw + ho;    
        dy = hs + hn;
 %        
        bc(k) = f(i+1,j+1); 
        bu(k) = bc(k);
       
%         if(k==3) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%erase
%             i,j,dxi(3),dxi(4),hs,hw,hn,ho,dx
%         end
        
        
 %  determine entry on the diagonal        
        wyhs = wy/hs;
        wxhw = wx/hw;
        wydy = wy/dy;
        wxdx = wx/dx;
 %        
        entryCc = eps2/(ho*hw) + eps2/(hn*hs)  + beta;
        entryCu = entryCc + wyhs + wxhw;
        scal    = abs(entryCc);
        scal    = 1; %unscaled matrix

 %  number of the neighbour in the north 
        entryNu = - eps2 / (hn*dy);        
        entryNc = entryNu + wydy; 
     %   if(k==20), [entryNu,entryNu/scal], end

        if i < n,
            NN = I(i+1,j);
            tmp = abs(entryNc);
  %          if (tmp > scal), scal = tmp; end;            
        else
            NN = 0;
            bc(k) = bc(k) - f(N+1,j+1)*entryNc; 
            bu(k) = bu(k) - f(N+1,j+1)*entryNu;            
        end;
        
%  number of the neighbour in the south      
        tmp = - eps2 / (hs*dy); 
        entrySc = tmp - wydy;
        entrySu = tmp - wyhs;
       %- eps2 / (hs*dy)- byhs
        if i > 1,
            NS = I(i-1,j);
            tmp = abs(entrySc);
   %         if (tmp > scal), scal = tmp; end;             
        else
            NS = 0;
            bc(k) = bc(k) - f(1,j+1)*entrySc;
            bu(k) = bu(k) - f(1,j+1)*entrySu;            
        end;

%  number of the neighbour in the west       
        tmp = - eps2 / (hw*dx);
        entryWc = tmp - wxdx;
        entryWu = tmp - wxhw;
%       
        if j > 1,
            NW = I(i,j-1);
            tmp = abs(entryWc);
    %        if (tmp > scal), scal = tmp; end;            
        else
            NW = 0;
            bc(k) = bc(k) - f(i+1,1)*entryWc;
            bu(k) = bu(k) - f(i+1,1)*entryWu;            
        end;

%  number of the neighbour in the east
        entryEu = - eps2 / (ho*dx);        
        entryEc = entryEu + wxdx;
 

        if j < n,
            NO = I(i,j+1);
            tmp = abs(entryEc);
     %       if (tmp > scal), scal = tmp; end;                    
        else
            NO = 0;
           bc(k) = bc(k) - f(i+1,N+1)*entryEc;
           bu(k) = bu(k) - f(i+1,N+1)*entryEu;            
        end;
        
        
%  scaling
        if NS > 0,  
            Ac(k,NS) = entrySc/scal; 
            Au(k,NS) = entrySu/scal;
        end;
%        
        if NW > 0,  
            Ac(k,NW) = entryWc/scal; 
            Au(k,NW) = entryWu/scal;
        end;
%        
        Ac(k,k) =  entryCc/scal;
        Au(k,k) =  entryCu/scal;
%        
        if NO > 0,  
            Ac(k,NO) = entryEc/scal; 
            Au(k,NO) = entryEu/scal;
        end;
%        
        if NN > 0,  
            Ac(k,NN) = entryNc/scal; 
            Au(k,NN) = entryNu/scal;
          
        end;
%        
        bc(k) = bc(k)/scal;
        bu(k) = bu(k)/scal;       
        
    %    if(i==4 && j==5), Au(k,NN), end

        
 end;   
   
   
   bu=bu';
   bc=bc';
   
   if disc==1
       A=Au;
       b=bu;
   elseif disc==2
       A=Ac;
       b=bc;
   end
   
   
end