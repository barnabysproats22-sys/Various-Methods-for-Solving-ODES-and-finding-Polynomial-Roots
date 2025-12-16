format longE
%=====Breif=Index============
%Functions for Qs 2.1 and 2.2
%Piecewise functions for Q2.3
%Question 2.1
%Question 2.2
%Question 2.3
%Methods to Approximate ODEs
%Evolve, Eocs and Newton
%Piecewise Functions for Q2.3
%Code to Produce Table
%============================

    %Functions for Qs 2.1 and 2.2
fexact = @(t) (1+(t*1.5*0.5))/(1+(t*0.5));
f = @(t,y) (1.5-y)^2;
df = @(t,y) -2*(1.5-y);

    %Piecewise functions for Q2.3
gexact = @(t) exp(t-(2/sqrt(2)));
gpwise = @(t,y) pwise(t,y);
Dgpwise = @(t,y) Dpwise(t,y);

    %Question 2.1
disp('Q2_1')
MakeTable(@backwardEuler,@(t,y) (1.5-y)^2,@(t,y) [-2*(1.5-y)],0,1,10,20,fexact,false)
    %Question 2.2
disp('Q2_2')
MakeTable(@CrkNich,@(t,y) (1.5-y)^2,@(t,y) [-2*(1.5-y)],0,1,10,20,fexact,false)
    %Question 2.3
disp('Q2_3 Forward Euler')
MakeTable(@forwardEuler, @pwise,@Dpwise, 0,1, 1, 20,gexact,true)
hold on
disp('Q2_3 Heun')
MakeTable (@Heun, @pwise,@Dpwise, 0,1, 1, 20,gexact,true)

xlabel('h');
ylabel('error');
legend({'Euler','Heun'});

    %Methods to Approximate ODEs
function y = forwardEuler (f,Df, t0,y0, h)           
      y = y0 + h*f(t0,y0);
end

function y = backwardEuler(f,Df, t0,y0, h)
    n = size(f(t0,y0));
    n = n(1,2);  
    F = @(d) d  - f(t0+h,y0 + h*d);
    DF = @(d) eye(n) - h*Df(t0+h,y0+h*d);
    [y,k] = newton (F, DF, y0, 1e-14,1e2); %epsillon and k?
    y = y0 + h*y; %Correct evolution?
end

function y = CrkNich(f,Df, t0,y0, h)
    n = size(f(t0,y0));
    n = n(1,2);  
    F = @(d) d  - 1/2*(f(t0,y0) + f(t0+h,y0 + h*d));
    DF = @(d) eye(n) - 0.5*h*Df(t0+h,y0+h*d);
    [y,k] = newton (F, DF, y0, 1e-14,1e2); %epsillon and k?
    y = y0 + h*y; %Correct evolution?
end

function y = Heun(f,Df, t0,y0, h)
    F = f(t0,y0);
    y = y0 + (h/2)*(F+f(t0+h,y0+h*F));
end

    %Evolve, Eocs and Newton

function [y,t] = evolve (phi, f,Df, t0,y0, T, N)   
      y = [y0;zeros(N-1,1)];                    % Set up for y
      % t = [t0;zeros(N-1,1)];
      t = t0;                                   % Initialise t
      h = T/N;                                  % Set up step
      for i = 2:N+1
         y(i) = phi(f,Df, t,y(i-1), h);
         % t(i) = t(i-1) + h;
         t = t+h;
      end
end

function eocs = computeEocs (herr)          
    m = size(herr);
    m = m(1,1);
    eocs = zeros(m-1,1);
    for i = 1:m-1
        eocs(i) = log(herr(i+1,2)/herr(i,2)) / (log(herr(i+1,1)/herr(i,1)));
    end
end

function [x, k] = newton (F, DF, x0, eps,K)          
      x = x0; %Initilise x
      k = 0;  %Initialise k
      while abs(F(x))>=eps
          DFeval = DF(x);  %Evaluate Df to allow division
          x = x - transpose(DFeval\transpose(F(x))); %Easier to invert by division. 
          k = k + 1;                          %Do trans as F row vctr but Deval normal          
          %At this point have x_(k+1) and (k+1)
          if k == K
              break
          end  
      end
end

    %Piecewise Functions for Q2.3

function y = pwise(t,y0)
if t>=(1/sqrt(2))
   x = y0;
else
   x = -y0;
end
y = x;
end

function y = Dpwise(t,y0)
if t<(1/sqrt(2))
   y = -1;
else
   y = 1;
end
end


    %Code to Produce Table

function T = MakeTable(phi, u,Du, t0,y0, T, N0,uexact,viz)
M = zeros(11,3);

for i = 0:10
    format short
    M(i+1,1) = i;            % i
    format longE
    M(i+1,2) = T/(N0*2^i);   %h_i
    y = evolve(phi, u,Du, t0,y0, T, N0*2^i); %generate all y
    M(i+1,3) = abs(uexact(T) - y(N0*2^i+1));   %Difference, abs error
end
herr = M(:,[2:3]);
noHeader = [M, [0;computeEocs(herr)]]; %Add Eocs
if viz == true
    loglog(noHeader(:,2),noHeader(:,3));
end
T = array2table(noHeader,'VariableNames',{'i','h_i','|Y(t_n)-y_N|','Eocs'});
end

