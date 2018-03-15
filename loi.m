function varargout=loi(varargin)
% Assigns input arguments
Mwvec = varargin{1};
RAKE  = varargin{2};
mu    = varargin{3};
DIP   = varargin{4};
Mo    = varargin{5};
Model = varargin{6};

A=zeros(8,numel(Mwvec));
D=zeros(8,numel(Mwvec));
Z=zeros(8,numel(Mwvec));
L=zeros(8,numel(Mwvec));
W=zeros(8,numel(Mwvec));

if nargin>8
    disp('error nargin must be <=8')
end

for k=Model
    if any (k==1) %Leonard 2010a
        %%parameter Leonard
        if (RAKE <=135 && RAKE >=45) % Dip-slip |reverse fault
            C1    = 17.5;                  % Dip-SLIP
            C2    = 3.8e-5;                % Dip-SLIP
            Beta  = 2/3;                     % Beta for Mwvec [5:7]
        end
        if (RAKE <=-45 && RAKE >=-135)   % normal fault
            C1    = 17.5;                  % Dip-SLIP
            C2    = 3.8e-5;                % Dip-SLIP
            Beta  = 2/3;                     % Beta for Mwvec [5:7]
        end
        if(RAKE >-45 && RAKE<=45) %Strike slip
            C1    = 15;                      % Strike-slip
            C2    = 3.7e-5;                  % Strike-slip
            Beta  = 2/3;                     % Beta for Mwvec [5:7]
        end
        if (RAKE >135) && (RAKE<=180)
            C1    = 15;                      % Strike-slip
            C2    = 3.7e-5;                  % Strike-slip
            Beta  = 2/3;                     % Beta for Mwvec [5:7]
        end
        if (RAKE <-135) && (RAKE>=-180) %Strike slip
            C1    = 15;                      % Strike-slip
            C2    = 3.7e-5;                  % Strike-slip
            Beta  = 2/3;                     % Beta for Mwvec [5:7]
        end
        
        %%% Define Model parameters
        %  Preferred values for C1 and C 2:
        % Data                          C1(m1/3)       C2x 10^5        sig (MPa)
        % Interplate dip-slip         17.5 (12-25)   3.8 (1.5-12)      3.0
        % Interplate strike-slip      15.0 (11-20)   3.7 (1.5-9.0)     3.0
        % SCR dip-slip                13.5 (11-17)   7.3 (5.0-10)      5.8
        
        % Earthquake properties and geometry
        L(1)=(Mo./(mu*C1^(3/2)*C2)).^(2/(3*(1+Beta))); % LENGTH : fault LENGTHa in the STRIKE direction (LENGTHa > 0) [m]
        W(1)=C1.*L(1).^Beta; % WIDTH  : fault WIDTHa in the DIP direction (WIDTHa > 0)[m]
        D(1)=C2*C1^(1/2)*L(1).^(0.5*(1+Beta)); % SLIP   : dislocation in RAKE direction [m]
        Z(1)=W(1).*sind(DIP)./2; %gorced to surfaced% DEPTH  : DEPTHa of the fault centroid (DEPTHa > 0) [m]
        % Assigns output arguments
    end
    if any(k==2) %Leonard 2010b
        %%parameter Leonard
        C1dip       = 17.5;                  % Dip-slip
        C1strikemin = 1;                  % Srike-slip 0-3.4 km
        C1strikemoy = 15;                 % Srike-slip 3.4-45 km
        C1strikemax = 19000;              % Srike-slip >45 km
        %C1SCR   = 13.5;                  % Stable Continental Region
        C2dip    = 3.8e-5;                % Dip-slip
        C2strike = 3.7e-5;                % Srike-slip
        %C2SCR   = 7.3;                   % Stable Continental Region
        Betamoy  = 2/3;                   % Beta for   Mwvec [5:7]
        Betamin  = 1;                     % Beta for   Mwvec [1:5]
        Betamax  = 0;                     % Beta for   Mwvec [7:10] if trilineare scaling law
        Beta=zeros(1,numel(Mwvec));
        
        %Beta definition depending on the type of fault and Magnitude
        if (RAKE <=135) && (RAKE >=45) % Dip-slip | reverse fault
            if   Mwvec <= 5.6    % Length for fimit
                Beta = Betamin;
                C1   = C1dip;
                C2   = C2dip;
                %disp ('here1')
            else
                Beta = Betamoy;
                C1   = C1dip;
                C2   = C2dip;
                %disp ('here2')
            end
        end
        
        if (RAKE <=-45) && (RAKE >=-135)  % Normal fault
            if  Mwvec <= 5.6   % Length for fimit
                Beta = Betamin;
                C1   = C1dip;
                C2   = C2dip;
                %disp ('here3')
            else
                Beta = Betamoy;
                C1   = C1dip;
                C2   = C2dip;
                %disp ('here4')
            end
        end
        if (RAKE >-45) && (RAKE<=45)% Beta for Strike-slip
            if  Mwvec< 5.1
                Beta = Betamin;
                C1   = C1strikemin;
                C2   = C2strike;
                %disp ('here5')
            elseif  (Mwvec>=5.1) && (Mwvec<7)
                Beta = Betamoy;
                C1   = C1strikemoy;
                C2   = C2strike;
                %disp ('here5a')
            else
                Beta = Betamax;
                C1   = C1strikemax;
                C2   = C2strike;
                %disp ('here6')
            end
        end
        if (RAKE >135) && (RAKE<=180)
            if  Mwvec< 5.1
                Beta = Betamin;
                C1   = C1strikemin;
                C2   = C2strike;
                %disp ('here7')
            elseif  (Mwvec>=5.1) && (Mwvec<7)
                Beta = Betamoy;
                C1   = C1strikemoy;
                C2   = C2strike;
                %disp ('here8')
            else
                Beta = Betamax;
                C1   = C1strikemax;
                C2   = C2strike;
                %disp ('here9')
            end
        end
        if (RAKE <-135) && (RAKE>=-180)
            if  Mwvec< 5.1
                Beta = Betamin;
                C1   = C1strikemin;
                C2   = C2strike;
                %disp ('here10')
            elseif  (Mwvec>=5.1) && ( Mwvec<7)
                Beta = Betamoy;
                C1   = C1strikemoy;
                C2   = C2strike;
                %disp ('here11')
            else
                Beta = Betamax;
                C1   = C1strikemax;
                C2   = C2strike;
                %disp ('here12')
            end
        end
        L(2)=(Mo./(mu*C1^(3/2)*C2)).^(2/(3*(1+Beta)));% LENGTH : fault length in the STRIKE direction (LENGTH > 0) [m]
        W(2)=C1.*L(2).^Beta;% WIDTH  : fault width in the DIP direction (WIDTH > 0) [m]
        D(2)=C2*C1^(1/2)*L(2).^(0.5*(1+Beta));% SLIP   : dislocation in RAKE direction [m]
        % Earthquake properties and geometry
        Z(2)=W(2).*sind(DIP)/2;% forced to surface% DEPTH  : depth of the fault centroid (DEPTH > 0) [m]
    end
if any (k==3)%Leonard 2010c
    %table 5 Leonard 2010
    if (RAKE <=135 && RAKE >=45)    % Dip-slip |reverse fault
        A (3) = (10^(( Mwvec-4)/1))*1e6;      % m²
        L (3) = 10^(( Mwvec-4.24)/1.67)*1e3;  % m
        W (3) = 10^(0.667*log10(L(3))+1.24); % m
        D (3) = 10^(0.833*log10 (L(3))-3.80);% m
    end
    if (RAKE <=-45 && RAKE >=-135)   % normal fault
        A (3) = (10^(( Mwvec-4)/1))*1e6;
        L (3) = 10^(( Mwvec-4.24)/1.67)*1e3;% m
        W (3) = 10^(0.667*log10(L(3))+1.24); % m
        D (3) = 10^(0.833*log10 (L(3))-3.80);% m
    end
    if(RAKE >-45 && RAKE<=45)
        if  Mwvec< 5.1% Strike-slip
            A (3) = (10^(( Mwvec-3.99)/1))*1e6;
            L (3) = (10^(( Mwvec-4.17)/1.67))*1e3; % m
            W (3) = 10^(0.667*log10(L(3))+1.18) ; % m
            D (3) = 10^(0.5*log10(A(3))-4.43);    % m
        elseif ( Mwvec>=5.1) && ( Mwvec<7)
            A (3) = (10^(( Mwvec-3.99)/1))*1e6;
            L (3) = (10^(( Mwvec-4.17)/1.67))*1e3; % m
            W (3) = 10^(0.667*log10(L(3))+1.18) ; % m
            D (3) = 10^(0.833*log10 (L(3))-3.84); % m
        else
            A (3) = (10^(( Mwvec-3.99)/1))*1e6;
            L (3) = (10^(( Mwvec-4.17)/1.67))*1e3; % m
            W (3) = 10^(0.667*log10(L(3))+1.18) ; % m
            D (3) = 10^(0.5*log10 (L(3))-2.29);   % m
        end
    end
    if (RAKE>135) && (RAKE<=180)
        if  Mwvec< 5.1% Strike-slip
            A (3) = (10^(( Mwvec-3.99)/1))*1e6;
            L (3) = (10^(( Mwvec-4.17)/1.67))*1e3; % m
            W (3) = 10^(0.667*log10(L(3))+1.18) ; % m
            %D (i) = 10^(1*log10 (L(i)));         % m
            D (3) = 10^(0.5*log10(A(3))-4.43);    % m
        elseif (Mwvec>=5.1) && (Mwvec<7)
            A (3) = (10^(( Mwvec-3.99)/1))*1e6;
            L (3) = (10^(( Mwvec-4.17)/1.67))*1e3; % m
            W (3) = 10^(0.667*log10(L(3))+1.18) ; % m
            D (3) = 10^(0.833*log10 (L(3))-3.84); % m
        else
            A (3) = (10^(( Mwvec-3.99)/1))*1e6;
            L (3) = (10^(( Mwvec-4.17)/1.67))*1e3; % m
            W (3) = 10^(0.667*log10(L(3))+1.18) ; % m
            D (3) = 10^(0.5*log10 (L(3))-2.29);   % m
        end
    end
    if (RAKE <-135) && (RAKE>=-180)
        if  Mwvec< 5.1% Strike-slip
            A (3) = (10^(( Mwvec-3.99)/1))*1e6;
            L (3) = (10^(( Mwvec-4.17)/1.67))*1e3; % m
            W (3) = 10^(0.667*log10(L(3))+1.18) ; % m
            D (3) = 10^(0.5*log10(A(3))-4.43);    % m
        elseif ( Mwvec>=5.1) && ( Mwvec<7)
            A (3) = (10^(( Mwvec-3.99)/1))*1e6;
            L (3) = (10^(( Mwvec-4.17)/1.67))*1e3; % m
            W (3) = 10^(0.667*log10(L(3))+1.18) ; % m
            D (3) = 10^(0.833*log10 (L(3))-3.84); % m
        else
            A (3) = (10^(( Mwvec-3.99)/1))*1e6;
            L (3) = (10^(( Mwvec-4.17)/1.67))*1e3; % m
            W (3) = 10^(0.667*log10(L(3))+1.18) ; % m
            D (3) = 10^(0.5*log10 (L(3))-2.29);   % m
        end
    end
    
    %     %SCR stable continental region  to define
    %     .
    %     .
    %     .
    %
    Z(3)=W(3).*sind(DIP)/2;%forced to surface
end
if any (k==4) %YenMa 2011a
        %Table 2 Yen Ma 2011
        A (4) = (10^(-13.79+ 0.87 *log10(Mo)))*1e6; % effective area m²
        L (4) = (10^(-7.46 + 0.47 *log10(Mo)))*1e3;  % effective length m
        W (4) = (10^(-6.30+0.40*log10(Mo)))*1e3;     % effective width m
        D (4) = (10^(-0.65+0.13*log10(Mo)))*1e-2;    % effective displacement m
        % Earthquake properties and geometry
        Z(4)=W(4).*sind(DIP)/2;
end
if any (k==5) %YenMa 2011b
    %Without Rutpure type consideration
    if Mwvec<6.6333 %Mo(i)<10^20 %Nm
        L(5)=(10^((1/2)*log10(Mo)- 8.08))*1e3; % m
        W(5)=(10^((1/2)*log10(Mo)- 8.08))*1e3; % do the assumption that W=L;
        D(5)=(10^(1.68))*1e-2;                     %m
    end
    if Mwvec>=6.6333 %Mo(i)>=10^20
        L(5)=(10^((1/3)*log10(Mo)- 4.84))*1e3;  % m
        W(5)=(10^((1/3)*log10(Mo)- 5.27))*1e3;  % m
        D(5)=(10^((1/3)*log10(Mo)- 4.37))*1e-2; % m
    end
    % Earthquake properties and geometry
    Z(5)=W(5).*sind(DIP)/2;
end
if any (k==6) %YenMa 2011c
    %TabL 2 Yen Ma 2011
    if (RAKE <=135) && (RAKE >=45) % Dip-slip |reverse fault
        A (6) = (10^(-12.45 + 0.80 *log10(Mo)))*1e6; % effective area m²
        L (6) = (10^(-6.66 + 0.42 *log10(Mo)))*1e3;  % effective Lngth m
        W (6) = (10^(-5.76+0.37*log10(Mo)))*1e3;     % effective width m
        D (6) = (10^(-2.01+0.20*log10(Mo)))*1e-2;    % effective displacement m
    end
    if (RAKE <=-45) && (RAKE >=-135) %Normal fault
        A (6) = (10^(-12.45 + 0.80 *log10(Mo)))*1e6; % effective area m²
        L (6) = (10^(-6.66 + 0.42 *log10(Mo)))*1e3;  % effective Lngth m
        W (6) = (10^(-5.76+0.37*log10(Mo)))*1e3;     % effective width m
        D (6) = (10^(-2.01+0.20*log10(Mo)))*1e-2;    % effective displacement m
    end
    if(RAKE >-45) && (RAKE<=45)% Strike-slip
        A (6) = (10^(-14.77 + 0.92*log10(Mo)))*1e6; % effective area m²
        L (6) = (10^(-8.11 + 0.50 *log10(Mo)))*1e3; % effective Lngth m
        W (6) = (10^(-6.67 + 0.42 *log10(Mo)))*1e3; % effective width m
        D (6) = (10^(0.36 + 0.08 *log10(Mo)))*1e-2; % effective Displacement m
    end
    if (RAKE >135) && (RAKE<=180)
        A (6) = (10^(-14.77 + 0.92*log10(Mo)))*1e6; % effective area m²
        L (6) = (10^(-8.11 + 0.50 *log10(Mo)))*1e3; % effective Lngth m
        W (6) = (10^(-6.67 + 0.42 *log10(Mo)))*1e3; % effective width m
        D (6) = (10^(0.36 + 0.08 *log10(Mo)))*1e-2; % effective Displacement m
    end
    if  (RAKE <-135) && (RAKE>=-180)
        A (6) = (10^(-14.77 + 0.92*log10(Mo)))*1e6; % effective area m²
        L (6) = (10^(-8.11 + 0.50 *log10(Mo)))*1e3; % effective Lngth m
        W (6) = (10^(-6.67 + 0.42 *log10(Mo)))*1e3; % effective width m
        D (6) = (10^(0.36 + 0.08 *log10(Mo)))*1e-2; % effective Displacement m
    end
    % Earthquake properties and geometry
    Z(6)=W(6).*sind(DIP)/2;%forced to surface
end
if any (k==7) %WC 1994a
    A (7) = (10^(Mwvec*0.91-3.49))*1e6;        % m²
    W (7) = (10^(Mwvec*0.32-1.01))*1e3 ;        % Rupture Width m
    L (7) = (10^(Mwvec*0.59-2.44))*1e3;         % Surface rupture Length km
    D (7) = (10^(-1.43+0.88*log10(L(7)*1e-3))); % Average displacement m
    % Earthquake properties and geometry
    Z(7)=W(7).*sind(DIP)/2;%forced tosurface
end
if any (k==8) %WC 1994b
    %%depending on the fault type
    if (RAKE <=135) && (RAKE >=45)  % Reverse
        A (8) = (10^(Mwvec*0.98-3.99))*1e6 ;                 % Rupture Area m²
        W (8) = (10^(Mwvec*0.41-1.61))*1e3;                 % Rupture Width m
        L (8) = (10^(Mwvec*0.63-2.86))*1e3;                 % Surface Rupture Length m
        D (8) = (10^(-0.60+0.31*log10(L(8)*1e-3)));         % Average Displacement m |table 2C| not consistent for reverse
    end
    if(RAKE <=-45) && (RAKE >=-135) % Normal
        A (8) = (10^(Mwvec*0.82-2.87))*1e6 ;
        W (8) =(10^(Mwvec*0.35-1.14))*1e3 ;
        L (8) = (10^(Mwvec*0.50-2.01))*1e3;
        D (8) = (10^(-1.99+1.24*log10(L(8)*1e-3)));         % Average Displacement m |table 2C
    end
    if (RAKE >-45) && (RAKE<=45)  % Strike-slip
        A (8) = (10^(Mwvec*0.90-3.42))*1e6 ;
        W (8) = (10^(Mwvec*0.27-0.76))*1e3 ;
        L (8) = (10^(Mwvec*0.74-3.55))*1e3;
        D (8) = (10^(-1.70+1.04*log10(L(8)*1e-3)));         % Average Displacement m |table 2C
    end
    if (RAKE >135) && (RAKE<=180)
        A (8) = (10^(Mwvec*0.90-3.42))*1e6 ;
        W (8) =(10^(Mwvec*0.27-0.76))*1e3 ;
        L (8) = (10^(Mwvec*0.74-3.55))*1e3;
        D (8) = (10^(-1.70+1.04*log10(L(8)*1e-3))); % Average displacement m         % Average Displacement m |table 2C
    end
    if (RAKE <-135) && (RAKE>=-180)
        A (8) = (10^(Mwvec*0.90-3.42))*1e6 ;
        W (8) =(10^(Mwvec*0.27-0.76))*1e3 ;
        L (8) = (10^(Mwvec*0.74-3.55))*1e3;
        D (8) = (10^(-1.70+1.04*log10(L(8)*1e-3)));         % Average Displacement m |table 2C
    end
    % Earthquake properties and geometry
    Z(8)=W(8).*sind(DIP)/2;%forced to surface
end
end
varargout {1} = D;
varargout {2} = Z ;
varargout {3} = L;
varargout {4} = W;
end
% SLIP   : dislocation in RAKE direction [m]
% DEPTH  : DEPTHa of the fault centroid (DEPTHa > 0) [m]
% LENGTH : fault LENGTHa in the STRIKE direction (LENGTHa > 0) [m]
% WIDTH  : fault WIDTHa in the DIP direction (WIDTHa > 0)[m]
% Assigns output arguments





