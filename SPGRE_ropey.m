P = phantom('Modified Shepp-Logan',80);
%%


RFPhase = 0;  %RF phase. radians. 

RFPhase0 = 117 * pi/180; %radians...recommended by zur et al.
%RFPhase0 = pi;
RFPhaseQuad = 0; %quadratic coefficient.
RFPhaseConst = 1; %constant coefficient. usually 2 for RF spoiling. 
RFPhaseLin  = 1; %linear coefficient.
RFOuterCoef = 1; %outer coeff... usually 1/2 for RF spoiling. 


Gmax = 40e-3; %T/mm 

%T1 = 1170e-3; %s
%T2 = 66e-3;   %s
T1 = 500e-3; %s
T2 = 500e-3; %s

TE = 3e-3; %s. 

%TR = 10e-3;
%TR = 10*T1; 

NPEy = 60;
NRO  = 60; 

lengthAlongPEy = 50; %mm.
lengthAlongRO  = 50; %mm.
%The actual object, P, has dimensions lengthAlongPEy x lengthAlongRO.

%IntElementsPerMMy = round(size(P,2)/lengthAlongPEy);
%IntElementsPerMMx = round(size(P,1)/lengthAlongRO);

IntElementsPerMMy = (size(P,2)/lengthAlongPEy);
IntElementsPerMMx = (size(P,1)/lengthAlongRO);

lengthPEy = (linspace(-lengthAlongPEy/2, lengthAlongPEy/2, IntElementsPerMMy * lengthAlongPEy))'; %mm
lengthRO = (linspace(-lengthAlongRO/2, lengthAlongRO/2, IntElementsPerMMx * lengthAlongRO))'; %mm

Object = zeros(size(P,1), size(P,2), 3);
Object(:, :, 3) = P;
%KSpace = zeros(NRO, NPEy, 3);

gamma = 2*pi*42.577478518*10^6; %rad/s/T. Gyromagnetic ratio.

BWpp = 50; %Hz/pixel.  ReadOut bandwidth per pixel. 
dwellTime = 1/(BWpp * NRO);
FOVro  = 50; %mm
FOVpey = 50; %mm
dKro = 1/FOVro;
dKpey = 1/FOVpey;

%unit check:
%(rad/T/s) * (s) * (T/mm) = rad/mm. 

%calculate the max time needed for the phase encoding gradient...
%assuming we ramp up to the max gradient amplitude ...
%you know what...fuck it.  because i'm lazy, i'll just make the PEy and 
%ROP have the same duration. 

%let's calculate the minimum TE.  we'll use that TE.  and of course... our 
%TR will be 2* TE. 
%oh...and we already have our 2D slice (the shepp-logan phantom) ... so
%we'll just use standard matrix rotation for the excitation. 
%we'll use the dwellTime as our increments.  because i'm too lazy to reraster. 

%first time point is for the excitation. 
%TEmin = dwellTime;
%then PEy and ROP:
%TEmin = TEmin + round(NRO/2) * dwellTime;
TEmin = NRO * dwellTime;
%then the readout's first half:
TEmin = TEmin + round(NRO/2) * dwellTime;

%sooooo....
if (TE < TEmin) 
    TE = TEmin;
end


%TR = 2*TE; %s
%TR = 10e-3;
TR = 10 * T1;

%this is a spoiled GRE scan.  
%spoiling will be applied along the Gro direction for now... will
%generalize this further later. 

%net gradient-induced intra voxel dephasing. 
NetPhiPerVOXro = 4 * pi; %rad/voxel
NetPhiPerVOXpe = 0 * pi; %rad/voxel
%NetPhiPerVOXsl = 0 * pi; %rad/voxel

%fadilali define the TimeToROBeginning and the TimeForDeph and RemTimeInTR.
%if we don't want any added net dephasing...

TimeToROBeginning = TE - round(NRO/2) * dwellTime; % ms the time to the RO's beginning. 
TimeForDeph =  round(NRO/2) * dwellTime; %ms the time for dephasing/ pe refocusing. 
TRMin = TimeToROBeginning + (NRO) * dwellTime + TimeForDeph ; %ms

if (TR < TRMin)
    TR = TRMin;
    disp('TR was set to TRMin because it was too small.')
end

RemTimeInTR = TR - (TRMin); %ms. The remaining time in the TR after the dephasing/pe refocusing gradient. 



Gro = dKro * 2 * pi * (BWpp * NRO) / gamma; %Readout gradient. T/mm. 

PhaseROPPMM = gamma * round(NRO/2) * dwellTime * -Gro; %phase/mm  caused by
                                                       %the readout 
                                                       %prephaser
                                                       
%for a B0 map. 
dFArr = zeros(size(P)); %Hz 
dFx = 0;
dFy = 0;

%{
%let's define an off-resonance gradient along the x direction: 
GxOffRes = 0 * Gro;
%GxOffRes = 0.2 * Gro;
%GxOffRes = 0.6 * Gro;
GyOffRes = 0 * Gro;

for n = 1 : size(P,1)
   for m = 1 : size (P,1)
       dFArr(n, m) = (gamma/(2*pi)) * (GxOffRes * squeeze(lengthRO(n,1)) + ...
           GyOffRes * squeeze(lengthPEy(m,1)) ) ;
   end
end
%}

%only at the edges
%dFArr = (P==1) * (-200);


%you'll use this as your prephasers' duration. 
%**************************************************************************
%for now we don't need the following calculation, but it's good to narrate 
%what's going on here.  We calculate the moments of the rop and ro gradients
%and show what the intravoxel dephasing will be a tthe end of applying both of them. 
%then an added lobe will be done to take that to 2*pi intravoxel dephasing. 

M0ro    =  gamma * Gro * (1/BWpp)  ; %rad/mm
M0rop   = -gamma * Gro * (1/BWpp)/2; %rad/mm
%M0rop = 0; %rad/mm
M0roTotPreDeph = M0rop + M0ro;              %rad/mm net after rop and ro. 

dPhiFOVro = M0roTotPreDeph * FOVro; %rad. intravoxel phase dispersion across FOV after rop and ro. 
dPhiVOXro = dPhiFOVro / NRO; %rad/voxel after rop and ro.  %should be pi rad/voxel. good. 

%"remaining" phase needed to have NetPhiPerVOXro intravoxel dephasing. 
dPhiRORem = NetPhiPerVOXro - dPhiVOXro; %rad/voxel caused solely by the dephasing gradient. "delta phi along readout remaining."
M0roGradDeph = dPhiRORem * NRO / FOVro; %rad/mm    caused solely by the dephasing gradient. 

M0roTotPostDeph = M0roGradDeph + M0roTotPreDeph; %rad/mm
%now M0roTotPostDeph * FOVro  should == (NRO* NetPhiPerVOXro)'s value.  good.   
%**************************************************************************

M0peTot = NetPhiPerVOXpe * NPEy / FOVro; %rad/mm. 
%this is the phase accrual you want after the PE refocused. usually zero. 
                                        

prepTRs = 1000;
MPrep = zeros(size(Object,1), size(Object,2), size(Object,3), prepTRs);
%MPrep(:, :, :, 1) = Object(:, :, :);

%df = 0;
%let's have off resonance varying by space an option. 
%ATR = zeros(size(P,1),size(P,2), 3, 3);
%BTR = zeros(size(P,1),size(P,2), 3, 1);

%{
for n = 1 : size(P,1)
    for m = 1: size(P,2)
        %off resonance slope factors along x or y. .. make them zero for now.
        dFx = 0;
        dFy = 0;
        
        df = dFx * (lengthRO(n,1)) + dFy * (lengthPEy(m,1));
        
        [ATR(n, m, :, :), BTR(n, m, :, :)] = freeprecess(TR,T1,T2,df);
    end
end
%}

%off resonance slope factors along x or y. .. make them zero for now.
  
flipAngle = pi/2;%flip angle. radians.  
%flipAngle = acos(exp(-TR/T1)); %ernst angle for spgre. 

%these arrays are for debugging purposes...to make sure we have the proper
%gradient moments.
M0roPerTRArr = zeros(prepTRs + NPEy, 1);
M0pePerTRArr = zeros(prepTRs + NPEy, 1);

RFPhaseArrPreph = zeros(prepTRs, 1);
%just to have each isochromat reach steady state. 
for k = 0 : prepTRs - 1
    %RFPhase = RFphase + RFphaseInc;
    RFPhase = RFOuterCoef * RFPhase0 * ( RFPhaseQuad * (k^2) + RFPhaseLin* k + RFPhaseConst * 2 );
    RFPhaseArrPreph(k+1, 1) = RFPhase;
    exciteMatrix = throt(flipAngle, RFPhase);
    
    %iterate through each spatial location. 
    for n = 1 : size(P,1) %ro axis
        for m = 1 : size(P,2) %pe axis. 
            
            %MPrep(n, m, :, p ) holds the magnetization after p excitations
            %for the location (n, m).  
            
            if (k+1 == 1)
                MPrep(n, m, :, k+1) = exciteMatrix * squeeze(Object(n, m, :)) ;
            else
                %df = dFx * (lengthRO(n,1)) + dFy * (lengthPEy(m,1)); %hz
                
                %M0Moments have units rad/mm.  length has units mm
                % M0Moments * mm /(2*pi * TR) then has units 1/s-> Hz. 
                dfGradMomentPerTR = (M0roTotPostDeph * squeeze(lengthRO(n, 1)) + ...
                    M0peTot * squeeze(lengthPEy(m,1))) / (2 * pi * TR); %hz
                
                df = dFArr(n, m) + dfGradMomentPerTR; %hz
                [ATR, BTR] = freeprecess(TR,T1,T2,df);
                MPrep(n, m, :, k+1) = exciteMatrix * (ATR * squeeze(MPrep(n, m, :, k)) + BTR * norm(squeeze(Object(n, m, :))) );
            end
            
            % so...
            % MPrep(x1, y1, :, k) holds the magnetization immediately after k
            % excitations at location (x1, y1). It does NOT pursue the decay afterwards. 
        end
    end
    M0roPerTRArr(k+1, 1) = M0roPerTRArr(k+1, 1) + M0roTotPostDeph;
    M0pePerTRArr(k+1, 1) = M0pePerTRArr(k+1, 1) + M0peTot;
    disp(k)
end
%%
 figure, plot(abs(squeeze(MPrep(40,40,1,: ) + 1i*MPrep(40,40,2,: ) ) ) )
 figure, imagesc(dFArr), colormap('gray')
%{

MPrepxy = squeeze(MPrep(:, :, 1, : ) + 1i * MPrep(:, :, 2, : ));

figure, plot(abs( squeeze(MPrepxy(40, 61, :)) ) ) 

figure,
imagesc(abs(squeeze(MPrepxy(:, :, size(MPrepxy,3) )))),
colormap('gray')

figure,
imagesc(abs(Object(:, :, 3))),
colormap('gray')
%}
%% now that we've reached steady state... let's do an imaging readout. 

MStart = squeeze(MPrep(:, :, :, size(MPrep,4)));

MPreph = zeros(size(MStart));
angleROP = M0rop * lengthRO;
%calculate the rotatoin about z due to the RO prephaser.
zrotROP = zeros(3, 3, length(lengthRO));
dGpey = dKpey * 2 * pi / (gamma * TimeToROBeginning ) ;
for n = 1 : length(lengthRO)
    zrotROP(:, :, n) = zrot(squeeze(angleROP(n)));
end

%calculate the rotatoin about z due to the RO dephser...
%angleRODephaser = M0roTotPostDeph * lengthRO; %this is actually wrong...
%M0roTotPostDeph is the sum of the RO's preph, RO, and RO deph gradients...
angleRODephaser = M0roGradDeph * lengthRO;
zrotRODephaser = zeros(3, 3, length(lengthRO));
for n = 1 : length(lengthRO)
    zrotRODephaser(:, :, n) = zrot(squeeze(angleRODephaser(n)));
end

zrotPEy = zeros(3, 3, length(lengthPEy));
peInt = 1; 
MReadOutData = zeros(size(Object,1), size(Object,2), 3, NRO, NPEy);
%remember...as we enter this loop...we cannot lose the steady state. 
RFPhaseArr = zeros(NPEy, 1);
for pestep = round(NPEy/2) : -1 : -1 * round(NPEy/2) + 1  %iterate through each phase encoding step
    %disp(peInt)
    disp(pestep)

    kRFPhaseInc = (prepTRs - 1) + peInt;
    RFPhase = RFOuterCoef * RFPhase0 * ( RFPhaseQuad * (kRFPhaseInc^2) + RFPhaseLin * kRFPhaseInc + RFPhaseConst * 2 );
    
    RFPhaseArr(peInt,1) = RFPhase; 
    exciteMatrix = throt(flipAngle, RFPhase);
    
    %first calculate the phase encoding table. 
    Gpey = pestep * dGpey; 
    disp(Gpey)
    
    %(round(NRO/2) * dwellTime) is just the time i designated for
    %the PE prephase. 
    PhasePEyPMM = gamma * Gpey * ( TimeToROBeginning) ; %rad/mm. 
    
    
    
    %untie test:  rad/T/s * T/mm * s = rad/mm. 
    %above is accurate:
    %(rad/T/s) * (T/mm) * (s) = rad/mm.
    
    anglePEy = squeeze(PhasePEyPMM) * lengthPEy; 
    
    %-----------------------------PREPHASERS-------------------------------
    %we are going to iterate through all space and change the magnetization
    %vector at each location according to the ROP, PEy, df phase
    %dispersions...and through T2 decay and T1 recovery. 
    for n = 1 : size(MStart, 1) 
        for m = 1 : size(MStart, 2)
            %df = squeeze(dFx * (squeeze(lengthRO(n,1))) + dFy * (squeeze(lengthPEy(m,1))));
            df = squeeze(dFArr(n, m));
            %{
            %debug
            %[APreph, BPreph] = freeprecess(round(NRO/2) * dwellTime, T1,T2,df);
            %}
            [APreph, BPreph] = freeprecess(TimeToROBeginning, T1,T2,df);
            
             zrotPEy(:, :, m ) = zrot(squeeze(anglePEy(m)));
             
             %calculate the magnetization vector at each location as a
             %result of the PEy gradient and the Readout prephaser. 
                                                       
 MPreph(n, m, : ) = APreph * squeeze(zrotPEy(:, :, m))*squeeze(zrotROP(:, :, n)) * squeeze(MStart(n, m, :)) + ...
     BPreph * norm(squeeze(Object(n, m, :)));
        end
    end
    
    
    M0roPerTRArr(kRFPhaseInc + 1, 1) = M0roPerTRArr(kRFPhaseInc + 1, 1) + M0rop;
    M0pePerTRArr(kRFPhaseInc + 1, 1) = M0pePerTRArr(kRFPhaseInc + 1, 1) + PhasePEyPMM;
    %-----------------------------PREPHASERS END---------------------------
    
    %debugging step
    %figure, imagesc(angle( squeeze(MPreph(:, :, 1)) + 1i*squeeze(MPreph(:, :, 2)) ) ) ,colormap('gray'), title(strcat('pe ', num2str(pestep), 'PhasePEyPMM ', num2str(PhasePEyPMM), ' Gpey ', num2str(Gpey)  ))

    %soo..above we iterated through each isochromat to prephase and phase
    %encode them. 
    
    
 %okay... so we now just pre encoded everything we need...now before the
 %readout.  we are going to be sampling at each dwell time...so we'll need
 %the appropriate decay matrix for each...at eac hlocation (because of off
 %resonance varying spatially).  of the format: 
 %           [Adt, Bdt] = freeprecess(dwellTime, T1,T2,df);
 
 %for each timepoint within the readout... we must go through each position
 %in space to calculate the magnetization rotation.  each rotation will
 %depend on the previous timepoint. 
   

    %debugging when commented out.
    %so now...we will iterate through each isochromat ... to record how its
    %magnetization precesses at each time point of the readout. 
    
    
    %debugging when commented out.
    %so now...we will iterate through each isochromat ... to record how its
    %magnetization precesses at each time point of the readout.  
    for n = 1 : size(MStart, 1) 
        for m = 1 : size(MStart, 2)
            for roInt = 1 : NRO
                if (roInt == 1)
                    MBefore =  squeeze(MPreph(n,m, :));
                else
                    MBefore = squeeze(MReadOutData(n, m, :, roInt - 1, peInt));
                end
                
                %df = dFx * squeeze(lengthRO(n,1)) + dFy * squeeze(lengthPEy(m,1));
                df = squeeze(dFArr(n, m));
               
                df = df + gamma * Gro * squeeze(lengthRO(n))/(2 * pi);
                %unit check:  (rad/(T*s)) * (T/mm) * (mm) * (1/rad) = 1/s.  
                [Adt, Bdt] = freeprecess(dwellTime, T1,T2,df);
                
                MReadOutData(n, m, :, roInt, peInt ) = Adt  * squeeze(MBefore) + Bdt * norm(squeeze(Object(n, m, :)));
            end
        end
    end
    
    M0roPerTRArr(kRFPhaseInc + 1, 1) = M0roPerTRArr(kRFPhaseInc + 1, 1) + gamma * Gro * NRO * dwellTime;
    M0pePerTRArr(kRFPhaseInc + 1, 1) = M0pePerTRArr(kRFPhaseInc + 1, 1) + 0;
    %unit test:   ( rad/(s*T) ) * ( T/mm ) * s = rad/mm. good.  
    
    %by this point all of our spins have gone through the readout... it is
    %now time go through the dephasing gradients and the rest of the TR... 
    %MTREnd = zeros(size(Object, 1), size(Object, 2), size(Object, 3));
     
    %debugging when commented out
    

     
     for n = 1 : size(MStart, 1) 
        for m = 1 : size(MStart, 2)

            
            %unite check:
            %(rad/mm ) / ( (rad/(s*T)) * s ) = T/mm. 
            
            %below is the gradient phase dispersion along the PEy direction
            %solely from the PE refocusing gradient.
            %we will have PhaseRePEyRefPMM.  It should be the exact
            %opposite of PhasePEyPMM to achieve refocusing. 
            %PhaseRePEyRefPMM = -PhasePEyPMM;
    
            GpeyRef = -PhasePEyPMM/(gamma * TimeForDeph);
            
            PhaseRePEyRefPMM = gamma * GpeyRef * TimeForDeph;
            %unit check:
            %( rad/(s*T) ) * (T/mm) * s = rad/mm.
            
            
            angleRePEy = ( PhaseRePEyRefPMM + M0peTot )* lengthPEy; 
            %phase encode rewinder. 
            zrotPEy(:, :, m ) = zrot(squeeze(angleRePEy(m)));
            
            %now we need the readout dephaser.  %remember above we
            %calculated the angle we want to net the readout dephaser.
            
            
 
            %below only involves the decay and inherent off-res during the 
            df = squeeze(dFArr(n, m));
            [ADeph, BDeph] = freeprecess(TimeForDeph, T1,T2,df);
            
            
            
            MDeph =  ...
                (  ...
                        ( ...
                                ADeph * squeeze(zrotPEy(:, :, m)) * ...
                                squeeze(zrotRODephaser(:, :, n)) * squeeze(MReadOutData(n, m, :, NRO, peInt))  ...
                        ) ...
                        + BDeph * norm(squeeze(Object(n, m, :)))  ...
                );
            
            [ARemTimeTR, BRemTimeTR] = freeprecess(RemTimeInTR, T1,T2,df);
            
            MStart(n, m, :) = exciteMatrix * ...
                ( ARemTimeTR * MDeph + BRemTimeTR * norm( squeeze(Object(n, m, :)) ) );

            
        end
     end
    
    M0roPerTRArr(kRFPhaseInc + 1, 1) = M0roPerTRArr(kRFPhaseInc + 1, 1) + M0roGradDeph;
    M0pePerTRArr(kRFPhaseInc + 1, 1) = M0pePerTRArr(kRFPhaseInc + 1, 1) + PhaseRePEyRefPMM + M0peTot;
    
    peInt = peInt + 1;
end

%for debugging purposes.  
RFPhaseArrTot = [RFPhaseArrPreph; RFPhaseArr];

%%

figure, plot(M0roPerTRArr*FOVro/NRO)
title('ro dir gradient-induced intravoxel phase dispersion (rad) each TR')
figure, plot(M0pePerTRArr*FOVpey/NPEy)
title('pe dir gradient-induced intravoxel phase dispersion (rad) each TR')

 %{
MPrephxy = squeeze(MPreph(:, :, 1 ) + 1i * MPreph(:, :, 2 ));


figure,
imagesc(abs(squeeze(MPrephxy(:, :)))),
colormap('gray')
%}
%%
% MReadOutData is 
%  Object's length along readout, object's length along pey, 
% [Mx, My, Mz]', each duration (segmented by dwellTime) along RO, for each
% PE step. 
%KSpace is a vector that holds the complex points during each readout. 
angleRemoveArr = zeros(NPEy, 1);
KSpace = zeros(NRO, NPEy);
for roInt = 1 : NRO
    for peInt = 1 : NPEy
        
        %sum up all of the isochromates for that PE step for that
        %particular point in the readout. 
        
        %just to demonstrate. 
        %rotRFPhaseConsideration = zrot(0); 
        %rotRFPhaseConsideration = zrot(-(RFPhaseArr(peInt) - RFPhaseArr(1))); 
        %rotRFPhaseConsideration = zrot(RFPhaseArr(peInt)); 
        
        kRFPhase = (prepTRs - 1) + peInt;
        angleRemove = RFOuterCoef * RFPhase0 * ( RFPhaseQuad * (kRFPhase^2) + RFPhaseLin * kRFPhase + RFPhaseConst * 2 );
        rotRFPhaseConsideration = zrot(-angleRemove);
        
        angleRemoveArr(peInt , 1) = angleRemove;
        
        for n = 1 : size(Object, 1)
            for m = 1 : size(Object, 2)
                %MTemp = rotRFPhaseConsideration * squeeze(MReadOutData(n, m, :, roInt, peInt)); 
                MTemp =  squeeze(MReadOutData(n, m, :, roInt, peInt)); 
                KSpace(roInt, peInt) = KSpace(roInt, peInt) + (squeeze(MTemp(1)) + ...
                    1i * squeeze(MTemp(2))) ;
            end
        end
        KSpace(roInt, peInt) = KSpace(roInt, peInt) * exp(-1i*angleRemove) ;% exp( 1i*( (kRFPhase +1)*RFPhase0 ) );
        
    end
    disp(roInt)
end

%%
figure,
imagesc(1:NRO, linspace(round(-NPEy/2), round(NPEy/2) - 1, NPEy), abs(KSpace)),
ylabel('phase encoding index')
xlabel(' ro index ')
colormap('gray')

%%
figure,
imagesc(abs(ifftnc(KSpace))),
colormap('gray')

%%
figure,
imagesc(angle(ifftnc(KSpace))),
colormap('gray')
%% real
figure,
imagesc(real(ifftnc(KSpace))),
colormap('gray')
%% imaginary
figure,
imagesc(imag(ifftnc(KSpace))),
colormap('gray')
