function integratedModel()
% intergratedModel models the cardiopulmonary system based on the papers
% "An integrated mathematical model of the human cardiopulmonary system: 
% model development" and "Interaction between carotid baroregulation
% and the pulsating heart: a mathematical model"
%
% Run the program and select one of the options listed. To change the time
% of the simulation or injection rate select "Set params" followed by
% "Additional model parameters"
%
% Once a model is selected wait for the button to turn from "running" to
% "complete" then click the "Show plots" button to view the plots.

%% Parameters
% Vasculature constants, basal conditions ---------------------------------

% sa = Systemic Arterial
sa.C = 0.28;
sa.Vu = 0;
sa.R = 0.06;
sa.L = 0.00022;

% sp = Splanchnic Arterial
sp.C = 2.05;
sp.Vu = 274.4;
sp.R = 3.307;

% ep = Extrasplanchnic Arterial
ep.C = 1.67;
ep.Vu = 336.6;
ep.R = 1.407;

% sv = Splanchnic Venous
sv.C = 61.11;
sv.Vu = 1121;
sv.R = 0.038;

% ep = Extrasplanchnic Venous
ev.C = 50;
ev.Vu = 1375;
ev.R = 0.016;

% pa = Pulmonary Arterial
pa.C = 0.76;
pa.Vu = 0;
pa.R = 0.023;
pa.L = 0.00018;

% pp = Pulmonary Peripheral
pp.C = 5.8;
pp.Vu = 123;
pp.R = 0.0894;

% pv = Pulmonary Venous
pv.C = 25.37;
pv.Vu = 120;
pv.R = 0.0056;

% Cardiac Parameters ------------------------------------------------------

% la = Left Atrium, lv = Left Ventricle
la.C = 19.23;
la.Vu = 25;
la.R = 0.0025;
lv.P0 = 1.5;
lv.Ke = 0.014;      % /ml; describes end-diastolic PV function of ventricle
lv.Vu = 16.77;
lv.Emax = 2.95;     % mmHg/ml; slope of end-diastolic relationship
lv.Kr = 0.000375;   % s/ml; parameter of dependence of ventricular resistance on
                    % isometric pressure

% ra = Right Atrium, rv = Right Ventricle
ra.C = 31.25;
ra.Vu = 25;
ra.R = 0.0025;
rv.P0 = 1.5;
rv.Ke = 0.011;
rv.Vu = 40.8;
rv.Emax = 1.75;
rv.Kr = 0.00014;

% Control Parameters ------------------------------------------------------

% Carotid sinus afferent pathway------
cs.Pn = 92;         % mmHg; intrasinus pressure at sigmoidal function central pt
cs.ka = 11.758;     % mmHg; slope of static function at central pt
cs.fmin = 2.52;     % spikes/s; minimum signaling freq
cs.fmax = 47.78;    % spikes/s; maximum signaling freq
cs.tauZ = 6.37;     % s; time constant (real pole)
cs.tauP = 2.076;    % s; time constant (real zero)

% Efferent sympathetic pathway------
fes.Inf = 2.10;     % spikes/s; all params describing shape of response curve
fes.o = 16.11;      % spikes/s
fes.Min = 2.66;     % spikes/s
fes.kes = 0.0675;   % s

% Efferent paraympathetic pathway------
fev.Inf = 6.3;           % spikes/s; all params describing shape of response curve
fev.o = 3.2;             % spikes/s
fev.Min = 25;            % spikes/s
fev.kev = 7.06;          % s

% Effectors-------
% end-systolic left ventricle elastance
emaxLv.G = 0.475;   % mmHg*s/ml/spikes; Gain
emaxLv.tau = 8;     % s; time constant of response
emaxLv.D = 2;       % s; time delay of response
emaxLv.o = 2.392;   % mmHg/ml

% end-systolic right ventricle elastance
emaxRv.G = 0.282;
emaxRv.tau = 8;
emaxRv.D = 2;
emaxRv.o = 1.412;

% splanchnic peripheral resistance
rSp.G = 0.695;
rSp.tau = 6;
rSp.D = 2;
rSp.o = 2.49;

% extrasplanchnic peripheral resistance
rEp.G = 0.53;
rEp.tau = 6;
rEp.D = 2;
rEp.o = 0.78;

% splanchnic venous unstressed volume
vUsv.G = -265.4;
vUsv.tau = 20;
vUsv.D = 5;
vUsv.o = 1435.4;

% extrasplanchnic venous unstressed volume
vUev.G = -132.5;
vUev.tau = 20;
vUev.D = 5;
vUev.o = 1537;

% sympathetic modulation of cardiac cycle period, T
sT.G = -0.13;       % s^s/spikes; HR period gain
sT.tau = 2;         % s; time constant
sT.D = 2;           % s; time delay

% parasympathetic modulation of cardiac cycle period, T
pT.G = 0.09;
pT.tau = 1.5;
pT.D = 0.2;

% Additional Model Parameters ---------------------------------------------

hr.Tbasal = 0.833;  % s; basal heart period
hr.Tsys0 = 0.5;     % s; period of systole
hr.ksys = 0.075;    % s^2; constant for determining hr T
hr.T = 60/72;       % s; initial heart rate

other.IR = 0;       % injection rate if positive; hemorrhage rate if negative
other.IRstart = 0;
other.IRend = 0;
other.tsim = 10;

% Parameters from integrated model paper ----------------------------------

% Table 3 lung mechanics parameters
vent.lC = 0.00127;      % larynx compliance, l/cmH2O
vent.tC = 0.00238;      % trachea
vent.bC = 0.0131;       % bronchi
vent.AC = 0.2;          % alveoli
vent.lVu = 34.4;        % larynx unstressed volume, ml
vent.tVu = 6.63;        % trachea
vent.bVu = 18.7;        % bronchi
vent.AVu = 1.263;       % alveoli
vent.mlR = 1.021;       % mouth->larynx resistance, cmH2O-s/l
vent.ltR = 0.3369;      % larynx->trachea
vent.tbR = 0.3063;      % trachea->bronchi
vent.bAR = 0.0817;      % bronchi->alveoli

vent.cwC = 0.2445;      % chest wall compliance
vent.RR0 = 12;          % respiration rate, basal
vent.IEratio = 0.6;     % inspiration/expiration ratio
vent.FRC = 2.4;         % functional residual capacity, liters
vent.PplEE = -5;        % cmH2O
vent.Pmus0 = -5;        % basal value Pmus, cmH2O
vent.tau = 1/5;         % ventilation time constant

% Table 4 Lung gas exchange parameters
lge.FIO2 = 0.210379;    % molar fraction O2 in ambient air
lge.FICO2 = 0.000421;   % molar fraction CO2 in ambient air
lge.K = 1.2103;         % BTPS -> STPD conversion factor
lge.Pb = 760;           % mmHg, barometric pressure
lge.Pws = 47;           % mmHg, water saturation
lge.CsatO2 = 9;         % mmol/l
lge.h1 = 0.3836;        % Hill coefficient O2
lge.alpha1 = 0.03198;   % dissociation constant, 1/mmHg
lge.beta1 = 0.008275;   % dissociation constant, 1/mmHg
lge.K1 = 14.99;         % dissociation constant, mmHg
lge.CsatCO2 = 86.11;    % mmol/l
lge.h2 = 1.819;         % Hill coefficient CO2
lge.alpha2 = 0.05591;   % dissociation constant, 1/mmHg
lge.beta2 = 0.03255;    % dissociation constant, 1/mmHg
lge.K2 = 194.4;         % dissociation constant, mmHg
lge.sh = 1.7;           % shunt fraction
lge.Hgb = 15;           % hemoglobin, g/dl

% Table 5 Tissue gas exchange parameters
tge.cV = 284;           % ml, coronary tissue volume
tge.bV = 1300;          % ml, brain tissue volume
tge.mV = 31200;         % ml, skeletal muscle tissue volume
tge.sV = 2673;          % ml, splanchnic tissue volume
tge.eV = 262;           % ml, extrasplanchnic tissue volume
tge.cO2 = 24;           % ml/min, coronary O2 consumption
tge.bO2 = 47.502;       % ml/min, brain O2 consumption
tge.mO2 = 51.6;         % ml/min, skeletal muscle O2 consumption
tge.sO2 = 108.419;      % ml/min, splanchnic O2 consumption
tge.eO2 = 14.683;       % ml/min, splanchnic O2 consumption
tge.cO2 = 20.16;        % ml/min, coronary CO2 production
tge.bO2 = 39.9017;      % ml/min, brain CO2 production
tge.mO2 = 43.344;       % ml/min, skeletal muscle CO2 production
tge.sCO2 = 91.072;      % ml/min, splanchnic CO2 production
tge.eCO2 = 12.3337;     % ml/min, splanchnic CO2 production
tge.tauLT = 18;         % s, blood transport delay
tge.tauVL = 10;         % s, blood transport delay

% Table 6 Respiratory control parameters
% Peripheral chemoreceptors
chemP.D = 7;            % s, delay
chemP.GA = 1310;        % cmH2O/spike/s, amplitude gain on Pmus
chemP.Gf = 0.8735;      % breaths/min/s, freq gain on Pmus
chemP.tauA = 83;        % s, time constant amplitude
chemP.tauF = 147.78;    % s, time constant freq
chemP.f = 3.7;          % spikes/s, chemmoreceptor firing rate
% Central chemoreceptors
chemC.D = 8;            % s, delay
chemC.GA = 850;         % cmH2O/spike/s, amplitude gain on Pmus
chemC.Gf = 0.9;         % breaths/min/s, freq gain on Pmus
chemC.tauA = 105;       % s, time constant amplitude
chemC.tauF = 400;       % s, time constant freq
chemC.pCO2 = 40;        % mmHg, setpoint

%% GUI
close all
mainFigure.h = figure('visible','off');
paramFigure.h = figure('visible','off');

initmainfigure();


    function initmainfigure()
        initialPosition = [0.1 0.1 0.3 0.6];
        mainFigure.h = makefigure(initialPosition,'Integrated Model');
        setbuttons()
        mainFigure.h.Visible = 'on';
    end

    function genFig = makefigure(initialPosition,name)
        % Creates a generic figure.

        genFig = figure();
        genFig.Name = name;
        genFig.Units = 'normalized';
        genFig.Position = initialPosition;
        genFig.MenuBar = 'none';
        genFig.NumberTitle = 'off';
        genFig.Visible = 'off';
    end

    function setbuttons()
        % create buttons for the figure
        buttonLabels = {'Set params','Run CV model','Run CV model (w/ control)',...
                        'Run resp model','Run full model','Show plots','Exit'};
        firstPosition = [0.05 0.82 0.9 0.12];
        callbacks = {@setparams,{@runmodel,2},{@runmodel,3},{@runmodel,4},...
                     {@runmodel,5},@showplots,'delete(gcf)'};
        mainFigure.btn = cell(1,length(buttonLabels));
        for iBtns = 1:length(buttonLabels)
            position = firstPosition - [0 (iBtns-1)*0.13 0 0];
            mainFigure.btn{iBtns} = makebuttons(mainFigure.h,position,buttonLabels{iBtns},callbacks{iBtns});
        end
    end

    function genBtn = makebuttons(fig,position,label,callback)
        % Creates a generic uicontrol button.

        genBtn = uicontrol();
        genBtn.Parent = fig;
        genBtn.String = label;
        genBtn.Units = 'normalized';
        genBtn.Position = position;
        genBtn.Style = 'pushbutton';
        genBtn.Callback = callback;
    end

%% Set parameters
paramExist = 0;
paramFigure.vasc = figure('Visible','Off');
paramFigure.cardiac = figure('Visible','Off');
paramFigure.control = figure('Visible','Off');
paramFigure.additional = figure('Visible','Off');
paramFigure.lungMech = figure('Visible','Off');
paramFigure.lungGas = figure('Visible','Off');
paramFigure.tissueGas = figure('Visible','Off');
paramFigure.respControl = figure('Visible','Off');

    function setparams(varargin)
        if ~paramExist || ~ishandle(paramFigure.h)
            initparamfigure();
            paramExist = 1;
        end

        paramFigure.h.Visible = 'On';
    end

    function initparamfigure()
        position = get(mainFigure.h,'Position');
        name = 'Set Parameters';
        paramFigure.h = makefigure(position,name);
        addbuttons();
    end

    function addbuttons()
        names = {'Vasculature constants',...
                 'Cardiac parameters',...
                 'Control parameters',...
                 'Additional model parameters',...
                 'Lung mechanics parameters',...
                 'Lung gas exchange parameters',...
                 'Tissue gas exchange parameters',...
                 'Respiratory control parameters',...
                 'Done'};

         callback = {@setvasculatur,@setcardiac,@setcontrol,...
                     @setadditional,@setlungmechanics,...
                     @setlunggasexchange,@settissuegasexchange,...
                     @setrespiratorycontrol,@hideparamfig};

         firstPosition = [0.05 0.87 0.9 0.1];
         for iName = 1:length(names)
             position = firstPosition - [0 (iName-1)*0.105 0 0];
             makebuttons(paramFigure.h,position,names{iName},callback{iName});
         end
    end

    function setvasculatur(varargin)
        if ~isfield(paramFigure,'vascexist') || ~ishandle(paramFigure.vasc)
            name = 'Vasculature control';
            paramFigure.vasc = makesubfig(name,sa,sp,ep,sv,ev,pa,pp,pv);
            paramFigure.vascexist = 1;
        else
            figure(paramFigure.vasc); % Brings to front.
        end

        paramFigure.vasc.Visible = 'On';
    end

    function setcardiac(varargin)
        doesExist = isfield(paramFigure,'cardiacexist') && ishandle(paramFigure.cardiac);
        if ~doesExist
            name = 'Cardiac parameters';
            paramFigure.cardiac = makesubfig(name,la,lv,ra,rv);
            paramFigure.cardiacexist = 1;
        else
            figure(paramFigure.cardiac);
        end

        paramFigure.cardiac.Visible = 'On';
    end

    function setcontrol(varargin)
        doesExist = isfield(paramFigure,'controlexist') && ishandle(paramFigure.control);
        if ~doesExist
            name = 'Control parameters';
            paramFigure.control = makesubfig(name,cs,fes,fev,...
                                         emaxLv,emaxRv,rSp,rEp,vUsv,vUev,...
                                         sT,pT);
            paramFigure.controlexist = 1;
        else
            figure(paramFigure.control);
        end

        paramFigure.control.Visible = 'On';
    end

    function setadditional(varargin)
        doesExist = isfield(paramFigure,'additionalexist') && ishandle(paramFigure.additional);
        if ~doesExist
            name = 'Additional model parameters';
            paramFigure.additional = makesubfig(name,hr,other);
            paramFigure.additionalexist = 1;
        else
            figure(paramFigure.additional);
        end

        paramFigure.additional.Visible = 'On';
    end

    function setlungmechanics(varargin)
        doesExist = isfield(paramFigure,'lungMechexist') && ishandle(paramFigure.lungMech);
        if ~doesExist
            name = 'Lung mechanics parameters';
            paramFigure.lungMech = makesubfig(name,vent);
            paramFigure.lungMechexist = 1;
        else
            figure(paramFigure.lungMech);
        end

        paramFigure.lungMech.Visible = 'On';
    end

    function setlunggasexchange(varargin)
        doesExist = isfield(paramFigure,'lgeexist') && ishandle(paramFigure.lungGas);
        if ~doesExist
            name = 'Lung gas exchange parameters';
            paramFigure.lungGas = makesubfig(name,lge);
            paramFigure.lgeexist = 1;
        else
            figure(paramFigure.lungGas);
        end

        paramFigure.lungGas.Visible = 'On';
    end

    function settissuegasexchange(varargin)
        doesExist = isfield(paramFigure,'tgeexist') && ishandle(paramFigure.tissueGas);
        if ~doesExist
            name = 'Tissue gas exchange parameters';
            paramFigure.tissueGas = makesubfig(name,tge);
            paramFigure.tgeexist = 1;
        else
            figure(paramFigure.tissueGas);
        end

        paramFigure.tissueGas.Visible = 'On';
    end

    function setrespiratorycontrol(varargin)
        doesExist = isfield(paramFigure,'respControlexist') && ishandle(paramFigure.respControl);
        if ~doesExist
            name = 'Respiratory control parameters';
            paramFigure.respControl = makesubfig(name,chemP,chemC);
            paramFigure.respControlexist = 1;
        else
            figure(paramFigure.respControl);
        end

        paramFigure.respControl.Visible = 'On';
    end

    function hideparamfig(varargin)
        paramFigure.h.Visible = 'Off';
    end

    function subfig = makesubfig(varargin)
        % first two inputs should be position and name rest should be
        % structures to change.

        position = [0.5 0.1 0.7 0.7];
        name = varargin{1};
        structs = cell(1,length(varargin)-1);
        structName = structs;

        for iVar = 2:length(varargin)
            structs{iVar-1} = varargin{iVar};
            structName{iVar-1} = inputname(iVar);
        end

        subfig = makefigure(position,name);
        fig = gcf;
        firstPosition = [0.01 0.95 0.1 0.03];
        DROP = [0 -0.04 0 0];
        SHIFTRIGHT = [0.07 0 0 0];
        SHIFTROW = [0.2 0 0 0];
        newPosition = firstPosition;

        for iStruct = 1:length(structs)
            currStruct = structs{iStruct};
            fields = fieldnames(currStruct);

            if pastBottom(newPosition(2),length(fields))
                firstPosition = firstPosition + SHIFTROW;
                newPosition = firstPosition;
            end

            maketxt(fig,newPosition,structName{iStruct},1);
            newPosition = newPosition + SHIFTRIGHT/2;

            for iParam = 1:length(fields)
                currParam = fields{iParam};
                newPosition = newPosition + DROP;
                maketxt(fig,newPosition,currParam,0);
                paramFigure.(structName{iStruct}).(currParam) = makeedt(fig,newPosition+SHIFTRIGHT,currStruct.(currParam));
            end

            newPosition = newPosition + DROP - SHIFTRIGHT/2;
        end

        lastBox = paramFigure.(structName{length(structs)}).(fields{end});
        resizefig(fig,lastBox);
        makebuttons(gcf,[0.55 0.02 0.4 0.07],'Set',{@setValue,structs,structName,fig});
        makebuttons(gcf,[0.1 0.02 0.4 0.07],'Cancel',{@cancel,fig});
    end

    function bool = pastBottom(currDepth,fieldLen)
        % Determines if the next set of variables would go beyond the bottom of
        % the window.

        bool = currDepth - 0.04*fieldLen < 0;
    end

    function genTxt = maketxt(fig,position,string,bold)
        % Creates a generic non-editable textbox. Bold is a boolean
        % indicating whether the label should be bolded or not.

        genTxt = annotation(fig,'textbox');
        genTxt.Units = 'normalized';
        genTxt.Interpreter = 'Tex';
        genTxt.FontSize = 8;
        genTxt.LineStyle = 'none';
        genTxt.Position = position;
        genTxt.String = string;

        if bold
            genTxt.FontWeight = 'bold';
        end
    end

    function genEdt = makeedt(fig,position,string)
        % Creates a generic user editable textbox.
        genEdt = uicontrol();
        genEdt.Parent = fig;
        genEdt.FontSize = 7;
        genEdt.Style = 'edit';
        genEdt.Units = 'normalized';
        genEdt.Position = position;
        genEdt.String = string;
    end

    function resizefig(fig,lastBox)
        % Change units to pixels so all objects don't change size with
        % figure.
        fig.Units = 'Pixels';
        objs = get(fig,'children');
        objs = [objs',findall(fig,'type','TextBox')'];

        for iUnits = 1:length(objs)
            objs(iUnits).Units = 'Pixels';
        end

        BUFFER = 10;
        figPosition = fig.Position;
        boxPosition = lastBox.Position;
        boxEdge = boxPosition(1) + boxPosition(3); % Edge with respect to screen instead of window.
        figWidth = boxEdge + BUFFER;
        figPosition(3) = figWidth;
        fig.Position = figPosition;
        fig.Units = 'Normalized';

        for iUnits = 1:length(objs)
        	objs(iUnits).Units = 'Normalized';
        end
    end

    function setValue(varargin)
        structs = varargin{3};
        structName = varargin{4};
        fig = varargin{5};

        for iStruct = 1:length(structs)
            currStruct = structs{iStruct};
            currName = structName{iStruct};
            fields = fieldnames(currStruct);

            for iParam = 1:length(fields)
                currParam = fields{iParam};
                newValue = str2double(paramFigure.(currName).(currParam).String);
                currStruct.(currParam) = newValue;
            end

            updateparams(currStruct,currName)
        end
        close(fig);
    end

    function updateparams(structure,name)
        % Updates global structure names since matlab makes a copy of the
        % structures when they're used as inputs.

        switch name
            case 'sa'
                sa = structure;
            case 'sp'
                sp = structure;
            case 'ep'
                ep = structure;
            case 'sv'
                sv = structure;
            case 'ev'
                ev = structure;
            case 'pa'
                pa = structure;
            case 'pp'
                pp = structure;
            case 'pv'
                pv = structure;
            case 'la'
                la = structure;
            case 'lv'
                lv = structure;
            case 'ra'
                ra = structure;
            case 'rv'
                rv = structure;
            case 'cs'
                cs = structure;
            case 'fes'
                fes = structure;
            case 'fev'
                fev = structure;
            case 'emaxLv'
                emaxLv = structure;
            case 'emaxRv'
                emaxRv = structure;
            case 'rSp'
                rSp = structure;
            case 'rEp'
                rEp = structure;
            case 'vUsv'
                vUsv = structure;
            case 'vUev'
                vUev = structure;
            case 'sT'
                sT = structure;
            case 'pT'
                pT = structure;
            case 'hr'
                hr = structure;
            case 'other'
                other = structure;
            case 'vent'
                vent = structure;
            case 'lge'
                lge = structure;
            case 'tge'
                tge = structure;
            case 'chemP'
                chemP = structure;
            case 'chemC'
                chemC = structure;
            otherwise
                return
        end

    end

    function cancel(varargin)
        fig = varargin{3};
        fig.Visible = 'Off';
    end

%% Model functions
global tout yout
global modelNum

    function runmodel(varargin)
        clear('tout','yout');
        modelNum = varargin{3};
        mainFigure.btn{modelNum}.String = 'Running';
        drawnow();
        tend = other.tsim; % time of simulation in seconds
        tspan = [0,tend];

        switch modelNum
            case 2
                runcvmodel(tspan);
            case 3
                runcvcontrolmodel(tspan);
            case 4
                runrespmodel(tspan);
            case 5
                runfullmodel(tspan);
            otherwise
                return
        end

        mainFigure.btn{modelNum}.String = 'Complete';
    end

    function runcvmodel(tspan)
        ic = createcvic();
        [tout,yout] = ode45(@(t,y)cvmodel(t,y),tspan,ic);
    end

    function ic = createcvic()
      % Creates initial condition (begining systole) vector for cardiovascular
      % model.

      ic = zeros(14,1);
      ic(1)  = 25;          % Ppa
      ic(2)  = 100;         % Fpa (20-8)/R_pa
      ic(3)  = 7;           % Ppp
      ic(4)  = 5;           % Ppv
      ic(5)  = 120;         % Psa
      ic(6)  = 100;         % Fsa (120-80)/R_sa
      ic(7)  = 90;          % Psp
      ic(8)  = 10;          % Pev
      ic(9)  = 5300;        % Vt
      ic(10) = 8;           % Pla
      ic(11) = 142;         % Vlv
      ic(12) = 0;           % u (start at beginning of systole.)
      ic(13) = 2;           % Pra
      ic(14) = rv.Vu;       % Vrv
    end

    function out = cvmodel(t,ic)
        % Just to keep everything consistent.
        out = cvfunctions(t,ic);
    end

global dDeltaVUev Fir
    function out = cvfunctions(t,ic)
        % Produces the system of equations for model in the paper "Interaction
        % between carotid baroregulation and the pulsating heart". For the
        % non-controlled functions (1-29).
        
        Ppa = ic(1); Fpa = ic(2); Ppp = ic(3); Ppv = ic(4);
        Psa = ic(5); Fsa = ic(6); Psp = ic(7); Pev = ic(8);
        Vt = ic(9); Pla = ic(10); Vlv = ic(11);
        Pra = ic(13); Vrv = ic(14);

        if length(ic) == 14         % For constant T; no baroregulation
            T = hr.T;
            spR = sp.R;
            epR = ep.R;
            evVu = ev.Vu;
            svVu = sv.Vu;
            dVuev = 0;

        else                        % For baroregulation
            deltaRSp = ic(18); deltaREp = ic(19); deltaVUsv = ic(20);
            deltaVUev = ic(21); deltaTs = ic(22); deltaTv = ic(23);

            T = deltaTs + deltaTv + hr.T;                              %Eqn(42)

            %Eqn 37
            spR = deltaRSp + rSp.o;
            epR = deltaREp + rEp.o;
            evVu = deltaVUev + vUev.o;
            svVu = deltaVUsv + vUsv.o;
            dVuev = dDeltaVUev;
        end

        [For,Fol,Fil,Fir] = ventriclefunctions(ic);

        Vu = sa.Vu + sp.Vu + ep.Vu + svVu + evVu + ra.Vu + ...       %Eqn(10)
              pa.Vu + pp.Vu + pv.Vu + la.Vu;

        sv.P = (Vt - sa.C*Psa - (sp.C + ep.C)*Psp - ...              %Eqn(9)
                ev.C*Pev - ra.C*Pra - Vrv - pa.C*Ppa - ...
                pp.C*Ppp - pv.C*Ppv - la.C*Pla - Vlv - Vu)/sv.C;

        % differintial eqns
        du = 1/T;
        dPpa = (For - Fpa)/pa.C;                                     %Eqn(1)
        dFpa = (Ppa - Ppp - pa.R*Fpa)/pa.L;                          %Eqn(2)
        dPpp = (Fpa - (Ppp - Ppv)/pp.R)/pp.C;                        %Eqn(3)
        dPpv = ((Ppp - Ppv)/pp.R - (Ppv - Pla)/pv.R)/pv.C;           %Eqn(4)
        dPsa = (Fol - Fsa)/sa.C;                                     %Eqn(5)
        dFsa = (Psa - Psp - sa.R*Fsa)/sa.L;                          %Eqn(6)

        dPsp = (Fsa - (Psp - sv.P)/spR - ...                         %Eqn(7)
               (Psp -Pev)/epR)/(sp.C + ep.C);

        dPev = ((Psp - Pev)/epR - (Pev - Pra)/ev.R - dVuev)/ev.C;    %Eqn(8)


        if t > other.IRstart && t < other.IRend
            dVt = other.IR;                                          %Eqn(11)
        else
            dVt = 0;
        end

        dPla = ((Ppv - Pla)/pv.R - Fil)/la.C;                        %Eqn(12)

        dVlv = Fil - Fol;                                            %Eqn(14)
        dPra = ((sv.P-Pra)/sv.R + (Pev-Pra)/ev.R - Fir)/ra.C;        %Eqn(23)

        dVrv = Fir - For;                                            %Eqn(25)

        out = [dPpa,dFpa,dPpp,dPpv,dPsa,dFsa,dPsp,...
               dPev,dVt,dPla,dVlv,du,dPra,dVrv]';
    end

    function [For,Fol,Fil,Fir,lvP] = ventriclefunctions(ic)
        Ppa = ic(1); Psa = ic(5); Pla = ic(10); Vlv = ic(11);
        u = ic(12); Pra = ic(13); Vrv = ic(14);

        if length(ic) == 14         % For constant T; no baroregulation
            T = hr.T;
            lvEmax = lv.Emax;
            rvEmax = rv.Emax;

        else                        % For baroregulation
            deltaEmaxLv = ic(16); deltaEmaxRv = ic(17);
            deltaTs = ic(22); deltaTv = ic(23);

            T = deltaTs + deltaTv + hr.T;                              %Eqn(42)

            %Eqn 37
            lvEmax = deltaEmaxLv + emaxLv.o;
            rvEmax = deltaEmaxRv + emaxRv.o;
        end

        Tsys = hr.Tsys0 - hr.ksys/T;                                 %Eqn(22)

        u = mod(u,1);
        if u <= Tsys/T;                                             %Eqn(19)
            phi = (sin(pi*T/Tsys*u))^2;
        else
            phi = 0;
        end

        Pmaxlv = phi*lvEmax*(Vlv - lv.Vu) + ...                      %Eqn(18)
                  (1 - phi)*lv.P0*(exp(lv.Ke*Vlv) - 1);

        Pmaxrv = phi*rvEmax*(Vrv - rv.Vu) + ...                      %Eqn(29)
                  (1 - phi)*rv.P0*(exp(rv.Ke*Vrv) - 1);

        lvR = lv.Kr*Pmaxlv;                                          %Eqn(16)
        rvR = rv.Kr*Pmaxrv;                                          %Eqn(27)

        if Pmaxrv <= Ppa                                             %Eqn(26)
          For = 0;
        else
          For = (Pmaxrv - Ppa)/rvR;
        end

        rvP = Pmaxrv - rvR*For;                                      %Eqn(28)

        if Pra <= rvP                                                %Eqn(24)
          Fir = 0;
        else
          Fir = (Pra - rvP)/ra.R;
        end

        if Pmaxlv <= Psa                                             %Eqn(15)
          Fol = 0;
        else
          Fol = (Pmaxlv - Psa)/lvR;
        end

        lvP = Pmaxlv - lvR*Fol;                                      %Eqn(17)
        
        if Pla <= lvP                                                %Eqn(13)
          Fil = 0;
        else
          Fil = (Pla - lvP)/la.R;
        end
    end

    function runcvcontrolmodel(tspan)
        ic = [createcvic();createcvcontrolic()];
        [tout,yout] = ode113(@(t,y)cvcontrolmodel(t,y),tspan,ic);
    end

    function ic = createcvcontrolic()
        ic = zeros(9,1);
        ic(1) = 9;          % Ptilda
        % All others should be zero.
    end

    function out = cvcontrolmodel(t,ic)
        outcontrol = cvcontrolfunctions(t,ic);
        outcv = cvfunctions(t,ic);
        out = [outcv;outcontrol];
    end

    function out = cvcontrolfunctions(t,ic)
        % Control functions (eqns 30-42) from "Interaction between carotid
        % and baroregulation and the pulsating heart". To be used with cv
        % functions.

        Psa = ic(5);

        Ptilda = ic(15); deltaEmaxLv = ic(16); deltaEmaxRv = ic(17);
        deltaRSp = ic(18); deltaREp = ic(19); deltaVUsv = ic(20);
        deltaVUev = ic(21); deltaTs = ic(22); deltaTv = ic(23);

        dPtilda = (Psa - Ptilda)/cs.tauP;                            %Eqn 30
        fcs = (cs.fmin + cs.fmax*exp((Ptilda - cs.Pn)/cs.ka))/...    %Eqn 31
              (1 + exp((Ptilda - cs.Pn)/cs.ka));

        % 1s used after variable name to differentiate between the structures.
        fcso = 25.15; % Central value occours when Ptilda = Pn; fcs = mean(fmin,fmax)
        fes1 = fes.Inf + (fes.o -  fes.Inf)*exp(-fes.kes*fcs);       %Eqn 33
        fev1 = (fev.o + fev.Inf*exp((fcs - fcso)/fev.kev))/...       %Eqn 34
              (1 + exp((fcs - fcso)/fev.kev));

        structs = {emaxLv,emaxRv,rSp,rEp,vUsv,vUev};
        names = {'emaxLv','emaxRv','rSp','rEp','vUsv','vUev'};
        for iSig = 1:length(structs)
            fesDelayed = getydelayed(names{iSig},t,fes1,structs{iSig}.D);
            if fes1 >= fes.Min                                       %Eqn 35
                sigma.(names{iSig}) = structs{iSig}.G* ...
                                        log(fesDelayed - fes.Min + 1);
            else
                sigma.(names{iSig}) = 0;
            end

        end

        % Eqn 36
        dDeltaEmaxLv = (-deltaEmaxLv + sigma.emaxLv)/emaxLv.tau;
        dDeltaEmaxRv = (-deltaEmaxRv + sigma.emaxRv)/emaxRv.tau;
        dDeltaRSp = (-deltaRSp + sigma.rSp)/rSp.tau;
        dDeltaREp = (-deltaREp + sigma.rEp)/rEp.tau;
        dDeltaVUsv = (-deltaVUsv + sigma.vUsv)/vUsv.tau;
        dDeltaVUev = (-deltaVUev + sigma.vUev)/vUev.tau; % Global so it can
                                                         % be used in
                                                         % cvfunctions

        fesDelayedTs = getydelayed('sT',t,fes1,sT.D);
        if fes1 >= fes.Min                                           %Eqn 38
            sigmaTs = sT.G*log(fesDelayedTs - fes.Min + 1);
        else
            sigmaTs = 0;
        end

        dDeltaTs = (-deltaTs + sigmaTs)/sT.tau;                      %Eqn 39

        fevDelayedTv = getydelayed('pT',t,fev1,pT.D);
        sigmaTv = pT.G*fevDelayedTv;                                 %Eqn 40
        dDeltaTv = (-deltaTv + sigmaTv)/pT.tau;                      %Eqn 41

        out = [dPtilda, dDeltaEmaxLv, dDeltaEmaxRv, dDeltaRSp, ...
               dDeltaREp, dDeltaVUsv, dDeltaVUev, dDeltaTs, dDeltaTv]';
    end

global tOld yOld
tOld.emaxLv = 0; tOld.emaxRv = 0; tOld.rSp = 0; tOld.rEp = 0;
tOld.vUsv = 0; tOld.vUev = 0; tOld.sT = 0; tOld.pT = 0;

yOld.emaxLv = fes.o; yOld.emaxRv = fes.o; yOld.rSp = fes.o;
yOld.rEp = fes.o; yOld.vUsv = fes.o; yOld.vUev = fes.o;
yOld.sT = fes.o; yOld.pT = fev.o;

    function yDelayed = getydelayed(name,t,y,delay)
        % Gets value of y from D seconds before current time, t.

        tOld.(name) = [tOld.(name),t];
        yOld.(name) = [yOld.(name),y];
        tDiff = t - tOld.(name);

        % Remove values for tDiff > D.
        tOld.(name) = tOld.(name)(tDiff < delay);
        yOld.(name) = yOld.(name)(tDiff < delay);

        % yOld(1) should be y at approximately t - D assuming small step size.
        yDelayed = yOld.(name)(1);
    end

    function runrespmodel(tspan)
        ic = createrespic();
        [tout,yout] = ode45(@(t,y)respmodel(t,y),tspan,ic);
    end

    function ic = createrespic()
        % Creates initial condition for respiratory system.
        ic = zeros(10,1);
        ic(6) = 0;    % Pl
        ic(7) = 0;    % Ptr
        ic(8) = 0;    % Pb
        ic(9) = 0;    % PA
        ic(10) = -5;  % Ppl
    end

    function out = respmodel(t,ic)
        out = lungmechanicsfunctions(t,ic);
    end

tOld.chemC = 0; tOld.ChemP = 0;
yOld.chemC = chemC.pCO2; yOld.chemP = chemP.f;

global PAco2
Pao = 0;

    function out = lungmechanicsfunctions(t,ic)
        % All equations are from "An integrated mathmatical model of the
        % human cardiopulmanary system" paper.

        Pmusminc = ic(1); RRc = ic(2); Pmusminp = ic(3); RRp = ic(4);
        Pl = ic(6); Ptr = ic(7); Pb = ic(8); PA = ic(9);

        if length(ic) > 1000 % For respiratory control. Not implemented.
            delayedPCo2 = getydelayed('chemC',t,PACo2,chemC.D);
            delayedFacp = getydelayed('chemP',t,facp,chemP.D);
            uc = delayedPCo2 - chemC.pCO2;
            up = delayedFacp - chemP.f;

            Pmusmin = vent.Pmus0 + Pmusminc + Pmusminp;              %Eqn A79
            RR = vent.RR0 + RRc + RRp;                               %Eqn A80

            dPmusminc = (-Pmusminc + chemC.GA*uc)/chemC.tauA;        %Eqn A81
            dRRc = (-RRc + chemC.Gf*uc)/chemC.tauF;                  %Eqn A82

            dPmusminp = (-Pmusminp + chemP.GA*up)/chemP.tauA;        %Eqn A8
            dRRp = (-RRp + chemP.Gf*up)/chemP.tauF;                  %Eqn A84
        else
            RR = vent.RR0;
            Pmusmin = vent.Pmus0;
            dPmusminc = 0;
            dPmusminp = 0;
            dRRc = 0;
            dRRp = 0;
        end

        TE = (60/RR)/(1+vent.IEratio);                               %Eqn 5
        TI = TE*vent.IEratio;
        T = TI + TE;
        vent.tau = TE/5;

        tCyclic = mod(t,T);
        if tCyclic < TI
            dPmus = -2*tCyclic*Pmusmin/(TI*TE) + Pmusmin*T/(TI*TE);
        else
            dPmus = -Pmusmin*exp(-(tCyclic-TI)/vent.tau)/((1-exp(-TE/vent.tau))*vent.tau);
        end

        dPl = ((Pao - Pl)/vent.mlR - (Pl - Ptr)/vent.ltR)/vent.lC;          %Eqn A30
        dPpl = ((Pl - Ptr)/vent.ltR)/vent.cwC + dPmus;                      %Eqn A34
        dPtr = ((Pl - Ptr)/vent.ltR - (Ptr - Pb)/vent.tbR)/vent.tC + dPpl;  %Eqn A31
        dPb = ((Ptr - Pb)/vent.tbR - (Pb - PA)/vent.bAR)/vent.bC + dPpl;    %Eqn A32
        dPA = ((Pb - PA)/vent.bAR)/vent.AC + dPpl;                          %Eqn A33

        out = [dPmusminc,dRRc,dPmusminp,dRRp,dPmus,dPl,dPtr,...
               dPb,dPA,dPpl]';
    end

    function runfullmodel(tspan)
        ic = createfullic();
        [tout,yout] = ode23tb(@(t,y)fullmodel(t,y),tspan,ic);
    end

    function ic = createfullic()
        cvIc = createcvic();                    % length = 14
        controlIc = createcvcontrolic();        % length = 9
        respIc = createrespic();                % length = 10
        geIc = creategasexchangeic();           % length = 14
        
        ic = [cvIc;controlIc;respIc;geIc];
    end

    function out = fullmodel(t,ic)
        cvOut = cvcontrolmodel(t,ic(1:23));
        geOut = gasexchangefunctions(t,ic);
        respOut = respmodel(t,ic(24:end));
        
        out = [cvOut;respOut;geOut];
    end

Patm = 760;
    function ic = creategasexchangeic()
        ic = zeros(14,1);
        ic(1) = 120/(Patm - lge.Pws);     % FDo2
        ic(2) = 40/(Patm - lge.Pws);      % FDco2
        ic(3) = 80/(Patm - lge.Pws);     % FAo2
        ic(4) = 40/(Patm - lge.Pws);      % FAco2
        ic(5) = 0.03;                       % Cepo2
        ic(6) = 0.06;                       % Cepco2
        ic(7) = 0.03;                       % Cspo2
        ic(8) = 0.06;                       % Cspco2
        ic(9) = 0.02;                       % Cevo2
        ic(10) = 0.065;                      % Cevco2
        ic(11) = 0.02;                     % Csvo2
        ic(12) = 0.065;                      % Csvco2
        ic(13) = 0.02;                       % Cvo2
        ic(14) = 0.065;                      % Cvco2
    end

tOld.Cvo2t = 0; tOld.Cvco2t = 0; tOld.Cao2t = 0; tOld.Caco2t = 0;
yOld.Cvo2t = 0.02; yOld.Cvco2t = 0.065; yOld.Cao2t = 0.03; yOld.Caco2t = 0.06;

    function out = gasexchangefunctions(t,ic)
        % All equations are from "An integrated mathmatical model of the
        % human cardiopulmanary system" paper.

        Qpa = ic(2); Psp = ic(7); Pep = Psp; Pev = ic(8);
        Pra = ic(13); Pl = ic(29); Ptr = ic(30); Pb = ic(31);
        PA = ic(32); Ppl = ic(33);

        FDo2 = ic(34); FDco2 = ic(35); FAo2 = ic(36); FAco2 = ic(37);
        Cepo2 = ic(38); Cepco2 = ic(39); Cspo2 = ic(40); Cspco2 = ic(41);
        Cevo2 = ic(42); Cevco2 = ic(43); Csvo2 = ic(44); Csvco2 = ic(45);
        Cvo2 = ic(46); Cvco2 = ic(47);

        Qep = (Pep - Pev)/ep.R;
        Qsp = (Psp - sv.P)/sp.R;
        Qev = (Pev - Pra)/ev.R;
        Qsv = (sv.P - Pra)/sv.R;

        % Lung-gas exchange.
        dV = 1000*(Pao - Pl)/vent.mlR;                                   %Eqn A35
        dVA = 1000*(Pb - PA)/vent.bAR;                                   %Eqn A36

        Vl = vent.lC*Pl*1000 + vent.lVu;                                 %Eqn A37
        Vtr = vent.tC*(Ptr - Ppl)*1000 + vent.tVu;                       %Eqn A38
        Vb = vent.bC*(Pb - Ppl)*1000 + vent.bVu;                         %Eqn A39
        VA = vent.AC*(PA - Ppl)*1000 + vent.AVu;                         %Eqn A40
        VD = Vl + Vtr + Vb;                                              %Eqn A41

        [Cao2,Caco2] = getatrialconcentrations(ic);
        Cppo2 = Cao2;                                                    %Eqn A54 minus shunt.
        Cppco2 = Caco2;                                                  %Eqn A55 minus shunt.

        Cvo2t = getydelayed('Cvo2t',t,Cvo2,tge.tauVL);
        Cvco2t = getydelayed('Cvco2t',t,Cvco2,tge.tauVL);

        dFDo2 = (heaviside(dV)*dV*(lge.FIO2 - FDo2) ...
                 + heaviside(-dV)*dVA*(FDo2 - FAo2))/VD;             %Eqn A42

        dFDco2 = (heaviside(dV)*dV*(lge.FICO2 - FDco2) ...
                  + heaviside(-dV)*dVA*(FDco2 - FAco2))/VD;          %Eqn A43

        dFAo2 = (heaviside(dV)*dVA*(FDo2 - FAo2) ...                 %Eqn A44
                 - lge.K*(Qpa*(Cppo2 - Cvo2t)))/VA;

        dFAco2 = (heaviside(dV)*dVA*(FDco2 - FAco2) ...              %Eqn A45
                  - lge.K*(Qpa*(Cppco2 - Cvco2t)))/VA;

        % tissue-gas exchange.
        Vep = ep.C*Pep + ep.Vu;
        Vsp = sp.C*Psp + sp.Vu;

        Cao2t = getydelayed('Cao2t',t,Cao2,tge.tauLT);
        Caco2t = getydelayed('Caco2t',t,Caco2,tge.tauLT);

        MV = 1/22.4;         % Reciprical of molar volume (mL/mmol).
        dCepo2 = (Qep*(Cao2t - Cepo2) - MV*tge.eO2/60)/(tge.eV + Vep);           %Eqn A63
        dCepco2 = (Qep*(Caco2t - Cepco2) + MV*tge.eCO2/60)/(tge.eV + Vep);       %Eqn A64
        dCspo2 = (Qsp*(Cao2t - Cspo2) - MV*tge.sO2/60)/(tge.sV + Vsp);           %Eqn A65
        dCspco2 = (Qsp*(Caco2t - Cspco2) + MV*tge.sCO2/60)/(tge.sV + Vsp);       %Eqn A66

        % venous-pool exchange.
        Vev = ev.C*Pev + ev.Vu;
        Vsv = sv.C*sv.P + sv.Vu;

        dCevo2 = Qep*(Cepo2 - Cevo2)/Vev;                                %Eqn A73
        dCevco2 = Qep*(Cepco2 - Cevco2)/Vev;                             %Eqn A74
        dCsvo2 = Qsp*(Cspo2 - Csvo2)/Vsv;                                %Eqn A75
        dCsvco2 = Qsp*(Cspco2 - Csvco2)/Vsv;                             %Eqn A76

        dCvo2 = (Qev*(Cevo2 - Cvo2) + Qsv*(Csvo2 - Cvo2))/(Vev + Vsv);           %Eqn A77
        dCvco2 = (Qev*(Cevco2 - Cvco2) + Qsv*(Csvco2 - Cvco2))/(Vev + Vsv);      %Eqn A78

        out = [dFDo2,dFDco2,dFAo2,dFAco2,dCepo2,dCepco2,...
               dCspo2,dCspco2,dCevo2,dCevco2,dCsvo2,dCsvco2,...
               dCvo2,dCvco2]';
    end

    function [Cao2,Caco2] = getatrialconcentrations(ic)
        FAo2 = ic(36); FAco2 = ic(37);
        
        PAo2 = FAo2*(Patm - lge.Pws);                                    %Eqn A52
        PAco2 = FAco2*(Patm - lge.Pws);                                  %Eqn A53
        Pppo2 = PAo2;                                                    %Eqn A50
        Pppco2 = PAco2;                                                  %Eqn A51

        Xppo2 = Pppo2*((1+lge.beta1*Pppo2)/((lge.K1)*(1+lge.alpha1*Pppo2)));         %Eqn 47
        Cppo2 = (lge.CsatO2/1000)*((Xppo2)^(1/lge.h1))/(1 + (Xppo2)^(1/lge.h1));     %Eqn 46
        Cao2 = Cppo2;                                                                %Eqn A54 minus shunt.
        Xppco2 = Pppco2*((1 + lge.beta2*Pppco2)/(lge.K2*(1 + lge.alpha2*Pppco2)));   %Eqn 49
        Cppco2 = (lge.CsatCO2/1000)*((Xppco2)^(1/lge.h2))/(1 + (Xppco2)^(1/lge.h2)); %Eqn 48
        Caco2 = Cppco2;                                                              %Eqn A55 minus shunt.
    end

    function y = heaviside(x)
        if x > 0
            y = 1;
        else
            y = 0;
        end
    end


%% show plots
    function showplots(varargin)
        switch modelNum;
            case 2
                cvplots();
            case 3
                cvcontrolplots();
            case 4
                respplots(yout);
            case 5
                fullplots();
            otherwise
                return
        end

        buttonLabels = {'Run CV model','Run CV model (w/ control)',...
                        'Run resp model','Run full model'};
        mainFigure.btn{modelNum}.String = buttonLabels{modelNum - 1};
    end

    function cvplots()
        Vlv = yout(:,11); Ppa = yout(:,1); Ppp = yout(:,3); 
        Pla = yout(:,10); Pra = yout(:,13); Psa = yout(:,5);
        
        lvP = zeros(length(tout),1);
        for i = 1:length(tout)
            [~,~,~,~,lvP(i)] = ventriclefunctions(yout(i,:));
        end
        
        figure('Name','CV')
        plot(Vlv,lvP)
        title('PV plot for left ventricle')
        ylabel('Pressure (mmHg)')
        xlabel('Volume (mL)')
        
        figure('Name','CV')
        subplot(321)
        plot(tout,Ppa)
        ylabel('P_{pa} (mmHg)')
        
        subplot(322)
        plot(tout,Ppp)
        ylabel('P_{pp} (mmHg)')
        
        subplot(323)
        plot(tout,Pla)
        ylabel('P_{la} (mmHg)')
        
        subplot(324)
        plot(tout,Pra)
        ylabel('P_{ra} (mmHg)')
        
        subplot(325)
        plot(tout,Vlv)
        ylabel('V_{lv} (mL)')
        xlabel('time (s)')
        
        subplot(326)
        plot(tout,Psa)
        ylabel('P_{sa} (mmHg)')
        xlabel('time (s)')
    end

    function cvcontrolplots()
        cvplots()

        figure('Name','CV Control')
        dTs = yout(:,22); dTv = yout(:,23);
        T = dTs + dTv + hr.T;
        plot(tout,[dTs,dTv,T])
        legend('{\Delta}Ts','{\Delta}Tv','Period')
        ylabel('Period (s)')
        xlabel('time (s)')
        title('Change in period due to baroreceptors')
    end

    function respplots(y)

        Pl = y(:,6); Ptr = y(:,7); Pb = y(:,8);
        PA = y(:,9); Ppl = y(:,10); Pmus = y(:,5);

        Vl = vent.lC.*Pl + vent.lVu/1000;                               %Eqn A37
        Vtr = vent.tC.*(Ptr - Ppl) + vent.tVu/1000;                     %Eqn A38
        Vb = vent.bC.*(Pb - Ppl) + vent.bVu/1000;                       %Eqn A39
        VA = vent.AC.*(PA - Ppl) + vent.AVu/1000;                       %Eqn A40
        VD = Vl + Vtr + Vb;                                             %Eqn A41
        VL = VD + VA;
        dV = (Pao - Pl)/vent.mlR;                                       %Eqn A35

        figure('Name','Respiratory')
        subplot(421)
        plot(tout,Pmus)
        ylabel('P_{mus} (cmH_2O)')

        subplot(423)
        plot(tout,Ppl)
        ylabel('P_{pl} (cmH_2O)')

        subplot(425)
        plot(tout,PA)
        ylabel('P_A (cmH_2O)')

        subplot(427)
        plot(tout,dV)
        ylabel('Airflow (L/s)')
        xlabel('time (s)')

        subplot(422)
        plot(tout,VL)
        ylabel('V_L (L)')

        subplot(424)
        plot(tout,VA)
        ylabel('V_A (L)')

        subplot(426)
        plot(tout,VD)
        ylabel('V_D (L)')
        xlabel('time (s)')

    end

    function gasexchangeplots(y)
        FDo2 = y(:,1); FDco2 = y(:,2); FAo2 = y(:,3); FAco2 = y(:,4);
        
        figure('Name','Gas Exchange')
        
        subplot(211)
        plot(tout,[FDo2,FAo2]*(Patm - lge.Pws))
        title('Partial pressure of O_2 and CO_2')
        ylabel('Partial pressure O_2 (mmHg)')
        legend('P_{DO_2}','P_{AO_2}')
        
        subplot(212)
        plot(tout,[FDco2,FAco2]*(Patm - lge.Pws))
        ylabel('Partial pressure CO_2 (mmHg)')
        xlabel('time (s)')
        legend('P_{DCO_2}','P_{ACO_2}')
    end

    function fullplots()
        cvcontrolplots()
        respplots(yout(:,24:33))
        gasexchangeplots(yout(:,34:end))
    end

end % main function end
