function OUT = oroscr(IN) %#codegen

fprintf('\nTime')

ZAxis = zeros(length(IN.ZAxis),1);
TAxis = zeros(length(IN.TAxis),1);
InitialField = complex(zeros(length(IN.InitialField),1));

ZAxis(:,1) = IN.ZAxis;
TAxis(:,1) = IN.TAxis;
INTERVALT = IN.INTERVALT;
INTERVALZ = IN.INTERVALZ;
Delta = IN.Delta;
I = IN.I;
Ne = IN.Ne;
InitialField(:,1) = IN.InitialField(:,1);

% dthdz0(:,1) = zeros(Ne,1);
% dthdz0(:,1) = Delta;

% ZBEG = 0;
% ZEND = 0.5;
% InitialField = zeros(length(ZAxis),1);
% IND1 = (ZAxis > ZBEG & ZAxis < ZEND);
% InitialField(IND1) = 0.001*sin(pi * (ZAxis(IND1) - ZBEG) / (ZEND - ZBEG)).^2;

% InitialFieldNormalized = InitialField/max(abs(InitialField));
if INTERVALZ > 1
    IZ = 0:INTERVALZ:length(ZAxis);
    IZ(1) = 1;
    SIZEZ = length(IZ);
else
    IZ = 1:INTERVALZ:length(ZAxis);
end

if INTERVALT > 1 && INTERVALZ > 1
    TAxisNew = zeros(fix(length(TAxis(:,1))/INTERVALT)+1,1);
    OUT.OUTJ = complex(zeros(SIZEZ,fix(length(TAxis(:,1))/INTERVALT)+1));
    OUT.OUTB = complex(zeros(SIZEZ,fix(length(TAxis(:,1))/INTERVALT)+1));
elseif INTERVALT ==1 && INTERVALZ > 1
    TAxisNew = zeros(fix(length(TAxis(:,1))/INTERVALT),1);
    OUT.OUTJ = complex(SIZEZ,fix(length(TAxis(:,1))/INTERVALT));
    OUT.OUTB = complex(SIZEZ,fix(length(TAxis(:,1))/INTERVALT));
elseif INTERVALT > 1 && INTERVALZ == 1
    TAxisNew = zeros(fix(length(TAxis(:,1))/INTERVALT)+1,1);
    OUT.OUTJ = complex(zeros(size(InitialField,1),fix(length(TAxis(:,1))/INTERVALT)+1));
    OUT.OUTB = complex(zeros(size(InitialField,1),fix(length(TAxis(:,1))/INTERVALT)+1));
else % (INTERVALT == 1) && (INTERVALZ == 1)
    TAxisNew = zeros(fix(length(TAxis(:,1))),1);
    OUT.OUTJ = complex(zeros(size(InitialField,1),length(TAxis(:,1))));
    OUT.OUTB = complex(zeros(size(InitialField,1),length(TAxis(:,1))));
end

Field = complex(zeros(size(InitialField,1),1));
OUT.ZAxis = zeros(length(IZ),1);
OUT.TAxis = zeros(length(TAxisNew),1);
J = complex(zeros(length(ZAxis),1));
Jp = complex(zeros(length(ZAxis),1));

% step = 1;
Field(:,1) = InitialField;
TAxisNew(1,1) = TAxis(1,1);

kpar2 = zeros(length(ZAxis),1);

DeltaZ = ZAxis(2) - ZAxis(1);
DeltaT = TAxis(2) - TAxis(1);

% SQR2 = sqrt(2);
SQR2M2 = 2.828427124746190;
SQR2D2 = 0.707106781186548;
SQRDT = sqrt(DeltaT);
SQRDZ = DeltaZ*DeltaZ;

% C0 = -1i;
C0 = 1i;
CR = 0;

OldOldSigmaNz = complex(0);
OldOldSigmaNzm1 = complex(0);

WNz = -((0.666666666666667*C0*DeltaZ/DeltaT + kpar2(end)*DeltaZ/3) - 1/DeltaZ);
WNzm1 = -((C0/3*DeltaZ/DeltaT + kpar2(end-1)*DeltaZ/6) + 1/DeltaZ);

% C2 = 1/sqrt(1i*pi);
C2 = 1/sqrt(-1i*pi);

N = length(ZAxis);
A = complex(zeros(N,1));
B = complex(zeros(N-1,1));
C = complex(zeros(N-1,1));
D = complex(zeros(N,1));

A(1) = 1;
A(2:end-1) = -2*(1 - DeltaZ/DeltaT*C0*DeltaZ - DeltaZ*kpar2(2:end-1)*DeltaZ/2);
A(end) = 1 + 1.333333333333333*C2*WNz*SQRDT;

B(1) = 0;
B(2:end) = 1;

C(1:end-1) = 1;
C(end) = 1.333333333333333*C2*WNzm1*SQRDT;

WR = complex(zeros(length(TAxis),1));
OldFNz = complex(zeros(length(TAxis)+1,1));
OldFNzm1 = complex(zeros(length(TAxis)+1,1));
% OldJNz = complex(zeros(length(TAxis)+1,1));
% OldJNzm1 = complex(zeros(length(TAxis)+1,1));


M = spdiags([[C; 0] A [0 ;B]], -1:1, N, N);


% hF= hField(length(ZAxis), 1);

IDX = @(j) (j + 2);

fmax = zeros(length(TAxis), 1);

% SHOW = 1;
% if SHOW == 1
%     [lhr, lha, hFig] = makeFig(ZAxis, TAxis);
% end

jout = 1;
OUT.OUTB(:, jout) = Field(IZ,1);

th0 = 2*pi*(0:Ne-1)/Ne;
theta = zeros(length(ZAxis), Ne);

steps = length(TAxis) - 1;

field = complex(zeros(length(Field),1));
tmp = complex(zeros(length(Field),1));
field_p = complex(zeros(length(Field),1));
rfield_p = complex(zeros(length(Field),1));
lfield_p = complex(zeros(length(Field),1));
NewfieldRaw = complex(zeros(length(Field),1));
rNewfieldRaw = complex(zeros(length(Field),1));
lNewfieldRaw = complex(zeros(length(Field),1));
IROldPart = complex(zeros(1,1));
OldOldJ = complex(zeros(length(Field),1));

timerVal = tic;

fprintf('\n');
% tmp = Field;
for step=1:steps
    
    field(:,1) = Field(:,1);
%     field = tmp;
    
%     if SHOW == 1
%         lHandleB.YData = abs(Field(:,1));
%         lHandleJ.YData = abs(J(:,1));
%         lhr.YData(1:step) = fmax(1:step);
%         lhr.XData(1:step) = TAxis(1:step);
%         lha.YData = abs(Field);
%         drawnow
%     end


%     J = zeros(length(ZAxis),1);
theta(:,:) = pendulumODE(field(:,1), ZAxis(:,1), th0(1,:), Delta);
J(:,1) = I * 2 / Ne * sum(exp(-1i*theta), 2);

if (step ~= 1) && (mod(step-1,INTERVALT) == 0)
    OUT.OUTJ(:,jout) = J(IZ,1);
    end
    
    OldFNz(IDX(step - 1)) =  field(end);
    OldFNzm1(IDX(step - 1)) = field(end - 1);
    
    OldSigmaNz = -(kpar2(end)/6 + C0/3/DeltaT) * OldFNz(IDX(step - 1)) ...
        + (C0/3/DeltaT - kpar2(end)/6) * OldFNz(IDX(step - 2)) ...
        + 0.166666666666667*(J(end) + OldOldJ(end)) - OldOldSigmaNz;
    OldSigmaNzm1 = -(kpar2(end - 1)/6 + C0/3/DeltaT) * OldFNzm1(IDX(step - 1)) ...
        + (C0/3/DeltaT - kpar2(end - 1)/6) * OldFNzm1(IDX(step - 2)) ...
        + 0.166666666666667*(J(end-1) + OldOldJ(end-1)) - OldOldSigmaNzm1;
    
    OldOldSigmaNz = OldSigmaNz;
    OldOldSigmaNzm1 = OldSigmaNzm1;
    
    WR_step_part = DeltaZ * ((C0 * 0.666666666666667 / DeltaT - kpar2(end) / 3) * field(end) ...
        + (C0 / 3 / DeltaT - kpar2(end - 1) / 6) * field(end - 1)...
        - (2 * OldSigmaNz + OldSigmaNzm1));
    WR(IDX(step)) = WR_step_part + DeltaZ * (0.166666666666667*(6*J(end) - 2*OldOldJ(end) + 3*J(end-1) - OldOldJ(end-1)));
    
    
    u = @(j) (WNzm1 * OldFNzm1(IDX(j)) + WNz * OldFNz(IDX(j)) + WR(IDX(j))).' .* exp(CR * DeltaT * (step - j));
    
    if step == 1
        IR = 0;
    elseif step == 2
        IR = 1.333333333333333 * SQRDT * (u(0)*(1 - SQR2D2) + u(1)*(SQR2M2 - 2.5));
    else
        j = 1:step-2;
        IR = 1.333333333333333 * SQRDT * (u(0)*((step - 1).^(1.5) - (step - 1.5)*sqrt(step))...
            + sum(u(j).*((step - j - 1).^(1.5) - 2*(step - j).^(1.5) + (step - j + 1).^(1.5)))...
            + u(step - 1)*(SQR2M2 - 2.5));
    end
    
    D(1) = 0;
    D_2_Nzm1_part = 2 * (1 + C0 * SQRDZ/DeltaT - SQRDZ * kpar2(2:end-1)/2) .* field(2:end - 1)...
        - (field(1:end - 2) + field(3:end));
    D(2:end - 1) = D_2_Nzm1_part + SQRDZ * (2*J(2:end - 1) - OldOldJ(2:end - 1));
    D_end_part = -C2*(IR + 0.666666666666667 * DeltaT / SQRDT * u(step - 1));
    D(end) = D_end_part - C2*(1.333333333333333 * WR(IDX(step)) * SQRDT);
    
%     if step == 3
%         tmp = SQRDZ * (2*J(2:end - 1) - OldOldJ(2:end - 1));
%         tmp=[real(tmp) imag(tmp)];
%         
%         WR(IDX(step))
%         %     tmp = ZAxis(:);
%         %     for i=1:Ne
%         %         tmp = [tmp theta(:,i)];
%         %     end
%         % tmp = th0(1,:)'
%         save('test.dat','tmp','-ascii')
%         pause
%     end    
    
    %     field_p = M \ D;
    rfield_p = rtridag(C,A,B,D);
    lfield_p = ltridag(C,A,B,D);
    field_p = (rfield_p + lfield_p)/2;
    
    
    %     Jp = zeros(length(ZAxis),1);
    theta(:,:) = pendulumODE(field_p(:,1), ZAxis(:,1), th0(1,:), Delta);
    Jp(:,1) = I * 2 / Ne * sum(exp(-1i*theta), 2);
    
    WR(IDX(step)) = WR_step_part + ...
        DeltaZ * (0.166666666666667 * (2 * Jp(end) + 2 * J(end) + Jp(end - 1) + J(end - 1)));
    
    D(2:end - 1) = D_2_Nzm1_part + SQRDZ * (Jp(2:end - 1) + J(2:end - 1));
    D(end) = D_end_part - C2*(1.333333333333333 * WR(IDX(step)) * SQRDT);
    
    
    %NewfieldRaw(:,1) = M \ D;
    rNewfieldRaw(:,1) = rtridag(C,A,B,D);
    lNewfieldRaw(:,1) = rtridag(C,A,B,D);
    NewfieldRaw = (rNewfieldRaw + lNewfieldRaw)/2;
    
    Field(:,1) = NewfieldRaw(:,1);
    fmax(step+1) = max(abs(Field(:,1)));
       
    OldOldJ(:,1) = J(:,1);
   
    k = step + 1;
    
    if mod(step,INTERVALT) == 0
        jout = jout + 1;
        OUT.OUTB(:, jout) = NewfieldRaw(IZ,1);
        TAxisNew(jout,1) = TAxis(k,1);
    end
    
    %     if app.Stop == 1
    %         OUT.ZAxis(:,1) = ZAxis(:,1);
    %         OUT.TAxis(:,1) = TAxisNew(:,1);
    %         app.Stop = jout;
    %         return
    %     end
    %     title(app.UIAxes, sprintf('B time=%8.5f',TAxis(k)));
    %     app.Time0Label.Text = sprintf('Time =%8.4f',TAxis(k));
%     fprintf(['\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...        
%         '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
%         '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
%         'Step = %8i   Time = %8.4f   Bmax = %15.10f   Jmax = %15.10f'], step, TAxis(k), fmax(k), max(abs(J(:,1))));

fprintf(['\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
         '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
         '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
         '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'...
         'Step = %8i   Time = %8.4f   Bmax = %15.10f   Jmax = %15.10f'],...
         step, TAxis(k), fmax(k), max(abs(J(:,1))));
end

OUT.OUTJ(:,jout) = J(IZ,1); % для последнего jout см. стр.212

fprintf("\n\n\n");

ExecutionTime = toc(timerVal);

hours = fix(ExecutionTime/3600);
minutes = fix((ExecutionTime - fix(hours*3600))/60);
seconds = ExecutionTime - hours*3600 - minutes*60;

fprintf("ExecitionTime = %8.4f [h]   %8.4f [m]   %8.4f [s]\n", hours, minutes, seconds);

OUT.ZAxis(:,1) = ZAxis(IZ,1);
OUT.TAxis(:,1) = TAxisNew(:,1);

fprintf(" \n\n");

% if SHOW == 1
%     close(hFig);
% end

end

