clc;


%% Variables
Ti = 32; % deg celsius
Tinf = 42; % deg celsius
m = 22*10^-3; % mass of the copper in kg
c = 385; %; %W/m-K ??? should be J/kg-K 385 -> copper; 1007 -> air
h_bar = 4; % guess, W/m-K turns out 4 is quite accurate instead of 15
%So basically: 4 is what I found to be a reasonable value
% 4 must be a good estimate because the flow rate is so small!
As = 0.791*0.838*0.0254^2*2; % surface area of the busbar; two times the Surface Area

%%

qbat = 0.4;
Tss = qbat/(h_bar*As)+Tinf;
tt = 2*m*c/(h_bar*As);

AccelDeccelValues = [5.5, 5.5, 5.5, 5.5, 5.5, 0, 0, 0, 0, 0, 0, 0, 4.5, 4.5, 4.5, 0, 0, 0, 0, 0, 5.5, 5.5, 5.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.5, 4.5, 4.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
vector = AccelDeccelValues;
i = 0;
j = 0;
k = 0;
Temp = [];
tempcounter = 0;
otimecounter = 0;
timevector = [];
laps = 0;
while laps <= 30
    while i < length(vector)
        if i > 0
                Ti_new = Temp(end);
        end
        
        if vector(i+1) > 0
            t_on_accel = abs(vector(i+1));
            timecounter = 0;
            while j < t_on_accel
                if i == 0 && laps == 0
                    Temp(tempcounter+1) = (Ti - Tss)*exp(-timecounter/tt)+Tss;
                else
                    Temp(tempcounter+1) = (Ti_new - Tss)*exp((-timecounter)/tt)+Tss;
                end
                timevector(tempcounter+1) = otimecounter;
                tempcounter = tempcounter + 1;
                otimecounter = otimecounter + 1;
                timecounter = timecounter + 1;
                j = j + 1;
            end
            j = 0;

        else
            t_on_deccel = abs(vector(i+1));
            timecounter = 0;
            while k < t_on_deccel
                if i == 0 && laps == 0
                    Temp(tempcounter+1) = (Ti - Tinf)*exp(-(h_bar*As)/(m*c)*timecounter) + Tinf;
                else
                    Temp(tempcounter+1) = (Ti_new-Tinf)*exp(-(h_bar*As)/(m*c)*(timecounter)) + Tinf;
                end
                timevector(tempcounter+1) = otimecounter;
                tempcounter = tempcounter + 1;
                otimecounter = otimecounter + 1;
                timecounter = timecounter + 1;
                k = k + 1;
            end
            k = 0;

        end

        i = i + 1;
    end
    i = 0;
    laps = laps + 1;
end

time = timevector/60;
figure(1)
plot(time, Temp);
xlabel('Time (min)');
ylabel('Temperature (C)');
