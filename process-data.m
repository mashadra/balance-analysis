%% paramters
COM_infile = "SRS/1.tsv"
COP_infile = "SRS/1_f_2.tsv";
COM_outfile = "SRS/trial1_COM.txt";
COP_outfile = "SRS/trial1_COP.txt";
COM_COP_diff_outfile = "SRS/trial1_DIFF.txt";
%% functions
function zerocrossings = getzerocross(rawdata)
   deriv2 = diff(rawdata, 2);
   zerocross = 0;
   for i = 2:length(deriv2)
       if deriv2(i) > 0 && deriv2(i-1) < 0
           zerocross = zerocross + 1;
       elseif deriv2(i) < 0 && deriv2(i-1) > 0
           zerocross = zerocross + 1;
       end
   end
   zerocrossings = zerocross;
end

function data = interpolation(rawdata, upsamplesize)
    data = [rawdata(1)];
    for i = 2:length(rawdata)
        p1 = rawdata(i-1);
        p2 = rawdata(i);

        slope = p2-p1;
        xin = p2-i*slope;

        for j = 1:upsamplesize
            data = [data (slope*((i-1)+(double(j))/upsamplesize)+xin)];
        end
    end
end
%% import COM data
t = readtable(COM_infile, "FileType","text",'Delimiter', '\t');
comf = table2array(t(1:end, "Frame"));
comx = table2array(t(1:end,"COMX"));
comy = table2array(t(1:end,"COMY"));
comz = table2array(t(1:end,"COMZ"));
%% COM math
% x com
meanx = mean(comx);
stdevx = std(comx);
rangex = range(comx);
iqrx = iqr(comx);
zerox = getzerocross(comx);
percentx = (double(zerox)/double(length(comf)))*100;
% y com
meany = mean(comy);
stdevy = std(comy);
rangey = range(comy);
iqry = iqr(comy);
zeroy = getzerocross(comy);
percenty = (double(zeroy)/double(length(comf)))*100;
% z com
meanz = mean(comz);
stdevz = std(comz);
rangez = range(comz);
iqrz = iqr(comz);
zeroz = getzerocross(comz);
percentz = (double(zeroz)/double(length(comf)))*100;
%% COM graphs
% 1D transient
figure(1)
subplot(3, 2, 1)
plot(comf, comx)
title('X Position vs Frame for COM')
xlabel('Frame')
ylabel('X (mm)')
subplot(3, 2, 3)
plot(comf, comy)
title('Y Position vs Frame for COM')
xlabel('Frame')
ylabel('Y (mm)')
subplot(3, 2, 5)
plot(comf, comz)
title('Z Position vs Frame for COM')
xlabel('Frame')
ylabel('Z (mm)')
subplot(3, 2, 2)
plot(comx, comy)

% 2D
title('X Position vs Y Position for COM')
xlabel('X (mm)')
ylabel('Y (mm)')
subplot(3, 2, 4)
plot(comx, comz)
title('X Position vs Z Position for COM')
xlabel('X (mm)')
ylabel('Z (mm)')
subplot(3, 2, 6)
plot(comy, comz)
title('Y Position vs Z Position for COM')
xlabel('Y (mm)')
ylabel('Z (mm)')
%% export COM data
Direction = {"x";"y";"z"};
Mean = [meanx;meany;meanz];
Stdev = [stdevx;stdevy;stdevz];
Rang = [rangex;rangey;rangez];
IQR = [iqrx;iqry;iqrz];
ZeroCross = [zerox;zeroy;zeroz];
ZeroCrossPercent = [percentx;percenty;percentz];
T2 = table(Direction,Mean,Stdev,Rang,IQR,ZeroCross,ZeroCrossPercent);
writetable(T2, COM_outfile);
%% import COP data
t = readtable(COP_infile, "FileType","text",'Delimiter', '\t');
top_label_t = ["FORCE_PLATE_CORNER_POSX_POSY_X",
"FORCE_PLATE_CORNER_POSX_POSY_Y",
"FORCE_PLATE_CORNER_POSX_POSY_Z",
"FORCE_PLATE_CORNER_NEGX_POSY_X",
"FORCE_PLATE_CORNER_NEGX_POSY_Y",
"FORCE_PLATE_CORNER_NEGX_POSY_Z",
"FORCE_PLATE_CORNER_NEGX_NEGY_X",
"FORCE_PLATE_CORNER_NEGX_NEGY_Y",
"FORCE_PLATE_CORNER_NEGX_NEGY_Z",
"FORCE_PLATE_CORNER_POSX_NEGY_X",
"FORCE_PLATE_CORNER_POSX_NEGY_Y",
"FORCE_PLATE_CORNER_POSX_NEGY_Z"
];
top_t = table2array(readtable(COP_infile, "FileType","text",'Delimiter', '\t', 'Range', 'B1:B26'));
f = table2array(t(1:end, "SAMPLE"));
x = table2array(t(1:end,"COP_X"));
y = table2array(t(1:end,"COP_Y"));
%% COP shift
xcorners = [];
x1 = 9 + find(top_label_t=="FORCE_PLATE_CORNER_POSX_POSY_X");
xcorners = [xcorners, top_t(x1)];
x2 = 9 + find(top_label_t=="FORCE_PLATE_CORNER_NEGX_POSY_X");
xcorners = [xcorners, top_t(x2)];
x3 = 9 + find(top_label_t=="FORCE_PLATE_CORNER_NEGX_NEGY_X");
xcorners = [xcorners, top_t(x3)];
x4 = 9 + find(top_label_t=="FORCE_PLATE_CORNER_POSX_NEGY_X");
xcorners = [xcorners, top_t(x4)];
xcenter = mean(xcorners)
ycorners = [];
y1 = 9 + find(top_label_t=="FORCE_PLATE_CORNER_POSX_POSY_Y");
ycorners = [ycorners, top_t(y1)];
y2 = 9 + find(top_label_t=="FFORCE_PLATE_CORNER_NEGX_POSY_Y");
ycorners = [ycorners, top_t(y2)];
y3 = 9 + find(top_label_t=="FORCE_PLATE_CORNER_NEGX_NEGY_Y");
ycorners = [ycorners, top_t(y3)];
y4 = 9 + find(top_label_t=="FORCE_PLATE_CORNER_POSX_NEGY_Y");
ycorners = [ycorners, top_t(y4)];
ycenter = mean(ycorners)
%% COP math
% x cop
meanx = mean(x);
stdevx = std(x);
rangex = range(x);
iqrx = iqr(x);
zerox = getzerocross(x);
percentx = (double(zerox)/double(length(f)))*100;
% y cop
meany = mean(y);
stdevy = std(y);
rangey = range(y);
iqry = iqr(y);
zeroy = getzerocross(y);
percenty = (double(zeroy)/double(length(f)))*100;
% fixing for offset
fixedx = x + xcenter;
fixedy = y + ycenter;
%% export COP data
Direction = {"x";"y"};
Mean = [meanx;meany];
Stdev = [stdevx;stdevy];
Rang = [rangex;rangey];
IQR = [iqrx;iqry];
ZeroCross = [zerox;zeroy];
ZeroCrossPercent = [percentx;percenty];
T = table(Direction,Mean,Stdev,Rang,IQR,ZeroCross,ZeroCrossPercent);
writetable(T, COP_outfile);
%% COP graphs

% 1D transient
figure(2)
subplot(3, 1, 1)
plot(f, x)
title('X Position vs Sample for COP')
xlabel('Sample')
ylabel('X (mm)')
subplot(3, 1, 2)
plot(f, y)
title('Y Position vs Sample for COP')
xlabel('Sample')
ylabel('Y (mm)')

% 2D
subplot(3, 1, 3)
plot(x, y)
title('Y Position vs X Position for COP')
xlabel('X (mm)')
ylabel('Y (mm)')
%% COM interpolation
interpolatedcomx = interpolation(comx, 4);
interpolatedcomy = interpolation(comy, 4);
%% COM - COP math
if length(interpolatedcomx) < length(f)
    final = length(interpolatedcomx);
else
    final = length(f);
end

% x com-cop
newx = interpolatedcomx(1:final) - fixedx(1:final)';
meanx = mean(newx);
stdevx = std(newx);
rangex = range(newx);
iqrx = iqr(newx);
% y com-cop
newy = interpolatedcomy(1:final) - fixedy(1:final)';
meany = mean(newy);
stdevy = std(newy);
rangey = range(newy);
iqry = iqr(newy);
%% export COM - COP data
Direction = {"x";"y"};
Mean = [meanx;meany];
Stdev = [stdevx;stdevy];
Rang = [rangex;rangey];
IQR = [iqrx;iqry];
T = table(Direction,Mean,Stdev,Rang,IQR);
writetable(T, COM_COP_diff_outfile);
%% COM and COP plots

% X position vs frame
figure (3)
subplot(2, 1, 1)
hold on
comxline = plot(1:final, interpolatedcomx(1:final));
copxline = plot(1:final, fixedx(1:final));
title('X Position vs Frame for Both')
xlabel('Frame')
ylabel('X (mm)')
legend([comxline copxline],["COM X", "COP X"]);
hold off
subplot(2, 1, 2)
hold on
comyline = plot(1:final, interpolatedcomy(1:final));
copyline = plot(1:final, fixedy(1:final));
title('Y Position vs Frame for Both')
xlabel('Frame')
ylabel('Y (mm)')
legend([comyline copyline],["COM Y", "COP Y"]);
hold off

% x vs y position
figure(4)
hold on
comyline = plot(interpolatedcomx(1:final), interpolatedcomy(1:final));
copyline = plot(fixedx(1:final), fixedy(1:final));
title('Y Position vs X Position for Both')
xlabel('X (mm)')
ylabel('Y (mm)')
legend([comyline copyline],["COM", "COP"]);
hold off

% acceleration plots
figure(5)

subplot(3, 2, 1)
plot(diff(comx, 2))
title('X Accelaration vs Frames for COM')
xlabel('Frame')
ylabel('X (mm/s^2)')

subplot(3, 2, 3)
plot(diff(comy, 2))
title('Y Accelaration vs Frames for COM')
xlabel('Frame')
ylabel('Y (mm/s^2)')

subplot(3, 2, 5)
plot(diff(comz, 2))
title('Z Accelaration vs Frames for COM')
xlabel('Frame')
ylabel('Z (mm/s^2)')

subplot(3, 2, 2)
plot(diff(x, 2))
title('X Accelaration vs Frames for COP')
xlabel('Frame')
ylabel('X (mm/s^2)')

subplot(3, 2, 4)
plot(diff(y, 2))
title('Y Accelaration vs Frames for COP')
xlabel('Frame')
ylabel('Y (mm/s^2)')

%% difference plots
figure(6)
subplot(2, 1, 1)
plot(newx)
title('X Difference vs Frames')
xlabel('Frame')
ylabel('X difference (mm)')

figure(6)
subplot(2, 1, 2)
plot(newy)
title('Y Difference vs Frames')
xlabel('Frame')
ylabel('Y difference (mm)')
