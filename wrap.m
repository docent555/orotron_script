IN.Ne = 16;
IN.Lz = 5;
IN.Tend = 1;
IN.dz = 0.0002;
IN.dt = 0.0005;


RES = orotron(16, 5, 1, 0, 0.0002, 0.0005);

plot(RES.ZAxis(:,1), abs(RES.OUTField(:,end)))

save('results.mat',"RES","-v7.3");