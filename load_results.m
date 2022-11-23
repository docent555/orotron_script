BRe = load('resultsBRe.dat');
BIm = load('resultsBIm.dat');
JRe = load('resultsJRe.dat');
JIm = load('resultsJIm.dat');
s.z = load('Z.txt');
s.t = load('Time.txt');

s.B=complex(BRe, BIm);
s.J=complex(JRe, JIm);

clear BRe BIm JRe JIm