OPENQASM 2.0;
include "qelib1.inc";
qreg q[20];
creg c[12];
u3(0,0,pi/2) q[0];
sx q[0];
u3(0,0,pi/2) q[0];
u3(0,0,pi/2) q[1];
sx q[1];
u3(0,0,pi/2) q[1];
u3(0,0,pi/2) q[2];
sx q[2];
u3(0,0,pi/2) q[2];
u3(0,0,pi/2) q[3];
sx q[3];
u3(0,0,pi/2) q[3];
u3(0,0,pi/2) q[4];
sx q[4];
u3(0,0,pi/2) q[4];
cx q[3],q[4];
cx q[4],q[3];
cx q[3],q[4];
u3(0,0,pi/2) q[5];
sx q[5];
u3(0,0,pi/2) q[5];
u3(0,0,pi/2) q[6];
sx q[6];
u3(0,0,pi/2) q[6];
u3(0,0,pi/2) q[7];
sx q[7];
u3(0,0,pi/2) q[7];
u3(0,0,pi/2) q[8];
sx q[8];
u3(0,0,pi/2) q[8];
u3(0,0,pi/2) q[9];
sx q[9];
u3(0,0,pi/2) q[9];
u3(0,0,pi/2) q[10];
sx q[10];
u3(0,0,pi/2) q[10];
u3(0,0,pi/2) q[11];
sx q[11];
u3(0,0,pi/2) q[11];
cx q[0],q[11];
cx q[11],q[0];
cx q[0],q[11];
cx q[11],q[2];
u3(pi/2,pi/2,3*pi/2) q[11];
u3(0,0,-1.0*g1) q[2];
u3(pi/2,-pi,pi/2) q[2];
cx q[2],q[11];
u3(pi/2,3*pi/2,0) q[11];
u3(pi/2,0,3*pi/2) q[2];
cx q[2],q[11];
u3(pi/2,pi,-pi/2) q[11];
cx q[11],q[0];
cx q[0],q[11];
cx q[11],q[0];
cx q[0],q[1];
u3(0,0,-1.0*g1) q[1];
cx q[0],q[1];
u3(pi/2,0,pi) q[2];
cx q[2],q[4];
u3(0,0,-1.0*g1) q[4];
cx q[2],q[4];
cx q[4],q[3];
cx q[3],q[4];
cx q[4],q[3];
cx q[2],q[4];
u3(0,0,-1.0*g1) q[4];
cx q[2],q[4];
cx q[2],q[9];
cx q[9],q[2];
cx q[2],q[9];
cx q[2],q[11];
cx q[11],q[2];
cx q[2],q[11];
cx q[9],q[5];
u3(0,0,-1.0*g1) q[5];
cx q[9],q[5];
cx q[6],q[14];
cx q[14],q[6];
cx q[6],q[14];
cx q[3],q[15];
cx q[15],q[3];
cx q[3],q[15];
cx q[15],q[13];
cx q[13],q[15];
cx q[15],q[13];
cx q[0],q[13];
u3(0,0,-1.0*g1) q[13];
cx q[0],q[13];
cx q[4],q[3];
cx q[3],q[4];
cx q[4],q[3];
cx q[9],q[16];
cx q[16],q[9];
cx q[9],q[16];
cx q[16],q[14];
u3(0,0,-1.0*g1) q[14];
cx q[16],q[14];
cx q[16],q[9];
cx q[9],q[16];
cx q[16],q[9];
cx q[9],q[7];
u3(0,0,-1.0*g1) q[7];
cx q[9],q[7];
cx q[5],q[18];
cx q[18],q[5];
cx q[5],q[18];
cx q[0],q[18];
u3(0,0,-1.0*g1) q[18];
cx q[0],q[18];
cx q[1],q[0];
cx q[0],q[1];
cx q[1],q[0];
cx q[1],q[7];
cx q[13],q[0];
u3(0,0,-1.0*g1) q[0];
cx q[13],q[0];
cx q[13],q[15];
cx q[15],q[13];
cx q[13],q[15];
cx q[0],q[13];
cx q[13],q[0];
cx q[0],q[13];
cx q[13],q[6];
cx q[15],q[3];
cx q[18],q[8];
u3(0,0,-1.0*g1) q[3];
cx q[15],q[3];
cx q[3],q[14];
cx q[14],q[3];
cx q[3],q[14];
cx q[6],q[13];
cx q[13],q[6];
cx q[14],q[6];
u3(0,0,-1.0*g1) q[6];
cx q[14],q[6];
cx q[7],q[1];
cx q[1],q[7];
cx q[7],q[3];
u3(0,0,-1.0*g1) q[3];
cx q[7],q[3];
cx q[15],q[3];
u3(0,0,-1.0*g1) q[3];
cx q[15],q[3];
cx q[14],q[3];
u3(0,0,-1.0*g1) q[3];
cx q[14],q[3];
cx q[7],q[1];
u3(0,0,-1.0*g1) q[1];
u3(pi/2,-pi,pi/2) q[1];
u3(pi/2,pi/2,3*pi/2) q[7];
cx q[1],q[7];
u3(pi/2,0,3*pi/2) q[1];
u3(pi/2,3*pi/2,0) q[7];
cx q[1],q[7];
u3(pi/2,0,pi) q[1];
u3(pi/2,pi,-pi/2) q[7];
cx q[8],q[18];
cx q[18],q[8];
cx q[8],q[6];
u3(0,0,-1.0*g1) q[6];
cx q[8],q[6];
cx q[6],q[15];
cx q[15],q[6];
cx q[6],q[15];
cx q[3],q[15];
u3(0,0,-1.0*g1) q[15];
cx q[3],q[15];
cx q[14],q[3];
cx q[3],q[14];
cx q[14],q[3];
cx q[3],q[7];
u3(0,0,-1.0*g1) q[7];
cx q[3],q[7];
cx q[3],q[14];
cx q[14],q[3];
cx q[3],q[14];
cx q[9],q[5];
cx q[5],q[9];
cx q[9],q[5];
cx q[5],q[18];
u3(0,0,-1.0*g1) q[18];
cx q[5],q[18];
cx q[1],q[18];
u3(0,0,-1.0*g1) q[18];
cx q[1],q[18];
cx q[5],q[11];
u3(0,0,-1.0*g1) q[11];
cx q[5],q[11];
cx q[11],q[0];
cx q[0],q[11];
cx q[11],q[0];
cx q[1],q[0];
u3(0,0,-1.0*g1) q[0];
cx q[1],q[0];
u3(b1,-pi/2,pi/2) q[1];
cx q[1],q[7];
cx q[5],q[17];
cx q[17],q[5];
cx q[5],q[17];
cx q[17],q[10];
u3(0,0,-1.0*g1) q[10];
cx q[17],q[10];
u3(b1,-pi/2,pi/2) q[17];
cx q[7],q[1];
cx q[1],q[7];
cx q[8],q[18];
cx q[18],q[8];
cx q[8],q[18];
cx q[18],q[1];
u3(0,0,-1.0*g1) q[1];
cx q[18],q[1];
cx q[6],q[8];
cx q[7],q[1];
cx q[1],q[7];
cx q[7],q[1];
cx q[3],q[7];
u3(0,0,-1.0*g1) q[7];
cx q[3],q[7];
cx q[15],q[3];
cx q[3],q[15];
cx q[15],q[3];
cx q[3],q[7];
cx q[7],q[3];
cx q[3],q[7];
u3(0,0,-1.0*g1) q[8];
cx q[6],q[8];
cx q[6],q[13];
cx q[13],q[6];
cx q[6],q[13];
cx q[13],q[0];
u3(0,0,-1.0*g1) q[0];
cx q[13],q[10];
u3(0,0,-1.0*g1) q[10];
cx q[13],q[10];
cx q[0],q[13];
cx q[13],q[0];
cx q[0],q[11];
cx q[11],q[0];
cx q[0],q[11];
cx q[11],q[2];
u3(0,0,-1.0*g1) q[2];
cx q[11],q[2];
u3(b1,-pi/2,pi/2) q[11];
cx q[8],q[6];
cx q[6],q[8];
cx q[8],q[6];
cx q[14],q[6];
cx q[18],q[8];
u3(0,0,-1.0*g1) q[6];
cx q[14],q[6];
cx q[14],q[10];
u3(0,0,-1.0*g1) q[10];
cx q[14],q[10];
cx q[14],q[16];
cx q[16],q[14];
cx q[14],q[16];
cx q[16],q[2];
u3(0,0,-1.0*g1) q[2];
cx q[16],q[2];
u3(b1,-pi/2,pi/2) q[16];
cx q[2],q[4];
cx q[4],q[2];
cx q[2],q[4];
cx q[8],q[18];
cx q[18],q[8];
cx q[8],q[6];
u3(0,0,-1.0*g1) q[6];
cx q[8],q[6];
cx q[6],q[13];
cx q[13],q[6];
cx q[6],q[13];
cx q[15],q[13];
u3(0,0,-1.0*g1) q[13];
cx q[8],q[6];
u3(0,0,-1.0*g1) q[6];
cx q[8],q[6];
cx q[15],q[6];
u3(0,0,-1.0*g1) q[6];
cx q[15],q[6];
cx q[13],q[15];
cx q[15],q[13];
cx q[13],q[10];
u3(0,0,-1.0*g1) q[10];
cx q[13],q[10];
u3(b1,-pi/2,pi/2) q[13];
cx q[3],q[15];
u3(0,0,-1.0*g1) q[15];
cx q[3],q[15];
cx q[6],q[14];
cx q[14],q[6];
cx q[6],q[14];
cx q[3],q[14];
u3(0,0,-1.0*g1) q[14];
cx q[3],q[14];
cx q[14],q[10];
cx q[10],q[14];
cx q[14],q[10];
cx q[3],q[14];
u3(0,0,-1.0*g1) q[14];
cx q[3],q[14];
cx q[3],q[4];
u3(0,0,-1.0*g1) q[4];
cx q[3],q[4];
u3(b1,-pi/2,pi/2) q[3];
cx q[3],q[15];
cx q[15],q[3];
cx q[3],q[15];
cx q[3],q[14];
u3(0,0,-1.0*g1) q[14];
cx q[3],q[14];
cx q[10],q[14];
u3(0,0,-1.0*g1) q[14];
cx q[10],q[14];
u3(b1,-pi/2,pi/2) q[3];
cx q[3],q[14];
cx q[14],q[3];
cx q[3],q[14];
u3(b1,-pi/2,pi/2) q[8];
cx q[10],q[19];
cx q[19],q[10];
cx q[10],q[19];
cx q[19],q[4];
u3(0,0,-1.0*g1) q[4];
cx q[19],q[4];
u3(b1,-pi/2,pi/2) q[19];
cx q[3],q[4];
u3(0,0,-1.0*g1) q[4];
cx q[3],q[4];
u3(b1,-pi/2,pi/2) q[3];
cx q[7],q[4];
u3(0,0,-1.0*g1) q[4];
cx q[7],q[4];
u3(b1,-pi/2,pi/2) q[4];
u3(b1,-pi/2,pi/2) q[7];
barrier q[18],q[17],q[19],q[1],q[4],q[16],q[6],q[0],q[13],q[9],q[2],q[7],q[11],q[14],q[8],q[5],q[3],q[15],q[10],q[12];
measure q[17] -> c[0];
measure q[7] -> c[1];
measure q[1] -> c[2];
measure q[11] -> c[3];
measure q[16] -> c[4];
measure q[8] -> c[5];
measure q[13] -> c[6];
measure q[15] -> c[7];
measure q[14] -> c[8];
measure q[19] -> c[9];
measure q[3] -> c[10];
measure q[4] -> c[11];