proof.all(False)
import sys
import random
import time
import importlib
import msidh

# Settings
settings = msidh.load_128()


p = settings.p
NA = settings.A
NB = settings.B
prime_A = settings.Af
prime_B = settings.Bf
E0 = settings.E0
PA = settings.PA
QA = settings.QA
PB = settings.PB
QB = settings.QB
G = settings.G

string_A = (NA).str(base=2);
string_B = (NB).str(base=2);

coefficients = [E0.a1(), E0.a2(), E0.a3(), E0.a4(), E0.a6()];
print("Coefficients of the curve: ", coefficients);

(x,) = G._first_ngens(1)
GX.<X> = G[];
print(PA, QA, PB, QB);


RA = PA - QA;

RB = PB - QB;

degree_bound = min(len(prime_A), len(prime_B));


T1 = time.time();
# Square root in F_{p^2}
def sqrFp2(x):
     [r, i] = x.polynomial().list();
     delta = r^2 + i^2;
     lamb = delta^((p - 3)/4);
     rho = delta * lamb;
     if rho^2 != delta:
         return 0;
     gamma = (r + rho)/2;
     if gamma == 0:
         gamma = -rho;
     mu = gamma^((p - 3)/4);
     sigma = gamma * mu;
     gamma_inv = sigma * mu^3;
     tau = sigma * gamma_inv;
     omega = (i/2)*tau;
     if sigma^2 == gamma:
         return sigma + omega * x;
     else:
         return omega - sigma*x;

# Point doubling
def xDBL(X1, Z1, A24, C24):
	t0 = X1 - Z1; t1 = X1 + Z1;
	t0 = t0^2;
	t1 = t1^2;
	Z1 = C24 * t0;
	X1 = t1 * Z1;
	t1 = t1 - t0;
	t0 = A24 * t1;
	Z1 = Z1 + t0;
	Z1 = Z1 * t1;
	return X1, Z1;

# Point doubling and point addition
def xDBLADD(X1, Z1, X2, Z2, X3, Z3, A24):
	t0 = X1 + Z1; t1 = X1 - Z1;
	X4 = t0^2;
	t2 = X2 - Z2; 
	X5 = X2 + Z2;
	t0 = t0 * t2;
	Z4 = t1^2;
	t1 = t1 * X5; t2 = X4 -Z4;
	X4 = X4 * Z4;
	X5 = A24 * t2;
	Z5 = t0 - t1;
	Z4 = X5 + Z4;
	X5 = t0 + t1;
	Z4 = Z4 * t2;
	Z5 = Z5^2;
	X5 = X5^2;
	Z5 = X3 * Z5;
	X5 = Z3 * X5;
	return X4, Z4, X5, Z5;
	

# Montgomery ladder
def Montgomery_ladder(ell, A, XP, ZP):
	X1 = XP;
	Z1 = ZP;
	A24 = A + G(2);
	C24 = G(4);
	X2, Z2 = xDBL(X1, Z1, A24, C24);
	ellbits = bin(ell);
	ellbits = tuple(ellbits);
	ellbits = ellbits[3:];
	l = len(ellbits);
	A24 = A24/C24;
	for j in range(l):
		if ellbits[j] == '1': 
			[X2, Z2, X1, Z1] = xDBLADD(X2, Z2, X1, Z1, XP, ZP, A24);
		else:
			[X1, Z1, X2, Z2] = xDBLADD(X1, Z1, X2, Z2, XP, ZP, A24);
	return X1, Z1, X2, Z2;

# Differential addition
def differential_addition_Montgomery_plus(X1, Z1, X2, Z2, X3, Z3, A0):
  if X1 == 0 or Z1 == 0 or [X2,Z2] == [0,0] or [X3,Z3] == [0,0]:
    return 0, 0;
  else:
    t0 = X1 + Z1; t1 = X1 - Z1;
    X4 = t0^2;
    t2 = X2 - Z2; 
    X5 = X2 + Z2;
    t0 = t0 * t2;
    Z4 = t1^2;
    t1 = t1 * X5; t2 = X4 - Z4;
    X4 = X4 * Z4;
    X5 = ((A0 + G(2))/G(4)) * t2;
    Z5 = t0 - t1;
    Z4 = X5 + Z4;
    X5 = t0 + t1;
    Z4 = Z4 * t2;
    Z5 = Z5^2;
    X5 = X5^2;
    Z5 = X3 * Z5;
    X5 = Z3 * X5;
    return X5,Z5;

# Point doubling
def double_point_Montgomery_plus(X2, Z2, A0):
  if Z2 == 0 or X2^3+A0*X2^2+X2 == 0:
    return 0, 0;
  else:
    return (X2^2-Z2^2)^2, 4*X2*Z2*(X2^2+A0*X2*Z2+Z2^2);

def recover_y3(X0, Y0, X1, Z1, X2, Z2, A):
	x = X0;
	y = Y0;
	t1 = x * Z1;
	t2 = X1 + t1;
	t3 = X1 - t1;
	t3 = t3^2;
	t3 = t3 * X2;
	t1 = 2 * Z1;
	t1 = t1 * A;
	t2 = t1 + t2;
	t4 = x * X1;
	t4 = t4 + Z1;
	t2 = t2 * t4;
	t1 = t1 * Z1;
	t2 = t2 - t1;
	t2 = t2 * Z2;
	Y = t2 - t3;
	t1 = 2 * y;
	t1 = t1 * Z1;
	t1 = t1 * Z2;
	X = X1 * t1;
	Z = Z1 * t1;
	if X == 0 and Y ==0 and Z == 0:
		Y = G(1);
	return X, Y, Z;

def step_zero_Montgomery_plus(X1, Z1, X2, Z2, X3, Z3, A0):
  return (X2^2-Z2^2)^2, 4*X2*Z2*(X2^2+A0*X2*Z2+Z2^2), 4*(X2*X3-Z2*Z3)^2*Z1, 4*(X2*Z3-Z2*X3)^2*X1;
  
def step_one_Montgomery_plus(X1, Z1, X3, Z3, X2, Z2, A0):
  return 4*(X2*X3-Z2*Z3)^2*Z1, 4*(X2*Z3-Z2*X3)^2*X1, (X2^2-Z2^2)^2, 4*X2*Z2*(X2^2+A0*X2*Z2+Z2^2);

# Scalar multiplication
def scalar_multiplication_Montgomery_plus(n, X1, Z1, A0):
  X2 = 1; Z2 = 0; X3 = X1; Z3 = Z1;
  nbits = [];
  while(n >= 1):
    nbits.append(n % 2);
    n = n // 2;
  nbits.reverse();
  if Z1 == 0:  
    print("Error");
  else:
    for i in range(len(nbits)):
      if nbits[i] == 0: 
        X2, Z2, X3, Z3 = step_zero_Montgomery_plus(X1, Z1, X2, Z2, X3, Z3, A0);
      else:
        X2, Z2, X3, Z3 = step_one_Montgomery_plus(X1, Z1, X2, Z2, X3, Z3, A0);
  return X2, Z2;
    
def F0(X1,X2,X3):
    return (X3*X1-X2)^2

def F1(X1,X2,X3,A0):
    return -2*((X1*X2+X3)*(X1*X3+X2)+2*A0*X1*X2*X3)

def F2(X1,X2,X3):
    return (X1*X2-X3)^2

# Compute I, J, K
def IJK(XP, ZP, l, A0):
  S = set(range(1,l-1,2));
  b = floor(sqrt(l - 1)/2);
  if b == 0:
    b2 = 0;
  else:
    b2 = floor((l - 1)/(4*b));
  I = set(2*b*(2*i + 1) for i in range(b2));
  J = set(2*j+1 for j in range(b));
  K = set(range(4*b*b2+1,l-1,2));
  xS = {s:scalar_multiplication_Montgomery_plus(s, XP, ZP, A0) for s in K.union(I).union(J)}
  return [xS, I, J, K];
  
 
# Compute hI(x) = \prod_{i \in I} (x - x_i)
def h_I(x1, I):
  hI = GX(prod(x1[i][1]*X-x1[i][0] for i in I));
  return hI;
 
# Compute hS(x)
def hS_for_iso(x1, alpha, A0, J, K, hI):
  EJ = GX(prod(F0(X,x1[j][0],x1[j][1])*alpha^2+F1(X,x1[j][0],x1[j][1],A0)*alpha+F2(X,x1[j][0],x1[j][1]) for j in J));
  R = hI.resultant(EJ);
  hK = prod(alpha*x1[s][1]-x1[s][0] for s in K);
  return hK*R;
  

# Large prime degree isogeny
def x_val_big(x1, l, X1, Z1, A0, J, K, hI):
  E_0J = GX(prod(F0(X,x1[j][0],x1[j][1])*Z1^2+F1(X,x1[j][0],x1[j][1],A0)*X1*Z1+F2(X,x1[j][0],x1[j][1])*X1^2 for j in J));
  E_1J = GX(prod(F0(X,x1[j][0],x1[j][1])*X1^2+F1(X,x1[j][0],x1[j][1],A0)*X1*Z1+F2(X,x1[j][0],x1[j][1])*Z1^2 for j in J));
  R0 = hI.resultant(E_0J);
  R1 = hI.resultant(E_1J);
  hK0 = prod(Z1*x1[s][1]- X1*x1[s][0] for s in K);
  hK1 = prod(X1*x1[s][1]- Z1*x1[s][0] for s in K);
  newX1 = X1 * (hK0 * R0)^2;
  newZ1 = Z1 * (hK1 * R1)^2;
  return newX1, newZ1;
  
def x_iso_big(x1, l, A0, J, K, hI):
  newA0 = 2*((A0 + 2)^l*hS_for_iso(x1, G(-1), A0, J, K, hI)^8 + (A0 - 2)^l*hS_for_iso(x1, G(1), A0, J, K, hI)^8)/((A0 + 2)^l*hS_for_iso(x1, G(-1), A0, J, K, hI)^8 - (A0 - 2)^l*hS_for_iso(x1, G(1), A0, J, K, hI)^8);
  return newA0;
 
# Small prime degree isogeny  
def x_val_small(XP, ZP, l, X1, Z1, A0):
  XQ = XP; ZQ = ZP;
  f1 = (X1*XQ - Z1*ZQ); f2 = (X1*ZQ - Z1*XQ);
  if l == 3:
    return X1 * f1^2, Z1 * f2^2;
  else:
    [XDBL, ZDBL] = double_point_Montgomery_plus(XQ, ZQ, A0);
    [XQ, ZQ] = differential_addition_Montgomery_plus(XDBL, ZDBL, XQ, ZQ, XQ, ZQ, A0);
    f1 *= (XQ*X1 - ZQ * Z1); f2 *= (ZQ*X1 - XQ * Z1);
    XPrev = XP; ZPrev = ZP;
    for i in range (3, (l - 1)//2 + 1):
        XTemp = XQ;
        ZTemp = ZQ;
        [XQ, ZQ] = differential_addition_Montgomery_plus(XQ, ZQ, XDBL, ZDBL, XPrev, ZPrev, A0);
        f1 *= (XQ*X1 - ZQ * Z1); f2 *= (ZQ*X1 - XQ*Z1);
        XPrev = XTemp; ZPrev = ZTemp;
  return X1 * f1^2, Z1 * f2^2;


def x_iso_small(XP, ZP, l, A0):
  XQ = XP; ZQ = ZP;
  pi = XQ/ZQ; sigma = (XQ^2 - ZQ^2)/(XQ*ZQ);
  if l == 3:
    return pi^2 * (A0 - 6*sigma);
  else:
    [XQ, ZQ] = double_point_Montgomery_plus(XQ, ZQ, A0);
    xQ = XQ/ZQ;
    pi *= xQ; sigma += xQ - 1/xQ;
    XPrev = XP; ZPrev = ZP;
    for i in range (3, (l - 1)//2 + 1):
        XTemp = XQ;
        ZTemp = ZQ;
        [XQ, ZQ] = differential_addition_Montgomery_plus(XQ, ZQ, XP, ZP, XPrev, ZPrev, A0);
        xQ = XQ/ZQ;
        pi *= xQ; sigma += xQ - 1/xQ;
        XPrev = XTemp; ZPrev = ZTemp;
  return pi^2 * (A0 - 6*sigma);

 
# Input: E0, P_B, Q_B, R_B, Q_A
# Output: \phi_A(Q_A), EA
def msidh_action_Alice(A0, XPush, private_key, degree_bound, X1, X2, Y1, Y2):
  S = []; 
  Ker = PA + private_key * QA;
  S.append([len(prime_A), [Ker[0], 1]]);
  i = 0;
  Aimage = A0;
  #APre.append(Aimage);
  while S != []:
    [h, [XP, ZP]] = S[0];
    S.pop(0);
    if h == 1 and S != []:
      ATemp = Aimage;
      if prime_A[len(prime_A) - i] >= degree_bound:
        [x1, I, J, K] = IJK(XP, ZP, prime_A[len(prime_A) - i], ATemp);
        hI = h_I(x1, I);
        Aimage = x_iso_big(x1, prime_A[len(prime_A) - i], ATemp, J, K, hI);
        S1 = [];
        while S != []:
          [h, [xTemp1, zTemp1]] = S[len(S)-1];
          [xTemp2, zTemp2] = [XP, ZP];
          S.pop();
          [XP, ZP] = x_val_big(x1,  prime_A[len(prime_A) - i], xTemp1 , zTemp1,  ATemp, J, K, hI);
          S1.insert(0, [h - 1, [XP,ZP]]);
        S = S1;
        for j in range(3):    
            [XPush[j][0], XPush[j][1]] = x_val_big(x1, prime_A[len(prime_A) - i], XPush[j][0], XPush[j][1] , ATemp, J, K, hI);  
      else:
        Aimage = x_iso_small(XP, ZP, prime_A[len(prime_A) - i], ATemp);
        S1 = [];
        while S != []:
          [h, [xTemp1, zTemp1]] = S[len(S)-1];
          [xTemp2, zTemp2] = [XP, ZP];
          S.pop();
          [XP, ZP] = x_val_small(XP, ZP, prime_A[len(prime_A) - i], xTemp1 , zTemp1, ATemp);
          S1.insert(0, [h - 1, [XP, ZP]]);
        S = S1;
        for j in range(3):
            [XPush[j][0], XPush[j][1]] = x_val_small(xTemp2, zTemp2, prime_A[len(prime_A) - i], XPush[j][0], XPush[j][1] , ATemp);
    elif h == 1 and S == []: 
      ATemp = Aimage;
      i = i + 1;
      if prime_A[len(prime_A) - i] >= degree_bound:
        [x1, I, J, K] = IJK(XP, ZP, prime_A[len(prime_A) - i], ATemp);
        hI = h_I(x1, I);
        Aimage = x_iso_big(x1, prime_A[len(prime_A) - i], ATemp, J, K, hI);
        for j in range(3):
            [XPush[j][0], XPush[j][1]] = x_val_big(x1, prime_A[len(prime_A) - i], XPush[j][0], XPush[j][1] , ATemp,  J, K, hI);
      else:
        Aimage = x_iso_small(XP, ZP, prime_A[len(prime_A) - i], ATemp);
        for j in range(3):
            [XPush[j][0], XPush[j][1]] = x_val_small(XP, ZP, prime_A[len(prime_A) - i], XPush[j][0], XPush[j][1] , ATemp); 
    elif s[i] > 0 and s[i] < h:
      S.insert(0, [h, [XP, ZP]]);
      for j in range(s[i]):
          [XP, ZP] = scalar_multiplication_Montgomery_plus(prime_A[j], XP, ZP, Aimage);
      S.insert(0, [h - s[i], [XP, ZP]]);
      i = i + 1;
  for i in range(3):
      XPush[i] = XPush[i][0]/XPush[i][1];
  return Aimage, XPush;
 
# Input: EB, U, V, SA(the result of the discrete logarithm computation)
# Output: EBA
def msidh_action_Alice_last(A0, private_key, degree_bound, U, V, C1, C2, C3, C4):
  S = [];
  P = (C1 + C2 * private_key)*U + (C3 + C4 * private_key)*V;
  S.append([286, [P[0], 1]]);
  i = 0;
  Aimage = A0;
  while S != []:
    [h, [XP,ZP]] = S[0];
    S.pop(0);
    if h == 1 and S != []:
      ATemp = Aimage;
      if prime_A[len(prime_A) - i] >= degree_bound:
        [x1, I, J, K] = IJK(XP, ZP, prime_A[len(prime_A) - i], ATemp);
        hI = h_I(x1, I);
        Aimage = x_iso_big(x1, prime_A[len(prime_A) - i], ATemp, J, K, hI);
        S1 = [];
        while S != []:
          [h, [xTemp1, zTemp1]] = S[len(S)-1];
          [xTemp2, zTemp2] = [XP, ZP];
          S.pop();
          [XP, ZP] = x_val_big(x1,  prime_A[len(prime_A) - i], xTemp1 , zTemp1,  ATemp,J, K, hI);
          S1.insert(0, [h - 1, [XP,ZP]]);
        S = S1;
      else:
        Aimage = x_iso_small(XP, ZP, prime_A[len(prime_A) - i], ATemp);
        S1 = [];
        while S != []:
          [h, [xTemp1, zTemp1]] = S[len(S)-1];
          [xTemp2, zTemp2] = [XP, ZP];
          S.pop();
          [XP, ZP] = x_val_small(XP, ZP, prime_A[len(prime_A) - i], xTemp1 , zTemp1, ATemp);
          S1.insert(0, [h - 1, [XP, ZP]]);
        S = S1;
    elif h == 1 and S == []: 
      ATemp = Aimage;
      i = i + 1;
      if prime_A[len(prime_A) - i] >= degree_bound:
        [x1, I, J, K] = IJK(XP, ZP, prime_A[len(prime_A) - i], ATemp);
        hI = h_I(x1, I);
        Aimage = x_iso_big(x1, prime_A[len(prime_A) - i], ATemp,J, K, hI);
      else:
        Aimage = x_iso_small(XP, ZP, prime_A[len(prime_A) - i], ATemp);
    elif s[i] > 0 and s[i] < h:
      S.insert(0, [h, [XP, ZP]]);
      for j in range(s[i]):
          [XP, ZP] = scalar_multiplication_Montgomery_plus(prime_A[j], XP, ZP, Aimage);
      S.insert(0, [h - s[i], [XP, ZP]]);
      i = i + 1;
  return Aimage;

# Input: E0, P_A, Q_A, R_A, Q_B
# Output: \phi_B(Q_B), EB
def msidh_action_Bob(A0, XPush, private_key, degree_bound, X1, X2, Y1, Y2):
  S = [];
  Ker = PB + private_key * QB;
  S.append([286, [Ker[0], 1]]);
  i = 0;
  Bimage = A0;
  while S != []:
    [h, [XP, ZP]] = S[0];
    S.pop(0);
    if h == 1 and S != []:
      BTemp = Bimage;
      if prime_B[len(prime_B) - i] >= degree_bound:
        [x1, I, J, K] = IJK(XP, ZP, prime_B[len(prime_B) - i], BTemp);
        hI = h_I(x1, I);
        Bimage = x_iso_big(x1, prime_B[len(prime_B) - i], BTemp, J, K, hI);
        S1 = [];
        while S != []:
          [h, [xTemp1, zTemp1]] = S[len(S)-1];
          [xTemp2, zTemp2] = [XP, ZP];
          S.pop();
          [XP, ZP] = x_val_big(x1,  prime_B[len(prime_B) - i], xTemp1 , zTemp1,  BTemp, J, K, hI);
          S1.insert(0, [h - 1, [XP,ZP]]);
        S = S1;
        for j in range(3):
            [XPush[j][0], XPush[j][1]] = x_val_big(x1, prime_B[len(prime_B) - i], XPush[j][0], XPush[j][1] , BTemp, J, K, hI);
      else:
        Bimage = x_iso_small(XP, ZP, prime_B[len(prime_B) - i], BTemp);
        S1 = [];
        while S != []:
          [h, [xTemp1, zTemp1]] = S[len(S)-1];
          [xTemp2, zTemp2] = [XP, ZP];
          S.pop();
          [XP, ZP] = x_val_small(XP, ZP, prime_B[len(prime_B) - i], xTemp1 , zTemp1, BTemp);
          S1.insert(0, [h - 1, [XP, ZP]]);
        S = S1;
        for j in range(3):
            [XPush[j][0], XPush[j][1]] = x_val_small(xTemp2, zTemp2, prime_B[len(prime_B) - i], XPush[j][0], XPush[j][1] , BTemp);
    elif h == 1 and S == []: 
      BTemp = Bimage;
      i = i + 1;
      if prime_B[len(prime_B) - i] >= degree_bound:
        [x1, I, J, K] = IJK(XP, ZP, prime_B[len(prime_B) - i], BTemp);
        hI = h_I(x1, I);
        Bimage = x_iso_big(x1, prime_B[len(prime_B) - i], BTemp, J, K, hI);
        for j in range(3):
            [XPush[j][0], XPush[j][1]] = x_val_big(x1, prime_B[len(prime_B) - i], XPush[j][0], XPush[j][1] , BTemp,  J, K, hI); 
      else:
        Bimage = x_iso_small(XP, ZP, prime_B[len(prime_B) - i], BTemp);
        for j in range(3):
            [XPush[j][0], XPush[j][1]] = x_val_small(XP, ZP, prime_B[len(prime_B) - i], XPush[j][0], XPush[j][1] , BTemp); 
    elif s[i] > 0 and s[i] < h:
      S.insert(0, [h, [XP, ZP]]);
      for j in range(s[i]):
          [XP, ZP] = scalar_multiplication_Montgomery_plus(prime_B[j], XP, ZP, Bimage);
      S.insert(0, [h - s[i], [XP, ZP]]);
      i = i + 1;
  for i in range(3):
      XPush[i] = XPush[i][0]/XPush[i][1];
  return Bimage, XPush;
  
# Input: EA, U, V, SB(the result of the discrete logarithm computation)
# Output: EAB 
def msidh_action_Bob_last(A0, private_key, degree_bound, U, V, C1, C2, C3, C4):
  S = [];
  P = (C1 + C2 * private_key)*U + (C3 + C4 * private_key)*V;
  [X4, Z4] = [P[0], P[2]];
  S.append([286, [X4, Z4]]);
  i = 0;
  Bimage = A0;
  while S != []:
    [h, [XP,ZP]] = S[0];
    S.pop(0);
    if h == 1 and S != []:
      BTemp = Bimage;
      if prime_B[len(prime_B) - i] >= degree_bound:
        [x1, I, J, K] = IJK(XP, ZP, prime_B[len(prime_B) - i], BTemp);
        hI = h_I(x1, I);
        Bimage = x_iso_big(x1, prime_B[len(prime_B) - i], BTemp, J, K, hI);
        S1 = [];
        while S != []:
          [h, [xTemp1, zTemp1]] = S[len(S)-1];
          [xTemp2, zTemp2] = [XP, ZP];
          S.pop();
          [XP, ZP] = x_val_big(x1,  prime_B[len(prime_B) - i], xTemp1 , zTemp1,  BTemp,J, K, hI);
          S1.insert(0, [h - 1, [XP,ZP]]);
        S = S1;
      else:
        Bimage = x_iso_small(XP, ZP, prime_B[len(prime_B) - i], BTemp);
        S1 = [];
        while S != []:
          [h, [xTemp1, zTemp1]] = S[len(S)-1];
          [xTemp2, zTemp2] = [XP, ZP];
          S.pop();
          [XP, ZP] = x_val_small(XP, ZP, prime_B[len(prime_B) - i], xTemp1 , zTemp1, BTemp);
          S1.insert(0, [h - 1, [XP, ZP]]);
        S = S1;
    elif h == 1 and S == []: 
      BTemp = Bimage;
      i = i + 1;
      if prime_B[len(prime_B) - i] >= degree_bound:
        [x1, I, J, K] = IJK(XP, ZP, prime_B[len(prime_B) - i], BTemp);
        hI = h_I(x1, I);
        Bimage = x_iso_big(x1, prime_B[len(prime_B) - i], BTemp,J, K, hI);
      else:
        Bimage = x_iso_small(XP, ZP, prime_B[len(prime_B) - i], BTemp);
    elif s[i] > 0 and s[i] < h:
      S.insert(0, [h, [XP, ZP]]);
      for j in range(s[i]):
          [XP, ZP] = scalar_multiplication_Montgomery_plus(prime_B[j], XP, ZP, Bimage);
      S.insert(0, [h - s[i], [XP, ZP]]);
      i = i + 1;
  return Bimage;

  
# Compute the j-invariant
def J(MA):
  return 256*(MA^2 - 3)^3 / (MA^2 - 4);

# Lucas sequences (Algorithm 1)
def Lucassequences(f,ell):
	v0 = G(2);
	tmp1 = f;
	v1 = tmp1;
	tmp2 = G(2);
	string = ell.str(base=2);
	Len = len(string);
	for j in range(Len):
		if string[j] == '1':
			v0 = v0 * v1 - tmp1;
			v1 = v1 * v1 - tmp2;
		else:
			v1 = v0 * v1 - tmp1;
			v0 = v0 * v0 - tmp2;
	return v0;

# Lucas sequences (Algorithm 1), but return v1 as well	
def Lucassequences2(f,ell):
	v0 = 2;
	tmp1 = f + f^-1;
	v1 = tmp1;
	tmp2 = G(2);
	string = ell.str(base=2);
	Len = len(string);
	for j in range(Len):
		if string[j] == '1':
			v0 = v0 * v1 - tmp1;
			v1 = v1 * v1 - tmp2;
		else:
			v1 = v0 * v1 - tmp1;
			v0 = v0 * v0 - tmp2;
	return v0, v1;
	
# Exponentiation using Lucas sequences (Algorithm 2)
def ELS(f,ell):
	if ell < 100000000000000000000000000:
		return f^ell;
	else:
		[tmp1,tmp2] = Lucassequences2(f,ell-1);
		tmp1 = tmp1/G(2); 
		tmp2 = tmp2/G(2);
		re = tmp2;
		f0 = (f.polynomial().list())[0];
		f1 = f - f0;		
		tmp2 = tmp2 * f0 - tmp1;
		tmp2 = tmp2* (f0^2-1)^(-1);
		im = tmp2 * f1;
		return re + im;

# elligator: randomly generate a point
def elligator(A, E):
  flag_2 = 0;
  u = 1 + x;
  while flag_2 == 0:
    r0 = random.randint(0, p);
    r1 = random.randint(0, p);
    r = r0 + r1 * x;
    tmp = u * r^2;
    cond1 = 1 + tmp;
    cond2l = A^2 * tmp;
    cond2r = cond1^2;
    if cond1 != 0 and cond2l != cond2r:
      v = -A/cond1;
      e = is_square(v * (v^2 + A*v + G(1)));
    if e == False:
       e = -1;
    else:
      e = 1;
    x = e * v - (1 - e)*(A/G(2));
    y = - e * sqrFp2(x * (x^2 + A * x + 1));
    flag_2 = 1;
  return E([x,y]);

# Batch cofactor exponentiation using trace (Algorithm 5)
def batch_cof_trace(fi,f,ind,n, L):
	if n == 1:
		fi.append(f);
		return 1;
	m=floor(n/2);
	u=1;
	v=1;
	for i in range(m):
		u=L[ind[i]]*u;
	for j in range(m,n):
		v=L[ind[j]]*v;
	left=Lucassequences(f,u);
	right=Lucassequences(f,v);
	batch_cof_trace(fi,left,[ind[i] for i in range(m,n)],n-m, L);
	batch_cof_trace(fi,right,[ind[i] for i in range(m)],m, L);

# Batch cofactor exponentiation in the cyclic group (Algorithm 7)
def batch_cof(fi,f,ind,n, L):
	if n == 1:
		fi.append(f);
		return 1;
	m=floor(n/2);
	u=1;
	v=1;
	for i in range(m):
		u=L[ind[i]]*u;
	for j in range(m,n):
		v=L[ind[j]]*v;
	left=ELS(f,u);
	right=ELS(f,v);
	batch_cof(fi,left,[ind[i] for i in range(m,n)],n-m, L);
	batch_cof(fi,right,[ind[i] for i in range(m)],m, L);

# Batch cofactor multiplication (Algorithm 4)
# Output:{Ui},the order of each point is ln,ln-1,...,l1.
def batch_cof_mul(Ui,U,ind,n,A, E, L):
    if n == 1:
        Ui.append(U);
        return 1;
    m=floor(n/2);
    u=1;
    v=1;
    for i in range(m):
        u=L[ind[i]]*u;
    for j in range(m,n):
        v=L[ind[j]]*v;
    R0X, R0Z, R1X, R1Z = Montgomery_ladder(u, A, U[0], U[2]);
    UX, UY, UZ = recover_y3(U[0], U[1], R0X, R0Z, R1X, R1Z, A);
    left = E([UX, UY, UZ]);
    R0X, R0Z, R1X, R1Z = Montgomery_ladder(v, A, U[0], U[2]);
    UX, UY, UZ = recover_y3(U[0], U[1], R0X, R0Z, R1X, R1Z, A);
    right = E([UX, UY, UZ]);
    batch_cof_mul(Ui,left,[ind[i] for i in range(m,n)],n-m, A, E, L);
    batch_cof_mul(Ui,right,[ind[i] for i in range(m)],m, A, E, L);

# Generate U 
def generate_U(Uj,A, ind, E, L, N):
	R=elligator(A, E); 
	R0X, R0Z, R1X, R1Z = Montgomery_ladder(4*N,A,R[0],R[2]);
	UX, UY, UZ = recover_y3(R[0], R[1], R0X, R0Z, R1X, R1Z, A);
	U = E([UX,UY,UZ]);
	n=len(ind);
	batch_cof_mul(Uj,U,ind,n,A, E, L);
	Iu=[];
	for i in range(n):
		if Uj[i] == E([0,1,0]):
			Iu.append(285-i);
	Seed = 1;
	while len(Iu) !=0:
		random.seed(int(Seed));
		R=elligator(A, E);
		R0X, R0Z, R1X, R1Z = Montgomery_ladder(4*N,A,R[0],R[2]);
		UX, UY, UZ = recover_y3(R[0], R[1], R0X, R0Z, R1X, R1Z, A);
		U1 = E([UX,UY,UZ]);
		Uj1=[];
		n=len(Iu);
		prime_prod = 1;
		for i in list(set(ind)-set(Iu)):
			prime_prod = prime_prod * L[i];
		R0X, R0Z, R1X, R1Z = Montgomery_ladder(prime_prod,A,U1[0],U1[2]);
		UX, UY, UZ = recover_y3(U1[0], U1[1], R0X, R0Z, R1X, R1Z, A);
		U1 = E([UX,UY,UZ]);
		batch_cof_mul(Uj1,U1,Iu,n,A, E, L);
		tmp=Iu;
		Iu=[];
		for i in range(n):
			if Uj1[i] != E([0,1,0]):
				U = U+Uj1[i];
				Uj[285-tmp[i]]=Uj1[i];
			else:
				Iu.append(tmp[i]);
				
		Seed += 1;
	return U;
       
# Generate V 
def generate_V(A,ind,U,Uj, E, primelist, N):
    V = elligator(A, E); 
    if N == NA:
      string_N = string_A;
    else:
      string_N = string_B;
    f = compressed_pairing(U,V,string_N,N, A);
    fj = [];
    n = len(ind);
    batch_cof_trace(fj,f,ind,n, primelist);
    Iv = [];
    for i in range(n):
        if fj[i] == 2:
            Iv.append(285-i);
    Seed = 1;
    while len(Iv) != 0:
        random.seed(int(Seed));
        V1 = elligator(A, E); 
        U1 = E([0,1,0]);
        L = 1;
        for i in Iv:
            U1 = U1 + Uj[285-i];
            L = L * primelist[i];
        n = len(Iv);
        string_L = (L).str(base=2);
        f1 = compressed_pairing(U1,V1,string_L,L, A);
        fj1 = [];
        batch_cof_trace(fj1,f1,Iv,n, primelist);
        flag = 0;
        for i in range(n):
          if fj[i] != 2:
            flag = 1;
        if flag == 1:
          prime_prod = 1;
          for k in list(set(ind)-set(Iv)):
            prime_prod = prime_prod * primelist[k];
          R0X, R0Z, R1X, R1Z = Montgomery_ladder(prime_prod,A,V1[0],V1[2]);
          VX, VY, VZ = recover_y3(V1[0], V1[1], R0X, R0Z, R1X, R1Z, A);
          V1 = E([VX,VY,VZ]);
          tmp = Iv;
          Iv = [];
          Vj1 = [];
          batch_cof_mul(Vj1,V1,tmp,n,A, E, primelist);
          for j in range(n):
              if fj1[j] == 2:
                  Iv.append(tmp[n-1-j]);
              else:
                  V = V+Vj1[j];
          Seed += 1;
    return V;

# Pairing
def miller_affine(P,Q,string,r, A):
    (Px,Py,Pz)=P;
    (Qx,Qy,Qz)=Q;
    f1=G(1);
    Rx=Px;
    Ry=Py;

    l=len(string);
    flag=0;
    for i in range(1,l):
        print(Ry)
        lamb=(3*Rx^2+2*A*Rx+1)/(2*Ry);

        f1 = f1^2*(Qy-Ry-lamb*(Qx-Rx));

        Nx = lamb^2-2*Rx-A;
        Ny = lamb*(Rx-Nx)-Ry;

        Rx = Nx;
        Ry = Ny;
        if Qx != Rx:
            f1 = f1/(Qx-Rx);
        else:
            flag=1;
        if string[i] == '1' and i != l-1:
            if Rx != Px:
                lamb_1 =(Ry-Py)/(Rx-Px);
            else:
                lamb_1 = 1;
                flag = 1;
            f1 = f1*(Qy-Py-lamb_1*(Qx-Px));
            Nx = lamb_1^2-(Rx+Px)-A;
            Ny = lamb_1*(Px-Nx)-Py;
            Rx = Nx;
            Ry = Ny;
            if Qx != Rx:
                f1 = f1/(Qx-Rx);
            else:
                flag=1;
    f1 = f1*(Qx-Rx);
    [f11,f12] = f1.polynomial().list();
    f1 = (f11 - f12 *x)/f1;
    f1 = ELS(f1,ceil((p+1)/r));
    if flag == 1:
        f1 = 1;
    return (f1);
    

# Compressed_pairing
def compressed_pairing(P,Q,string,r, A):
    (Px,Py,Pz)=P;
    (Qx,Qy,Qz)=Q;
    f1=G(1);
    Rx=Px;
    Ry=Py;

    l=len(string);
    flag=0;
    for i in range(1,l):
        lamb=(3*Rx^2+2*A*Rx+1)/(2*Ry);

        f1 = f1^2*(Qy-Ry-lamb*(Qx-Rx));

        Nx = lamb^2-2*Rx-A;
        Ny = lamb*(Rx-Nx)-Ry;

        Rx = Nx;
        Ry = Ny;
        if Qx != Rx:
            f1 = f1/(Qx-Rx);
        else:
            flag=1;
        if string[i] == '1' and i != l-1:
            if Rx != Px:
                lamb_1 =(Ry-Py)/(Rx-Px);
            else:
                lamb_1 = 1;
                flag = 1;
            f1 = f1*(Qy-Py-lamb_1*(Qx-Px));
            Nx = lamb_1^2-(Rx+Px)-A;
            Ny = lamb_1*(Px-Nx)-Py;
            Rx = Nx;
            Ry = Ny;
            if Qx != Rx:
                f1 = f1/(Qx-Rx);
            else:
                flag=1;
    f1 = f1*(Qx-Rx);
    [f11, f12] = f1.polynomial().list();
    f1 = (f11 - f12 * x)/f1;
    [f11, f12] = f1.polynomial().list();
    f1 = 2 * f11;
    f1 = Lucassequences(f1, ceil((p+1)/r));
    if flag == 1:
        f1 = 2;
    return (f1);

# Chinese Remainder Theorem
def chinese_remainder(s, L):
	if len(s) == 0:
		return 1;
	while len(s) > 1:
		label = 0;
		if mod(len(s),2) == 1:
			tmps = s[len(s)-1];
			tmpL = L[len(s)-1];
			label = 1;
		s = [crt(s[i*2], s[i*2 + 1], L[i*2], L[i*2 + 1]) for i in range(floor(len(s)/2))]; 
		L = [prod(L[i*2 : (i + 1)*2])  for i in range(floor(len(L)/2))];
		if label == 1:
			s.append(tmps);
			L.append(tmpL);
	return s[0];
  

# New discrete logarithm computation (Algorithm 9)
def DLC_new(ind, hh0, hh1, hh2, hh3, hh4, Primes):
  n = 286; 
  ans = zero_vector(4*n);
  hh = [];
  for i in range(4*n):
    hh.append(G(0));
  for i in range(n):
    hh[i] = hh1[n-1-i]; 
    hh[n+i] = hh2[n-1-i];
    hh[2*n+i] = hh3[n-1-i];
    hh[3*n+i] = hh4[n-1-i];
  for i in range(n):
    tmp = hh0[n-1-i];
    hhtmp = zero_vector(8);
    label4 = [0,0,0,0];
    for j in range(4):
      if hh[i+j*n] == 1:
        label4[j] = 1;
      else:
        [hhtmp[2*j], hhtmp[2*j+1]] = hh[i+j*n].polynomial().list();
        hhtmp[2*j] = 2 * hhtmp[2*j];
    [tmp1,tmp2] = hh0[n-1-i].polynomial().list();
    tmp0 = 2 * tmp1; v0 = 2; v1 = tmp0; tmp3 = tmp0^2 - G(4);
    j = 0; v2 = v1^2 - 2;
    while label4 != [1,1,1,1]:
      for k in range(4):
        if hhtmp[2*k] == v1:
          label4[k] = 1;
          if j != 0:	
            if hhtmp[2*k+1] * tmp3 == tmp2 * (v2 - v0):
              ans[k*n+i] = j+1; label4[k] = 1;
            else:
              ans[k*n+i] = -j-1; label4[k] = 1; 
          else:
            if hhtmp[2*k+1] == tmp2:
              ans[k*n+i] = 1; label4[k] = 1;
            else:
              ans[k*n+i] = -1; label4[k] = 1;
      j = j + 1; 
      tmp = v2;
      v2 = tmp * tmp0 - v1;
      v0 = v1;
      v1 = tmp;
  S = zero_vector(3); label = 0;
  s1 = []; s2 = []; s3 = [];
  for i in range(n):
    s1.append(1); s2.append(1); s3.append(1);
  for i in range(n):
    if ans[3*n+i] == 0:
      tmp4 = mod(ans[2*n+i], Primes[i]);
      tmp4 = tmp4^-1;
      s2[i] = int(-tmp4 * ans[n+i]); s3[i] = int(tmp4 * ans[i]);
    else:
      tmp4 = mod(ans[3*n+i], Primes[i]);
      tmp4 = tmp4^-1;
      s1[i] = int(-tmp4 * ans[2*n+i]); s2[i] = int(tmp4 * ans[n+i]); s3[i] = int(-tmp4 * ans[i]);
      label = label + 2^i; 
  S[0] = chinese_remainder(s1, Primes);
  S[1] = chinese_remainder(s2, Primes);
  S[2] = chinese_remainder(s3, Primes);
  return S, label;


s = [];

for i in range(1,len(prime_A)):
    s.append(i);

ind = [];

for i in range(len(prime_A)):
    ind.append(i);
   
s.reverse();

Apush = G(E0.a2());

random.seed(int(0)); 
skA = random.randint(1, NA - 1);
random.seed(int(0));
skB = random.randint(1, NB - 1);

h0 = miller_affine(PA, QA, string_A, NA, Apush);
hh0 = [];
ind = [];
for i in range(len(prime_A)):
    ind.append(i);
batch_cof(hh0,h0,ind,len(ind), prime_A);
h5 = miller_affine(PB, QB, string_B, NB, Apush);
hh5 = [];
ind = [];
for i in range(len(prime_B)):
    ind.append(i);
batch_cof(hh5,h5,ind,len(ind), prime_B);

T2 = time.time();

XPushA = [[PB[0], G(1)], [QB[0], G(1)], [RB[0], G(1)]];

[A1, XPushA] = msidh_action_Alice(Apush, XPushA, skA, degree_bound, QA[0], PA[0], QA[1], PA[1]);

xphi_APB =  XPushA[0];

xphi_AQB = XPushA[1];

xphi_ARB = XPushA[2];

EA = EllipticCurve(G,[0,A1,0,1,0]);
phi_APB = EA([xphi_APB, sqrFp2(xphi_APB^3+A1*xphi_APB^2+xphi_APB)]);
phi_AQB = EA([xphi_AQB, sqrFp2(xphi_AQB^3+A1*xphi_AQB^2+xphi_AQB)]);
if (phi_APB[0] + phi_AQB[0] + xphi_ARB + A1) * (phi_APB[0] - phi_AQB[0])^2 == (phi_APB[1] - phi_AQB[1])^2:
  phi_AQB = - phi_AQB; 

T3 = time.time();

print("Alice's isogeny computation:", T3-T2);

XPushB = [[PA[0], G(1)], [QA[0], G(1)], [RA[0], G(1)]];

[B1, XPushB] = msidh_action_Bob(Apush, XPushB, skB, degree_bound, QB[0], PB[0], QB[1], PB[1]);

xphi_BPA = XPushB[0]; 

xphi_BQA = XPushB[1]; 

xphi_BRA = XPushB[2]; 


EB = EllipticCurve(G,[0,B1,0,1,0]);
phi_BPA = EB([xphi_BPA, sqrFp2(xphi_BPA^3+B1*xphi_BPA^2+xphi_BPA)]);
phi_BQA = EB([xphi_BQA, sqrFp2(xphi_BQA^3+B1*xphi_BQA^2+xphi_BQA)]);
if (phi_BPA[0] + phi_BQA[0] + xphi_BRA + B1) * (phi_BPA[0] - phi_BQA[0])^2 == (phi_BPA[1] - phi_BQA[1])^2:
  phi_BQA = - phi_BQA; 

T4 = time.time();

print("Bob's isogeny computation:", T4-T3);

UAi=[];
random.seed(int(2022));
UA=generate_U(UAi,B1,ind, EB, prime_A, NB);
random.seed(int(2023));
VA=generate_V(B1,ind,UA,UAi, EB, prime_A, NA);
VA = 4*NB*VA;

T5 = time.time();
print("Bob's torsion basis generation:", T5-T4);

h1 = miller_affine(phi_BQA,UA,string_A,NA,B1);
h2 = miller_affine(phi_BQA,VA,string_A,NA,B1);
h3 = miller_affine(phi_BPA,UA,string_A,NA,B1);
h4 = miller_affine(phi_BPA,VA,string_A,NA,B1);

T6 = time.time();
  
print("Bob's pairing computation:", T6 - T5);
  
hh1 = [];
hh2 = [];
hh3 = [];
hh4 = [];
batch_cof(hh1,h1,ind,len(ind), prime_A);
batch_cof(hh2,h2,ind,len(ind), prime_A);
batch_cof(hh3,h3,ind,len(ind), prime_A);
batch_cof(hh4,h4,ind,len(ind), prime_A);

[SA, lA] = DLC_new(ind, hh0, hh1, hh2, hh3, hh4, prime_A);

T7 = time.time();

print("Bob's discrete logarithm computation", T7 - T6);
 
UBi=[];
random.seed(int(2022));
UB=generate_U(UBi,A1,ind, EA, prime_B, NA);
random.seed(int(2023));
VB=generate_V(A1,ind,UB,UBi, EA, prime_B, NB);
VB=4*NA*VB;

T8 = time.time();

print("Alice's torsion basis generation:", T8 - T7);

h6 = miller_affine(phi_AQB,UB,string_B,NB,A1);
h7 = miller_affine(phi_AQB,VB,string_B,NB,A1);
h8 = miller_affine(phi_APB,UB,string_B,NB,A1);
h9 = miller_affine(phi_APB,VB,string_B,NB,A1);


T9 = time.time();
print("Alice's pairing computation:", T9 - T8);

hh6 = [];
hh7 = [];
hh8 = [];
hh9 = [];
batch_cof(hh6,h6,ind,len(ind),prime_B);
batch_cof(hh7,h7,ind,len(ind), prime_B);
batch_cof(hh8,h8,ind,len(ind), prime_B);
batch_cof(hh9,h9,ind,len(ind), prime_B);

[SB, lB] = DLC_new(ind, hh5, hh6, hh7, hh8, hh9, prime_B);
    
T10 = time.time();
print("Alice's discrete logarithm computation", T10 - T9);

# Alice generates the same torsion basis with Bob.
 
UAi=[];
random.seed(int(2022));
UA=generate_U(UAi,B1,ind, EB, prime_A, NB);
random.seed(int(2023));
VA=generate_V(B1,ind,UA,UAi, EB, prime_A, NA);
VA = 4*NB*VA;

ss4 = [];
for i in range(286):
  if mod(lA, 2) == 1:
    ss4.append(1);
  else:
    ss4.append(0);
  lA = floor(lA/2);
SS4 = chinese_remainder(ss4, prime_A);
[SS3, SS2, SS1] = SA;
AB  = msidh_action_Alice_last(B1, skA, degree_bound, UA, VA, SS4, SS2, SS3, SS1);        
T11 = time.time();

print("Alice's Key agreement", T11 - T10);

# Bob generates the same torsion basis with Alice.

UBi=[];
random.seed(int(2022));
UB=generate_U(UBi,A1,ind, EA, prime_B, NA);
random.seed(int(2023));
VB=generate_V(A1,ind,UB,UBi, EA, prime_B, NB);
VB=4*NA*VB;

ss4 = [];
for i in range(286):
  if mod(lB, 2) == 1:
    ss4.append(1);
  else:
    ss4.append(0);
  lB = floor(lB/2);
SS4 = chinese_remainder(ss4, prime_B);
[SS3, SS2, SS1] = SB;       
BA  = msidh_action_Bob_last(A1, skB, degree_bound, UB, VB, SS4, SS2, SS3, SS1);        

T12 = time.time();
print("Bob's Key agreement", T12 - T11);

if J(AB) == J(BA):
  print("The key exchange is successful.");
else:
  print("The key exchange is unsuccessful.");




  

















