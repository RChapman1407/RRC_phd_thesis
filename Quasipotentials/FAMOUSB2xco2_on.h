//three box model parameters FAMOUS B 2xCO2

#define VN 0.3683e+17
#define VT 0.5418e+17
#define VS 0.6097e+17
#define VIP 1.4860e+17
#define VB 9.9250e+17

#define SN0 0.034912
#define ST0 0.035425
#define SS0 0.034427
#define SIP0 0.034668
#define SB0 0.034538

#define FN 0.486e+6
#define FT -0.997e+6
#define FS 1.265e+6
#define FIP -0.754e+6

#define Alpha 0.12
#define Beta 790.0
#define So 0.0350
#define TS 7.919
#define To 3.870
#define eta 33.264e+6

#define KN 1.762e+6
#define KS 1.872e+6
#define KIP 99.977e+6
#define Lambda 1.62e+7
#define Gamma 0.36
#define Mu 22.0e-8

// parameter smoothing Heaviside
#define Xi 1.0e+2

// domain to scan and initial point

#define XL 35.0-3.0
#define XH 35.0+4.0
#define YL 35.0-0.0
#define YH 35.0+10.0

// QP from here

#define XINIT 35.3722
#define YINIT 36.3670

// minimum action path to here

#define Xshoot 33.9;
#define Yshoot 43.3;

#define SALT (VN*SN0 + VT*ST0 + VS*SS0 + VIP*SIP0 + VB*SB0)

#define Y 3.15e7 //seconds per year
