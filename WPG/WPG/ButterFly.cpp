#include <stdio.h>
//#define _USE_MATH_DEFINES
//#include <math.h>
#include <iostream>

#include "LinAlgebra.h"
#include "stLUdecomp.h"
#include "stGaussJ.h"

const int	 No_Stride = 3;    //  No_Stide >= 1
const double T_s = 0.002;
const double T_dsp = 0.3;
const double T_ssp = 1.2;
const double T_step = T_dsp + T_ssp;
const double T_stride = 2.0 * T_step;
const double Gravity = 9.8;
const double l_G = 0.65;
const double omega_s = sqrt(Gravity / l_G);
const double omega_d = sqrt(Gravity / (0.9*l_G));


template <unsigned dim>
void SetBias(const double & a, const double & b, const double & c, const double & d, CstVector<dim> & ret);

int main()
{
	FILE *fp;
	double time = 0, t = 0, tt = 0;
	double p_ydot0, p_ydot1, p_ydot2, p_ydot3;
	double y_zmp, y_G, ydot_G, yddot_G = 0.;
	double z_zmp = 0., z_G, zdot_G, zddot_G;
	double C_D, S_D, C_S, S_S, C_DS, S_DS;
	unsigned step = 0;
	unsigned N_MAX = unsigned(No_Stride * T_stride / T_s);

	CstMatrix<2, 2> A_d, A_s, A_ds, E2;
	CstMatrix<2, 2> Tmp1, Inv;
	CstVector<2> b1, b2, b3, b4;
	CstVector<2> b_s, b_d, b_ds(0.);
	CstVector<2> y_IC0(0.), y_IC1(0.), y_IC2(0.), y_IC3(0.);
	CstVector<2> z_IC0(0.), z_IC1(0.);

	printf("Omega_s * T_ssp = %lf\n", omega_s*T_ssp);
	printf("Omega_d * T_dsp = %lf\n", omega_d*T_dsp);

	E2.setIdentity();

	C_D = cosh(omega_d*T_dsp);
	S_D = sinh(omega_d*T_dsp);
	C_S = cosh(omega_s*T_ssp);
	S_S = sinh(omega_s*T_ssp);
	C_DS = cosh(omega_d*T_dsp + omega_s*T_ssp);
	S_DS = sinh(omega_d*T_dsp + omega_s*T_ssp);

	A_d[0][0] = C_D;				A_d[0][1] = S_D / omega_d;
	A_d[1][0] = omega_d * S_D;		A_d[1][1] = C_D;

	A_s[0][0] = C_S;				A_s[0][1] = S_S / omega_s;
	A_s[1][0] = omega_s * S_S;		A_s[1][1] = C_S;

	A_ds = A_d * A_s;

	fopen_s(&fp, "output.txt", "wt");

	//	Y-Axis ZMP Setting
	const double p_y0 =  0.1;
	const double p_y1 =  0.1;
	const double p_y2 = -0.1;
	const double p_y3 = -0.1;
	const double p_y4 = p_y0;

	p_ydot0 = (p_y1 - p_y0) / T_ssp;
	p_ydot1 = (p_y2 - p_y1) / T_dsp;
	p_ydot2 = (p_y3 - p_y2) / T_ssp;
	p_ydot3 = (p_y4 - p_y3) / T_dsp;

	//	For Y-Axis CoM Motion
	SetBias(omega_s, T_ssp, p_ydot0, p_y0, b1);
	SetBias(omega_d, T_dsp, p_ydot1, p_y1, b2);
	SetBias(omega_s, T_ssp, p_ydot2, p_y2, b3);
	SetBias(omega_d, T_dsp, p_ydot3, p_y3, b4);

	//	y_IC0를 구해야 ~~
	Tmp1 = E2 - A_d * A_s * A_d * A_s;
	disp(Tmp1, "Tmp1 Matrix :");

	LU_dcmp<2> LU(Tmp1);

	Inv = LU.InvMat();
	printf("Determinant = %lf\t %lf\n", LU.det(), Tmp1[0][0]*Tmp1[1][1] - Tmp1[0][1]*Tmp1[1][0]);

	disp(Inv*Tmp1, "Verrify !");

	y_IC0 = Inv * (A_d*A_s*(A_d*b1 + b2) + A_d * b3 + b4);
	y_IC1 = A_s * y_IC0 + b1;
	y_IC2 = A_d * y_IC1 + b2;
	y_IC3 = A_s * y_IC2 + b3;

	//	For Z-axis CoM
	double p_s = Gravity / SQR(omega_s);
	double p_d = Gravity / SQR(omega_d);
	double pp, v_s, v_d, z_max = 0.;
	double gamma;

	v_s = 0;// (p_d - p_s + 0.025) / T_ssp;
	v_d = 0;// (p_s - p_d) / T_dsp;

	SetBias(omega_s, T_ssp, v_s, p_s, b1);
	SetBias(omega_d, T_dsp, v_d, p_d, b2);

	b3 = A_d * b1 + b2;

	Tmp1 = E2 - A_ds;
	Inv = Tmp1;
	Inverse_GJ(Inv);

	z_IC0 = Inv * b3;
	z_IC1 = A_s * z_IC0 + b1;

	do {
		time = step * T_s;
		t = time - T_stride * (int)(time / T_stride);

		if (t <= T_ssp) {
			tt = t;
			y_zmp = p_y0 + p_ydot0 * t;
			double c_t = cosh(omega_s*tt);
			double s_t = sinh(omega_s*tt);

			y_G = (y_IC0[0] - p_y0) * c_t + (y_IC0[1] - p_ydot0) * s_t / omega_s + y_zmp;
			ydot_G = omega_s * (y_IC0[0] - p_y0) * s_t + (y_IC0[1] - p_ydot0) * c_t + p_ydot0;

			z_G = (z_IC0[0] - p_s) * c_t + (z_IC0[1] - v_s) * s_t / omega_s + p_s + v_s * tt;
			zdot_G = omega_s * (z_IC0[0] - p_s) * s_t + (z_IC0[1] - v_s) * c_t + v_s;
			zddot_G = SQR(omega_s) * (z_IC0[0] - p_s) * c_t + omega_s * (z_IC0[1] - v_s) * s_t;
			pp = p_s + v_s * tt;
		}
		else if (t > T_ssp && t <= T_step) {
			tt = t - T_ssp;
			y_zmp = p_y1 + p_ydot1 * tt;
			double c_t = cosh(omega_d*tt);
			double s_t = sinh(omega_d*tt);

			y_G = (y_IC1[0] - p_y1) * c_t + (y_IC1[1] - p_ydot1) * s_t / omega_d + y_zmp;
			ydot_G = omega_d * (y_IC1[0] - p_y1) * s_t + (y_IC1[1] - p_ydot1) * c_t + p_ydot1;

			z_G = (z_IC1[0] - p_d) * c_t + (z_IC1[1] - v_d) * s_t / omega_d + p_d + v_d * tt;
			zdot_G = omega_d * (z_IC1[0] - p_d) * s_t + (z_IC1[1] - v_d) * c_t + v_d;
			zddot_G = SQR(omega_d)*(z_IC1[0] - p_d) * c_t + omega_d * (z_IC1[1] - v_d) * s_t;
			pp = p_d + v_d * tt;
		}
		else if (t > T_step && t <= T_step + T_ssp) {
			tt = t - T_step;
			y_zmp = p_y2 + p_ydot2 * tt;
			double c_t = cosh(omega_s*tt);
			double s_t = sinh(omega_s*tt);

			y_G = (y_IC2[0] - p_y2) * c_t + (y_IC2[1] - p_ydot2) * s_t / omega_s + y_zmp;
			ydot_G = omega_s * (y_IC2[0] - p_y2) * s_t + (y_IC2[1] - p_ydot2) * c_t + p_ydot2;

			z_G = (z_IC0[0] - p_s) * c_t + (z_IC0[1] - v_s) * s_t / omega_s + p_s + v_s * tt;
			zdot_G = omega_s * (z_IC0[0] - p_s) * s_t + (z_IC0[1] - v_s) * c_t + v_s;
			zddot_G = SQR(omega_s) * (z_IC0[0] - p_s) * c_t + omega_s * (z_IC0[1] - v_s) * s_t;
			pp = p_s + v_s * tt;
		}
		else if (t > T_step + T_ssp && t <= T_stride) {
			tt = t - T_step - T_ssp;
			y_zmp = p_y3 + p_ydot3 * tt;
			double c_t = cosh(omega_d*tt);
			double s_t = sinh(omega_d*tt);

			y_G = (y_IC3[0] - p_y3) * c_t + (y_IC3[1] - p_ydot3) * s_t / omega_d + y_zmp;
			ydot_G = omega_d * (y_IC3[0] - p_y3) * s_t + (y_IC3[1] - p_ydot3) * c_t + p_ydot3;

			z_G = (z_IC1[0] - p_d) * c_t + (z_IC1[1] - v_d) * s_t / omega_d + p_d + v_d * tt;
			zdot_G = omega_d * (z_IC1[0] - p_d) * s_t + (z_IC1[1] - v_d) * c_t + v_d;
			zddot_G = SQR(omega_d)*(z_IC1[0] - p_d) * c_t + omega_d * (z_IC1[1] - v_d) * s_t;
			pp = p_d + v_d * tt;
		}

		z_max = MAX(z_G, z_max);
		gamma = (zddot_G + Gravity) / (z_G - z_zmp);

		fprintf_s(fp, "%lf\t %lf\t %lf\t %lf\t %lf\t", time, y_zmp, y_G, ydot_G, yddot_G);
		fprintf_s(fp, "%lf\t %lf\t %lf\t %lf\t %lf\n", pp, (z_G - z_max + l_G), zdot_G, zddot_G, gamma);

	} while (++step <= N_MAX);

	fclose(fp);

	printf("N-MAX = %d\n", N_MAX);

	WaitForEnterKey();
}


template<unsigned dim>
void SetBias(const double & omega, const double & time, const double & vel, const double & zmp, CstVector<dim> & bias)
{
	double C_T = cosh(omega*time);
	double S_T = sinh(omega*time);

	bias[0] = vel * (time - S_T / omega) + zmp * (1 - C_T);
	bias[1] = vel * (1 - C_T) - omega * zmp * S_T;
}
