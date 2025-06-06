//
//	3D_Pattern_LTV.cpp : 이 파일에는 'main' 함수가 포함됩니다. 거기서 프로그램 실행이 시작되고 종료됩니다.
//
#include <stdio.h>
#include "MathCore/QuadProg.h"
#include <chrono>
using namespace std;

//#define _SAGITTAL_ONLY

#define dT					double(0.05)	//	Sampling time [s]
#define TIME_DSP			double(0.2)		//	Time duration for Double Support Phase : dT의 정수배
#define TIME_SSP			double(0.8)		//	Time duration for Single Support Phase : dT의 정수배
#define TIME_S2W			double(1.0)		//	Time duration for Stance to Walk : ~ TIME_PREVIEW
#define TIME_PREVIEW		double(1.5)		//	Preview Time : Approx 1.5 [s]

#define M_G					double(50.0)	//	Total mass [kg]
#define CoM_H				double(0.65)	//	CoM - ZMP : 현재는 constant
#define GRAVITY				double(9.81)	//	중력가속도
#define FOOT_WIDTH			double(0.08)	//	발 너비 ~ Y-axis ZMP boundary (Frontal plane)
#define FOOT_LENGTH			double(0.20)	//	발 길이 ~ X-axis ZMP boundary (Sagittal plane)
#define TOE_OUT				double(0.02)	//	Toe out distance [m]
#define HEEL2TOE_ZMP		double(0.08)	//	CoP/ZMP moving distance from heel to toe [m]
#define STEP_WIDTH			double(0.20)	//	양발 중심간의 거리 [m]
#define STEP_LENGTH			double(0.5)		//	스텝 길이 [m]
#define NO_STEPS			unsigned(20)	//	Number of steps >= 0
#define FIRST_FOOT			int(1)			//	첫 발 : 왼발 = -1, 오른발 = 1
#define WALKING_DIRECTION	int(1)			//	보행 방향 : 전진 = 1, 후진 = -1

#define BOUND_X_ZMP			double(0.5*FOOT_LENGTH)
#define BOUND_Y_ZMP			double(0.5*FOOT_WIDTH)

#ifdef _SAGITTAL_ONLY
#define INPUT_DOF			unsigned(2)
#else
#define INPUT_DOF			unsigned(4)		//	Input DoF : Currently (u_x, v_y, u_y, v_x)
#endif

#define STATE_LIN			unsigned(3)		//	State dimension for Linear Motion
#define STATE_ROT			unsigned(2)		//	State dimension for Rotational Motion

#define ALPHA_U				double(0.001)	//	Weight for input : u = CoM jerk
#define	ALPHA_V				double(0.001)	//	Weight for input : v = Hddot_G
#define ALPHA_ZMP			double(1000.)	//	Weight for ZMP : ZMP - ref_ZMP
#define ALPHA_POS			double(10.)		//	Weight for CoM position : CoM - ref_CoM
#define ALPHA_VEL			double(0.1)		//	Weight for CoM velocity : CoM_dot - ref_CoM_dot
#define ALPHA_THETA			double(10.0)	//	Weight for Angular momentum : H_G
#define ALPHA_OMEGA			double(1.0)		//	Weight for Rate of change of angular momentum : Hdot_G

constexpr auto DIM_SSP = unsigned(TIME_SSP / dT + 0.001);
constexpr auto DIM_DSP = unsigned(TIME_DSP / dT + 0.001);
constexpr auto DIM_S2W = unsigned(TIME_S2W / dT + 0.001);
constexpr auto DIM_STEP = DIM_SSP + DIM_DSP;
constexpr auto DIM_PREVIEW = unsigned(TIME_PREVIEW / dT + 0.001);	//	Preview length : Approx. 1.6 / dT
constexpr auto DIM_TOTAL = 2 * DIM_S2W + (NO_STEPS + 1)*DIM_STEP + DIM_DSP;// +DIM_PREVIEW;
constexpr auto MAX_ARRAY = DIM_TOTAL + DIM_PREVIEW;

constexpr auto DIM_UNKNOWN = INPUT_DOF * DIM_PREVIEW;
constexpr auto ZMP_CONSTRAINT = 4 * DIM_PREVIEW;
constexpr auto FRICTION_CONSTRAINT = 4 * DIM_PREVIEW;
constexpr auto DIM_CONSTRAINT = ZMP_CONSTRAINT + FRICTION_CONSTRAINT;


CstMatrix <DIM_PREVIEW, STATE_LIN> P_ps, P_vs, P_as, P_zs;
CstMatrix <DIM_PREVIEW, STATE_ROT> H_ps, H_vs;
CstMatrix <DIM_PREVIEW, DIM_PREVIEW> P_pu, P_vu, P_au, P_zu, H_pu, H_vu;


struct ZMP_TRJ {
	double X_comp[MAX_ARRAY];
	double Y_comp[MAX_ARRAY];
	double Z_comp[MAX_ARRAY];
	double X_r[DIM_PREVIEW];
	double Y_r[DIM_PREVIEW];
	double X, Y;
};

struct VerticalCoM {
	double pos[MAX_ARRAY];
	double vel[MAX_ARRAY];
	double acc[MAX_ARRAY];
};

struct FOOT_POS {
	double X_pos[MAX_ARRAY];
	double Y_pos[MAX_ARRAY];
	double X_f[DIM_PREVIEW];
	double Y_f[DIM_PREVIEW];
};

template <unsigned nDim1, unsigned nDim2>
void ComputeDiscreteModel(CstMatrix<nDim1, nDim1>& A_lin, CstMatrix<nDim2, nDim2>& A_rot, CstVector<nDim1>& B_lin, CstVector<nDim2>& B_rot);

void ComputePredictionModel();
void Set_ZMP_TotalTime(ZMP_TRJ & ref_zmp, FOOT_POS & foot_pos);
void Set_CoM_Z_TotalTime(double * z, double * zdot, double * zddot);

void Get_Preview_ZMP_FootPos(const unsigned indx, ZMP_TRJ & zmp, FOOT_POS & foot);

template <unsigned Dim>
void Get_Preview_CoM_Z(const unsigned indx, const VerticalCoM & Z_G, const double * z_o,
				CstVector<Dim> & z_acc, CstMatrix<Dim, Dim> & Omega, CstMatrix<Dim, Dim> & Theta);


int main()
{
	FILE *fout;
	double time(0.);
	double u_x(0.), v_y(0.), u_y(0.), v_x(0.);
	double gamma(0.), eta(0.);
	double C_f = 0.5;

	unsigned indx = 0;

	VerticalCoM z_CoM = { 0. };
	ZMP_TRJ Ref_ZMP = { 0. };
	FOOT_POS FootPos = { 0. };

	CstMatrix<STATE_LIN, STATE_LIN> A_lin;
	CstMatrix<STATE_ROT, STATE_ROT> A_rot;
	CstVector<STATE_LIN> B_lin, c_x, c_y, tmp_c_x, tmp_c_y;
	CstVector<STATE_ROT> B_rot, zeta_x, zeta_y, tmp_zeta_x, tmp_zeta_y;

	CstMatrix<DIM_PREVIEW, DIM_PREVIEW> C_THETA, S_THETA;
	CstMatrix<DIM_PREVIEW, DIM_PREVIEW> Omega, Theta, E_n;
	CstMatrix<DIM_PREVIEW, DIM_UNKNOWN> S_ux, S_uy, S_vx, S_vy, U_x, U_y, Tmp1, Tmp2;
	CstMatrix<DIM_PREVIEW, DIM_UNKNOWN> A_pos_X, A_vel_X, A_zmp_X, A_theta_X, A_omega_X;
	CstMatrix<DIM_PREVIEW, DIM_UNKNOWN> A_pos_Y, A_vel_Y, A_zmp_Y, A_theta_Y, A_omega_Y;
	CstVector<DIM_PREVIEW> bias_pos_X, bias_vel_X, bias_zmp_X, bias_theta_X, bias_omega_X;
	CstVector<DIM_PREVIEW> bias_pos_Y, bias_vel_Y, bias_zmp_Y, bias_theta_Y, bias_omega_Y;
	CstVector<DIM_PREVIEW> zmpXupper(BOUND_X_ZMP), zmpXlower(BOUND_X_ZMP), zmpYupper(BOUND_Y_ZMP), zmpYlower(BOUND_X_ZMP);
	CstVector<DIM_PREVIEW> tmp_X, tmp_Y, temp1, temp2, temp3, temp4, z_acc;

	CstMatrix<ZMP_CONSTRAINT, DIM_UNKNOWN> K_zmp;
	CstVector<ZMP_CONSTRAINT> d_zmp;
	CstMatrix<FRICTION_CONSTRAINT, DIM_UNKNOWN> K_fric;
	CstVector<FRICTION_CONSTRAINT> d_fric;

	CstMatrix<DIM_UNKNOWN, DIM_UNKNOWN> Q_mat, G_mat;
	CstVector<DIM_UNKNOWN> g_vec, u_vec;
	CstMatrix<DIM_CONSTRAINT, DIM_UNKNOWN> CI_mat;
	CstVector<DIM_CONSTRAINT> ci_vec;

	printf("SSP Dimension = %d\t SSP Time = %lf[s]\n", int(DIM_SSP), DIM_SSP*dT);
	printf("DSP Dimension = %d\t DSP Time = %lf[s]\n", int(DIM_DSP), DIM_DSP*dT);
	printf("S2W Dimension = %d\t S2W Time = %lf[s]\n", int(DIM_S2W), DIM_S2W*dT);
	printf("Total Dimension = %d\t Total Exec Time = %lf[s]\n", int(DIM_TOTAL), DIM_TOTAL*dT);
	printf("Maximum Array = %d\n", int(MAX_ARRAY));

	errno_t err = fopen_s(&fout, "Output.txt", "w");


	ComputeDiscreteModel(A_lin, A_rot, B_lin, B_rot);
	ComputePredictionModel();

	Set_ZMP_TotalTime(Ref_ZMP, FootPos);
	Set_CoM_Z_TotalTime(z_CoM.pos, z_CoM.vel, z_CoM.acc);

	E_n.setIdentity();
	C_THETA.setIdentity();
	S_THETA.setZero();

	//	Inputs : u = [ u_x ; v_y ; u_y ; v_x ]
	S_ux.setBlockMatrix(0, 0 * DIM_PREVIEW, E_n);
	S_vy.setBlockMatrix(0, 1 * DIM_PREVIEW, E_n);
	S_uy.setBlockMatrix(0, 2 * DIM_PREVIEW, E_n);
	S_vx.setBlockMatrix(0, 3 * DIM_PREVIEW, E_n);

	//	Sagittal Plane Cost
	A_pos_X = P_pu * S_ux;
	A_vel_X = P_vu * S_ux;
	A_theta_Y = H_pu * S_vy;
	A_omega_Y = H_vu * S_vy;

	//	Frontal Plane Cost
	A_pos_Y = P_pu * S_uy;
	A_vel_Y = P_vu * S_uy;
	A_theta_X = H_pu * S_vx;
	A_omega_X = H_vu * S_vx;

	//===== Compute Time Invariant Part of Quadratic Terms in Cost Function =====

	//	Sagittal Plane
	Q_mat = ALPHA_U * Transpose(S_ux)*S_ux + ALPHA_V * Transpose(S_vy)*S_vy;
	Q_mat += ALPHA_POS * Transpose(A_pos_X)*A_pos_X + ALPHA_VEL * Transpose(A_vel_X)*A_vel_X;
	Q_mat += ALPHA_THETA * Transpose(A_theta_Y)*A_theta_Y + ALPHA_OMEGA * Transpose(A_omega_Y)*A_omega_Y;

	//	Frontal Plane
	Q_mat += ALPHA_U * Transpose(S_uy)*S_uy + ALPHA_V * Transpose(S_vx)*S_vx;
	Q_mat += ALPHA_POS * Transpose(A_pos_Y)*A_pos_Y + ALPHA_VEL * Transpose(A_vel_Y)*A_vel_Y;
	Q_mat += ALPHA_THETA * Transpose(A_theta_X)*A_theta_X + ALPHA_OMEGA * Transpose(A_omega_X)*A_omega_X;

	//	Inequlity Constraint
	d_zmp.setVector(0, zmpXupper);
	d_zmp.setVector(DIM_PREVIEW, zmpXlower);
	d_zmp.setVector(2 * DIM_PREVIEW, zmpYupper);
	d_zmp.setVector(3 * DIM_PREVIEW, zmpYlower);


	do {
		//	Time Stamp Routine
		auto start_time = chrono::steady_clock::now();

		time = double(indx) * dT;

		gamma = (z_CoM.pos[indx] - Ref_ZMP.Z_comp[indx]) / (z_CoM.acc[indx] + GRAVITY);
		eta = 1.0 / (M_G*(z_CoM.acc[indx] + GRAVITY));

		Get_Preview_ZMP_FootPos(indx, Ref_ZMP, FootPos);
		Get_Preview_CoM_Z(indx, z_CoM, Ref_ZMP.Z_comp, z_acc, Omega, Theta);

		Ref_ZMP.X = c_x[0] - gamma * c_x[2] - eta * zeta_y[1];
		Ref_ZMP.Y = c_y[0] - gamma * c_y[2] + eta * zeta_x[1];

		//	Cost Function에 있는 Time Varying Matrix Term
		P_zs = P_ps - Omega * P_as;
		P_zu = P_pu - Omega * P_au;
		U_x = P_zu * S_ux - Theta * H_vu*S_vy;
		U_y = P_zu * S_uy + Theta * H_vu*S_vx;

		G_mat = Q_mat + ALPHA_ZMP * Transpose(U_x)*U_x + ALPHA_ZMP * Transpose(U_y)*U_y;

		//	Cost Function에 있는 Time Varying Linear Term
		bias_pos_X = P_ps * c_x - Ref_ZMP.X_r;
		bias_vel_X = P_vs * c_x;
		bias_theta_Y = H_ps * zeta_y;
		bias_omega_Y = H_vs * zeta_y;

		bias_pos_Y = P_ps * c_y - Ref_ZMP.Y_r;
		bias_vel_Y = P_vs * c_y;
		bias_theta_X = H_ps * zeta_x;
		bias_omega_X = H_vs * zeta_x;

		tmp_X = P_zs * c_x - Theta * bias_omega_Y;
		tmp_Y = P_zs * c_y + Theta * bias_omega_X;

		bias_zmp_X = tmp_X - Ref_ZMP.X_r;
		bias_zmp_Y = tmp_Y - Ref_ZMP.Y_r;

		g_vec = ALPHA_POS * bias_pos_X*A_pos_X + ALPHA_VEL * bias_vel_X*A_vel_X;
		g_vec += ALPHA_THETA * bias_theta_Y*A_theta_Y + ALPHA_OMEGA * bias_omega_Y*A_omega_Y;
		g_vec += ALPHA_ZMP * bias_zmp_X*U_x;

		g_vec += ALPHA_POS * bias_pos_Y*A_pos_Y + ALPHA_VEL * bias_vel_Y*A_vel_Y;
		g_vec += ALPHA_THETA * bias_theta_X*A_theta_X + ALPHA_OMEGA * bias_omega_X*A_omega_X;
		g_vec += ALPHA_ZMP * bias_zmp_Y*U_y;


		//	Set Inequlity Constraint : ZMP
		Tmp1 = C_THETA * U_x + S_THETA * U_y;
		Tmp2 = S_THETA * U_x - C_THETA * U_y;

		temp1 = C_THETA * (tmp_X - FootPos.X_f) + S_THETA * (tmp_Y - FootPos.Y_f);
		temp2 = S_THETA * (tmp_X - FootPos.X_f) - C_THETA * (tmp_Y - FootPos.Y_f);

		K_zmp.setBlockMatrix(0, 0, -Tmp1);
		K_zmp.setBlockMatrix(DIM_PREVIEW, 0, Tmp1);
		K_zmp.setBlockMatrix(2 * DIM_PREVIEW, 0,  Tmp2);
		K_zmp.setBlockMatrix(3 * DIM_PREVIEW, 0, -Tmp2);

		d_zmp.setVector(0,               zmpXupper - temp1);
		d_zmp.setVector(DIM_PREVIEW,     zmpXlower + temp1);
		d_zmp.setVector(2 * DIM_PREVIEW, zmpYupper + temp2);
		d_zmp.setVector(3 * DIM_PREVIEW, zmpYlower - temp2);

		//	Set Inequality Constraint : Friction
		Tmp1 = P_au * S_ux;
		Tmp2 = P_au * S_uy;

		temp1 = C_f * z_acc - P_as * c_x;
		temp2 = C_f * z_acc + P_as * c_x;
		temp3 = C_f * z_acc - P_as * c_y;
		temp4 = C_f * z_acc + P_as * c_y;

		K_fric.setBlockMatrix(0, 0, -Tmp1);
		K_fric.setBlockMatrix(DIM_PREVIEW, 0, Tmp1);
		K_fric.setBlockMatrix(2 * DIM_PREVIEW, 0, -Tmp2);
		K_fric.setBlockMatrix(3 * DIM_PREVIEW, 0,  Tmp2);

		d_fric.setVector(0, temp1);
		d_fric.setVector(DIM_PREVIEW, temp2);
		d_fric.setVector(2 * DIM_PREVIEW, temp3);
		d_fric.setVector(3 * DIM_PREVIEW, temp4);

		//	Set Inequality Constraint
		CI_mat.setBlockMatrix(0, 0, K_zmp);
		CI_mat.setBlockMatrix(ZMP_CONSTRAINT, 0, K_fric);
		ci_vec.setVector(0, d_zmp);
		ci_vec.setVector(ZMP_CONSTRAINT, d_fric);

		//	QP Solver : CASE-2 (LTV + Inequality Constraint)
		QuadProgPP::Solve_QuadProg(G_mat, g_vec, CI_mat, ci_vec, u_vec);

		u_x = u_vec[0];
		v_y = u_vec[DIM_PREVIEW];
		u_y = u_vec[2 * DIM_PREVIEW];
		v_x = u_vec[3 * DIM_PREVIEW];

		tmp_c_x = A_lin * c_x + B_lin * u_x;
		tmp_c_y = A_lin * c_y + B_lin * u_y;
		tmp_zeta_x = A_rot * zeta_x + B_rot * v_x;
		tmp_zeta_y = A_rot * zeta_y + B_rot * v_y;

		chrono::duration<double> elapsed = chrono::steady_clock::now() - start_time;

		fprintf_s(fout, "%lf\t %lf\t %lf\t %lf\t %lf\t", time, Ref_ZMP.X_comp[indx], Ref_ZMP.Y_comp[indx], Ref_ZMP.X, Ref_ZMP.Y);
		fprintf_s(fout, "%lf\t %lf\t %lf\t", c_x[0], c_x[1], c_x[2]);
		fprintf_s(fout, "%lf\t %lf\t %lf\t", c_y[0], c_y[1], c_y[2]);
		fprintf_s(fout, "%lf\t %lf\t %lf\t", z_CoM.pos[indx], z_CoM.vel[indx], z_CoM.acc[indx]);
		fprintf_s(fout, "%lf\t %lf\t %lf\t %lf\t", zeta_x[0], zeta_x[1], zeta_y[0], zeta_y[1]);
		fprintf_s(fout, "%lf\t %lf\t", FootPos.X_pos[indx], FootPos.Y_pos[indx]);
		fprintf_s(fout, "%lf\n", elapsed.count()*1000.);

		//	Update States : Linear & Angular momentum...
		c_x = tmp_c_x;
		c_y = tmp_c_y;
		zeta_x = tmp_zeta_x;
		zeta_y = tmp_zeta_y;

	} while (++indx <= DIM_TOTAL);

	fclose(fout);

	WaitForEnterKey();
}

//====================== End of Main() =========================


//
//	Compute Discrete Model for Linear/Angular Parts
//
template <unsigned nDim1, unsigned nDim2>
void ComputeDiscreteModel(CstMatrix<nDim1, nDim1>& A_lin, CstMatrix<nDim2, nDim2>& A_rot, CstVector<nDim1>& B_lin, CstVector<nDim2>& B_rot)
{
	//	Linear momentum (CoM Motion)
	A_lin[0][0] = 1.;				A_lin[0][1] = dT;				A_lin[0][2] = SQR(dT) / 2.;
	A_lin[1][0] = 0.;				A_lin[1][1] = 1.;				A_lin[1][2] = dT;
	A_lin[2][0] = 0.;				A_lin[2][1] = 0.;				A_lin[2][2] = 1.;
	B_lin[0] = CUBE(dT) / 6.;		B_lin[1] = SQR(dT) / 2.;		B_lin[2] = dT;

	//	Change of the rate of angular momentum
	A_rot[0][0] = 1.;				A_rot[0][1] = dT;
	A_rot[1][0] = 0.;				A_rot[1][1] = 1.;
	B_rot[0] = SQR(dT) / 2.;		B_rot[1] = dT;
}


//
//	Compute Prediction Model for Linear/Angular Parts
//
void ComputePredictionModel()
{
	unsigned	i, j;
	double		N;

	P_pu.setZero();
	P_vu.setZero();
	P_au.setZero();

	H_pu.setZero();
	H_vu.setZero();

	for (i = 0; i < DIM_PREVIEW; i++) {
		N = double(i) + 1.0;

		P_ps[i][0] = 1.0;		P_ps[i][1] = N * dT;		P_ps[i][2] = SQR(N * dT) / 2.0;
		P_vs[i][0] = 0.0;		P_vs[i][1] = 1.0;			P_vs[i][2] = N * dT;
		P_as[i][0] = 0.0;		P_as[i][1] = 0.0;			P_as[i][2] = 1.0;

		H_ps[i][0] = 1.0;		H_ps[i][1] = N * dT;
		H_vs[i][0] = 0.0;		H_vs[i][1] = 1.0;

		for (j = 0; j <= i; j++) {
			N = i - j;
			P_pu[i][j] = (3. * N * N + 3. * N + 1.)*CUBE(dT) / 6.0;
			P_vu[i][j] = (2. * N + 1.)*SQR(dT) / 2.0;
			P_au[i][j] = dT;

			H_pu[i][j] = (2. * N + 1.)*SQR(dT) / 2.0;
			H_vu[i][j] = dT;
		}
	}
}

//
//	Set Full Length ZMP Profile
//
void Set_ZMP_TotalTime(ZMP_TRJ & zmp, FOOT_POS & foot)
{
	double sign;
	double Foot_Info = 0.5 * FIRST_FOOT * STEP_WIDTH;
	unsigned i, j, jj;
	unsigned nStep;
	unsigned DIM_OFFSET = DIM_S2W + DIM_STEP;

	//	CASE-0 : 대부분 0~ .. 
	for (i = 0; i < DIM_S2W; i++) {
		zmp.X_comp[i] = zmp.Y_comp[i] = zmp.Z_comp[i] = 0.0;
		foot.X_pos[i] = foot.Y_pos[i] = 0.;
	}

	//	CASE-1 : 초기 Stand to Walking Phase
	for (i = DIM_S2W; i < DIM_OFFSET; i++) {
		j = i - DIM_S2W;

		if (j < DIM_DSP) {
			zmp.X_comp[i] = 0.0;
			zmp.Y_comp[i] = 0.5 * FIRST_FOOT * (STEP_WIDTH - TOE_OUT) * j / DIM_DSP;
			zmp.Z_comp[i] = 0.0;

			foot.X_pos[i] = 0.0;
			foot.Y_pos[i] = 0.5 * FIRST_FOOT * STEP_WIDTH * j / DIM_DSP;
		}
		else {
			jj = j - DIM_DSP;
			zmp.X_comp[i] = 0.5 * WALKING_DIRECTION * HEEL2TOE_ZMP * double(jj) / DIM_SSP;
			zmp.Y_comp[i] = FIRST_FOOT * (0.5 * (STEP_WIDTH - TOE_OUT) + TOE_OUT * jj / DIM_SSP);
			zmp.Z_comp[i] = 0.0;

			foot.X_pos[i] = 0.0;
			foot.Y_pos[i] = 0.5 * FIRST_FOOT * STEP_WIDTH;
		}
	}

	//	CASE-2 : Steady Walking Phase
	for (i = DIM_OFFSET; i < (DIM_S2W + NO_STEPS * DIM_STEP + DIM_DSP); i++) {
		j = i - DIM_OFFSET;

		nStep = j / DIM_STEP;			//	몇 번째 step ?
		sign = (nStep & 1) ? 1. : -1.;	//	홀/짝 구분; 홀수 = 1, 짝수 = -1

		jj = j - nStep * DIM_STEP;

		if (jj < DIM_DSP) {
			zmp.X_comp[i] = WALKING_DIRECTION * (0.5*HEEL2TOE_ZMP + nStep * STEP_LENGTH + (STEP_LENGTH - HEEL2TOE_ZMP)*jj / DIM_DSP);
			zmp.Y_comp[i] = sign * FIRST_FOOT * (STEP_WIDTH * double(jj) / DIM_DSP - 0.5 * (STEP_WIDTH + TOE_OUT));
			zmp.Z_comp[i] = 0.0;

			foot.X_pos[i] = WALKING_DIRECTION * (nStep * STEP_LENGTH + STEP_LENGTH * double(jj) / DIM_DSP);
			foot.Y_pos[i] = sign * FIRST_FOOT * STEP_WIDTH * (double(jj) / DIM_DSP - 0.5);
		}
		else {
			zmp.X_comp[i] = WALKING_DIRECTION * ((nStep + 1.0)*STEP_LENGTH - 0.5*HEEL2TOE_ZMP + HEEL2TOE_ZMP * (jj - DIM_DSP) / DIM_SSP);
			zmp.Y_comp[i] = sign * FIRST_FOOT*(0.5*(STEP_WIDTH - TOE_OUT) + (jj - DIM_DSP) * TOE_OUT / DIM_SSP);
			zmp.Z_comp[i] = 0.0;

			foot.X_pos[i] = WALKING_DIRECTION * (nStep + 1.0) * STEP_LENGTH;
			foot.Y_pos[i] = 0.5 * sign * FIRST_FOOT * STEP_WIDTH;
		}
	}

	//	CASE-3 : Walking to Stand Phase
	sign = (NO_STEPS & 1) ? -1. : 1.;
	for (i = (DIM_S2W + NO_STEPS * DIM_STEP + DIM_DSP); i < MAX_ARRAY; i++) {
		j = i - (DIM_S2W + NO_STEPS * DIM_STEP + DIM_DSP);

		if (j < DIM_SSP) {
			zmp.X_comp[i] = WALKING_DIRECTION * (NO_STEPS * STEP_LENGTH - 0.5*HEEL2TOE_ZMP + 0.5*HEEL2TOE_ZMP * j / DIM_SSP);
			zmp.Y_comp[i] = sign * FIRST_FOOT*(0.5*(STEP_WIDTH - TOE_OUT) + j * TOE_OUT / DIM_SSP);
			zmp.Z_comp[i] = 0.0;

			foot.X_pos[i] = WALKING_DIRECTION * double(NO_STEPS) * STEP_LENGTH;
			foot.Y_pos[i] = 0.5 * sign * FIRST_FOOT * STEP_WIDTH;
		}
		else if (j >= DIM_SSP && j < DIM_STEP) {
			jj = j - DIM_SSP;
			zmp.X_comp[i] = WALKING_DIRECTION * double(NO_STEPS) * STEP_LENGTH;
			zmp.Y_comp[i] = 0.5 * sign * FIRST_FOOT * (STEP_WIDTH + TOE_OUT) * (1.0 - double(jj) / DIM_DSP);
			zmp.Z_comp[i] = 0.0;

			foot.X_pos[i] = zmp.X_comp[i];
			foot.Y_pos[i] = 0.5 * sign * FIRST_FOOT * STEP_WIDTH * (1.0 - double(jj) / DIM_DSP);
		}
		else {
			zmp.X_comp[i] = WALKING_DIRECTION * double(NO_STEPS) * STEP_LENGTH;
			zmp.Y_comp[i] = 0.0;
			zmp.Z_comp[i] = 0.0;

			foot.X_pos[i] = zmp.X_comp[i];
			foot.Y_pos[i] = 0.0;
		}
	}
}


constexpr auto H_MAX = 0.7;
constexpr auto AMP = 0.04;
void Set_CoM_Z_TotalTime(double * pos, double * vel, double * acc)
{
	register unsigned i;
	double local_t(0.);
	double omega_n = 2. * M_PI / (TIME_DSP + TIME_SSP);
	unsigned DIM_SHIFT = unsigned(DIM_S2W + DIM_DSP + DIM_SSP / 2);

	for (i = 0; i < DIM_SHIFT; i++) {
		pos[i] = H_MAX;
		vel[i] = 0.;
		acc[i] = 0.;
	}

	for (i = DIM_SHIFT; i < NO_STEPS*DIM_STEP + DIM_SHIFT; i++) {
		local_t = (i - DIM_SHIFT)*dT;
		pos[i] = (H_MAX - 0.5*AMP) + 0.5*AMP*cos(omega_n*local_t);
		vel[i] = -0.5*AMP*omega_n*sin(omega_n*local_t);
		acc[i] = -0.5*AMP*SQR(omega_n)*cos(omega_n*local_t);
	}

	for (i = NO_STEPS*DIM_STEP + DIM_SHIFT; i < MAX_ARRAY; i++) {
		pos[i] = H_MAX;
		vel[i] = 0.;
		acc[i] = 0.;
	}
}


void Get_Preview_ZMP_FootPos(const unsigned k, ZMP_TRJ & zmp, FOOT_POS & p_f)
{
	memcpy(&zmp.X_r[0], &zmp.X_comp[k + 1], sizeof(double)*(DIM_PREVIEW));
	memcpy(&zmp.Y_r[0], &zmp.Y_comp[k + 1], sizeof(double)*(DIM_PREVIEW));

	memcpy(&p_f.X_f[0], &p_f.X_pos[k + 1], sizeof(double)*(DIM_PREVIEW));
	memcpy(&p_f.Y_f[0], &p_f.Y_pos[k + 1], sizeof(double)*(DIM_PREVIEW));
}



template<unsigned Dim>
void Get_Preview_CoM_Z(const unsigned k, const VerticalCoM & Z_G, const double * z_o,
		CstVector<Dim> & z_acc, CstMatrix<Dim, Dim>& Omega, CstMatrix<Dim, Dim>& Theta)
{
	unsigned i, j;

	for (i = 0; i < DIM_PREVIEW; i++) {
		j = k + i + 1;
		z_acc[i] = Z_G.acc[j] + GRAVITY;

		Omega[i][i] = (Z_G.pos[j] - z_o[j]) / z_acc[i];
		Theta[i][i] = 1.0 / (M_G*z_acc[i]);
	}
}
