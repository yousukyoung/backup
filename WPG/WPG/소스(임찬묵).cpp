//
//	LMPC_WalkingPattern.cpp : �� ���Ͽ��� 'main' �Լ��� ���Ե˴ϴ�. �ű⼭ ���α׷� ������ ���۵ǰ� ����˴ϴ�.
//
#include "pch.h"
#include <math.h>
#include "QuadProg.h"
#include <chrono>
using namespace std::chrono;


#define dT					((double)0.05)		//	Sampling time
#define DOF					((unsigned)2)
#define STATE_DIM			((unsigned)3)		//	State dimension
#define STATE_ANG			((unsigned)2)		//	State angle

#define	TIME_DSP			((double)0.2)		//	Time duration for Double Support Phase : ����� dT
#define	TIME_SSP			((double)0.8)		//	Time duration for Single Support Phase : dT�� ������
#define TIME_S2W			TIME_SSP         	//	Time duration for Stance to Walk : ����� TIME_SSP
#define	NO_STEPS			((unsigned)20)		//	Number of steps >= 0
#define	FIRST_FOOT			((int)1)			//	ù �� : �޹� = -1, ������ = 1
#define WALKING_DIRECTION	((int)1)			//	���� ���� : ���� = 1, ���� = -1
#define TOE_OUT             ((double)10)        //  toe out ����

#define CoM_H               ((double)0.65)      //	CoM - ZMP : ����� constant
#define G_GRAVITY			((double)9.81)		//	�߷°��ӵ�
#define ALPHA_m             G_GRAVITY * 0.5
#define ALPHA_M             G_GRAVITY * 1.5
#define BETA_m              ((double)0.5)
#define BETA_M              ((double)0.8)

#define m_G                 ((double)50)        //	���� : �� 50kg?
#define GAMMA               (BETA_m / ALPHA_M + BETA_M / ALPHA_m) / 2          //  ���� ���
#define ETA                 (1 / m_G * ALPHA_M + 1 / m_G * ALPHA_m) / 2        //  ��Ÿ ���

#define FOOT_WIDTH			((double)0.1)		//	�ߺ��� �ʺ�
#define FOOT_LENGTH			((double)0.2)		//	�� ����
#define STEP_LENGTH			((double)0.65)		//	���� ����
#define STEP_WIDTH			((double)0.2)		//	��� �߽ɰ��� �Ÿ�

#define THETA_STEP_START    ((int)3)            //  ���° �������� ȸ���� �����ϴ°� (THETA_STEP_START <= NO_STEP)
#define THETA_STEP          ((int)0)            //  ȸ�� �ϴ� ������ ���� (THETA_STEP_START + THETA_STEP < NO_STEP)
#define THETA		    	((double)15)		//	�߸��� ȸ�� ����

#define ALPHA1				((double)0.0001) 	    //	Weight for input : u = CoM Jerk (X_three_dot)
#define	ALPHA2				((double)0.0001)    	//	Weight for input : u = CoM Jerk (Y_three_dot)
#define	ALPHA3				((double)0.001)  	  //	Weight for input : u = CoM Jerk (Z_three_dot)
#define ALPHA4				((double)0.0001)      //	Weight for angular momentum : H_x_two_dot
#define	ALPHA5				((double)0.0001)      //	Weight for angular momentum : H_y_two_dot
#define	ALPHA6				((double)0.0001)      //	Weight for angular momentum : H_z_two_dot
#define ALPHA7				((double)1000)  	//	Weight for x_ZMP : ZMP - ref_ZMP
#define ALPHA8				((double)1000)  	//	Weight for y_ZMP : ZMP - ref_ZMP
#define ALPHA9				((double)0)       	//	Weight for x_CoM : CoM - ref_CoM
#define	ALPHA10				((double)0)         //	Weight for y_CoM : CoM - ref_CoM
#define ALPHA11		     	((double)0)     	//	Weight for x_CoM_velocity : CoM_velocity - ref_CoM_velocity
#define	ALPHA12				((double)0)         //	Weight for y_CoM_velocity : CoM_velocity - ref_CoM_velocity
#define ALPHA13				((double)100)     	//	Weight for angular momentum : H_x
#define	ALPHA14				((double)100)     	//	Weight for angular momentum : H_y
#define ALPHA15				((double)10)     	//	Weight for angular momentum : H_x_dot
#define	ALPHA16				((double)10)	        //	Weight for angular momentum : H_y_dot
#define	ALPHA17				((double)50)	        //	Weight for height :  Z_p - Z_o - h


#define PREVIEW_LENGTH		((unsigned)30)		//	Preview length

#define	DIM_MAX				((unsigned)700)

#define FRICTION            ((double)0.1)     //  Low friction value
#define FRIC_START_TIME     ((double)8)         //  Low friction star time :' �Է½ð� + (PREVIEW_LENGTH * dT)'���� �����
#define FRIC_TIME_LENGTH    ((double)0)         //  Low friction time length

template <unsigned nRows, unsigned nCols>
void ComputeDiscreteModel(CstMatrix<nRows, nCols>& A_d, CstVector<nRows>& B_d, CstVector<nCols>& G_d);
template <unsigned nRows, unsigned nCols>
void ComputeDiscreteModel_2(CstMatrix<nRows, nCols>& Alpha_d, CstVector<nRows>& Delta_d);

void SetPredictiveSystemModel();

void Set_ZMP_Reference(const unsigned nSSP, const unsigned nDSP, const unsigned nEX, const unsigned nTotal, double* xZMP_ref, double* yZMP_ref, double* zZMP_ref, double* x_foot, double* y_foot);
void Get_ZMP_Reference(const unsigned index, const unsigned nDim, double* ZMP_ref, CstVector<PREVIEW_LENGTH>& preview_ZMP_ref);

void Predict_friction(double* C_f, double fric_start_time, double fric_time_length, double friction, const unsigned DIM_TOTAL);
void Get_Friction(const unsigned k, const unsigned DIM_TOTAL, double* friction, CstVector<PREVIEW_LENGTH>& Friction);

void Set_z_Bound(const unsigned DIM_SSP, const unsigned DIM_DSP, const unsigned DIM_EXTRA, const unsigned DIM_TOTAL, double* z_bound);
void Get_z_bound(const unsigned k, const unsigned DIM_TOTAL, double* z_Bound, CstVector<PREVIEW_LENGTH>& z_bound);

CstMatrix<PREVIEW_LENGTH, STATE_DIM> P_ps, P_vs, P_zs, P_as;
CstMatrix<PREVIEW_LENGTH, STATE_ANG> H_ps, H_vs;
CstMatrix<PREVIEW_LENGTH, PREVIEW_LENGTH> P_pu, P_vu, P_zu, P_au, H_pu, H_vu;

double gravity = G_GRAVITY;
double height = CoM_H;

int main()
{
	FILE* fp;
	FILE* fp_2;
	unsigned indx = 0;
	double time = 0, x_zmp = 0, y_zmp = 0;
	double xZMP_ref[DIM_MAX], yZMP_ref[DIM_MAX], zZMP_ref[DIM_MAX];
	double x_foot[4 * (NO_STEPS + 3)], y_foot[4 * (NO_STEPS + 3)];
	double C_fric[DIM_MAX];
	double z_Bound[DIM_MAX];

	CstMatrix<PREVIEW_LENGTH, PREVIEW_LENGTH > E_N, C_THETA, S_THETA, C_f_matrix;
	CstMatrix<PREVIEW_LENGTH * 4, PREVIEW_LENGTH * 2> THETA_matrix;
	CstMatrix<PREVIEW_LENGTH * 2, PREVIEW_LENGTH * 6> D_MM, D_mM, D_Mm, D_mm;
	CstMatrix<PREVIEW_LENGTH * 4, PREVIEW_LENGTH> K_z;
	CstMatrix<PREVIEW_LENGTH * 4, PREVIEW_LENGTH * 6> K_fric;
	CstMatrix<PREVIEW_LENGTH * 16, PREVIEW_LENGTH * 6> K_zmp;
	CstMatrix<PREVIEW_LENGTH * 6, PREVIEW_LENGTH * 6> E_Mat, Q_Mat, Q_Mat_x, Q_Mat_y, A_x, B_x, C_x, A_y, B_y, C_y;
	CstMatrix<PREVIEW_LENGTH * 24, PREVIEW_LENGTH * 6> CI;
	CstMatrix<PREVIEW_LENGTH * 6, PREVIEW_LENGTH * 24> Transpose_CI;
	CstMatrix<PREVIEW_LENGTH, PREVIEW_LENGTH * 6> S_ux, S_uy, S_uz, S_vx, S_vy, S_vz, U_x, U_y;

	CstVector<PREVIEW_LENGTH> preview_ZMP_ref_x, preview_ZMP_ref_y, preview_ZMP_ref_z, preview_CoM_ref_x, preview_CoM_ref_y, preview_CoM_velocity_ref_x, preview_CoM_velocity_ref_y;
	CstVector<PREVIEW_LENGTH> g_c, u_x, u_y, u_z, nu_x, nu_y, nu_z, C_f, alpha_M, alpha_m, beta_M, beta_m;
	CstVector<PREVIEW_LENGTH > l_M, l_m, w_M, w_m;
	CstVector<PREVIEW_LENGTH * 2> ci0_d_zmp_mm, ci0_d_zmp_mM, ci0_d_zmp_Mm, ci0_d_zmp_MM;
	CstVector<PREVIEW_LENGTH * 4> d_fric, d_z, Bound;
	CstVector<PREVIEW_LENGTH * 16> d_zmp;
	CstVector<PREVIEW_LENGTH * 6> g0, solution, E_x, F_x, G_x, E_y, F_y, G_y, I_z;
	CstVector<PREVIEW_LENGTH * 24> ci0;

	CstVector<PREVIEW_LENGTH> z_bound;

	//CstVector<PREVIEW_LENGTH> tmp_state_y, tmp_state_x, tmp_velocity_y, tmp_velocity_x, tmp_acceleration_y, tmp_acceleration_x;

	CstMatrix<STATE_DIM, STATE_DIM> A_d;
	CstVector<STATE_DIM> c_x, c_y, c_z, B_d, G_d, tmp_state_y, tmp_state_x, tmp_state_z;


	CstMatrix<STATE_ANG, STATE_ANG> Alpha_d;
	CstVector<STATE_ANG> Delta_d, zeta_x, zeta_y, zeta_z, tmp_angle_x, tmp_angle_y, tmp_angle_z;


	const unsigned DIM_SSP = (unsigned)round(TIME_SSP / dT);
	const unsigned DIM_DSP = (unsigned)round(TIME_DSP / dT);
	const unsigned DIM_EXTRA = (unsigned)round(3 * TIME_S2W / dT);
	const unsigned DIM_TOTAL = DIM_EXTRA + (NO_STEPS + 1) * (DIM_SSP + DIM_DSP) + DIM_DSP + 3 * PREVIEW_LENGTH;



	if (DIM_TOTAL > DIM_MAX) {
		ERROR("Total Dimension > Maximum Dimension !!");
	}

	printf("SSP Dimension = %d\t SSP Time = %lf[s]\n", DIM_SSP, DIM_SSP * dT);
	printf("DSP Dimension = %d\t DSP Time = %lf[s]\n", DIM_DSP, DIM_DSP * dT);
	printf("EXTRA Dimension = %d\t EXTRA Time = %lf[s]\n", DIM_EXTRA, DIM_EXTRA * dT);
	printf("Total Dimension = %d\t Total Exec Time = %lf[s]\n", DIM_TOTAL, DIM_TOTAL * dT);
	//cin.get();

	fopen_s(&fp, "Test9.txt", "w");
	fopen_s(&fp_2, "Foot.txt", "w");

	double yZMP_bound = FOOT_WIDTH / 2;
	double xZMP_bound = FOOT_LENGTH / 2;

	c_z[0] = CoM_H;
	//c_z[1] = 0;
	//c_z[2] = 0;


	Set_ZMP_Reference(DIM_SSP, DIM_DSP, DIM_EXTRA, DIM_TOTAL, xZMP_ref, yZMP_ref, zZMP_ref, x_foot, y_foot);


	SetPredictiveSystemModel();
	ComputeDiscreteModel_2(Alpha_d, Delta_d);

	Predict_friction(C_fric, FRIC_START_TIME, FRIC_TIME_LENGTH, FRICTION, DIM_TOTAL);

	//Set_z_Bound(DIM_SSP, DIM_DSP, DIM_EXTRA, DIM_TOTAL, z_Bound);

	do {
		auto start_time = high_resolution_clock::now();

		time = indx * dT;

		ComputeDiscreteModel(A_d, B_d, G_d);
		Get_ZMP_Reference(indx, DIM_TOTAL, yZMP_ref, preview_ZMP_ref_y);
		Get_ZMP_Reference(indx, DIM_TOTAL, xZMP_ref, preview_ZMP_ref_x);
		Get_ZMP_Reference(indx, DIM_TOTAL, zZMP_ref, preview_ZMP_ref_z);

		//Get_z_bound(indx, DIM_TOTAL, z_Bound, z_bound);

		Get_Friction(indx, DIM_TOTAL, C_fric, C_f);


		E_N.setDiagonal(1);

		E_Mat.setBlockMatrix(0, 0, ALPHA1 * E_N);
		E_Mat.setBlockMatrix(PREVIEW_LENGTH, PREVIEW_LENGTH, ALPHA2 * E_N);
		E_Mat.setBlockMatrix(PREVIEW_LENGTH * 2, PREVIEW_LENGTH * 2, ALPHA3 * E_N);
		E_Mat.setBlockMatrix(PREVIEW_LENGTH * 3, PREVIEW_LENGTH * 3, ALPHA4 * E_N);
		E_Mat.setBlockMatrix(PREVIEW_LENGTH * 4, PREVIEW_LENGTH * 4, ALPHA5 * E_N);
		E_Mat.setBlockMatrix(PREVIEW_LENGTH * 5, PREVIEW_LENGTH * 5, ALPHA6 * E_N);

		S_ux.setBlockMatrix(0, 0, E_N);
		S_uy.setBlockMatrix(0, PREVIEW_LENGTH, E_N);
		S_uz.setBlockMatrix(0, PREVIEW_LENGTH * 2, E_N);
		S_vx.setBlockMatrix(0, PREVIEW_LENGTH * 3, E_N);
		S_vy.setBlockMatrix(0, PREVIEW_LENGTH * 4, E_N);
		S_vz.setBlockMatrix(0, PREVIEW_LENGTH * 5, E_N);

		P_zs = P_ps - GAMMA * P_as;
		P_zu = P_pu - GAMMA * P_au;

		U_x = P_zu * S_ux - ETA * H_vu * S_vy;
		U_y = P_zu * S_uy + ETA * H_vu * S_vx;

		//	Q_Mat = E_Mat + ALPHA2 * Transpose_P_zu * P_zu;
		A_x = ALPHA7 * Transpose(U_x) * U_x;
		B_x = Transpose(S_ux) * (ALPHA11 * Transpose(P_vu) * P_vu + ALPHA9 * Transpose(P_pu) * P_pu) * S_ux;
		C_x = Transpose(S_vy) * (ALPHA14 * Transpose(H_pu) * H_pu + ALPHA16 * Transpose(H_vu) * H_vu) * S_vy;

		A_y = ALPHA8 * Transpose(U_y) * U_y;
		B_y = Transpose(S_uy) * (ALPHA12 * Transpose(P_vu) * P_vu + ALPHA10 * Transpose(P_pu) * P_pu) * S_uy;
		C_y = Transpose(S_vx) * (ALPHA13 * Transpose(H_pu) * H_pu + ALPHA15 * Transpose(H_vu) * H_vu) * S_vx;

		Q_Mat_x = E_Mat + A_x + B_x + C_x;
		Q_Mat_y = E_Mat + A_y + B_y + C_y;
		Q_Mat = E_Mat + A_x + B_x + C_x + A_y + B_y + C_y + ALPHA17 * Transpose(S_uz) * Transpose(P_pu) * P_pu * S_uz;
		//	LTI ������ ��� : �ѹ��� Ǯ�� ��!
		//QuadProg::Cholesky_Decomposition<PREVIEW_LENGTH * 4>(Q_Mat_x);
		//QuadProg::Cholesky_Decomposition<PREVIEW_LENGTH * 4>(Q_Mat_y);
		QuadProg::Cholesky_Decomposition<PREVIEW_LENGTH * 6>(Q_Mat);


		E_x = ALPHA7 * Transpose(U_x) * (P_zs * c_x - ETA * H_vs * zeta_y - preview_ZMP_ref_x);
		F_x = Transpose(S_ux) * (ALPHA11 * Transpose(P_vu) * (P_vs * c_x) + ALPHA9 * Transpose(P_pu) * (P_ps * c_x - preview_ZMP_ref_x));
		G_x = Transpose(S_vy) * (ALPHA14 * Transpose(H_pu) * (H_ps * zeta_y) + ALPHA16 * Transpose(H_vu) * (H_vs * zeta_y));

		E_y = ALPHA8 * Transpose(U_y) * (P_zs * c_y + ETA * H_vs * zeta_x - preview_ZMP_ref_y);
		F_y = Transpose(S_uy) * (ALPHA12 * Transpose(P_vu) * (P_vs * c_y) + ALPHA10 * Transpose(P_pu) * (P_ps * c_y - preview_ZMP_ref_y));
		G_y = Transpose(S_vx) * (ALPHA13 * Transpose(H_pu) * (H_ps * zeta_x) + ALPHA15 * Transpose(H_vu) * (H_vs * zeta_x));

		I_z = ALPHA17 * Transpose(S_uz) * Transpose(P_pu) * (P_ps * c_z - preview_ZMP_ref_z);

		//g0_x = E_x + F_x + G_x;
		//g0_y = E_y + F_y + G_y;
		g0 = E_x + F_x + G_x + E_y + F_y + G_y + I_z;


		//	Set Inequality Constraints

		C_THETA.setDiagonal(cos(THETA * M_PI / 180));
		S_THETA.setDiagonal(sin(THETA * M_PI / 180));

		alpha_M.setFill(ALPHA_M);
		alpha_m.setFill(ALPHA_m);
		beta_M.setFill(BETA_M);
		beta_m.setFill(BETA_m);


		THETA_matrix.setBlockMatrix(0, 0, C_THETA);
		THETA_matrix.setBlockMatrix(0, PREVIEW_LENGTH, S_THETA);
		THETA_matrix.setBlockMatrix(PREVIEW_LENGTH, 0, -C_THETA);
		THETA_matrix.setBlockMatrix(PREVIEW_LENGTH, PREVIEW_LENGTH, -S_THETA);
		THETA_matrix.setBlockMatrix(PREVIEW_LENGTH * 2, 0, -S_THETA);
		THETA_matrix.setBlockMatrix(PREVIEW_LENGTH * 2, PREVIEW_LENGTH, C_THETA);
		THETA_matrix.setBlockMatrix(PREVIEW_LENGTH * 3, 0, S_THETA);
		THETA_matrix.setBlockMatrix(PREVIEW_LENGTH * 3, PREVIEW_LENGTH, -C_THETA);

		D_mm.setBlockMatrix(0, 0, (P_pu - (BETA_m / ALPHA_M) * P_au) * S_ux - (1 / (m_G * ALPHA_M)) * H_vu * S_vy);
		D_mm.setBlockMatrix(PREVIEW_LENGTH, 0, (P_pu - (BETA_m / ALPHA_M) * P_au)* S_uy + (1 / (m_G * ALPHA_M)) * H_vu * S_vx);

		D_mM.setBlockMatrix(0, 0, (P_pu - (BETA_m / ALPHA_M) * P_au)* S_ux - (1 / (m_G * ALPHA_m)) * H_vu * S_vy);
		D_mM.setBlockMatrix(PREVIEW_LENGTH, 0, (P_pu - (BETA_m / ALPHA_M) * P_au)* S_uy + (1 / (m_G * ALPHA_m)) * H_vu * S_vx);

		D_Mm.setBlockMatrix(0, 0, (P_pu - (BETA_M / ALPHA_m) * P_au)* S_ux - (1 / (m_G * ALPHA_M)) * H_vu * S_vy);
		D_Mm.setBlockMatrix(PREVIEW_LENGTH, 0, (P_pu - (BETA_M / ALPHA_m) * P_au)* S_uy + (1 / (m_G * ALPHA_M)) * H_vu * S_vx);

		D_MM.setBlockMatrix(0, 0, (P_pu - (BETA_M / ALPHA_m) * P_au)* S_ux - (1 / (m_G * ALPHA_m)) * H_vu * S_vy);
		D_MM.setBlockMatrix(PREVIEW_LENGTH, 0, (P_pu - (BETA_M / ALPHA_m) * P_au)* S_uy + (1 / (m_G * ALPHA_m)) * H_vu * S_vx);


		K_zmp.setBlockMatrix(0, 0, THETA_matrix * D_mm);
		K_zmp.setBlockMatrix(PREVIEW_LENGTH * 4, 0, THETA_matrix* D_mM);
		K_zmp.setBlockMatrix(PREVIEW_LENGTH * 8, 0, THETA_matrix* D_Mm);
		K_zmp.setBlockMatrix(PREVIEW_LENGTH * 12, 0, THETA_matrix* D_MM);


		
		for (int i = 0; i < PREVIEW_LENGTH; i++) {
			C_f_matrix[i][i] = C_f[i];
		}

		K_fric.setBlockMatrix(0, 0, P_au * (S_ux - C_f_matrix * S_uz));
		K_fric.setBlockMatrix(PREVIEW_LENGTH, 0, -P_au * (S_ux + C_f_matrix * S_uz));
		K_fric.setBlockMatrix(PREVIEW_LENGTH * 2, 0, P_au * (S_uy - C_f_matrix * S_uz));
		K_fric.setBlockMatrix(PREVIEW_LENGTH * 3, 0, -P_au * (S_uy + C_f_matrix * S_uz));

		K_z.setBlockMatrix(0, 0, P_au);
		K_z.setBlockMatrix(PREVIEW_LENGTH, 0, -P_au);
		K_z.setBlockMatrix(PREVIEW_LENGTH * 2, 0, P_pu);
		K_z.setBlockMatrix(PREVIEW_LENGTH * 3, 0, -P_pu);

		CI.setBlockMatrix(0, 0, K_zmp);
		CI.setBlockMatrix(PREVIEW_LENGTH * 16, 0, K_fric);
		CI.setBlockMatrix(PREVIEW_LENGTH * 20, 0, K_z * S_uz);
		Transpose_CI = -Transpose(CI);


		l_M.setFill(xZMP_bound);
		l_m.setFill(xZMP_bound);
		w_M.setFill(yZMP_bound);
		w_m.setFill(yZMP_bound);

		Bound.setVector(0, l_M);
		Bound.setVector(PREVIEW_LENGTH, l_m);
		Bound.setVector(PREVIEW_LENGTH * 2, w_M);
		Bound.setVector(PREVIEW_LENGTH * 3, w_m);

		ci0_d_zmp_mm.setVector(0, preview_ZMP_ref_x - (P_ps - (BETA_m / ALPHA_M) * P_as) * c_x + (1 / (m_G * ALPHA_M)) * H_vs * zeta_y);
		ci0_d_zmp_mm.setVector(PREVIEW_LENGTH, preview_ZMP_ref_y - (P_ps - (BETA_m / ALPHA_M) * P_as) * c_y - (1 / (m_G * ALPHA_M)) * H_vs * zeta_x);
		ci0_d_zmp_mM.setVector(0, preview_ZMP_ref_x - (P_ps - (BETA_m / ALPHA_M) * P_as) * c_x + (1 / (m_G * ALPHA_m)) * H_vs * zeta_y);
		ci0_d_zmp_mM.setVector(PREVIEW_LENGTH, preview_ZMP_ref_y - (P_ps - (BETA_m / ALPHA_M) * P_as) * c_y - (1 / (m_G * ALPHA_m)) * H_vs * zeta_x);
		ci0_d_zmp_Mm.setVector(0, preview_ZMP_ref_x - (P_ps - (BETA_M / ALPHA_m) * P_as) * c_x + (1 / (m_G * ALPHA_M)) * H_vs * zeta_y);
		ci0_d_zmp_Mm.setVector(PREVIEW_LENGTH, preview_ZMP_ref_y - (P_ps - (BETA_M / ALPHA_m) * P_as) * c_y - (1 / (m_G * ALPHA_M)) * H_vs * zeta_x);
		ci0_d_zmp_MM.setVector(0, preview_ZMP_ref_x - (P_ps - (BETA_M / ALPHA_m) * P_as) * c_x + (1 / (m_G * ALPHA_m)) * H_vs * zeta_y);
		ci0_d_zmp_MM.setVector(PREVIEW_LENGTH, preview_ZMP_ref_y - (P_ps - (BETA_M / ALPHA_m) * P_as) * c_y - (1 / (m_G * ALPHA_m)) * H_vs * zeta_x);


		d_zmp.setVector(0, Bound + THETA_matrix * ci0_d_zmp_mm);
		d_zmp.setVector(PREVIEW_LENGTH * 4, Bound + THETA_matrix * ci0_d_zmp_mM);
		d_zmp.setVector(PREVIEW_LENGTH * 8, Bound + THETA_matrix * ci0_d_zmp_Mm);
		d_zmp.setVector(PREVIEW_LENGTH * 12, Bound + THETA_matrix * ci0_d_zmp_MM);

		g_c.setFill(G_GRAVITY);

		d_fric.setVector(0, C_f_matrix* (g_c + P_as * c_z) - P_as * c_x);
		d_fric.setVector(PREVIEW_LENGTH, C_f_matrix* (g_c + P_as * c_z) + P_as * c_x);
		d_fric.setVector(PREVIEW_LENGTH * 2, C_f_matrix* (g_c + P_as * c_z) - P_as * c_y);
		d_fric.setVector(PREVIEW_LENGTH * 3, C_f_matrix* (g_c + P_as * c_z) + P_as * c_y);

		d_z.setVector(0, alpha_M - P_as * c_z - g_c);
		d_z.setVector(PREVIEW_LENGTH, -alpha_m + P_as * c_z + g_c);
		d_z.setVector(PREVIEW_LENGTH * 2, beta_M - P_ps * c_z);
		d_z.setVector(PREVIEW_LENGTH * 3, -beta_m + P_ps * c_z);


		ci0.setVector(0, d_zmp);
		ci0.setVector(PREVIEW_LENGTH * 16, d_fric);
		ci0.setVector(PREVIEW_LENGTH * 20, d_z);



		//	Solve Quadratic Programming
		QuadProg::Solve_QuadProg<6 * PREVIEW_LENGTH, 0, 24 * PREVIEW_LENGTH>(Q_Mat, g0, Transpose_CI, ci0, solution);

		//	Compute state & ZMP

		for (int i = 0; i < PREVIEW_LENGTH; i++) {
			u_x[i] = solution[i];
			u_y[i] = solution[PREVIEW_LENGTH + i];
			u_z[i] = solution[PREVIEW_LENGTH * 2 + i];

			nu_x[i] = solution[PREVIEW_LENGTH * 3 + i];
			nu_y[i] = solution[PREVIEW_LENGTH * 4 + i];
			nu_z[i] = solution[PREVIEW_LENGTH * 5 + i];
		}

		tmp_state_x = A_d * c_x + B_d * u_x[0];
		tmp_state_y = A_d * c_y + B_d * u_y[0];
		tmp_state_z = A_d * c_z + B_d * u_z[0];


		tmp_angle_y = Alpha_d * zeta_y + Delta_d * nu_y[0];
		tmp_angle_x = Alpha_d * zeta_x + Delta_d * nu_x[0];

		x_zmp = G_d * c_x - zeta_y[1] / (m_G * G_GRAVITY);
		y_zmp = G_d * c_y + zeta_x[1] / (m_G * G_GRAVITY);


		auto end_time = high_resolution_clock::now();
		auto del_time = end_time - start_time;

		double elapsed = ((double)duration_cast<microseconds>(del_time).count()) / 1000;


		//		cout.precision(3);
		//		cout << "Time : " << std::fixed << elapsed << " [us] to run.\n";
		//		cout << "Time : " << std::fixed << (double)std::chrono::duration_cast<std::chrono::nanoseconds>(time).count() / 1000. << " [us] to run.\n";

		fprintf_s(fp, "%lf\t %lf\t %lf\t %lf\t %lf\t", time, xZMP_ref[indx], yZMP_ref[indx], x_zmp, y_zmp);
		fprintf_s(fp, "%lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t", c_x[0], c_y[0], c_x[1], c_y[1], c_x[2], c_y[2], zeta_x[0], zeta_y[0], zeta_x[1], zeta_y[1], ((double)elapsed));
		fprintf_s(fp, "%lf\t %lf\t %lf\t %lf\n", c_z[0], c_z[1], c_z[2], preview_ZMP_ref_z[0]);

		//	Update system state

		c_x = tmp_state_x;
		c_y = tmp_state_y;
		c_z = tmp_state_z;

		zeta_y = tmp_angle_y;
		zeta_x = tmp_angle_x;

		gravity = c_z[2] + G_GRAVITY;
		height = c_z[0];

	} while (++indx <= DIM_TOTAL);




	fclose(fp);

	if (NO_STEPS % 2 == 1) {
		for (int f = 0; f < 2 + NO_STEPS / 2; f++) {
			for (int g = 0; g < 4; g++) {
				fprintf_s(fp_2, "%lf\t %lf\t %lf\t %lf\n", x_foot[8 * f + g], y_foot[8 * f + g], x_foot[(8 * f) + 4 + g], y_foot[(8 * f) + 4 + g]);
			}
			fprintf_s(fp_2, "%lf\t %lf\t %lf\t %lf\n\n", x_foot[8 * f], y_foot[8 * f], x_foot[(8 * f) + 4], y_foot[(8 * f) + 4]);

		}
	}
	else {
		for (int f = 0; f < 1 + NO_STEPS / 2; f++) {
			for (int g = 0; g < 4; g++) {
				fprintf_s(fp_2, "%lf\t %lf\t %lf\t %lf\n", x_foot[8 * f + g], y_foot[8 * f + g], x_foot[(8 * f) + 4 + g], y_foot[(8 * f) + 4 + g]);
			}
			fprintf_s(fp_2, "%lf\t %lf\t %lf\t %lf\n\n", x_foot[8 * f], y_foot[8 * f], x_foot[(8 * f) + 4], y_foot[(8 * f) + 4]);

		}
		for (int h = 0; h < 4; h++) {
			if (FIRST_FOOT || 1) {
				fprintf_s(fp_2, "\t \t");
			}
			fprintf_s(fp_2, "%lf\t %lf\n", x_foot[4 * (NO_STEPS + 2) + h], y_foot[4 * (NO_STEPS + 2) + h]);
		}
		fprintf_s(fp_2, "%lf\t %lf\n\n", x_foot[4 * (NO_STEPS + 2)], y_foot[4 * (NO_STEPS + 2)]);
	} x_foot[4 * (NO_STEPS + 3)], y_foot[4 * (NO_STEPS + 3)];

	fclose(fp_2);

	return 0;
}


void Predict_friction(double* C_fric, double fric_start_time, double fric_time_length, double friction, unsigned DIM_TOTAL) {

	double time = fric_start_time / dT;
	double length = fric_time_length / dT;

	for (int i = 0; i <= DIM_TOTAL; i++) {
		if (i > time && i <= time + length) {
			C_fric[i] = friction;
		}
		else {
			C_fric[i] = 0.6;
		}
	}
}

//	Discrete System Model
template <unsigned nRows, unsigned nCols>
void ComputeDiscreteModel(CstMatrix<nRows, nCols>& A_d, CstVector<nRows>& B_d, CstVector<nCols>& G_d)
{
	A_d[0][0] = 1;
	A_d[0][1] = dT;
	A_d[0][2] = dT * dT / 2;
	A_d[1][0] = 0;
	A_d[1][1] = 1;
	A_d[1][2] = dT;
	A_d[2][0] = 0;
	A_d[2][1] = 0;
	A_d[2][2] = 1;

	B_d[0] = dT * dT * dT / 6;
	B_d[1] = dT * dT / 2;
	B_d[2] = dT;

	G_d[0] = 1;
	G_d[1] = 0;
	G_d[2] = -height / gravity;
}

template <unsigned nRows, unsigned nCols>
void ComputeDiscreteModel_2(CstMatrix<nRows, nCols>& Alpha_d, CstVector<nRows>& Delta_d)
{
	Alpha_d[0][0] = 1;
	Alpha_d[0][1] = dT;
	Alpha_d[1][0] = 0;
	Alpha_d[1][1] = 1;

	Delta_d[0] = dT * dT / 2;
	Delta_d[1] = dT;
}


//	Compute Predictive Models
void SetPredictiveSystemModel()
{
	unsigned	i, j;
	double		N;

	P_pu.setZero();
	P_vu.setZero();
	P_zu.setZero();
	P_as.setZero();

	H_ps.setZero();
	H_pu.setZero();
	H_vs.setZero();
	H_vu.setZero();

	for (i = 0; i < PREVIEW_LENGTH; i++) {
		N = (double)i + 1;

		P_ps[i][0] = 1;		P_ps[i][1] = N * dT;	P_ps[i][2] = SQR(N * dT) / 2;
		P_vs[i][0] = 0;		P_vs[i][1] = 1;			P_vs[i][2] = N * dT;
		P_as[i][0] = 0;		P_as[i][1] = 0;			P_as[i][2] = 1;

		H_ps[i][0] = 1;     H_ps[i][1] = N * dT;
		H_vs[i][0] = 0;     H_vs[i][1] = 1;

		for (j = 0; j <= i; j++) {
			N = i - j + 1;
			P_pu[i][j] = (3 * N * N - 3 * N + 1) * CUBE(dT) / 6;
			P_vu[i][j] = (2 * N - 1) * SQR(dT) / 2;
			P_au[i][j] = dT;

			H_pu[i][j] = (2 * N - 1) * SQR(dT) / 2;
			H_vu[i][j] = dT;
		}
	}
}

void Set_z_Bound(const unsigned DIM_SSP, const unsigned DIM_DSP, const unsigned DIM_EXTRA, const unsigned DIM_TOTAL, double* z_Bound){
	for (int i = 0; i <= DIM_TOTAL; i++) {
		const unsigned DIM_STEP = DIM_SSP + DIM_DSP;
		if (i >= DIM_EXTRA && i < (NO_STEPS)* DIM_STEP + DIM_EXTRA) {
			if ((i - DIM_EXTRA - 1) % DIM_STEP - DIM_SSP == 0) {
				z_Bound[i] = BETA_M;
			}
			else {
				z_Bound[i] = BETA_M;
			}

		}
		else {
			z_Bound[i] = BETA_M;
		}
	}
}

//	Set Reference ZMP Trajectory
void Set_ZMP_Reference(const unsigned DIM_SSP, const unsigned DIM_DSP, const unsigned DIM_EXTRA, const unsigned DIM_TOTAL,
	double* xZMP_ref, double* yZMP_ref, double* zZMP_ref, double* x_foot, double* y_foot)
{
	unsigned i;
	int j;
	unsigned k = 0;
	double l_s, l_d;
	double theta = 0, theta_2 = 0;
	unsigned Quotient, residual;
	double	foot_info = FIRST_FOOT * STEP_WIDTH / 2;
	const unsigned DIM_STEP = (DIM_DSP + DIM_SSP);

	if (STEP_LENGTH >= FOOT_LENGTH * DIM_STEP / DIM_SSP) {
		l_s = FOOT_LENGTH;
		l_d = STEP_LENGTH - l_s;
	}
	else {
		l_s = STEP_LENGTH * DIM_SSP / DIM_STEP;
		l_d = STEP_LENGTH - l_s;
	}

	for (int i = 0; i <= DIM_TOTAL; i++) {

		j = i - DIM_EXTRA - 1;
		Quotient = j / DIM_STEP;
		residual = Quotient % 2;

		if (i <= DIM_EXTRA) {
			Quotient = 0;
			xZMP_ref[i] = yZMP_ref[i] = 0;
			zZMP_ref[i] = CoM_H;

		}

		else if (i > DIM_EXTRA && i <= DIM_SSP + DIM_EXTRA) {
			xZMP_ref[i] = ((j + 1) / (double)DIM_SSP) * (l_s / 2);
			yZMP_ref[i] = pow(-1, (double)residual) * (foot_info + FIRST_FOOT * l_s * sin(TOE_OUT * M_PI / 180) * ((j + 1) / (double)DIM_SSP));

			if (i <= DIM_EXTRA + DIM_SSP / 2) {
				zZMP_ref[i] = CoM_H + 0.05 * sin(M_PI * j / (double)DIM_SSP);
			}
			else {
				zZMP_ref[i] = CoM_H - 0.05 * cos(2 * M_PI * j / (double)DIM_SSP);
			}
		}
		else if (i > DIM_SSP + DIM_EXTRA && i <= DIM_STEP + DIM_EXTRA) {
			xZMP_ref[i] = xZMP_ref[DIM_SSP + DIM_EXTRA] + ((j - DIM_SSP + 1) / (double)DIM_DSP) * l_d;
			yZMP_ref[i] = yZMP_ref[DIM_SSP + DIM_EXTRA] - pow(-1, (double)residual) * FIRST_FOOT * (STEP_WIDTH + l_s * sin(TOE_OUT * M_PI / 180)) * ((j - DIM_SSP + 1) / (double)DIM_DSP);
			zZMP_ref[i] = CoM_H - 0.05 * cos(2 * M_PI * j / (double)DIM_SSP);

			if (THETA_STEP_START == 1 && THETA_STEP >= 1) {
				xZMP_ref[i] = xZMP_ref[i] + (l_s / 2) * (1 - cos((THETA * ((j % DIM_STEP) - DIM_SSP + 1) / (double)DIM_DSP) * M_PI / 180));
				yZMP_ref[i] = yZMP_ref[i] - (l_s / 2) * sin((THETA * ((j % DIM_STEP) - DIM_SSP + 1) / (double)DIM_DSP) * M_PI / 180);
			}
		}

		else if (i > DIM_STEP + DIM_EXTRA && i <= (NO_STEPS)* DIM_STEP + DIM_EXTRA) {
			if (j % DIM_STEP == 0) {
				k++;
			}

			if (i > THETA_STEP_START * DIM_STEP + DIM_EXTRA && i <= (THETA_STEP_START + THETA_STEP) * DIM_STEP + DIM_EXTRA && j % DIM_STEP == 0) {
				theta = theta + THETA;
			}
			if (i > (THETA_STEP_START - 1) * DIM_STEP + DIM_EXTRA && i <= (THETA_STEP_START + THETA_STEP - 1) * DIM_STEP + DIM_EXTRA && j % DIM_STEP == 0) {
				theta_2 = theta + THETA;

			}


			if (i > k * DIM_STEP + DIM_EXTRA && i <= k * (DIM_STEP) + DIM_SSP + DIM_EXTRA) {
				xZMP_ref[i] = xZMP_ref[k * DIM_STEP + DIM_EXTRA] + ((j % DIM_STEP + 1) / (double)DIM_SSP) * l_s * cos(theta * M_PI / 180);
				yZMP_ref[i] = yZMP_ref[k * DIM_STEP + DIM_EXTRA] + ((j % DIM_STEP + 1) / (double)DIM_SSP) * (FIRST_FOOT * pow(-1, (double)residual) * l_s * sin(TOE_OUT * M_PI /180) + l_s * sin(theta * M_PI / 180));
				if ((j % DIM_STEP) < DIM_SSP / 4) {
					zZMP_ref[i] = CoM_H - 0.05 * cos(2 * M_PI * (j % DIM_STEP + DIM_DSP) / (double)DIM_SSP);
				}
				else {
					zZMP_ref[i] = CoM_H + 0.05 * cos(2 * M_PI * ((j % DIM_STEP) - DIM_SSP / 4) / ((double)DIM_SSP * 1.5));
				}
			}
			else {
				xZMP_ref[i] = xZMP_ref[k * (DIM_STEP)+DIM_SSP + DIM_EXTRA] + (((j % DIM_STEP) - DIM_SSP + 1) / (double)DIM_DSP) * (l_d * cos(theta * M_PI / 180) + FIRST_FOOT * pow(-1, (double)residual) * STEP_WIDTH * sin(theta * 2 * M_PI / 360));
				yZMP_ref[i] = yZMP_ref[k * (DIM_STEP)+DIM_SSP + DIM_EXTRA] - FIRST_FOOT * pow(-1, (double)residual) * (((j % DIM_STEP) - DIM_SSP + 1) / (double)DIM_DSP) * ((STEP_WIDTH * cos(theta * M_PI / 180) - FIRST_FOOT * pow(-1, (double)residual) * l_d * sin(theta * M_PI / 180)) + l_s * sin(TOE_OUT * M_PI / 180));

				zZMP_ref[i] = CoM_H + 0.05 * cos(2 * M_PI * ((j % DIM_STEP) - DIM_SSP / 4) / ((double)DIM_SSP * 1.5));

				xZMP_ref[i] = xZMP_ref[i] + (l_s / 2) * (cos(theta * M_PI / 180) - cos((theta + (theta_2 - theta) * (double)((j % DIM_STEP) - DIM_SSP + 1) / (double)DIM_DSP) * M_PI / 180));
			    yZMP_ref[i] = yZMP_ref[i] - (l_s / 2) * (-sin(theta * M_PI / 180) + sin((theta + (theta_2 - theta) * (double)((j % DIM_STEP) - DIM_SSP + 1) / (double)DIM_DSP) * M_PI / 180));
			}
		}
		else if (i > (NO_STEPS)* DIM_STEP + DIM_EXTRA && i <= DIM_SSP + (NO_STEPS)* DIM_STEP + DIM_EXTRA) {
			xZMP_ref[i] = xZMP_ref[(NO_STEPS)* DIM_STEP + DIM_EXTRA] + (((j % DIM_STEP) + 1) / (double)DIM_SSP) * (l_s / 2) * cos(theta * 2 * M_PI / 360);
			yZMP_ref[i] = yZMP_ref[(NO_STEPS)* DIM_STEP + DIM_EXTRA] + (((j % DIM_STEP) + 1) / (double)DIM_SSP) * ((l_s / 2) * sin(theta * 2 * M_PI / 360) + FIRST_FOOT * pow(-1, (double)residual) * l_s * sin(TOE_OUT * M_PI / 180));
			if ((j % DIM_STEP) < DIM_SSP / 4) {
				zZMP_ref[i] = CoM_H - 0.05 * cos(2 * M_PI * (j % DIM_STEP + DIM_DSP) / (double)DIM_SSP);
			}
			else {
				zZMP_ref[i] = (CoM_H + 0.025) + 0.025 * cos(2 * M_PI * ((j % DIM_STEP) - DIM_SSP / 4) / ((double)DIM_SSP * 1.5));
			}
		}
		else {
			xZMP_ref[i] = xZMP_ref[DIM_SSP + (NO_STEPS) * (DIM_STEP)+DIM_EXTRA] + foot_info * pow(-1, (double)NO_STEPS) * sin(theta * 2 * M_PI / 360);
			yZMP_ref[i] = yZMP_ref[DIM_SSP + (NO_STEPS) * (DIM_STEP)+DIM_EXTRA] - pow(-1, (double)NO_STEPS) * (foot_info * cos(theta * 2 * M_PI / 360) + FIRST_FOOT * l_s * sin(TOE_OUT * M_PI / 180));
			zZMP_ref[i] = CoM_H;
		}
	}

	int b = 0;

	x_foot[0] = x_foot[1] = x_foot[4] = x_foot[5] = -FOOT_LENGTH / 2;
	x_foot[2] = x_foot[3] = x_foot[6] = x_foot[7] = FOOT_LENGTH / 2;
	y_foot[0] = y_foot[3] = -FIRST_FOOT * (STEP_WIDTH + FOOT_WIDTH) / 2;
	y_foot[1] = y_foot[2] = -FIRST_FOOT * (STEP_WIDTH - FOOT_WIDTH) / 2;
	y_foot[4] = y_foot[7] = FIRST_FOOT * (STEP_WIDTH - FOOT_WIDTH) / 2;
	y_foot[5] = y_foot[6] = FIRST_FOOT * (STEP_WIDTH + FOOT_WIDTH) / 2;

	for (int a = 1; a < NO_STEPS; a++) {
		x_foot[4 + 4 * a] = xZMP_ref[a * DIM_STEP + DIM_EXTRA] - FOOT_WIDTH * sin(b * THETA * M_PI / 180) / 2;
		x_foot[5 + 4 * a] = xZMP_ref[a * DIM_STEP + DIM_EXTRA] + FOOT_WIDTH * sin(b * THETA * M_PI / 180) / 2;
		x_foot[6 + 4 * a] = xZMP_ref[a * DIM_STEP + DIM_SSP + DIM_EXTRA] + FOOT_WIDTH * sin(b * THETA * M_PI / 180) / 2;
		x_foot[7 + 4 * a] = xZMP_ref[a * DIM_STEP + DIM_SSP + DIM_EXTRA] - FOOT_WIDTH * sin(b * THETA * M_PI / 180) / 2;
		if (l_s < FOOT_LENGTH) {
			x_foot[4 + 4 * a] = x_foot[4 + 4 * a] - (FOOT_LENGTH - l_s) / 2;
			x_foot[5 + 4 * a] = x_foot[5 + 4 * a] - (FOOT_LENGTH - l_s) / 2;
			x_foot[6 + 4 * a] = x_foot[6 + 4 * a] + (FOOT_LENGTH - l_s) / 2;
			x_foot[7 + 4 * a] = x_foot[7 + 4 * a] + (FOOT_LENGTH - l_s) / 2;
		}
		y_foot[4 + 4 * a] = yZMP_ref[a * DIM_STEP + DIM_EXTRA] + FOOT_WIDTH * cos(b * THETA * M_PI / 180) / 2;
		y_foot[5 + 4 * a] = yZMP_ref[a * DIM_STEP + DIM_EXTRA] - FOOT_WIDTH * cos(b * THETA * M_PI / 180) / 2;
		y_foot[6 + 4 * a] = yZMP_ref[a * DIM_STEP + DIM_SSP + DIM_EXTRA] - FOOT_WIDTH * cos(b * THETA * M_PI / 180) / 2;
		y_foot[7 + 4 * a] = yZMP_ref[a * DIM_STEP + DIM_SSP + DIM_EXTRA] + FOOT_WIDTH * cos(b * THETA * M_PI / 180) / 2;

		if (a >= THETA_STEP_START - 1 && a < THETA_STEP_START + THETA_STEP - 1) {
			b++;
		}
	}
	x_foot[4 + 4 * NO_STEPS] = xZMP_ref[NO_STEPS * DIM_STEP + DIM_EXTRA] - FOOT_WIDTH * sin(THETA_STEP * THETA * M_PI / 180) / 2;
	x_foot[5 + 4 * NO_STEPS] = xZMP_ref[NO_STEPS * DIM_STEP + DIM_EXTRA] + FOOT_WIDTH * sin(THETA_STEP * THETA * M_PI / 180) / 2;
	x_foot[6 + 4 * NO_STEPS] = xZMP_ref[NO_STEPS * DIM_STEP + DIM_SSP + DIM_EXTRA] + FOOT_WIDTH * sin(THETA_STEP * THETA * M_PI / 180) / 2 + FOOT_LENGTH * cos(THETA_STEP * THETA * M_PI / 180) / 2;
	x_foot[7 + 4 * NO_STEPS] = xZMP_ref[NO_STEPS * DIM_STEP + DIM_SSP + DIM_EXTRA] - FOOT_WIDTH * sin(THETA_STEP * THETA * M_PI / 180) / 2 + FOOT_LENGTH * cos(THETA_STEP * THETA * M_PI / 180) / 2;
	y_foot[4 + 4 * NO_STEPS] = yZMP_ref[NO_STEPS * DIM_STEP + DIM_EXTRA] + FOOT_WIDTH * cos(THETA_STEP * THETA * M_PI / 180) / 2;
	y_foot[5 + 4 * NO_STEPS] = yZMP_ref[NO_STEPS * DIM_STEP + DIM_EXTRA] - FOOT_WIDTH * cos(THETA_STEP * THETA * M_PI / 180) / 2;
	y_foot[6 + 4 * NO_STEPS] = yZMP_ref[NO_STEPS * DIM_STEP  + DIM_EXTRA] - FOOT_WIDTH * cos(THETA_STEP * THETA * M_PI / 180) / 2 + FOOT_LENGTH * sin(THETA_STEP * THETA * M_PI / 180) / 2;
	y_foot[7 + 4 * NO_STEPS] = yZMP_ref[NO_STEPS * DIM_STEP  + DIM_EXTRA] + FOOT_WIDTH * cos(THETA_STEP * THETA * M_PI / 180) / 2 + FOOT_LENGTH * sin(THETA_STEP * THETA * M_PI / 180) / 2;

	x_foot[4 * NO_STEPS + 8] = x_foot[4 + 4 * NO_STEPS] + FIRST_FOOT * pow(-1, (double)NO_STEPS) * STEP_WIDTH * sin(THETA_STEP * THETA * M_PI / 180);
	x_foot[4 * NO_STEPS + 9] = x_foot[5 + 4 * NO_STEPS] + FIRST_FOOT * pow(-1, (double)NO_STEPS) * STEP_WIDTH * sin(THETA_STEP * THETA * M_PI / 180);
	x_foot[4 * NO_STEPS + 10] = x_foot[6 + 4 * NO_STEPS] + FIRST_FOOT * pow(-1, (double)NO_STEPS) * STEP_WIDTH * sin(THETA_STEP * THETA * M_PI / 180);
	x_foot[4 * NO_STEPS + 11] = x_foot[7 + 4 * NO_STEPS] + FIRST_FOOT * pow(-1, (double)NO_STEPS) * STEP_WIDTH * sin(THETA_STEP * THETA * M_PI / 180);
	y_foot[4 * NO_STEPS + 8] = y_foot[4 + 4 * NO_STEPS] - FIRST_FOOT * pow(-1, (double)NO_STEPS) * STEP_WIDTH * cos(THETA_STEP * THETA * M_PI / 180);
	y_foot[4 * NO_STEPS + 9] = y_foot[5 + 4 * NO_STEPS] - FIRST_FOOT * pow(-1, (double)NO_STEPS) * STEP_WIDTH * cos(THETA_STEP * THETA * M_PI / 180);
	y_foot[4 * NO_STEPS + 10] = y_foot[6 + 4 * NO_STEPS] - FIRST_FOOT * pow(-1, (double)NO_STEPS) * STEP_WIDTH * cos(THETA_STEP * THETA * M_PI / 180);
	y_foot[4 * NO_STEPS + 11] = y_foot[7 + 4 * NO_STEPS] - FIRST_FOOT * pow(-1, (double)NO_STEPS) * STEP_WIDTH * cos(THETA_STEP * THETA * M_PI / 180);
}


void Get_ZMP_Reference(const unsigned k, const unsigned DIM_TOTAL, double* ref_zmp, CstVector<PREVIEW_LENGTH>& Zmp_ref_vec)
{
	if (k + PREVIEW_LENGTH <= DIM_TOTAL) {
		memcpy(&Zmp_ref_vec[0], &ref_zmp[k + 1], sizeof(double) * (PREVIEW_LENGTH));
	}
}

void Get_Friction(const unsigned k, const unsigned DIM_TOTAL, double* friction, CstVector<PREVIEW_LENGTH>& Friction)
{
	if (k + PREVIEW_LENGTH <= DIM_TOTAL) {
		memcpy(&Friction[0], &friction[k + 1], sizeof(double) * (PREVIEW_LENGTH));
	}
}

void Get_z_bound(const unsigned k, const unsigned DIM_TOTAL, double* z_Bound, CstVector<PREVIEW_LENGTH>& z_bound)
{
	if (k + PREVIEW_LENGTH <= DIM_TOTAL) {
		memcpy(&z_bound[0], &z_Bound[k + 1], sizeof(double) * (PREVIEW_LENGTH));
	}
}