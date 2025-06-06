//
//	2D_Pattern_LTI.cpp : 이 파일에는 'main' 함수가 포함됩니다. 거기서 프로그램 실행이 시작되고 종료됩니다.
//
#include "QuadProg.h"
#include <chrono>

#define dT					double(0.05)	//	Sampling time [s]
#define	TIME_DSP			dT				//	Time duration for Double Support Phase : 현재는 dT
#define TIME_SSP			double(0.95)	//	Time duration for Single Support Phase : dT의 정수배
#define TIME_S2W			double(1.0)		//	Time duration for Stance to Walk : ~ TIME_PREVIEW
#define TIME_PREVIEW		double(1.6)		//	Preview Time : Approx 1.5 [s]

#define DOF					unsigned(2)
#define STATE_DIM			unsigned(3)		//	State dimension

#define	NO_STEPS			unsigned(10)		//	Number of steps >= 0
#define	FIRST_FOOT			int(-1)				//	첫 발 : 왼발 = -1, 오른발 = 1
#define WALKING_DIRECTION	int(1)				//	보행 방향 : 전진 = 1, 후진 = -1

#define CoM_H				((double)0.65)		//	CoM - ZMP : 현재는 constant
#define GRAVITY				((double)9.81)		//	중력가속도
#define FOOT_WIDTH			((double)0.08)		//	발볼의 너비
#define STEP_LENGTH			((double)0.2)		//	스텝 길이
#define STEP_WIDTH			((double)0.2)		//	양발 중심간의 거리

#define ALPHA1				double(0.001)		//	Weight for input : u = CoM Jerk
#define	ALPHA2				double(1000)		//	Weight for ZMP : ZMP - ref_ZMP
#define ALPHA3				double(0)			//	Weight for CoM velocity : CoM_dot - ref_CoM_dot
#define ALPHA4				double(0)			//	Weight for CoM position : CoM - ref_CoM

constexpr auto DIM_SSP = unsigned(TIME_SSP / dT + 0.001);
constexpr auto DIM_DSP = unsigned(TIME_DSP / dT + 0.001);
constexpr auto DIM_S2W = unsigned(TIME_S2W / dT + 0.001);
constexpr auto PREVIEW_DIM = unsigned(TIME_PREVIEW / dT);	//	Preview length : Approx. 1.6 / dT
constexpr auto DIM_TOTAL = 2 * DIM_S2W + (NO_STEPS + 1)*(DIM_SSP + DIM_DSP) + DIM_DSP;
constexpr auto MAX_ARRAY = DIM_TOTAL + PREVIEW_DIM;

template <unsigned nRows, unsigned nCols>
void ComputeDiscreteModel(CstMatrix<nRows, nCols>& A_d, CstVector<nRows>& B_d, CstVector<nCols>& G_d);

void SetPredictiveSystemModel();

void Set_ZMP_Reference(double *xZMP_ref, double *yZMP_ref);
void Get_ZMP_Reference(const unsigned index, const unsigned nDim, double *yZMP_ref, CstVector<PREVIEW_DIM> & preview_ZMP_ref);

CstMatrix<PREVIEW_DIM, STATE_DIM> P_ps, P_vs, P_zs;
CstMatrix<PREVIEW_DIM, PREVIEW_DIM> P_pu, P_vu, P_zu;


int main()
{
	FILE* fp;
	unsigned indx = 0;
	double time = 0, x_zmp = 0, y_zmp = 0;
	double xZMP_ref[MAX_ARRAY] = { 0 }, yZMP_ref[MAX_ARRAY] = { 0 };

	CstMatrix<PREVIEW_DIM, PREVIEW_DIM> Q_Mat, E_Mat, ALPHA_trP_zu;
	CstMatrix<PREVIEW_DIM * 2, PREVIEW_DIM> CI;
	CstVector<PREVIEW_DIM> u_x, u_y;
	CstVector<PREVIEW_DIM> g0, preview_ZMP_ref, ci0_upper, ci0_lower, foot_width_array;
	CstVector<PREVIEW_DIM * 2> ci0;

	CstMatrix<STATE_DIM, STATE_DIM> A_d;
	CstVector<STATE_DIM> c_y, tmp_state, B_d, G_d;

	printf("SSP Dimension = %d\t SSP Time = %lf[s]\n", DIM_SSP, DIM_SSP*dT);
	printf("DSP Dimension = %d\t DSP Time = %lf[s]\n", DIM_DSP, DIM_DSP*dT);
	printf("S2W Dimension = %d\t S2W Time = %lf[s]\n", DIM_S2W, DIM_S2W*dT);
	printf("Total Dimension = %d\t Total Exec Time = %lf[s]\n", DIM_TOTAL, DIM_TOTAL*dT);

	fopen_s(&fp, "Traj.txt", "wt");

	ComputeDiscreteModel(A_d, B_d, G_d);
	SetPredictiveSystemModel();

	Set_ZMP_Reference(xZMP_ref, yZMP_ref);

	//	Set Quadratic cost
	E_Mat.setIdentity();
	ALPHA_trP_zu = ALPHA2 * Transpose(P_zu);
	Q_Mat = ALPHA1 * E_Mat + ALPHA_trP_zu * P_zu;

	//
	//	LTI QP + Inequality Constraints 문제의 경우 :
	//	--> 'G' matrix의 trace & Cholesky_Decomposition() 한번만 풀면 됨!
	//
	double trace = 0.;
	for (int i = 0; i < PREVIEW_DIM; i++) {
		trace += Q_Mat[i][i];
	}
	QuadProgPP::QPCholesky_dcmp(Q_Mat);

	//	Set Inequality Constraints
	CI.setBlockMatrix(0, 0, -P_zu);
	CI.setBlockMatrix(PREVIEW_DIM, 0, P_zu);

	double yZMP_bound = FOOT_WIDTH / 2;

	do {
		//	Time Stamp Routine
		auto start_time = std::chrono::steady_clock::now();

		time = indx * dT;

		Get_ZMP_Reference(indx, DIM_TOTAL, yZMP_ref, preview_ZMP_ref);

		//	Set Cost function matrix & vector
		g0 = ALPHA_trP_zu * (P_zs * c_y - preview_ZMP_ref);

		//	Set inequality constraints
		ci0_upper = yZMP_bound + preview_ZMP_ref - P_zs * c_y;
		ci0_lower = yZMP_bound - preview_ZMP_ref + P_zs * c_y;

		memcpy(&ci0[0], &ci0_upper[0], sizeof(double)*PREVIEW_DIM);
		memcpy(&ci0[PREVIEW_DIM], &ci0_lower[0], sizeof(double)*PREVIEW_DIM);

		//	Solve Quadratic Programming : CASE-3 (LTI QP + Inequality Constraints)
		QuadProgPP::Solve_QuadProg(trace, Q_Mat, g0, CI, ci0, u_y);

		//	Compute state & ZMP
		tmp_state = A_d * c_y + B_d * u_y[0];
		y_zmp = G_d * c_y;

		std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time;
		//		auto del_time = std::chrono::steady_clock::now() - start_time;
		//		double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(del_time).count() / 1000.;

		//		cout.precision(3);
		//		cout << "Time : " << std::fixed << elapsed << " [us] to run.\n";
		//		cout << "Time : " << std::fixed << (double)std::chrono::duration_cast<std::chrono::nanoseconds>(time).count() / 1000. << " [us] to run.\n";

		fprintf_s(fp, "%lf\t %lf\t %lf\t %lf\t %lf\t", time, xZMP_ref[indx], yZMP_ref[indx], x_zmp, y_zmp);
		fprintf_s(fp, "%lf\t %lf\t %lf\t %lf\n", c_y[0], c_y[1], c_y[2], elapsed.count()*1000.);

		//	Update system state
		c_y = tmp_state;

	} while (++indx < DIM_TOTAL);

	fclose(fp);

	WaitForEnterKey();

	return 0;
}


//	Discrete System Model
template <unsigned nRows, unsigned nCols>
void ComputeDiscreteModel(CstMatrix<nRows, nCols> & A_d, CstVector<nRows> & B_d, CstVector<nCols> & G_d)
{
	A_d[0][0] = 1;				A_d[0][1] = dT;				A_d[0][2] = SQR(dT) / 2;
	A_d[1][0] = 0;				A_d[1][1] = 1;				A_d[1][2] = dT;
	A_d[2][0] = 0;				A_d[2][1] = 0;				A_d[2][2] = 1;

	B_d[0] = CUBE(dT) / 6;		B_d[1] = SQR(dT) / 2;		B_d[2] = dT;

	G_d[0] = 1.;				G_d[1] = 0.;				G_d[2] = -CoM_H / GRAVITY;
}


//	Compute Predictive Models
void SetPredictiveSystemModel()
{
	unsigned	i, j;
	double		N;

	P_pu.setZero();
	P_vu.setZero();
	P_zu.setZero();

	for (i = 0; i < PREVIEW_DIM; i++) {
		N = (double)i + 1;

		P_ps[i][0] = 1;		P_ps[i][1] = N * dT;	P_ps[i][2] = SQR(N * dT) / 2;
		P_vs[i][0] = 0;		P_vs[i][1] = 1;			P_vs[i][2] = N * dT;
		P_zs[i][0] = 1;		P_zs[i][1] = N * dT;	P_zs[i][2] = P_ps[i][2] - CoM_H / GRAVITY;

		for (j = 0; j <= i; j++) {
			N = i - j;
			P_pu[i][j] = (3 * N * N + 3 * N + 1)*CUBE(dT) / 6;
			P_vu[i][j] = (2 * N + 1)*SQR(dT) / 2;

			P_zu[i][j] = P_pu[i][j] - CoM_H * dT / GRAVITY;
		}
	}
}


//	Set Reference ZMP Trajectory
void Set_ZMP_Reference(double *xZMP_ref, double *yZMP_ref)
{
	unsigned i, j;
	unsigned Quotient, residual;
	unsigned DIM_STEP = DIM_DSP + DIM_SSP;
	double	foot_info = FIRST_FOOT * STEP_WIDTH / 2;

	for (i = 0; i < MAX_ARRAY; i++) {
//		if (i <= DIM_S2W) {
//			xZMP_ref[i] = yZMP_ref[i] = 0;
//		}
		if (i > DIM_S2W && i <= DIM_STEP + DIM_S2W) {
			j = i - DIM_S2W - 1;
			Quotient = j / DIM_STEP;
			residual = Quotient % 2;

			xZMP_ref[i] = 0;
			yZMP_ref[i] = pow(-1, (double)residual)*foot_info;
		}
		else if (i > DIM_STEP + DIM_S2W && i <= (NO_STEPS + 1)*DIM_STEP + DIM_S2W) {
			j = i - DIM_S2W - 1;
			Quotient = j / DIM_STEP;
			residual = Quotient % 2;

			xZMP_ref[i] = Quotient * STEP_LENGTH;
			yZMP_ref[i] = pow(-1, (double)residual)*foot_info;
		}
		else {
			xZMP_ref[i] = xZMP_ref[(NO_STEPS + 1)*DIM_STEP + DIM_S2W];
			yZMP_ref[i] = 0;
		}
	}
}


//
void Get_ZMP_Reference(const unsigned k, const unsigned DIM_TOTAL, double *ref_zmp, CstVector<PREVIEW_DIM> & Zmp_ref_vec)
{
	if (k + PREVIEW_DIM < MAX_ARRAY) {
		memcpy(&Zmp_ref_vec[0], &ref_zmp[k + 1], sizeof(double)*(PREVIEW_DIM));
	}
}
