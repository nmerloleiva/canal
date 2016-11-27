/*
**	Universidad de Buenos Aires
**	Facultad de Ingeniería
**	75.12 Análisis Numérico I
**	Trabajo Práctico 2
**	Curso 3
**	26/11/2016
**
**	Merlo Leiva Nahuel
**	Padrón 92115
*/

#define _USE_MATH_DEFINES // for C++  
#include <cmath>  
#include <memory>
#include <limits.h>
#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <iomanip>
#include <iostream>
#include <locale> 

// Precision type abstraction.
typedef double REAL;
#define REAL_MAX std::numeric_limits<REAL>::max()

#define PANAMA_COCOLI_LOCKS_LATITUDE 8.985470

struct REAL_PAIR
{
	REAL first;
	REAL second;
};

const REAL_PAIR Kv_table[] = 
{
	{ 0, 99999 },
	{ 0.01, 12000 },
	{ 0.04, 900 },
	{ 0.1, 105 },
	{ 0.2, 33 },
	{ 0.3, 12.5 },
	{ 0.4, 6.0 },
	{ 0.5, 3.2 },
	{ 0.6, 1.6 },
	{ 0.7, 1.0 },
	{ 0.8, 0.56 },
	{ 0.9, 0.28 },
	{ 1, 0.11 },
};

const auto Kv_table_size = sizeof(Kv_table) / sizeof(Kv_table[0]);

struct SOLVE_CONTEXT
{
	REAL tap = 0;
	REAL t = 0;
	REAL k = 0;
	REAL H_1_n = 0;
	REAL H_2_n = 0;
	REAL Q_n = 0;
	REAL H_1_n_1 = 0;
	REAL H_2_n_1 = 0;
	REAL Q_n_1 = 0;

	REAL const_1 = 0;
	REAL const_2 = 0;
	REAL const_3 = 0;

	REAL var_1 = 0;
	REAL var_2 = 0;

	REAL H_2_max = 0;
	REAL H_2_delta = 0;
	REAL H_2_osc_min = 0;
	REAL H_2_osc_max = 0;
	REAL H_1_truncation_error = 0;
	REAL H_2_truncation_error = 0;
};

struct CANAL_CONTEXT
{
	REAL omega_1 = 0;
	REAL omega_2 = 0;
	REAL g = 0;
	REAL A = 0;
	REAL L = 0;
	REAL D_e = 0;
	REAL f = 0;
	REAL K_e = 0;
};

struct SIMULATION_STEP
{
	REAL time = 0;
	REAL H_1 = 0;
	REAL H_2 = 0;
};

struct SIMULATION_RESULT
{
	std::vector<SIMULATION_STEP> steps;
	REAL k = 0;
	REAL H_2_max = 0;
	
	REAL GetTotalTime()
	{
		return steps.back().time;
	}

	void AddStep(REAL time, REAL H_1, REAL H_2)
	{
		SIMULATION_STEP step;
		step.time = time;
		step.H_1 = H_1;
		step.H_2 = H_2;
		steps.push_back(step);

		if (H_2 > H_2_max)
		{
			H_2_max = H_2;
		}
	}

	void Reset(REAL k)
	{
		this->steps.clear();
		this->k = k;
		this->H_2_max = 0;
	}
};

struct SIMULATION_CONTEXT
{
	SIMULATION_RESULT result;
	CANAL_CONTEXT canalContext;
	REAL openTime = 0;
	REAL oscillationTime = 0;
	REAL closeTime = 0;
	REAL maxOscillationAmplitude = 0;
	REAL startStep = 0;
	REAL maxTruncationError = 0;
	REAL startHeight1 = 0;
	REAL startHeight2 = 0;
};

REAL ComputeKv(REAL openPercent)
{
	for (size_t i = 0; i < Kv_table_size; i++)
	{
		if (openPercent == Kv_table[i].first)
		{
			return Kv_table[i].second;
		}
		else if ((i + 1) < Kv_table_size &&
			openPercent < Kv_table[i + 1].first)
		{
			REAL m = (Kv_table[i + 1].second - Kv_table[i].second) / Kv_table[i + 1].first - Kv_table[i].first;
			REAL b = Kv_table[i].second;
			return m * (openPercent - Kv_table[i].first) + b;
		}
	}

	// openPercent > 1: return minimum Kv
	return Kv_table[Kv_table_size - 1].second;
}

void CanalSolveRK2Step(CANAL_CONTEXT& c, SOLVE_CONTEXT& s)
{
	REAL q_1_H_1 = s.k * ((-2 * s.Q_n) / c.omega_1);
	REAL q_1_H_2 = s.k * ((2 * s.Q_n) / c.omega_2);
	REAL q_1_Q = s.k * (
		s.const_1 * (s.H_1_n - s.H_2_n) -
		s.const_2 * s.Q_n * abs(s.Q_n) -
		s.const_3 * s.Q_n * abs(s.Q_n) -
		s.var_1 * s.Q_n * abs(s.Q_n));
	REAL q_2_H_1 = s.k * ((-2 * (s.Q_n + q_1_Q)) / c.omega_1);
	REAL q_2_H_2 = s.k * ((2 * (s.Q_n + q_1_Q)) / c.omega_2);
	REAL q_2_Q = s.k * (
		s.const_1 * (s.H_1_n + q_1_H_1 - (s.H_2_n + q_1_H_2)) -
		s.const_2 * (s.Q_n + q_1_Q) * abs((s.Q_n + q_1_Q)) -
		s.const_3 * (s.Q_n + q_1_Q)* abs((s.Q_n + q_1_Q)) -
		s.var_2 * (s.Q_n + q_1_Q) * abs((s.Q_n + q_1_Q)));
	s.H_1_n_1 = s.H_1_n + 0.5 * (q_1_H_1 + q_2_H_1);
	s.H_2_n_1 = s.H_2_n + 0.5 * (q_1_H_2 + q_2_H_2);
	s.Q_n_1 = s.Q_n + 0.5 * (q_1_Q + q_2_Q);
}

REAL ComputeGravityFromLatitude(REAL l)
{
	return 9.780327 * (1.0 + 0.0053024 * pow(sin(l), 2) - 0.0000058 * pow(sin(2 * l), 2));
}

//std::fstream file;
//std::string fileName = "openLockSimulation.csv";
//file.open(fileName, std::ios_base::out);
//file << "Simulation start k=" << solveContext.k << "\n";

//if (floor(solveContext.t) == solveContext.t)
//{
//	file << solveContext.t << ";" << solveContext.H_1_n << ";" << solveContext.H_2_n << ";" << "\n";
//}

//file << "Simulation abort:" << " truncation error exceeded limit!" << " H1 error=" << truncationErrorH_1 << " H2 error=" << truncationErrorH_2 << "\n" << "\n";
//file << "Simulation completed!" << "\n";
//file.close();

void RunCanalSimulation(SIMULATION_CONTEXT& simulation)
{
	CANAL_CONTEXT canalContext = simulation.canalContext;

	SOLVE_CONTEXT solveContext;
	solveContext.const_1 = (canalContext.g * canalContext.A) / canalContext.L;
	solveContext.const_2 = canalContext.f / (2.0 * canalContext.D_e * canalContext.A);
	solveContext.const_3 = canalContext.K_e / (2.0 * canalContext.A * canalContext.L);
	solveContext.tap = simulation.openTime; // s
	solveContext.k = simulation.startStep; // s
	
	simulation.result.Reset(solveContext.k);

	bool truncationErrorBelowMax = false;
	while (!truncationErrorBelowMax)
	{
		solveContext.t = 0; // s
		solveContext.H_1_n = simulation.startHeight1; // H_1_0
		solveContext.H_2_n = simulation.startHeight2; // H_2_0
		solveContext.Q_n = 0.0; // Q_0
		solveContext.H_2_delta = solveContext.H_2_n_1 - solveContext.H_2_n; // m

		truncationErrorBelowMax = true;
		bool openTimeHit = solveContext.t >= solveContext.tap;
		bool oscillationHit = solveContext.H_1_n < solveContext.H_2_n;
		bool stableOscillationHit = abs(solveContext.H_1_n - solveContext.H_2_n) < simulation.maxOscillationAmplitude;
		while (truncationErrorBelowMax &&
			(!openTimeHit ||
			!oscillationHit ||
			!stableOscillationHit))
		{
			if (solveContext.t <= solveContext.tap)
			{
				REAL Kv_n = ComputeKv(solveContext.t / solveContext.tap);
				REAL Kv_n_1 = ComputeKv((solveContext.t + solveContext.k) / solveContext.tap);
				solveContext.var_1 = Kv_n / (2.0 * canalContext.A * canalContext.L);
				solveContext.var_2 = Kv_n_1 / (2.0 * canalContext.A * canalContext.L);
			}

			CanalSolveRK2Step(canalContext, solveContext);

			simulation.result.AddStep(solveContext.t, solveContext.H_1_n, solveContext.H_2_n);

			REAL truncationErrorH_1 = abs(solveContext.H_1_n_1 - solveContext.H_1_n);
			REAL truncationErrorH_2 = abs(solveContext.H_2_n_1 - solveContext.H_2_n);
			if (truncationErrorH_1 > simulation.maxTruncationError ||
				truncationErrorH_2 > simulation.maxTruncationError)
			{
				truncationErrorBelowMax = false;
				solveContext.k /= 2;
				simulation.result.Reset(solveContext.k);
			}
			else
			{
				truncationErrorBelowMax = true;

				if (!openTimeHit)
				{
					openTimeHit = solveContext.t >= solveContext.tap;
				}
				else
				{
					if (solveContext.H_2_delta < 0 &&
						solveContext.H_2_n_1 > solveContext.H_2_n)
					{
						solveContext.H_2_osc_min = solveContext.H_2_n;
					}
					else if (solveContext.H_2_delta > 0 &&
						solveContext.H_2_n_1 < solveContext.H_2_n)
					{
						solveContext.H_2_osc_max = solveContext.H_2_n;
					}

					if (!oscillationHit)
					{
						oscillationHit = solveContext.H_2_osc_min > 0 && solveContext.H_2_osc_max > 0;
					}
					else
					{
						stableOscillationHit = abs(solveContext.H_2_osc_max - solveContext.H_2_osc_min) < simulation.maxOscillationAmplitude;
					}
				}

				solveContext.H_2_delta = solveContext.H_2_n_1 - solveContext.H_2_n;
				solveContext.H_1_n = solveContext.H_1_n_1;
				solveContext.H_2_n = solveContext.H_2_n_1;
				solveContext.Q_n = solveContext.Q_n_1;
				solveContext.t += solveContext.k;
			}
		}
	}
}

int main()
{
	printf("Merlo Leiva Nahuel\n");
	printf("Padrón 92115\n");

	CANAL_CONTEXT canalContext;
	canalContext.omega_1 = 28000; // m^2
	canalContext.omega_2 = 28000; // m^2
	canalContext.g = ComputeGravityFromLatitude(PANAMA_COCOLI_LOCKS_LATITUDE); // m/s^2 
	canalContext.A = 8.30 * 6.50; // m
	canalContext.L = 600; // m
	canalContext.D_e = sqrt((4 * canalContext.A) / M_PI); // m
	canalContext.f = 0.02;
	canalContext.K_e = 0;

	SIMULATION_CONTEXT simulation1;
	simulation1.canalContext = canalContext;
	simulation1.maxTruncationError = 0.01;
	simulation1.maxOscillationAmplitude = 0.1;
	simulation1.openTime = 240;
	simulation1.startHeight1 = 18;
	simulation1.startHeight2 = 0;
	simulation1.startStep = 1;
	RunCanalSimulation(simulation1);

	return 0;
}