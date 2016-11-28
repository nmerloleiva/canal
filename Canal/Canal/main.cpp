/*
**	Universidad de Buenos Aires
**	Facultad de Ingeniería
**	75.12 Análisis Numérico I
**	Trabajo Práctico 2
**	Curso 3
**	30/11/2016
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
	REAL var_1_2 = 0;
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
	REAL Q = 0;
};

struct TEST_STEP
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
	REAL Q_max = 0;
	
	REAL GetTotalTime()
	{
		return steps.back().time;
	}

	void AddStep(REAL time, REAL H_1, REAL H_2, REAL Q)
	{
		SIMULATION_STEP step;
		step.time = time;
		step.H_1 = H_1;
		step.H_2 = H_2;
		step.Q = Q;
		steps.push_back(step);

		if (H_2 > H_2_max)
		{
			H_2_max = H_2;
		}
		if (abs(Q) > Q_max)
		{
			Q_max = abs(Q);
		}
	}

	void Reset(REAL k)
	{
		this->steps.clear();
		this->k = k;
		this->H_2_max = 0;
		this->Q_max = 0;
	}
};

typedef void(SOLVE_METHOD)(CANAL_CONTEXT&, SOLVE_CONTEXT&);

struct SIMULATION_CONTEXT
{
	SIMULATION_RESULT result;
	CANAL_CONTEXT canalContext;
	REAL openTime = 0;
	REAL closeStartTime = 0;
	REAL closeTime = 0;
	REAL timeLimit = 0;
	REAL maxOscillationAmplitude = 0;
	REAL startStep = 0;
	REAL maxTruncationError = 0;
	REAL startHeight1 = 0;
	REAL startHeight2 = 0;
	SOLVE_METHOD* solveMethod = nullptr;
};

REAL ComputeKv(REAL openPercent)
{
	if (openPercent > 1)
	{
		return Kv_table[Kv_table_size - 1].second;
	}

	if (openPercent < 0)
	{
		return Kv_table[0].second;
	}

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

	return 0;
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

void CanalSolveRK4Step(CANAL_CONTEXT& c, SOLVE_CONTEXT& s)
{
	REAL q_1_H_1 = s.k * ((-2 * s.Q_n) / c.omega_1);
	REAL q_1_H_2 = s.k * ((2 * s.Q_n) / c.omega_2);
	REAL q_1_Q = s.k * (
		s.const_1 * (s.H_1_n - s.H_2_n) -
		s.const_2 * s.Q_n * abs(s.Q_n) -
		s.const_3 * s.Q_n * abs(s.Q_n) -
		s.var_1 * s.Q_n * abs(s.Q_n));

	REAL q_2_H_1 = s.k * ((-2 * (s.Q_n + 0.5 * q_1_Q)) / c.omega_1);
	REAL q_2_H_2 = s.k * ((2 * (s.Q_n + 0.5 * q_1_Q)) / c.omega_2);
	REAL q_2_Q = s.k * (
		s.const_1 * (s.H_1_n + 0.5 * q_1_H_1 - (s.H_2_n + 0.5 * q_1_H_2)) -
		s.const_2 * (s.Q_n + 0.5 * q_1_Q) * abs((s.Q_n + 0.5 * q_1_Q)) -
		s.const_3 * (s.Q_n + 0.5 * q_1_Q)* abs((s.Q_n + 0.5 * q_1_Q)) -
		s.var_1_2 * (s.Q_n + 0.5 * q_1_Q) * abs((s.Q_n + 0.5 * q_1_Q)));

	REAL q_3_H_1 = s.k * ((-2 * (s.Q_n + 0.5 * q_2_Q)) / c.omega_1);
	REAL q_3_H_2 = s.k * ((2 * (s.Q_n + 0.5 * q_2_Q)) / c.omega_2);
	REAL q_3_Q = s.k * (
		s.const_1 * (s.H_1_n + 0.5 * q_2_H_1 - (s.H_2_n + 0.5 * q_2_H_2)) -
		s.const_2 * (s.Q_n + 0.5 * q_2_Q) * abs((s.Q_n + 0.5 * q_2_Q)) -
		s.const_3 * (s.Q_n + 0.5 * q_2_Q)* abs((s.Q_n + 0.5 * q_2_Q)) -
		s.var_1_2 * (s.Q_n + 0.5 * q_2_Q) * abs((s.Q_n + 0.5 * q_2_Q)));

	REAL q_4_H_1 = s.k * ((-2 * (s.Q_n + q_3_Q)) / c.omega_1);
	REAL q_4_H_2 = s.k * ((2 * (s.Q_n + q_3_Q)) / c.omega_2);
	REAL q_4_Q = s.k * (
		s.const_1 * (s.H_1_n + q_3_H_1 - (s.H_2_n + q_3_H_2)) -
		s.const_2 * (s.Q_n + q_3_Q) * abs((s.Q_n + q_3_Q)) -
		s.const_3 * (s.Q_n + q_3_Q)* abs((s.Q_n + q_3_Q)) -
		s.var_2 * (s.Q_n + q_3_Q) * abs((s.Q_n + q_3_Q)));

	s.H_1_n_1 = s.H_1_n + (q_1_H_1 + 2 * q_2_H_1 + 2 * q_3_H_1 + q_4_H_1) / 6.0;
	s.H_2_n_1 = s.H_2_n + (q_1_H_2 + 2 * q_2_H_2 + 2 * q_3_H_2 + q_4_H_2) / 6.0;
	s.Q_n_1 = s.Q_n + (q_1_Q + 2 * q_2_Q + 2 * q_3_Q + q_4_Q) / 6.0;
}

REAL ComputeGravityFromLatitude(REAL l)
{
	return 9.780327 * (1.0 + 0.0053024 * pow(sin(l), 2) - 0.0000058 * pow(sin(2 * l), 2));
}

void RunCanalSimulation(SIMULATION_CONTEXT& simulation)
{
	CANAL_CONTEXT canalContext = simulation.canalContext;

	SOLVE_CONTEXT solveContext;
	solveContext.const_1 = (canalContext.g * canalContext.A) / canalContext.L;
	solveContext.const_2 = canalContext.f / (2.0 * canalContext.D_e * canalContext.A);
	solveContext.const_3 = canalContext.K_e / (2.0 * canalContext.A * canalContext.L);
	solveContext.k = simulation.startStep; // s
	
	simulation.result.Reset(solveContext.k);

	bool truncationErrorBelowMax = false;
	while (!truncationErrorBelowMax)
	{
		solveContext.t = 0; // s
		solveContext.H_1_n = simulation.startHeight1; // H_1_0
		solveContext.H_2_n = simulation.startHeight2; // H_2_0
		solveContext.H_1_n_1 = 0; 
		solveContext.H_2_n_1 = 0;
		solveContext.Q_n = 0.0; // Q_0
		solveContext.H_2_delta = solveContext.H_2_n_1 - solveContext.H_2_n; // m

		truncationErrorBelowMax = true;
		bool openTimeHit = solveContext.t >= simulation.openTime;
		bool oscillationHit = solveContext.H_1_n < solveContext.H_2_n;
		bool stableOscillationHit = abs(solveContext.H_1_n - solveContext.H_2_n) < simulation.maxOscillationAmplitude;
		bool timeLimitHit = simulation.timeLimit ? solveContext.t >= simulation.timeLimit : false;
		while (!timeLimitHit &&
			truncationErrorBelowMax &&
			(!openTimeHit ||
			!oscillationHit ||
			!stableOscillationHit))
		{
			if (solveContext.t <= simulation.openTime)
			{
				REAL Kv_n = ComputeKv(solveContext.t / simulation.openTime);
				REAL Kv_n_1_2 = ComputeKv((solveContext.t + (solveContext.k / 2)) / simulation.openTime);
				REAL Kv_n_1 = ComputeKv((solveContext.t + solveContext.k) / simulation.openTime);
				solveContext.var_1 = Kv_n / (2.0 * canalContext.A * canalContext.L);
				solveContext.var_1_2 = Kv_n_1_2 / (2.0 * canalContext.A * canalContext.L);
				solveContext.var_2 = Kv_n_1 / (2.0 * canalContext.A * canalContext.L);
			}
			else if (solveContext.t >= simulation.closeStartTime &&
				solveContext.t <= (simulation.closeStartTime + simulation.closeTime))
			{
				REAL closeTimeLeft = simulation.closeTime - (solveContext.t - simulation.closeStartTime);
				REAL Kv_n = ComputeKv(closeTimeLeft / simulation.closeTime);
				REAL Kv_n_1_2 = ComputeKv((closeTimeLeft - (solveContext.k / 2)) / simulation.closeTime);
				REAL Kv_n_1 = ComputeKv((closeTimeLeft - solveContext.k) / simulation.closeTime);
				solveContext.var_1 = Kv_n / (2.0 * canalContext.A * canalContext.L);
				solveContext.var_1_2 = Kv_n_1_2 / (2.0 * canalContext.A * canalContext.L);
				solveContext.var_2 = Kv_n_1 / (2.0 * canalContext.A * canalContext.L);
			}

			simulation.solveMethod(canalContext, solveContext);

			simulation.result.AddStep(solveContext.t, solveContext.H_1_n, solveContext.H_2_n, solveContext.Q_n);

			if (simulation.maxTruncationError)
			{
				REAL truncationErrorH_1 = abs(solveContext.H_1_n_1 - solveContext.H_1_n);
				REAL truncationErrorH_2 = abs(solveContext.H_2_n_1 - solveContext.H_2_n);
				if (truncationErrorH_1 > simulation.maxTruncationError ||
					truncationErrorH_2 > simulation.maxTruncationError)
				{
					truncationErrorBelowMax = false;
					solveContext.k /= 10;
					simulation.result.Reset(solveContext.k);
				}
				else
				{
					truncationErrorBelowMax = true;
				}
			}
			else
			{
				truncationErrorBelowMax = true;
			}
	
			if (truncationErrorBelowMax)
			{
				if (!openTimeHit)
				{
					openTimeHit = solveContext.t >= simulation.openTime;
				}
				else
				{
					timeLimitHit = simulation.timeLimit ? solveContext.t >= simulation.timeLimit : false;

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

				// Move to next step
				solveContext.H_2_delta = solveContext.H_2_n_1 - solveContext.H_2_n;
				solveContext.H_1_n = solveContext.H_1_n_1;
				solveContext.H_2_n = solveContext.H_2_n_1;
				solveContext.Q_n = solveContext.Q_n_1;
				solveContext.t += solveContext.k;
			}
		}
	}
}

void SplitCSV(std::vector<std::string>& values, const std::string& line)
{
	std::string buff;
	for each (auto c in line)
	{
		if (c == ';')
		{
			values.push_back(buff);
			buff.clear();
		}
		else
		{
			buff += c;
		}
	}
	values.push_back(buff);
}

void RunKeAdjustment(SIMULATION_CONTEXT& simulation, std::string testFileName)
{
	std::fstream file;
	std::vector<TEST_STEP> testSteps;
	file.open(testFileName, std::ios_base::in);
	std::string line;
	if (std::getline(file, line))
	{
		std::vector<std::string> values;
		SplitCSV(values, line);
		simulation.openTime = std::stod(values.at(0));
		simulation.closeStartTime = std::stod(values.at(1));
		simulation.closeTime = std::stod(values.at(2));
	}
	while (std::getline(file, line))
	{
		std::vector<std::string> values;
		SplitCSV(values, line);
		TEST_STEP step;
		step.time = std::stod(values.at(0));
		step.H_1 = std::stod(values.at(1));
		step.H_2 = std::stod(values.at(2));
		testSteps.push_back(step);
	}
	file.close();

	simulation.timeLimit = testSteps.back().time;
	simulation.startHeight1 = testSteps.front().H_1;
	simulation.startHeight2 = testSteps.front().H_2;

	bool keAdjusted = false;
	simulation.canalContext.K_e = 0;
	REAL lastMaxStepError = 0;
	while (!keAdjusted)
	{
		RunCanalSimulation(simulation);

		REAL maxStepError = 0;
		size_t simulationStepIndex = 0;
		for (size_t testStepIndex = 0; testStepIndex < testSteps.size(); testStepIndex++)
		{
			TEST_STEP testStep = testSteps.at(testStepIndex);
			while (simulation.result.steps.at(simulationStepIndex).time < testStep.time)
			{
				simulationStepIndex++;
			}
			REAL stepError = abs(testStep.H_2 - simulation.result.steps.at(simulationStepIndex).H_2);
			if (stepError > maxStepError)
			{
				maxStepError = stepError;
			}
		}

		keAdjusted = simulation.canalContext.K_e > 0 && maxStepError > lastMaxStepError;
		simulation.canalContext.K_e += 0.001;
		lastMaxStepError = maxStepError;
	}
}

void SaveResult(SIMULATION_CONTEXT& context, std::string fileName)
{
	std::fstream file;
	file.open(fileName, std::ios_base::out);
	file << "k;" << context.result.k << "\n";
	file << "f;" << context.canalContext.f << "\n";
	file << "ke;" << context.canalContext.K_e << "\n";
	file << "max H2;" << context.result.H_2_max << "\n";
	file << "max Q;" << context.result.Q_max << "\n";
	REAL oscillationt = 0;
	if (context.closeStartTime)
	{
		oscillationt = context.closeStartTime - context.openTime;
	}
	else
	{
		oscillationt = context.result.GetTotalTime() - context.openTime;
	}
	file << "oscillation t;" << oscillationt << "\n";
	file << "\ntime;Q;H1;H2\n";
	for each (auto step in context.result.steps)
	{
		if (floor(step.time) == step.time)
		{
			file << step.time << ";" << step.Q << ";" << step.H_1 << ";" << step.H_2 << "\n";
		}
	}
	file.close();
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

	printf("Run simulation 0...\n");
	SIMULATION_CONTEXT simulation0;
	simulation0.canalContext = canalContext;
	simulation0.maxTruncationError = 0.01;
	simulation0.maxOscillationAmplitude = 0.1;
	simulation0.openTime = 240;
	simulation0.startHeight1 = 18;
	simulation0.startHeight2 = 0;
	simulation0.startStep = 1;
	simulation0.solveMethod = CanalSolveRK2Step;
	RunCanalSimulation(simulation0);
	SaveResult(simulation0, "simulation0.result.csv");

	printf("Run simulation 1...\n");
	SIMULATION_CONTEXT simulation1;
	simulation1.canalContext = canalContext;
	simulation1.startStep = simulation0.result.k;
	simulation1.solveMethod = CanalSolveRK2Step;
	RunKeAdjustment(simulation1, "simulation1.csv");
	SaveResult(simulation1, "simulation1.result.csv");

	printf("Run simulation 2...\n");
	SIMULATION_CONTEXT simulation2;
	simulation2.canalContext = canalContext;
	simulation2.startStep = simulation0.result.k;
	simulation2.solveMethod = CanalSolveRK2Step;
	RunKeAdjustment(simulation2, "simulation2.csv");
	SaveResult(simulation2, "simulation2.result.csv");

	printf("Run simulation 3...\n");
	SIMULATION_CONTEXT simulation3;
	simulation3.canalContext = canalContext;
	simulation3.startStep = simulation0.result.k;
	simulation3.solveMethod = CanalSolveRK2Step;
	RunKeAdjustment(simulation3, "simulation3.csv");
	SaveResult(simulation3, "simulation3.result.csv");

	printf("Run simulation 4...\n");
	SIMULATION_CONTEXT simulation4;
	simulation4.canalContext = canalContext;
	simulation4.startStep = simulation0.result.k;
	simulation4.solveMethod = CanalSolveRK2Step;
	RunKeAdjustment(simulation4, "simulation4.csv");
	SaveResult(simulation4, "simulation4.result.csv");

	printf("Run simulation 5...\n");
	SIMULATION_CONTEXT simulation5;
	simulation5.canalContext = canalContext;
	simulation5.maxTruncationError = 0.01;
	simulation5.maxOscillationAmplitude = 0.1;
	simulation5.openTime = 240;
	simulation5.startHeight1 = 18;
	simulation5.startHeight2 = 0;
	simulation5.startStep = 1;
	simulation5.solveMethod = CanalSolveRK4Step;
	RunCanalSimulation(simulation5);
	SaveResult(simulation5, "simulation5.result.csv");

	return 0;
}