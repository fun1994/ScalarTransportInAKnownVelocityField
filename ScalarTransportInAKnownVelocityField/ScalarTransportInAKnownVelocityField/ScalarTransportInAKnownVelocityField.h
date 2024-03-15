#pragma once
#include "Data.h"

class ScalarTransportInAKnownVelocityField {
	double rho;
	double Gamma;
	int Nx;
	int Ny;
	double dx;
	double dy;
	std::string convection;
	double tol;
	double omega;
	void initialize(Data& data, std::vector<std::vector<double>>& A_W, std::vector<std::vector<double>>& A_S, std::vector<std::vector<double>>& A_P, std::vector<std::vector<double>>& A_N, std::vector<std::vector<double>>& A_E, std::vector<std::vector<double>>& Q_P, std::vector<std::vector<double>>& phi) {
		data.x.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			data.x[i] = (i + 0.5) * dx;
		}
		data.y.resize(Ny);
		for (int i = 0; i < Ny; i++) {
			data.y[i] = (i + 0.5) * dy;
		}
		data.phi.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			data.phi[i].resize(Ny);
		}
		A_W.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			A_W[i].resize(Ny);
		}
		A_S.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			A_S[i].resize(Ny);
		}
		A_P.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			A_P[i].resize(Ny);
		}
		A_N.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			A_N[i].resize(Ny);
		}
		A_E.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			A_E[i].resize(Ny);
		}
		Q_P.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			Q_P[i].resize(Ny);
		}
		phi.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			phi[i].resize(Ny);
		}
	}
	void discretize(Data& data, std::vector<std::vector<double>>& A_W, std::vector<std::vector<double>>& A_S, std::vector<std::vector<double>>& A_P, std::vector<std::vector<double>>& A_N, std::vector<std::vector<double>>& A_E, std::vector<std::vector<double>>& Q_P) {
		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {
				if (i == Nx - 1) {
					A_P[i][j] += rho * dy;
				}
				else {
					if (convection == "UDS") {
						A_P[i][j] += rho * (data.x[i] + data.x[i + 1]) * dy / 2;
					}
					else if (convection == "CDS") {
						A_P[i][j] += rho * (data.x[i] + data.x[i + 1]) * dy / 4;
						A_E[i][j] += rho * (data.x[i] + data.x[i + 1]) * dy / 4;
					}
				}
				if (i > 0) {
					if (convection == "UDS") {
						A_W[i][j] += -rho * (data.x[i - 1] + data.x[i]) * dy / 2;
					}
					else if (convection == "CDS") {
						A_W[i][j] += -rho * (data.x[i - 1] + data.x[i]) * dy / 4;
						A_P[i][j] += -rho * (data.x[i - 1] + data.x[i]) * dy / 4;
					}
				}
				if (j < Ny - 1) {
					if (convection == "UDS") {
						A_N[i][j] += -rho * (data.y[j] + data.y[j + 1]) * dx / 2;
					}
					else if (convection == "CDS") {
						A_P[i][j] += -rho * (data.y[j] + data.y[j + 1]) * dx / 4;
						A_N[i][j] += -rho * (data.y[j] + data.y[j + 1]) * dx / 4;
					}
				}
				if (j > 0) {
					if (convection == "UDS") {
						A_P[i][j] += rho * (data.y[j - 1] + data.y[j]) * dx / 2;
					}
					else if (convection == "CDS") {
						A_S[i][j] += rho * (data.y[j - 1] + data.y[j]) * dx / 4;
						A_P[i][j] += rho * (data.y[j - 1] + data.y[j]) * dx / 4;
					}
				}
				if (i < Nx - 1) {
					A_P[i][j] += Gamma * dy / dx;
					A_E[i][j] += -Gamma * dy / dx;
				}
				if (i == 0) {
					A_P[i][j] += 2 * Gamma * dy / dx;
					Q_P[i][j] += 2 * Gamma * dy / dx * (1 - data.y[j]);
				}
				else {
					A_W[i][j] += -Gamma * dy / dx;
					A_P[i][j] += Gamma * dy / dx;
				}
				if (j == Ny - 1) {
					A_P[i][j] += 2 * Gamma * dx / dy;
				}
				else {
					A_P[i][j] += Gamma * dx / dy;
					A_N[i][j] += -Gamma * dx / dy;
				}
				if (j > 0) {
					A_S[i][j] += -Gamma * dx / dy;
					A_P[i][j] += Gamma * dx / dy;
				}
			}
		}
	}
	void SOR(Data& data, std::vector<std::vector<double>>& A_W, std::vector<std::vector<double>>& A_S, std::vector<std::vector<double>>& A_P, std::vector<std::vector<double>>& A_N, std::vector<std::vector<double>>& A_E, std::vector<std::vector<double>>& Q_P, std::vector<std::vector<double>>& phi) {
		while (true) {
			for (int i = 0; i < Nx; i++) {
				for (int j = 0; j < Ny; j++) {
					phi[i][j] = Q_P[i][j] - A_P[i][j] * data.phi[i][j];
					if (i < Nx - 1) {
						phi[i][j] -= A_E[i][j] * data.phi[i + 1][j];
					}
					if (i > 0) {
						phi[i][j] -= A_W[i][j] * phi[i - 1][j];
					}
					if (j < Ny - 1) {
						phi[i][j] -= A_N[i][j] * data.phi[i][j + 1];
					}
					if (j > 0) {
						phi[i][j] -= A_S[i][j] * phi[i][j - 1];
					}
					phi[i][j] *= omega / A_P[i][j];
					phi[i][j] += data.phi[i][j];
				}
			}
			double res = residual(data, phi);
			data.phi = phi;
			if (res < tol) {
				break;
			}
		}
	}
	double residual(Data& data, std::vector<std::vector<double>>& phi) {
		double res = 0.0;
		for (int i = 0; i < Nx; i++) {
			for (int j = 0; j < Ny; j++) {
				res = res < abs(data.phi[i][j] - phi[i][j]) ? abs(data.phi[i][j] - phi[i][j]) : res;
			}
		}
		return res;
	}
public:
	ScalarTransportInAKnownVelocityField(double rho, double Gamma, int Nx, int Ny, std::string convection, double tol, double omega) :rho(rho), Gamma(Gamma), Nx(Nx), Ny(Ny), dx(1.0 / Nx), dy(1.0 / Ny), convection(convection), tol(tol), omega(omega) {}
	void solve(Data& data) {
		std::vector<std::vector<double>> A_W, A_S, A_P, A_N, A_E, Q_P, phi;
		initialize(data, A_W, A_S, A_P, A_N, A_E, Q_P, phi);
		discretize(data, A_W, A_S, A_P, A_N, A_E, Q_P);
		SOR(data, A_W, A_S, A_P, A_N, A_E, Q_P, phi);
	}
};
