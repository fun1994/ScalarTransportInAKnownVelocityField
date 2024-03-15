#include "ScalarTransportInAKnownVelocityField.h"

void test(std::string Gamma, int N, std::string convection, double omega) {
	std::cout << "Gamma=" << Gamma << " " << "N=" << N << " " << "convection=" << convection << std::endl;
	ScalarTransportInAKnownVelocityField STIAKVF(1.0, std::stod(Gamma), N, N, convection, 1e-8, omega);
	Data data;
	STIAKVF.solve(data);
	data.save("Gamma=" + Gamma + "/N=" + std::to_string(N), convection);
}

void test(std::string Gamma, std::string convection, double omega) {
	test(Gamma, 40, convection, omega);
	test(Gamma, 80, convection, omega);
	test(Gamma, 160, convection, omega);
	test(Gamma, 320, convection, omega);
}

void test() {
	test("0.01", "UDS", 1.0);
	test("0.01", "CDS", 1.0);
	test("0.001", "UDS", 1.0);
	test("0.001", "CDS", 0.1);
}

int main() {
	test();
	return 0;
}
