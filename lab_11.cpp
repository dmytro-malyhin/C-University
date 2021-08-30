
#include <iostream>
#include <math.h>
#include <fstream>
#include <fstream>

const double m = 0.35; 


class InputData
{
public:
	int radius;
	int	height;
	double pressure;

	InputData(int radius, int height, double pressure)
	{
		this->height = height;
		this->pressure = pressure;
		this->radius = radius;
	}
};

class Fixed
{
public:
	virtual double tangentMoment(double p, double m, double R, double r) {
		return p;
	}

	virtual double radialMoment(double p, double m, double R, double r) {
		return p;
	}

	virtual double deflection(double p, double R, double r, double v, double E, double height, double m) {
		return p;
	}

	double tangentStraining(double height, double z, double p, double m, double R, double r) {
		return tangentMoment(p, m, R, r)*z / (pow(height, 3) / 12);
	}

	double radialStraining(double height, double z, double p, double m, double R, double r) {
		return radialMoment(p, m, R, r)*z / (pow(height, 3) / 12);
	}

	double cylindricalStiffness(double E, double height, double m) {
		return (E*pow(height, 3))*(12 * (1 - pow(m, 2)));
	}

};

class HardFixed : public Fixed
{
public:

	double tangentMoment(double p, double m, double R, double r)
	{
		return (p / 16)*((1 + m)*pow(R, 2) - (3 + m)*pow(r, 2));
	}

	double radialMoment(double p, double m, double R, double r)
	{
		return (p / 16)*((1 + m)*pow(R, 2) - (1 + 3 * m)*pow(r, 2));
	}

	double deflection(double p, double R, double r, double v, double E, double height, double m) {
		return (p * pow((pow(R, 2) - pow(r, 2)), 2)) / (64 * cylindricalStiffness(E, height, m));
	}

};

class HingeFixed : public Fixed
{
public:

	double tangentMoment(double p, double m, double R, double r)
	{
		return (((3 + m)*p) / 16)*(pow(R, 2) - pow(r, 2));
	}

	double radialMoment(double p, double m, double R, double r)
	{
		return (p / 16)*((3 + m)*pow(R, 2) - (1 + 3 * m)*pow(r, 2));
	}

	double deflection(double p, double R, double r, double v, double E, double height, double m) {
		return (p*pow(r, 4)) / (64 * cylindricalStiffness(E, height, m)) -
			(((3 + v) / (1 + v)) * ((p*pow(R, 2) * pow(r, 2)) / (32 * cylindricalStiffness(E, height, m)))) +
			(((5 + v) / (1 + v)) * ((p*pow(R, 4)) / (64 * cylindricalStiffness(E, height, m))));
	}

};



int main()
{
	double r, p, m, R, h, E;
	double points = 11;
	std::ifstream input("D:\\Games\\Lab11\\input.txt");
	std::ofstream output("D:\\Games\\Lab11\\output.txt");
	input >> r >> p >> m >> R >> h >> E;
	input.close();
	double z = R / 2;
	double delta = (R - r) / ((int)points);
	double max = -1;
	double between = r;
	Fixed fixed;
	HardFixed hardFixed;
	HingeFixed hingeFixed;
	std::cout << "hingeFixed:	\n";
	std::cout << "r" << "		tangentStraining" << "		radialStraining" << "			deflection				MAX\n";
	for (int i = 0; i < points; i++) {
		if (hingeFixed.deflection(p, R, r, m, E, h, m) > max) {
			max = hingeFixed.deflection(p, R, r, m, E, h, m);
		}
		std::cout << r << "		" << hingeFixed.tangentStraining(h, z, p, m, R, r) << "			" << hingeFixed.radialStraining(h, z, p, m, R, r) << "			" << hingeFixed.deflection(p, R, r, m, E, h, m) << "				" << max << std::endl;
		output << r << "		" << hingeFixed.tangentStraining(h, z, p, m, R, r) << "			" << hingeFixed.radialStraining(h, z, p, m, R, r) << "			" << hingeFixed.deflection(p, R, r, m, E, h, m) << "				" << max << std::endl;
		r += delta;
	}
	output << "\n\n";
	r = between;
	max = -1;
	std::cout << "\nhardFixed:	\n";
	std::cout << "r" << "		tangentStraining" << "		radialStraining" << "			deflection				MAX\n";
	for (int i = 0; i < points; i++) {
		if (hardFixed.deflection(p, R, r, m, E, h, m) > max) {
			max = hardFixed.deflection(p, R, r, m, E, h, m);
		}
		std::cout << r << "		" << hardFixed.tangentStraining(h, z, p, m, R, r) << "			" << hardFixed.radialStraining(h, z, p, m, R, r) << "			" << hardFixed.deflection(p, R, r, m, E, h, m) << "			" << max << std::endl;
		output << r << "		" << hardFixed.tangentStraining(h, z, p, m, R, r) << "			" << hardFixed.radialStraining(h, z, p, m, R, r) << "			" << hardFixed.deflection(p, R, r, m, E, h, m) << "			" << max << std::endl;
		r += delta;
	}
//	std::cout << "\n" << half_flate <<"\n";
//	std::cout << r << "	" << p << "	" << m << "	" << R << "	" << h << "	" << E;
	output.close();
	std::cout << "Hello World!\n";
	system("pause");
	return 200;
}



