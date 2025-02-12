﻿// ConsoleApplication2.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"
#include <iostream>
#include <math.h>

const double r1 = 0.015;	//mm
const double r2 = 0.020;	//mm
const double p1 = 15 * pow(10, 6);	// MPa pressure number 1
const double p2 = 0;	// MPa pressure number 2
const double elasticity = 100000 * pow(10, 6);
const double coefPua = 0.35;
const double points = 11;

double rad_tension(double r) {
	double delta = ((pow(r1, 2)*p1) - (pow(r2, 2)*p2)) / (r2*r2 - r1 * r1);
	double delta1 = ((pow(r1, 2) * pow(r2, 2) * (p1 - p2)) / (r2*r2 - r1 * r1));
	double delta2 = (1 / (r*r));
	return delta - (delta1 * delta2);
}

double tan_tension(double r) {
	double delta = ((pow(r1, 2)*p1) - (pow(r2, 2)*p2)) / (r2*r2 - r1 * r1);
	double delta1 = ((pow(r1, 2) * pow(r2, 2) * (p1 - p2)) / (r2*r2 - r1 * r1));
	double delta2 = (1 / (r*r));
	return delta + (delta1 * delta2);
}

double moving(double r) {
	double res1 = (1 - coefPua) / elasticity;
	double res2 = (1 + coefPua) / elasticity;
	double delta = ((pow(r1, 2)*p1) - (pow(r2, 2)*p2)) / (r2*r2 - r1 * r1);
	double delta1 = ((pow(r1, 2) * pow(r2, 2) * (p1 - p2)) / (r2*r2 - r1 * r1));
	double delta2 = (1 / (r*r));
	return res1 * delta + (delta1 * delta2) * res2;
}

double equivalent_tension(double r, double ten1, double ten2) {
	double res = ten2 * ten2 - ten1 * ten2 + ten1 * ten1;
	double result = 1 / (pow(2, 0.5)) * pow(res, 0.5);
	return result;
}


int main() {
	std::cout << "Calculating for external radius:\n";
	std::cout << "	Radial tension is : " << rad_tension(r2) << std::endl;
	std::cout << "	Tangential tension is : " << tan_tension(r2) << std::endl;
	std::cout << "	Moving is : " << moving(r2) << std::endl;
	std::cout << "Calculating for internal radius:\n";
	std::cout << "	Radial tension is : " << rad_tension(r1) << std::endl;
	std::cout << "	Tangential tension is : " << tan_tension(r1) << std::endl;
	std::cout << "	Moving is : " << moving(r1) << std::endl;
	std::cout << "==============================================================================================================" << std::endl;
	std::cout << "Radius		Radial tension		Tangential tension		Moving			Constant" << std::endl;
	double radius = r1;
	for (int i = 0; i < points; i++) {
		if (radius == 0.015) {
			std::cout << radius << "		" << rad_tension(radius) << "		" << tan_tension(radius) << "			" << moving(radius) << "		" << rad_tension(radius) + tan_tension(radius) << std::endl;
		}
		if (radius == 0.020) {
			std::cout << radius << "		" << rad_tension(radius) << "			" << tan_tension(radius) << "			" << moving(radius) << "		" << rad_tension(radius) + tan_tension(radius) << std::endl;
		}
		if (radius != 0.020 && radius != 0.015) {
			std::cout << radius << "	" << rad_tension(radius) << "		" << tan_tension(radius) << "			" << moving(radius) << "		" << rad_tension(radius) + tan_tension(radius) << std::endl;
		}
		radius += (0.005 / 9);
	}
}
