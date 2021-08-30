#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
using namespace std;

double r1;	//mm
double r2;	//mm
double p1;	// MPa pressure number 1
double p2;	// MPa pressure number 2
double elasticity;
double coefPua;
double points;

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

double equivalent_tension(double r, double ten1, double ten2) {
	double res = ten2 * ten2 - ten1 * ten2 + ten1 * ten1;
	double result = 1 / (pow(2, 0.5)) * pow(res, 0.5);
	return result;
}


int main() {
	ifstream input;
	ofstream output;
	input.open("C:\\Users\\tumur\\Desktop\\ТДДТ\\Лабороторная 7\\lab007\\Debug\\data.txt", std::ifstream::in);
	input >> r1 >> r2 >> p1 >> p2 >> elasticity >> coefPua >> points;
	input.close();
	output.open("C:\\Users\\tumur\\Desktop\\ТДДТ\\Лабороторная 7\\lab007\\Debug\\result.txt");
	double Pc, delta, k1, k2, Rc, R;
	Rc = (r1 + r2) / 2;
	k1 = r1 / Rc;
	k2 = Rc / r2;
	delta = (r2 - r1)*0.005;

	Pc = (delta*elasticity / (2 * Rc))*(((1 - pow(k1, 2))*(1 - pow(k2, 2))) / ((1 + pow(k1, 2))*(1 - pow(k2, 2)) + (1 + pow(k2, 2))*(1 - pow(k1, 2))));
	//cout << endl << "Pc=   " << Pc << endl;
	int nn = points;
	int n = (nn/2) + 1;
	double deltaR1 = (Rc - r1) / (n - 1);
	double deltaR2 = (r2 - Rc) / (n - 1);
	//double deltaR = (R2 - R1) / (n - 1);
	double* SigRvnutr = new double[nn];
	double* SigTvnutr = new double[nn] ;
	double* SigEkvvnutr = new double[nn];


	double max1 = 0;
	cout << fixed << left << setw(20) << "R:" << right << setw(20) << "SigRvnutr:" << right << setw(20) << "SigTvnutr:" << right << setw(20) << "SigEkvvnutr:" << right << setw(20) << "max1:" << endl;

	for (int i = 0; i < nn; i++)
	{
		R = (r1 + i * deltaR1);
		//r2 = Rc;
		//double P2 = Pc;

		SigRvnutr[i] = rad_tension(R);
		SigTvnutr[i] = tan_tension(R);
		SigEkvvnutr[i] = equivalent_tension(R, SigTvnutr[i], SigRvnutr[i]);


		if (SigEkvvnutr[i] > max1)
		{
			max1 = SigEkvvnutr[i];
		}
		cout << fixed << left << setw(20) << R << right << setw(20) << SigRvnutr[i] << right << setw(20) << SigTvnutr[i] << right << setw(20) << SigEkvvnutr[i] << right << setw(20) << max1 << endl;
		output << fixed << left << setw(20) << R << right << setw(20) << SigRvnutr[i] << right << setw(20) << SigTvnutr[i] << right << setw(20) << SigEkvvnutr[i] << right << setw(20) << max1 << endl;

		//cout << R << "	" << SigRvnutr[i] << "	" << SigTvnutr[i] << "	" << SigEkvvnutr[i] << "	" << max1 << endl;
		//output << R << "	" << SigRvnutr[i] << "	" << SigTvnutr[i] << "	" << SigEkvvnutr[i] << "	" << max1 << endl;


	}

	double* SigRzovn = new double[nn];
	double* SigTzovn = new double[nn];
	double* SigEkvzovn = new double[nn];

	double max2 = 0;

	cout << fixed << left << setw(20) << "R:" << right << setw(20) << "SigRzovn:" << right << setw(20) << "SigTzovn:" << right << setw(20) << "SigEkvzovn:" << right << setw(20) << "max2:" << endl;

	r1 = 0.015;
	r2 = 0.020;
	Rc = (r1 + r2) / 2;
	//cout << endl << "Rc: " << Rc << endl;


	for (int f = 0; f < nn; f++)
	{
		r1 = Rc;
		R = (r1 + f * deltaR2);
		p1 = Pc;
		double P2 = 0;

		SigRzovn[f] = rad_tension(R);
		SigTzovn[f] = tan_tension(R);
		SigEkvzovn[f] = equivalent_tension(R, SigRzovn[f], SigTzovn[f]);

		if (SigEkvzovn[f] > max2)
		{
			max2 = SigEkvzovn[f];
		}

		cout << fixed << left << setw(20) << R << right << setw(20) << SigRzovn[f] << right << setw(20) << SigTzovn[f] << right << setw(20) << SigEkvzovn[f] << right << setw(20) << max2 << endl;
		output << fixed << left << setw(20) << R << right << setw(20) << SigRzovn[f] << right << setw(20) << SigTzovn[f] << right << setw(20) << SigEkvzovn[f] << right << setw(20) << max2 << endl;

	}



	output.close();
	system("PAUSE");
	return 200;
}