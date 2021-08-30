#include "pch.h"
#include <iostream>
#include "math.h"
#include <fstream>
#include <iomanip>	
using namespace std;


int main()
{
	double r, rc, vr1, vr2, r1, r2, E, k1, k2, d, pc, p1, p2, vp1;
	double maxsi1 = 1000, maxsi2 = 1001, max = 0;
	int n, n1;
//	cout << "Vvedit kilkist tochok" << endl;
//	cin >> n;
	ifstream inp;
	inp.open("D:\\Games\\lab9\\data.txt");
	inp >> vr1 >> vr2 >> E >> vp1 >> n;
	inp.close();
	ofstream out;
	out.open("D:\\Games\\lab9\\result.txt");
	double rs = (vr2 - vr1) / (n - 1);
	d = 0.035;
	n1 = n / 2 + 1;
	rc = vr1 + (vr2 - vr1) / 2;
	k1 = vr1 / rc;
	k2 = rc / vr2;
	double *sir = new double[n];
	double *sit = new double[n];
	double *siec = new double[n];
	double *sir1 = new double[n];
	double *sit1 = new double[n];
	double *siec1 = new double[n];
	double *sirs = new double[n];
	double *sits = new double[n];
	double *siecs = new double[n];
	while (maxsi1 < maxsi2)
	{
		maxsi2 = maxsi1;
		pc = d * E*(1 - pow(k1, 2)*(1 - pow(k2, 2))) / 2 / rc / ((1 + pow(k1, 2))*(1 - pow(k2, 2)) + (1 + pow(k2, 2))*(1 - pow(k1, 2)));
		r1 = vr1;
		r2 = rc;
		p1 = 0;
		p2 = pc;
		for (int i = 0; i < n1; i++)
		{
			r = r1 + rs * i;
			sir[i] = (pow(r1, 2)*p1 - pow(r2, 2)*p2) / (pow(r2, 2) - pow(r1, 2)) - pow(r1, 2)*pow(r2, 2)*(p1 - p2) / (pow(r2, 2) - pow(r1, 2)) / pow(r, 2);
			sit[i] = (pow(r1, 2)*p1 - pow(r2, 2)*p2) / (pow(r2, 2) - pow(r1, 2)) + pow(r1, 2)*pow(r2, 2)*(p1 - p2) / (pow(r2, 2) - pow(r1, 2)) / pow(r, 2);
			siec[i] = sqrt(pow(sit[i], 2) - sit[i] * sir[i] + pow(sit[i], 2)) / sqrt(2);

		}
		r1 = rc;
		r2 = vr2;
		p1 = pc;
		p2 = 0;
		for (int i = 0; i < n1; i++)
		{
			r = r1 + rs * i;
			sir[i + n1 - 1] = (pow(r1, 2)*p1 - pow(r2, 2)*p2) / (pow(r2, 2) - pow(r1, 2)) - pow(r1, 2)*pow(r2, 2)*(p1 - p2) / (pow(r2, 2) - pow(r1, 2)) / pow(r, 2);
			sit[i + n1 - 1] = (pow(r1, 2)*p1 - pow(r2, 2)*p2) / (pow(r2, 2) - pow(r1, 2)) + pow(r1, 2)*pow(r2, 2)*(p1 - p2) / (pow(r2, 2) - pow(r1, 2)) / pow(r, 2);
			siec[i + n1 - 1] = sqrt(pow(sit[i + n1 - 1], 2) - sit[i + n1 - 1] * sir[i + n1 - 1] + pow(sit[i + n1 - 1], 2)) / sqrt(2);

		}
		r1 = vr1;
		r2 = vr2;
		p1 = vp1;
		p2 = 0;
		for (int i = 0; i < n; i++)
		{
			r = r1 + rs * i;
			sir1[i] = (pow(r1, 2)*p1 - pow(r2, 2)*p2) / (pow(r2, 2) - pow(r1, 2)) - pow(r1, 2)*pow(r2, 2)*(p1 - p2) / (pow(r2, 2) - pow(r1, 2)) / pow(r, 2);
			sit1[i] = (pow(r1, 2)*p1 - pow(r2, 2)*p2) / (pow(r2, 2) - pow(r1, 2)) + pow(r1, 2)*pow(r2, 2)*(p1 - p2) / (pow(r2, 2) - pow(r1, 2)) / pow(r, 2);
			siec1[i] = sqrt(pow(sit1[i], 2) - sit1[i] * sir1[i] + pow(sit1[i], 2)) / sqrt(2);

		}
		max = 0;
		for (int i = 0; i < n; i++)
		{
			sirs[i] = sir[i] + sir1[i];
			sits[i] = sit[i] + sit1[i];
			siecs[i] = sqrt(pow(sits[i], 2) - sits[i] * sirs[i] + pow(sits[i], 2)) / sqrt(2);
			if (siecs[i] > max)
			{
				max = siecs[i];
			}
			maxsi1 = max;
		}


		out.precision(4);
		out << left << setw(15) << d << left << setw(15) << maxsi1 << endl;
		cout << left << setw(15) << d << left << setw(15) << maxsi1 << endl;
		d = d + 0.000001;


		if (maxsi1 >= maxsi2)
		{
			out << "Napruzhenna vid dii vnytrishnogo tusky" << endl;
			cout << "Napruzhenna vid dii vnytrishnogo tusky" << endl;
			for (int i = 0; i < n; i++)
			{
				r = r1 + rs * i;
				out.precision(4);
				out << left << setw(15) << r << left << setw(15) << sir1[i] << left << setw(15) << sit1[i] << left << setw(15) << siec1[i] << endl;
				cout << left << setw(15) << r << left << setw(15) << sir1[i] << left << setw(15) << sit1[i] << left << setw(15) << siec1[i] << endl;
			}
			out << "Sumarni napruzhenna" << endl;
			cout << "Sumarni napruzhenna" << endl;
			for (int i = 0; i < n; i++)
			{
				out.precision(4);
				out << left << setw(15) << sirs[i] << left << setw(15) << sits[i] << left << setw(15) << siecs[i] << endl;
				cout << left << setw(15) << sirs[i] << left << setw(15) << sits[i] << left << setw(15) << siecs[i] << endl;
			}
		}
	}
	out.close();
	delete[] sir;
	delete[] sit;
	delete[] siec;
	delete[] sir1;
	delete[] sit1;
	delete[] siec1;
	delete[] sirs;
	delete[] sits;
	delete[] siecs;
}
