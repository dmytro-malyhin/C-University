//#include "stdafx.h"
#include <iostream>
#include "math.h"
#include <fstream>
#include <iomanip>

using namespace std;

int main()
{
	int n = 0, n1;
	//cout << "How much dots?" << endl;
	//cin >> n;//const int n=13
	
	double R1, R2, P1, P2, ModUng, Pua, Ro, PIns, POut;
	ifstream inp;
	inp.open("D:\\C++\\TNN_DIMA\\EightLab\\lab08\\Debug\\data.txt");
	inp >> R1 >> R2 >> PIns >> ModUng >> Pua >> n;
	inp.close();
	double *MasRad = new double[n];
	double *MasTan = new double[n];
	double *MasRad1 = new double[n];
	double *MasTan1 = new double[n];
	double *MasSumEcv = new double[n];
	double *MasSumRad = new double[n];
	double *MasSumTan = new double[n];
	ofstream out;
	out.open("D:\\C++\\TNN_DIMA\\EightLab\\lab08\\Debug\\result.txt");
	n1 = n / 2 + 1;
	double DR = (R2 - R1) / (n - 1);
	double Rc = R1 + (R2 - R1) / 2;
	double R;
	out.precision(4);
	out << left << setw(20) << "R" << left << setw(20) << "MasSumRad[i]" << left << setw(20) << "MasSumTan[i]" << left << setw(20) << "MasSumEcv[i]" << left << setw(20) << "MasRad[i]" << left << setw(20) << "MasTan[i]" << left << setw(20) << "MasRad1[i]" << left << setw(20) << "MasTan1[i]" << endl;
	Ro = (R2 - R1) / 100;
	double MaxEcvMove1 = 0;
	double MaxEcvMove2 = 0;
	double MaxEcvMove3 = 0;
	double MaxEcvMove4 = 0;
	double k1 = R1 / Rc;
	double k2 = Rc / R2;
	double Pc = Ro * ModUng / 2 / Rc * ((1 - pow(k1, 2))*(1 - pow(k2, 2))) / ((1 + pow(k1, 2))*(1 - pow(k2, 2)) + (1 + pow(k2, 2))*(1 - pow(k1, 2)));
	double RIns, ROut;
	RIns = R1;
	ROut = Rc;
	P1 = 0;
	P2 = Pc;
	for (int i = 0; i < n1; i++)
	{
		R = R1 + DR * i;
		MasRad[i] = ((pow(RIns, 2)*P1 - pow(ROut, 2)*P2) / (pow(ROut, 2) - pow(RIns, 2))) - (pow(RIns, 2)*pow(ROut, 2)*(P1 - P2) / (pow(ROut, 2) - pow(RIns, 2))*(1 / pow(R, 2)));
		MasTan[i] = ((pow(RIns, 2)*P1 - pow(ROut, 2)*P2) / (pow(ROut, 2) - pow(RIns, 2))) + (pow(RIns, 2)*pow(ROut, 2)*(P1 - P2) / (pow(ROut, 2) - pow(RIns, 2))*(1 / pow(R, 2)));
	}
	RIns = Rc;
	ROut = R2;
	P1 = Pc;
	P2 = 0;
	for (int i = 0; i < (n1); i++)
	{
		R = R1 + DR * (n1 - 1) + DR * i;
		MasRad[i + n1 - 1] = ((pow(RIns, 2)*P1 - pow(ROut, 2)*P2) / (pow(ROut, 2) - pow(RIns, 2))) - (pow(RIns, 2)*pow(ROut, 2)*(P1 - P2) / (pow(ROut, 2) - pow(RIns, 2))*(1 / pow(R, 2)));
		MasTan[i + n1 - 1] = ((pow(RIns, 2)*P1 - pow(ROut, 2)*P2) / (pow(ROut, 2) - pow(RIns, 2))) + (pow(RIns, 2)*pow(ROut, 2)*(P1 - P2) / (pow(ROut, 2) - pow(RIns, 2))*(1 / pow(R, 2)));
	}
	RIns = R1;
	ROut = R2;
	P1 = PIns;
	P2 = 0;
	for (int i = 0; i < n; i++)
	{
		R = R1 + DR * i;
		MasRad1[i] = ((pow(RIns, 2)*P1 - pow(ROut, 2)*P2) / (pow(ROut, 2) - pow(RIns, 2))) - (pow(RIns, 2)*pow(ROut, 2)*(P1 - P2) / (pow(ROut, 2) - pow(RIns, 2))*(1 / pow(R, 2)));
		MasTan1[i] = ((pow(RIns, 2)*P1 - pow(ROut, 2)*P2) / (pow(ROut, 2) - pow(RIns, 2))) + (pow(RIns, 2)*pow(ROut, 2)*(P1 - P2) / (pow(ROut, 2) - pow(RIns, 2))*(1 / pow(R, 2)));
		MasSumEcv[i] = 1 / sqrt(2.) * sqrt(pow((MasTan[i] + MasTan1[i]), 2) - (MasTan[i] + MasTan1[i]) * (MasRad[i] + MasRad1[i]) + pow((MasRad[i] + MasRad1[i]), 2));
		MasSumRad[i] = MasRad[i] + MasRad1[i];
		MasSumTan[i] = MasTan[i] + MasTan1[i];
		if (MasSumEcv[i] > MaxEcvMove3)
		{
			MaxEcvMove3 = MasSumEcv[i];
		}
		out << left << setw(20) << R << left << setw(20) << MasSumRad[i] << left << setw(20) << MasSumTan[i] << left << setw(20) << MasSumEcv[i] << left << setw(20) << MasRad[i] << left << setw(20) << MasTan[i] << left << setw(20) << MasRad1[i] << left << setw(20) << MasTan1[i] << endl;

	}
	for (int i = 0; i < (n1); i++)
	{
		MasSumEcv[i + n1 - 1] = 1 / sqrt(2.) * sqrt(pow((MasTan[i] + MasTan1[i]), 2) - (MasTan[i] + MasTan1[i]) * (MasRad[i] + MasRad1[i]) + pow((MasRad[i] + MasRad1[i]), 2));
		if (MasSumEcv[i + n1 - 1] > MaxEcvMove4)
		{
			MaxEcvMove4 = MasSumEcv[i + n1 - 1];
		}

	}
	if (MaxEcvMove3 > MaxEcvMove4)
	{
		out << "ћаксимальные напр€жени€ на внутренний стороне" << endl;
	}
	else
	{
		out << "ћаксимальные напр€жени€ на внешней стороне" << endl;
	}

	out.close();
	return 0;

}
