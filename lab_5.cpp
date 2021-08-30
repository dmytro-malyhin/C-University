#include <iostream>
#include <math.h>
#include <fstream>
#define L (p1 - p2 + ((3+coefPua)/8)*specmassmat*4*10*frequency*frequency*(r1*r1 - r2*r2) - deltaT*alpha*elasticity/3)/(1/(r1*r1) - 1/(r2*r2)) 	// ������� ����������� L
#define K p1 - L/(r1*r1) + ((3+coefPua)/8)*specmassmat*4*10*frequency*frequency*(r1*r1) + (deltaT*alpha*elasticity*r1)/(3*(r2-r1))				// ������� ����������� K
// ����� ��������� K � L, ����������� ��������������� ��������� ������� ��������� ������������ �������, ���������������� ��� ������� ...
// ... � ���������� � ��� �������������� �� ������� ������ �������������� ��������� (������� ������������� ������ ����������)
using namespace std;

double points;	                                     // ���������� ��� ������������� �����
const double r1 = 0.075;                             // ������ R1
const double r2 = 0.350;                             // ������ R2
const double p1 = 0;	                             // ��� ������������ �� r1
const double p2 = 6 * pow(10, 6);                    // ��� ������������ �� r2
const double elasticity = 0.115 * pow(10,6);         // ����������� ������������
const double deltaT = 175;                           // T* - ������ ���������� T1 - T2	
const double alpha = 18 * pow(10, -6);               // ������������� ����������� 
const double coefPua = 0.35;                         // ����������� ��������
//const double points = 11;	                         // ����� �����
const double frequency = 50;                         // ������� ��/���
const double specmassmat = 8700;                     // �������� ����� ��������� ��/�3


double getL (){	                                     // �������-��������� ��� �������� ������������ �������� ������������ L 
	return (p1 - p2 + ((3+coefPua)/8)*specmassmat*4*10*frequency*frequency*(r1*r1 - r2*r2) - deltaT*alpha*elasticity/3)/(1/(r1*r1) - 1/(r2*r2));
}
double getK (){	                                     // �������-��������� ��� �������� ������������ �������� ������������ K
	return p1 - L/(r1*r1) + ((3+coefPua)/8)*specmassmat*4*10*frequency*frequency*(r1*r1) + (deltaT*alpha*elasticity*r1)/(3*(r2-r1));
}

double rad_tension (double r){	                     // ������� ��� ���������� ���������� �������������
	return K + L/(r*r) - ((3+coefPua)/8)*specmassmat*4*10*frequency*frequency*(r*r) - (deltaT*alpha*elasticity*r)/(3*(r2-r1));
}

double tan_tension (double r){	                     // ������� ��� ���������� �������������� �������������
	return K - L/(r*r) - ((1+3*coefPua)/8)*specmassmat*4*10*frequency*frequency*(r*r) - (2*deltaT*alpha*elasticity*r)/(3*(r2-r1));
}

double equivalent_tension (double ten1, double ten2){      // ������� ��� ���������� ������������� ������������� �� 4�� ������� ���������
	double res = ten2 * ten2 - ten1 * ten2 + ten1 * ten1;
	double result = 1/(pow(2, 0.5)) * pow(res, 0.5);
	return result;
}


int main(){
	std::ofstream out("output.txt");                 // ������������� �������� ��������� ������
	cout << " Please, type number of your points: ";
	cin >> points;                                   // ���� ����� �����
	double rad[(int)points], radt[(int)points], tant[(int)points], equiv[(int)points];	// ������������� �������� ��� �������� ������
	std::cout << "	Calculating of K constant:	(using macros)" << "		" << K << std::endl;;//	����� � ������� ��. � ��������� ������
	std::cout << "	Calculating of K constant:	(using function)" << "	" << getK() << std::endl;//	����� ... ��. � ��������� �������
	std::cout << "	Calculating of L constant:	(using macros)" << "		" << L << std::endl;//	����� ... ��. L ��������� ������
	std::cout << "	Calculating of L constant:	(using function)" << "	" << getL() << std::endl;//	����� ... ��. L ��������� �������
	std::cout << "==============================================================================================================" << std::endl;
	std::cout << "Radius		Radial tension		Tangential tension		Equivalent tension" << std::endl;
	double radius = r1;	                             // ������������� �������� ��� ������ � ��������
	for (int i = 0; i < points; i++){                // � ���������� ������ ���������� ����� ������� ��� ���������� ��������� ������� ������ ����� | ���� ��� �������� ����� �������
		rad[i] = radius;                             // ������������� i-�� �������� ������� rad[]
		radt[i] = rad_tension(radius);	             // ������������� i-�� �������� ������� radt[]
		tant[i] = tan_tension(radius);	             // ������������� i-�� �������� ������� tant[]
		equiv[i] = equivalent_tension(radt[i], tant[i]);	// ������������� i-�� �������� ������� equiv[]
		cout << rad[i] << "	" << radt[i] << "	" << tant[i] << "	" << equiv[i] << endl;	   // ����� � ������� ��������������� ��������� �������, ��� ������������� ������ �������� �������  + �������
		out << rad[i] << " " << radt[i] << " " << tant[i] << " " << equiv[i] << endl;         // ... � ���� ...
		radius+=((r2-r1)/(points-1));	             // ���������������� ���������� radius �� ��������, ������� �� ���������� �� deltaR
	}
	out.close();                                     // �������� ��������� ������
}
