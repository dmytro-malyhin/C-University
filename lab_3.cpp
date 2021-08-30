#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

const double alpha = 17.06 * pow(10, -6); 		// alpha for temperature 100 C
const double T1 = 100;
const double T2 = 0;

//const double r1 = 0.015;	//mm
//const double r2 = 0.020;	//mm
//const double p1 = 15 * pow(10, 6);	// MPa pressure number 1
//const double p2 = 0;	// MPa pressure number 2
//const double elasticity = 100000 * pow(10,6);
//const double coefPua = 0.35;
//const double points = 11;

/**********************************************UTILITY*FUNCTIONS*********************************************/
double T_class (double r, double r1, double r2){						// T(r) for classical variant
	double res = (T1 - T2) * (r2 - r)/(r2 - r1);						
	return res;															
}

double T_log (double r, double r1, double r2){							// T(r) for logarithmical variant
	double res = ((T1 - T2)/log(r2/r1)) * log(r2/r);					
	return res;
}

double module_class (double elasticity, double coefPua, double r1, double r2){		// first (class.) module. Was made only for convinient calculatings
	double res = (elasticity*alpha*(T1 - T2))/(3*(1-coefPua) * (r2 - r1));
	return res;
}

double module_log (double elasticity, double coefPua, double r1, double r2){		// first (log.) module. Was made only for convinient calculatings.
	double res = (-1)*(elasticity*alpha*(T1 - T2))/(2*(1-coefPua) * log(r2/r1));
	return res;
}
/*******************************************ENDING OF UTILITY FUNCTIONS************************************************/
//
//
/**********************************BEGINING OF MAIN FUNCTIONS (CLASSICAL VARIANTS)*****************************************/
double rad_tension_class (double r, double r1, double r2, double elasticity, double coefPua){		
	double delta = module_class(elasticity, coefPua, r1, r2);
	double delta1 = (r - (pow(r1, 3))/(r*r) - (1 - (r1*r1)/(r*r))*(pow(r2, 3) - pow(r1, 3))/(r2*r2 - r1*r1));
	//cout << endl << delta1 << endl << endl;
	return delta * delta1;
}

double tan_tension_class (double r, double r1, double r2, double elasticity, double coefPua){
	double delta = module_class(elasticity, coefPua, r1, r2);
	double delta1 = (2*r + (pow(r1, 3))/(r*r) - (1 + (r1*r1)/(r*r))*(pow(r2, 3) - pow(r1, 3))/(r2*r2 - r1*r1));
	//cout << endl << delta1 << endl << endl;
	return delta * delta1;
}

double axis_tension_class (double r, double r1, double r2, double elasticity, double coefPua){
	double delta = module_class(elasticity, coefPua, r1, r2);
	double delta1 = (3*r - (2*(pow(r2, 3) - pow(r1, 3)))/(r2*r2 - r1*r1));
	//cout << endl << delta1 << endl << endl;
	return delta * delta1;
}
/**********************************ENDING OF MAIN FUNCTIONS (CLASSICAL VARIANTS)*****************************************/
//
//
//
/*****************************************BEGINING OF MAIN FUNCTIONS (LOG VARIANTS)*************************************/
double rad_tension_log (double r, double r1, double r2, double elasticity, double coefPua){		
	double delta = module_log(elasticity, coefPua, r1, r2);
	double delta1 = (log(r2/r) + ((r1*r1)/(r2*r2 - r1*r1)) * (1 - (r2*r2)/(r*r)) * log(r2/r1));
	//cout << endl << delta1 << endl << endl;
	return delta * delta1;
}

double tan_tension_log (double r, double r1, double r2, double elasticity, double coefPua){
	double delta = module_log(elasticity, coefPua, r1, r2);
	double delta1 = (1 - log(r2/r) - ((r1*r1)/(r2*r2 - r1*r1)) * (1 + (r2*r2)/(r*r)) * log(r2/r1));
	//cout << endl << delta1 << endl << endl;
	return delta * delta1;
}

double axis_tension_log (double r, double r1, double r2, double elasticity, double coefPua){
	double delta = module_log(elasticity, coefPua, r1, r2);
	double delta1 = (1 - 2*log(r2/r) - (2*(r1*r1)/(r2*r2 - r1*r1)) * log(r2/r1));
	//cout << endl << delta1 << endl << endl;
	return delta * delta1;
}
/**********************************ENDING OF MAIN FUNCTIONS (LOG VARIANTS)*****************************************/

//double moving (double r, double r1, double r2, double p1, double p2, double coefPua, double elasticity){
//	double res1 = (1 - coefPua)/elasticity;
//	double res2 = (1 + coefPua)/elasticity;
//	double delta = ((pow(r1, 2)*p1) - (pow(r2, 2)*p2))/(r2*r2-r1*r1);
//	double delta1 = ((pow(r1, 2) * pow(r2, 2) * (p1 - p2))/(r2*r2-r1*r1));
//	double delta2 = (1/(r*r));
//	return res1 * delta + (delta1 * delta2) * res2;
//}
//
//double equivalent_tension (double ten1, double ten2){
//	double res = ten2 * ten2 - ten1 * ten2 + ten1 * ten1;
//	double result = 1/(pow(2, 0.5)) * pow(res, 0.5);
//	return result;
//}


int main(){
	double r1, r2, p1, p2, elasticity, coefPua, points;
	double element[7];
	std::ifstream in("input.txt");
	std::ofstream out("output.txt");
	for (int i = 0; i < 7; i++){
		in >> element[i];
	}
//	in >> r1;
//	in >> r2;
//	in >> p1;
//	in >> p2;
//	in >> elasticity;
//	in >> coefPua;
//	in >> points;
	in.close();
	std::cout << "Calculating for class variant:\n";
	std::cout << "	Radial tension is : " << rad_tension_class((element[0]+element[1])/2, element[0], element[1], element[4], element[5]) << std::endl;
	std::cout << "	Tangential tension is : " << tan_tension_class((element[0]+element[1])/2, element[0], element[1], element[4], element[5]) << std::endl;
	std::cout << "	Axis tension is : " << axis_tension_class((element[0]+element[1])/2, element[0], element[1], element[4], element[5]) << std::endl;
	std::cout << "Calculating for log variant:\n";
	std::cout << "	Radial tension is : " << rad_tension_log((element[0]+element[1])/2, element[0], element[1], element[4], element[5]) << std::endl;
	std::cout << "	Tangential tension is : " << tan_tension_log((element[0]+element[1])/2, element[0], element[1], element[4], element[5]) << std::endl;
	std::cout << "	Axis tension is : " << axis_tension_log((element[0]+element[1])/2, element[0], element[1], element[4], element[5]) << std::endl;
	out << "========================================================================================================================================================================================================" << std::endl;
	out << "Radius			Rad tension			Tang tension			Axis tension			Log Rad tension			Log Tang tension		Log Axis tension" << std::endl;
	double radius = element[0];
	double equi = 0;
	for (int i = 0; i < element[6]; i++){
//		if (radius!=0.020 && radius!=0.015){
//			out << radius <<"		" << rad_tension_class(radius, r1, r2, elasticity, coefPua) << "			"<<tan_tension_class(radius, r1, r2, elasticity, coefPua)<<"				"<<axis_tension_class(radius, r1, r2, elasticity, coefPua)<<"		"<<rad_tension_log(radius, r1, r2, elasticity, coefPua) << "			" << tan_tension_log(radius, r1, r2, elasticity, coefPua) << "			" << axis_tension_class(radius, r1, r2, elasticity, coefPua) << std::endl;
//		}
		if(radius==0.020){
			out << radius <<"			" << rad_tension_class(radius, element[0], element[1], element[4], element[5]) << "				"<<tan_tension_class(radius, element[0], element[1], element[4], element[5])<<"			"<<axis_tension_class(radius, element[0], element[1], element[4], element[5])<<"			"<<rad_tension_log(radius, element[0], element[1], element[4], element[5]) <<"				" << tan_tension_log(radius, element[0], element[1], element[4], element[5]) << "			" << axis_tension_class(radius, element[0], element[1], element[4], element[5]) << std::endl;
		}
		else if(radius==0.015){
			out << radius <<"			" << rad_tension_class(radius, element[0], element[1], element[4], element[5]) << "				"<<tan_tension_class(radius, element[0], element[1], element[4], element[5])<<"			"<<axis_tension_class(radius, element[0], element[1], element[4], element[5])<<"			"<<rad_tension_log(radius, element[0], element[1], element[4], element[5]) <<"				" << tan_tension_log(radius, element[0], element[1], element[4], element[5]) << "			" << axis_tension_class(radius, element[0], element[1], element[4], element[5]) << std::endl;
		}

//		if (equivalent_tension (rad_tension(radius, r1, r2, p1, p2), tan_tension(radius, r1, r2, p1, p2)) > equi){
//			equi = equivalent_tension (rad_tension(radius, r1, r2, p1, p2), tan_tension(radius, r1, r2, p1, p2));
//		}
		if (radius != 0.020 && radius != 0.015){
			out << radius <<"		" << rad_tension_class(radius, element[0], element[1], element[4], element[5]) << "			"<<tan_tension_class(radius, element[0], element[1], element[4], element[5])<<"			"<<axis_tension_class(radius, element[0], element[1], element[4], element[5])<<"			"<<rad_tension_log(radius, element[0], element[1], element[4], element[5]) <<"			" << tan_tension_log(radius, element[0], element[1], element[4], element[5]) << "			" << axis_tension_class(radius, element[0], element[1], element[4], element[5]) << std::endl;
		}
		
		radius+=(0.005/9);
	}
//	out << std::endl << std::endl << "Max equivalent tension using 4th theory of strength is " << equi;
	out.close();
}
