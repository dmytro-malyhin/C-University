#include <iostream>
#include <math.h>
#include <fstream>

//const double r1 = 0.015;	//mm
//const double r2 = 0.020;	//mm
//const double p1 = 15 * pow(10, 6);	// MPa pressure number 1
//const double p2 = 0;	// MPa pressure number 2
//const double elasticity = 100000 * pow(10,6);
//const double coefPua = 0.35;
//const double points = 11;

double rad_tension (double r, double r1, double r2, double p1, double p2){
	double delta = ((pow(r1, 2)*p1) - (pow(r2, 2)*p2))/(r2*r2-r1*r1);
	double delta1 = ((pow(r1, 2) * pow(r2, 2) * (p1 - p2))/(r2*r2-r1*r1));
	double delta2 = (1/(r*r));
	return delta - (delta1 * delta2);
}

double tan_tension (double r, double r1, double r2, double p1, double p2){
	double delta = ((pow(r1, 2)*p1) - (pow(r2, 2)*p2))/(r2*r2-r1*r1);
	double delta1 = ((pow(r1, 2) * pow(r2, 2) * (p1 - p2))/(r2*r2-r1*r1));
	double delta2 = (1/(r*r));
	return delta + (delta1 * delta2);
}

double moving (double r, double r1, double r2, double p1, double p2, double coefPua, double elasticity){
	double res1 = (1 - coefPua)/elasticity;
	double res2 = (1 + coefPua)/elasticity;
	double delta = ((pow(r1, 2)*p1) - (pow(r2, 2)*p2))/(r2*r2-r1*r1);
	double delta1 = ((pow(r1, 2) * pow(r2, 2) * (p1 - p2))/(r2*r2-r1*r1));
	double delta2 = (1/(r*r));
	return res1 * delta + (delta1 * delta2) * res2;
}

double equivalent_tension (double ten1, double ten2){
	double res = ten2 * ten2 - ten1 * ten2 + ten1 * ten1;
	double result = 1/(pow(2, 0.5)) * pow(res, 0.5);
	return result;
}


int main(){
	double r1, r2, p1, p2, elasticity, coefPua, points;
	std::ifstream in("input.txt");
	std::ofstream out("output.txt");
	in >> r1;
	in >> r2;
	in >> p1;
	in >> p2;
	in >> elasticity;
	in >> coefPua;
	in >> points;
	in.close();
//	std::cout << "Calculating for external radius:\n";
//	std::cout << "	Radial tension is : " << rad_tension(r2, r1, r2, p1, p2) << std::endl;
//	std::cout << "	Tangential tension is : " << tan_tension(r2, r1, r2, p1, p2) << std::endl;
//	std::cout << "	Moving is : " << moving(r2, r1, r2, p1, p2, coefPua, points) << std::endl;
//	std::cout << "Calculating for internal radius:\n";
//	std::cout << "	Radial tension is : " << rad_tension(r1, r1, r2, p1, p2) << std::endl;
//	std::cout << "	Tangential tension is : " << tan_tension(r1, r1, r2, p1, p2) << std::endl;
//	std::cout << "	Moving is : " << moving(r1, r1, r2, p1, p2, coefPua, points) << std::endl;
	out << "====================================================================================================================================================" << std::endl;
	out << "Radius			Radial tension			Tangential tension			Moving			Constant			Equivalent tension" << std::endl;
	double radius = r1;
	double equi = 0;
	for (int i = 0; i < points; i++){
		if (radius!=0.020 && radius!=0.015){
			out << radius <<"		" << rad_tension(radius, r1, r2, p1, p2) << "			"<<tan_tension(radius, r1, r2, p1, p2)<<"				"<<moving(radius, r1, r2, p1, p2, coefPua, elasticity)<<"		"<<rad_tension(radius, r1, r2, p1, p2)+tan_tension(radius, r1, r2, p1, p2) << "			" << equivalent_tension (rad_tension(radius, r1, r2, p1, p2), tan_tension(radius, r1, r2, p1, p2)) <<std::endl;
		}
		else if(radius==0.020){
			out << radius <<"			" << rad_tension(radius, r1, r2, p1, p2) << "				"<<tan_tension(radius, r1, r2, p1, p2)<<"				"<<moving(radius, r1, r2, p1, p2, coefPua, elasticity)<<"		"<<rad_tension(radius, r1, r2, p1, p2)+tan_tension(radius, r1, r2, p1, p2) <<"			" << equivalent_tension (rad_tension(radius, r1, r2, p1, p2), tan_tension(radius, r1, r2, p1, p2)) <<std::endl;
		}
		else if(radius==0.015){
			out << radius <<"			" << rad_tension(radius, r1, r2, p1, p2) << "			"<<tan_tension(radius, r1, r2, p1, p2)<<"				"<<moving(radius, r1, r2, p1, p2, coefPua, elasticity)<<"		"<<rad_tension(radius, r1, r2, p1, p2)+tan_tension(radius, r1, r2, p1, p2) << "			" << equivalent_tension (rad_tension(radius, r1, r2, p1, p2), tan_tension(radius, r1, r2, p1, p2)) <<std::endl;
		}
		if (equivalent_tension (rad_tension(radius, r1, r2, p1, p2), tan_tension(radius, r1, r2, p1, p2)) > equi){
			equi = equivalent_tension (rad_tension(radius, r1, r2, p1, p2), tan_tension(radius, r1, r2, p1, p2));
		}
		radius+=(0.005/9);
	}
	out << std::endl << std::endl << "Max equivalent tension using 4th theory of strength is " << equi;
	out.close();
}
