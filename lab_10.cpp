
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
using namespace std;

double r1;  //mm
double r2;  //mm
double p1;  // MPa pressure number 1
double p2;  // MPa pressure number 2
double elasticity;
double coefPua;
double points;

double rad_tension(double r, double p1, double p2) {                          // function for calculation of radial tension
	double delta = ((pow(r1, 2)*p1) - (pow(r2, 2)*p2)) / (r2*r2 - r1 * r1);
	double delta1 = ((pow(r1, 2) * pow(r2, 2) * (p1 - p2)) / (r2*r2 - r1 * r1));
	double delta2 = (1 / (r*r));
	return delta - (delta1 * delta2);
}

double tan_tension(double r, double p1, double p2) {                          // function for calculation of tangential tension
	double delta = ((pow(r1, 2)*p1) - (pow(r2, 2)*p2)) / (r2*r2 - r1 * r1);
	double delta1 = ((pow(r1, 2) * pow(r2, 2) * (p1 - p2)) / (r2*r2 - r1 * r1));
	double delta2 = (1 / (r*r));
	return delta + (delta1 * delta2);
}

double equivalent_tension(double r, double ten1, double ten2) {          // function for calculation of equivalent tension
	double res = ten2 * ten2 - ten1 * ten2 + ten1 * ten1;
	double result = 1 / (pow(2, 0.5)) * pow(res, 0.5);
	return result;
}

double contact_presure(double sigma, double elasticity, double rc, double k1, double k2) {      // function for calculation of contact presure
	return ((sigma*elasticity) / (2 * rc)) * (1 - pow(k1, 2))*(1 - pow(k2, 2)) / ((1 + pow(k1, 2))*(1 - pow(k2, 2)) + (1 + pow(k2, 2))*(1 - pow(k1, 2)));
}


int main() {
	ifstream input;            // initialization of input file stream
	ofstream output;          // the same, but of "output"
	input.open("D:\\C++\\TNN_DIMA\\TenLab\\lab10\\Debug\\data.txt", std::ifstream::in);        // openning file strean for file by this path
	input >> r1 >> r2 >> p1 >> elasticity >> coefPua >> points;                  // inputting of elements
	input.close();                                            // closing of input file stream
	output.open("D:\\C++\\TNN_DIMA\\TenLab\\lab10\\Debug\\result.txt");                // openning output file stream for file by this path
	double delta, k1, k2, Rc, R;                                    // initialization of temperaly variables:  
	R = 0;


	Rc = r1 + (r2 - r1) / 2;                                          //  middle value of radius
//	int Ro = (r2 - r1) / 100;
	k1 = r1 / Rc;                                            //  k1 - coef for calculation contact presure
	k2 = Rc / r2;                                            //  k2 - also coef for calculation contact presure   
	double con_pres = contact_presure(coefPua, elasticity, Rc, k1, k2);
	double* SigRvnutr = new double[(int)points];                            //  initialization of arrays for INSIDE
	double* SigTvnutr = new double[(int)points];                          //  <===^
	double* SigEkvvnutr = new double[(int)points];                          //  <=======^

	double max1 = 0;                                          // variable for first MAX
//s	double R;
	std::cout << "  R:  presure_rad_inside:  presure_tan_inside:  presure_equiv_inside:  max1:" << endl;
	bool parameter = true;
	SigRvnutr[0] = rad_tension(r1, 0, con_pres);                                  // saving result of calculation in each array for INSIDE side cylinder
	SigTvnutr[0] = tan_tension(r1, 0, con_pres);                                  // ^
	SigEkvvnutr[0] = equivalent_tension(r1, SigTvnutr[0], SigRvnutr[0]);
	max1 = SigEkvvnutr[0];
//	cout << "  " << r1 << "    " << SigRvnutr[0] << "    " << SigTvnutr[0] << "    " << SigEkvvnutr[0] << "    " << max1 << endl;
//	output << "  " << r1 << "    " << SigRvnutr[0] << "    " << SigTvnutr[0] << "    " << SigEkvvnutr[0] << "    " << max1 << endl;
	int i = 1;
	double con_pres_rvnutr = con_pres;
	while (parameter)
	{
		if (R == (Rc - (Rc - r1) / (points))) break;
		con_pres_rvnutr = contact_presure((coefPua + i * 0.005), elasticity, Rc, k1, k2);
		parameter = false;
		R = (r1 + i * (Rc - r1) / (points));                              // changing value of "coming" radius
		SigRvnutr[i] = rad_tension(R, 0, con_pres_rvnutr);                                  // saving result of calculation in each array for INSIDE side cylinder
		SigTvnutr[i] = tan_tension(R, 0, con_pres_rvnutr);                                  // ^
		SigEkvvnutr[i] = equivalent_tension(R, SigTvnutr[i], SigRvnutr[i]);                          // |

		if (SigEkvvnutr[i] < max1)                      // condition iterator for reinitialization of MAX value for equivalent_tension INSIDE
		{
			parameter = true;
			max1 = SigEkvvnutr[i];                      // reinitialization
		}

		std::cout << "  " << R << "  " << SigRvnutr[i] << "    " << SigTvnutr[i] << "    " << SigEkvvnutr[i] << "    " << max1 << endl;
		output << "  " << R << "  " << SigRvnutr[i] << "    " << SigTvnutr[i] << "    " << SigEkvvnutr[i] << "    " << max1 << endl;
		i += 1;
	}
	//  cout << "Maximum value for contact presure from inside is: " << max1 << endl;
	//  output << "max1:          " << max1 << endl;                            // outputting value of MAX INSIDE

	double* SigRzovn = new double[(int)points];                                //  initialization of arrays for OUTSIDE
	double* SigTzovn = new double[(int)points];                                // <======^
	double* SigEkvzovn = new double[(int)points+1];                              // <==========^

	double max2 = 0;                                              // initialization variable for outside equivalent_tension MAX value

	std::cout << "\n\n  R:  presure_rad_outside:  presure_tan_outside:  presure_equiv_outside:  max2:" << endl;
	parameter = true;
	SigRzovn[0] = rad_tension(r1, contact_presure(coefPua, elasticity, Rc, k1, k2), 0);                                      // saving result of calculation in each array for OUTSIDE side cylinder
	SigTzovn[0] = tan_tension(r1, contact_presure(coefPua, elasticity, Rc, k1, k2), 0);                                      // ^
	SigEkvzovn[0] = equivalent_tension(R, SigRzovn[0], SigTzovn[0]);
	max2 = SigEkvzovn[0];
//	std::cout << "  " << r1 << "    " << SigRzovn[0] << "    " << SigTzovn[0] << "    " << SigEkvzovn[0] << "    " << max2 << endl;
//	output << "  " << r1 << "    " << SigRzovn[0] << "    " << SigTzovn[0] << "    " << SigEkvzovn[0] << "    " << max2 << endl;
	i = 1;
//	R = 0;
	double con_pres_rvzovn = con_pres;
	while (parameter)
	{
		if (R >= (r2 - (r2 - Rc) / (points))) break;
		R = (Rc + i * (r2 - Rc) / (points - 1));                                    // changing value of "coming" radius
	//	std::cout << std::endl << R << std::endl;
		con_pres_rvzovn = contact_presure((coefPua + i * 0.005), elasticity, Rc, k1, k2);
		SigRzovn[i] = rad_tension(R, con_pres_rvzovn, 0);                                      // saving result of calculation in each array for OUTSIDE side cylinder
		SigTzovn[i] = tan_tension(R, con_pres_rvzovn, 0);                                      // ^
		SigEkvzovn[i] = equivalent_tension(R, SigRzovn[i], SigTzovn[i]);                          // |


			if (SigEkvzovn[i] < max2)                    // condition iterator for reinitialization of MAX value for equivalent_tension OUTSIDE
			{
				max2 = SigEkvzovn[i];                    // reinitialization
			}

		std::cout << "  " << R << "    " << SigRzovn[i] << "    " << SigTzovn[i] << "    " << SigEkvzovn[i] << "    " << max2 << endl;
		output << "  " << R << "    " << SigRzovn[i] << "    " << SigTzovn[i] << "    " << SigEkvzovn[i] << "    " << max2 << endl;
		i += 1;
	}

	//  cout << "Maximum value for contact presure from outside is: " << max2 << endl;
	//  output << "max2:          " << max2 << endl;                            // outputting value of MAX OUTSIDE

	double* SigR_without_counter = new double[(int)points];                                //  initialization of arrays for cylinder without using of counter pression
	double* SigT_withour_counter = new double[(int)points];                                //  <======^
	double* SigEkv_withour_counter = new double[(int)points];                              //  <==========^

	std::cout << endl << "\n\n        " << "R:" << "        " << "tension_rad_without_cp:" << "        " << "tension_tan_without_cp:" << "        " << "tension_equiv_without_cp:" << endl;
	for (int i = 0; i < (int)points; i++) {
		R = r1 + i * (r2 - r1) / (points - 1);
		SigR_without_counter[i] = rad_tension(R, p1, p2);
		SigT_withour_counter[i] = tan_tension(R, p1, p2);
		SigEkv_withour_counter[i] = equivalent_tension(R, SigT_withour_counter[i], SigR_without_counter[i]);

		std::cout << "        " << R << "          " << SigR_without_counter[i] << "          " << SigT_withour_counter[i] << "            " << SigEkv_withour_counter[i] << endl;
		output << "        " << R << "          " << SigR_without_counter[i] << "          " << SigT_withour_counter[i] << "            " << SigEkv_withour_counter[i] << endl;

	}

	double* SigR_sum = new double[(int)points];                                //  initialization of arrays for cylinder sum of with using of counter pression + without
	double* SigT_sum = new double[(int)points];                                //  <======^
	double* SigEkv_sum = new double[(int)points];

	std::cout << endl << "        " << "R:" << "        " << "tension_rad_sum:" << "        " << "tension_tan_sum:" << "        " << "tension_equiv_sum:" << endl;
	for (int i = 0; i < (int)points / 2; i++)
	{
		R = (r1 + i * (r2 - r1) / (int)(points));                              // changing value of "coming" radius
		if (SigRvnutr[i] > SigRzovn[i])    SigR_sum[i] = SigR_without_counter[i] + SigRvnutr[i];
		else SigR_sum[i] = SigR_without_counter[i] + SigRzovn[i];
		if (SigTvnutr[i] > SigTzovn[i])    SigT_sum[i] = SigT_withour_counter[i] + SigTvnutr[i];
		else SigR_sum[i] = SigT_withour_counter[i] + SigTzovn[i];
		if (SigEkvvnutr[i] > SigEkvzovn[i])    SigEkv_sum[i] = SigEkv_withour_counter[i] + SigEkvvnutr[i];
		else SigR_sum[i] = SigEkv_withour_counter[i] + SigEkvzovn[i];
		// |
		std::cout << "        " << R << "          " << SigR_sum[i] << "          " << SigT_sum[i] << "            " << SigEkv_sum[i] << endl;
		output << "        " << R << "          " << SigR_sum[i] << "          " << SigT_sum[i] << "            " << SigEkv_sum[i] << endl;


	}


	output.close();                                                // closing of output file stream
	std::system("PAUSE");
	return 200;
}