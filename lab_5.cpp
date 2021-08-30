#include <iostream>
#include <math.h>
#include <fstream>
#define L (p1 - p2 + ((3+coefPua)/8)*specmassmat*4*10*frequency*frequency*(r1*r1 - r2*r2) - deltaT*alpha*elasticity/3)/(1/(r1*r1) - 1/(r2*r2)) 	// функция определения L
#define K p1 - L/(r1*r1) + ((3+coefPua)/8)*specmassmat*4*10*frequency*frequency*(r1*r1) + (deltaT*alpha*elasticity*r1)/(3*(r2-r1))				// функция определения K
// здесь константы K и L, вычисляемые непосредственно используя входные параметры технического задания, инициализированы как макросы ...
// ... и просчитаны в них соответственно из свойств данной синтаксической структуры (которая предусмотрена языком разработки)
using namespace std;

double points;	                                     // переменная для инициализации точек
const double r1 = 0.075;                             // радиус R1
const double r2 = 0.350;                             // радиус R2
const double p1 = 0;	                             // МПа напряжённость на r1
const double p2 = 6 * pow(10, 6);                    // МПа напряжённость на r2
const double elasticity = 0.115 * pow(10,6);         // коэффициент эластичности
const double deltaT = 175;                           // T* - дельта температур T1 - T2	
const double alpha = 18 * pow(10, -6);               // температурный коэффициент 
const double coefPua = 0.35;                         // коэффициент Пуассона
//const double points = 11;	                         // число точек
const double frequency = 50;                         // частота об/сек
const double specmassmat = 8700;                     // удельная масса материала кг/м3


double getL (){	                                     // функция-поддержки для проверки правильности просчёта коэффициента L 
	return (p1 - p2 + ((3+coefPua)/8)*specmassmat*4*10*frequency*frequency*(r1*r1 - r2*r2) - deltaT*alpha*elasticity/3)/(1/(r1*r1) - 1/(r2*r2));
}
double getK (){	                                     // функция-поддержки для проверки правильности просчёта коэффициента K
	return p1 - L/(r1*r1) + ((3+coefPua)/8)*specmassmat*4*10*frequency*frequency*(r1*r1) + (deltaT*alpha*elasticity*r1)/(3*(r2-r1));
}

double rad_tension (double r){	                     // функция для вычисления радиальных напряжённостей
	return K + L/(r*r) - ((3+coefPua)/8)*specmassmat*4*10*frequency*frequency*(r*r) - (deltaT*alpha*elasticity*r)/(3*(r2-r1));
}

double tan_tension (double r){	                     // функция для вычисления тангенциальных напряжённостей
	return K - L/(r*r) - ((1+3*coefPua)/8)*specmassmat*4*10*frequency*frequency*(r*r) - (2*deltaT*alpha*elasticity*r)/(3*(r2-r1));
}

double equivalent_tension (double ten1, double ten2){      // функция для вычисления эквивалентных напряжённостей за 4ой теорией прочности
	double res = ten2 * ten2 - ten1 * ten2 + ten1 * ten1;
	double result = 1/(pow(2, 0.5)) * pow(res, 0.5);
	return result;
}


int main(){
	std::ofstream out("output.txt");                 // инициализация открытия файлового потока
	cout << " Please, type number of your points: ";
	cin >> points;                                   // ввод числа точек
	double rad[(int)points], radt[(int)points], tant[(int)points], equiv[(int)points];	// инициализация массивов для хранения данных
	std::cout << "	Calculating of K constant:	(using macros)" << "		" << K << std::endl;;//	вывод в консоль эл. К используя макрос
	std::cout << "	Calculating of K constant:	(using function)" << "	" << getK() << std::endl;//	вывод ... эл. К используя функцию
	std::cout << "	Calculating of L constant:	(using macros)" << "		" << L << std::endl;//	вывод ... эл. L используя макрос
	std::cout << "	Calculating of L constant:	(using function)" << "	" << getL() << std::endl;//	вывод ... эл. L используя функцию
	std::cout << "==============================================================================================================" << std::endl;
	std::cout << "Radius		Radial tension		Tangential tension		Equivalent tension" << std::endl;
	double radius = r1;	                             // инициализация пеменной для работы с радиусом
	for (int i = 0; i < points; i++){                // в дальнейшем данная переменная будет служить для сохранения изменений радиуса внутри цикла | цикл для перехода между точками
		rad[i] = radius;                             // инициализация i-го элемента массива rad[]
		radt[i] = rad_tension(radius);	             // инициализация i-го элемента массива radt[]
		tant[i] = tan_tension(radius);	             // инициализация i-го элемента массива tant[]
		equiv[i] = equivalent_tension(radt[i], tant[i]);	// инициализация i-го элемента массива equiv[]
		cout << rad[i] << "	" << radt[i] << "	" << tant[i] << "	" << equiv[i] << endl;	   // вывод в консоль соответствующих элементов массива, что соответствуют данной итерации массива  + радиусу
		out << rad[i] << " " << radt[i] << " " << tant[i] << " " << equiv[i] << endl;         // ... в файл ...
		radius+=((r2-r1)/(points-1));	             // переинциализация переменной radius на значение, большее за предыдущее на deltaR
	}
	out.close();                                     // закрытие файлового потока
}
