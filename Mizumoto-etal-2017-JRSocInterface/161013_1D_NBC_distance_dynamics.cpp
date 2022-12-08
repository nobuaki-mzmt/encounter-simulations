// 1D_sim_pair_annihilation
// version for individual distances
// 161002_Mizumoto update

// libraries
#include "stdafx.h"
#include <iostream>	// cout; cin
#include <sstream>	// to_string
#include <fstream>	// writing
using namespace std;

// for random sampling
#include <random>
std::uniform_real_distribution<double> distribution(0.0, 1.0) ;
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
auto rnd = [&]{ return distribution(engine) ; } ;

//// parameters ////
#define n1 1	// num of male
#define n2 1	// num of female
int n = (n1 + n2);
#define steps 1001	// num of steps

// individual parameters
double sM = 1.0;
double sF = 1.0;

double Mloc, Mangle;
double Floc, Fangle;

// Levy
double Mnokori, Fnokori;

// result
int sepdis;
std::string FileName;

// function
// power low
double U, OBJ;
double rbeki(double Myu){
	U = 1 - rnd();
	OBJ = 1 * pow(U, 1/(1-Myu));
	return(OBJ);
}

int main(){
	int Malemyu, Femalemyu;
	double dis;
	int OK, iter;
	char bs=8;

	// input data
	cout << "Enter the Distance" << endl;
	cin >> dis;

	cout << "Enter the myu values" << endl;
	cout << "0:1.1, 9:2.0, 19:3.0" << endl;
	cout << "male" << endl;
	cin >> Malemyu;
	cout << "female" << endl;
	cin >> Femalemyu;

	cout << "Enter the number of replications" << endl;
	cin >> iter;

	cout << "Input data : " << dis << ", " << endl;

	cout << "Are you OK?" << endl;
	cout << "Yes: 1, No: 2" << endl;
	cin >> OK;

	if(OK == 2){ return(2); }

	// output data
	int h,i,j,p,time,cum;
	double R;
	vector<int> res(steps);
	FileName = "sumres_" ;	
	FileName += "1D_NBC_ddcalic_n";
	FileName += std::to_string(n) ;
	FileName += "_M";
	FileName += std::to_string(Malemyu) ;
	FileName += "_F";
	FileName += std::to_string(Femalemyu) ;
	FileName += "_rep";
	FileName += std::to_string(iter) ;
	FileName += "_d";
	FileName += std::to_string(int(dis)) ;
	FileName += ".csv";
	std::ofstream ofs( FileName );

	ofs << "male, female, rep, t5, t10, t20, t50, t100, t200, t500, t1000";
	ofs << endl;

	// simulations
	for(j=0; j<iter; j++){
		ofs << 1+(Malemyu+1)*0.1 << "," << 1+(Femalemyu+1)*0.1 << "," << j;
		sepdis = dis;
		cum = 1;

		// initialization
		Mloc = 0;
		Floc = dis;

		R = rnd() - 0.5;
		Mangle = R/abs(R);
		Mnokori = rbeki(1+(Malemyu+1)*0.1);

		R = rnd() - 0.5;
		Fangle = R/abs(R);
		Fnokori = rbeki(1+(Femalemyu+1)*0.1);

		// movement
		for(time = 1; time < steps+1; time ++){
			Mloc += Mangle;
			Mnokori -= 1;
			if(Mnokori <= 0){
				Mloc += Mangle * Mnokori;
				R = rnd() - 0.5;
				Mangle = R/abs(R);
				Mloc -= Mangle * Mnokori;
				Mnokori += rbeki(1+(Malemyu+1)*0.1);
			}

			Floc += Fangle;
			Fnokori -= 1;
			if(Fnokori <= 0){
				Floc += Fangle * Fnokori;
				R = rnd() - 0.5;
				Fangle = R/abs(R);
				Floc -= Fangle * Fnokori;
				Fnokori += rbeki(1+(Femalemyu+1)*0.1);
			}

			sepdis = Floc - Mloc;
			if( Mloc >= Floc ){
				sepdis = 0;
				break;
			}

			switch(time){
				case 5:
				case 10:
				case 20:
				case 50:
				case 100:
				case 200:
				case 500:
				case 1000:
					ofs << "," << sepdis;
					cum ++;
					break;
			}
		}
		if(cum < 9){
			for(h=cum; h<9; h++){
				ofs << "," << sepdis;
			}
		}
		ofs << endl;
		printf("male: %.1f, female: %.1f, rep: %d\r", 1.1+Malemyu*0.1, 1.1+Femalemyu*0.1, j);
	}
}

