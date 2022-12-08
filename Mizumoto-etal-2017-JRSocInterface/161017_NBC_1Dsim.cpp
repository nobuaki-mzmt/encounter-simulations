// 1D_sim_pair_annihilation
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

// individual parameters
double sM = 1.0;
double sF = 1.0;

double Mloc, Mangle;
double Floc, Fangle;

// Levy
double Mnokori, Fnokori;
double Malemyu, Femalemyu;
#define MinMyu 1.1
#define MaxMyu 3
#define MyuRange 0.1

// result
int endtime;
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

	double dis;
	int iter, steps, OK;
	char bs=8;

	// input data
	cout << "Enter the Distance" << endl;
	cin >> dis;
	
	cout << "Enter the num of rep" << endl;
	cin >> iter;

	cout << "Enter the num of steps" << endl;
	cin >> steps;


	cout << "Input data : " << dis << ", " << iter << endl;

	cout << "Are you OK?" << endl;
	cout << "Yes: 1, No: 2" << endl;
	cin >> OK;

	if(OK == 2){ return(2); }

	// output data
	int h,i,j,p,time;
	double R;
	vector<int> res(steps);
	FileName = "sumres_" ;	
	FileName += "SA_Levy";
	FileName += "_NBC";
	FileName += "_n";
	FileName += std::to_string(n) ;
	FileName += "_rep";
	FileName += std::to_string(iter) ;
	FileName += "_d";
	FileName += std::to_string(int(dis)) ;
	FileName += ".csv";
	std::ofstream ofs( FileName );

	ofs << "male, female";
	for(i = 0; i<steps; i++){
		ofs << "," << i+1;
	}
	ofs << endl;
	
	// simulations
	for(Malemyu = MinMyu; Malemyu < (MaxMyu+0.1); Malemyu+=MyuRange){
		for(Femalemyu = MinMyu; Femalemyu < (MaxMyu+0.1); Femalemyu+=MyuRange){
			if(Malemyu >= Femalemyu){
				fill(res.begin(), res.end(), 0);
				for(j=0; j<iter; j++){
					endtime = steps+1;
					// initialization
					Mloc = 0;
					Floc = dis;

					R = rnd() - 0.5;
					Mangle = R/abs(R);
					Mnokori = rbeki(Malemyu);

					R = rnd() - 0.5;
					Fangle = R/abs(R);
					Fnokori = rbeki(Femalemyu);

					// movement
					for(time = 0; time < steps; time ++){
						Mloc += Mangle;
						Mnokori -= 1;
						if(Mnokori <= 0){
							Mloc += Mangle * Mnokori;
							R = rnd() - 0.5;
							Mangle = R/abs(R);
							Mloc -= Mangle * Mnokori;
							Mnokori += rbeki(Malemyu);
						}

						Floc += Fangle;
						Fnokori -= 1;
						if(Fnokori <= 0){
							Floc += Fangle * Fnokori;
							R = rnd() - 0.5;
							Fangle = R/abs(R);
							Floc -= Fangle * Fnokori;
							Fnokori += rbeki(Femalemyu);
						}

						if( Floc - Mloc < 1 ){
							endtime = time+1;
							break;
						}else if( Floc - Mloc > (steps-time)/2 ){
							break;
						}
					}

					for(h=endtime; h<steps+1; h++){
						res[h-1] ++; 
					}
					printf("male: %.1f, female: %.1f, rep: %d\r", Malemyu, Femalemyu, j);
				}
				// output data
				std::ofstream ofs( FileName, std::ios::out | std::ios::app );
				ofs << Malemyu << "," << Femalemyu;
				for(p = 0; p < steps; p++){
					ofs<< "," << res[p];
				}
				ofs << endl;
			}
		}
	}
}

