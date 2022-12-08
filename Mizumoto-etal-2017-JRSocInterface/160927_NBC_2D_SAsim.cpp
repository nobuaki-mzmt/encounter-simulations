// Sim. of Evolution of sex different movement patterns
// Mizumoto Nobuaki 160914

// libraries
#include "stdafx.h"
#include <fstream>	
#include <regex>	// std::regex_replace
#include <iostream>	// cout; cin
#include <numeric>	// iota
#include <sstream>
using namespace std;

// for random sampling
#include <random>
std::uniform_real_distribution<double> distribution(0.0, 1.0) ;
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
auto rnd = [&]{ return distribution(engine) ; } ;

#define pi 3.1415926535
#define nipi 3.1415926535*2

//// parameters ////
#define moment 0.1	// increment of a step
#define dis 10

#define steps 1501	// num of steps
#define nsteps 15010	// steps/moment

// individual parameters
double sM = 1.0;
double sF = 1.0;

// result
std::string FileName;
double sumMmyu, sumFmyu;
int endtime;

// power law
double U, OBJ;
double rbeki(double Myu){
	U = 1 - rnd();
	OBJ = 1 * pow(U, 1/(1-Myu));
	return(OBJ);
}

vector<string> split(string& input, char delimiter){
    istringstream stream(input);
    string field;
    vector<string> result;
    while (getline(stream, field, delimiter)) {
        result.push_back(field);
    }
    return result;
}

// simulations
int main(){

	int simtype, iter, n, n1, n2, prerep;
	char bs=8;
	double mindis, karidis;

	// setting
	int OK;

	cout << "Simulations for Sensitive analysis" << endl;
	cout << "Select the Population size" << endl;
	cout << "100, 144, 196, 256, 324, 400, 484, 576, 676, 784, 900, 1024, 1156, 1296, 1444, 1600, 1764, 1936, 2116, 2304, 2500" << endl;

	cin >> n;
	n1 = n2 = n/2;
	int nroot = sqrt(n);

	cout << "create new: 1, update: 2" << endl;
	cin >> simtype;

	cout << "Enter the number of replications" << endl;
	cin >> iter;

	cout << "input : " << n << " " << simtype << " " << iter << endl;

	cout << "Are you OK?" << endl;  
	cout << "Yes: 1, No: 2" << endl; 
	cin >> OK;


	ifstream ifs1(FileName);
	ifstream ifs2(FileName);
	string line;
	int countraw1 = 0;
	int countraw2 = 0;
	while (getline(ifs1, line)) {
		countraw1 += 1;
	}
	while (getline(ifs2, line)) {
		countraw2 += 1;
		if(countraw1 == countraw2){
			vector<string> strvec = split(line, ',');
			prerep = stoi(strvec.at(3))+1;
		}
	}

	if(OK == 2){ return(2); }

	vector<double> MlocX(n1), FlocX(n2);
	vector<double> MlocY(n1), FlocY(n2);
	vector<double> Mangle(n1), Fangle(n2);
	vector<int> Mstatus(n1), Fstatus(n2);

	// Levy
	int setMmyu, setFmyu;
	vector<int> Mmyu(n1), Fmyu(n2);
	vector<double> Mnokori(n1), Fnokori(n2);
	vector<int> id(n);

	int pairs;

	vector<double> Mmindis(n1);

	int i,p,q,z;

	FileName = "res_SA_Levy_n";
	FileName += std::to_string(n) ;
	//FileName += "_d";
	//FileName += std::to_string(dis) ;
	FileName = std::regex_replace(FileName, std::regex("\\."), "");
	FileName += ".csv";

	if(simtype == 1){
		std::ofstream ofs( FileName );
		ofs << "male, female, rep, t20, t50, t100, t200, t500, t1000, t1500";
		ofs << std::endl;
	}
	std::ofstream ofs( FileName, std::ios::out | std::ios::app );
	for(z = 0; z < iter; z++){
		for(setMmyu = 0; setMmyu < 20; setMmyu++){
			for(setFmyu = 0; setFmyu < 20; setFmyu++){
				if(setFmyu <= setMmyu){ 

					ofs << 1+(setMmyu+1)*0.1 << "," << 1+(setFmyu+1)*0.1 << "," << prerep+z;

					iota(id.begin(),id.end(), 0);
					shuffle(id.begin(), id.end(), mt19937());

					// new world
					fill(Mmyu.begin(), Mmyu.end(), setMmyu);
					fill(Fmyu.begin(), Fmyu.end(), setFmyu);
					fill(Mstatus.begin(), Mstatus.end(), 1);
					fill(Fstatus.begin(), Fstatus.end(), 1);
					fill(Mmindis.begin(), Mmindis.end(), 0);

					for(i = 0; i < n1; i++){
						MlocX[i] =  dis * (id[i] % nroot);
						MlocY[i] =  dis * (id[i] / nroot);
						FlocX[i] =  dis * (id[i+n1] % nroot);
						FlocY[i] =  dis * (id[i+n1] / nroot);				

						Mangle[i] = rnd();
						Mangle[i] *= nipi;
						Mnokori[i] = rbeki(1.1 + 0.1*Mmyu.at(i));

						Fangle[i] = rnd();
						Fangle[i] *= nipi;
						Fnokori[i] = rbeki(1.1 + 0.1*Fmyu.at(i));
					}

					pairs = 0;
					for(int time = 0; time < nsteps; time ++){
						for(i = 0; i < n1; i ++){
							if(Mstatus[i]>0){
								MlocX[i] += cos(Mangle[i]) * moment;
								MlocY[i] += sin(Mangle[i]) * moment;
								Mnokori[i] -= moment;
								if(Mnokori[i] <= 0){
									MlocX[i] += cos(Mangle[i]) * Mnokori[i];
									MlocY[i] += sin(Mangle[i]) * Mnokori[i];
									Mangle[i] = rnd();
									Mangle[i] *= nipi;
									MlocX[i] -= cos(Mangle[i]) * Mnokori[i];
									MlocY[i] -= sin(Mangle[i]) * Mnokori[i];
									Mnokori[i] += rbeki(1.1 + 0.1*Mmyu.at(i));
								}
								Mmindis[i] -= 2*moment;
							}
						}
						for(i = 0; i < n1; i ++){
							if(Fstatus[i]>0){
								FlocX[i] += cos(Fangle[i]) * moment;
								FlocY[i] += sin(Fangle[i]) * moment;
								Fnokori[i] -= moment;
								if(Fnokori[i] <= 0){
									FlocX[i] += cos(Fangle[i]) * Fnokori[i];
									FlocY[i] += sin(Fangle[i]) * Fnokori[i];
									Fangle[i] = rnd();
									Fangle[i] *= nipi;
									FlocX[i] -= cos(Fangle[i]) * Fnokori[i];
									FlocY[i] -= sin(Fangle[i]) * Fnokori[i];
									Fnokori[i] += rbeki(1.1 + 0.1*Fmyu.at(i));
								}
							}
						}

						for(p=0; p<n1; p++){
							if(Mstatus[p]>0){
								if(Mmindis[p] < 1){
									mindis = 99999;
									for(q=0; q<n2; q++){
										if(Fstatus[q]>0){
											karidis = sqrt(((MlocX[p]-FlocX[q])*(MlocX[p]-FlocX[q])) + ((MlocY[p]-FlocY[q])*(MlocY[p]-FlocY[q])));
											if(karidis < mindis){
												mindis = karidis;
												if( mindis < 1 ){
													Mstatus[p] = 0;
													Fstatus[q] = 0;
													pairs += 1;
													break;
												}
											}
										}										
									}
									Mmindis[p] = mindis;
									if(mindis > steps-time/10){
										Mstatus[p] = 0;
									}
								}
							}
						}
						switch(time){
							case 20:
							case 50:
							case 100:
							case 200:
							case 500:
							case 1000:
							case 1500:
								ofs << "," << pairs;
								break;
						}
					}
					ofs << endl;
					printf("male:%.1f, female:%.1f, rep:%5d\r", 1.1+setMmyu*0.1, 1.1+setFmyu*0.1,prerep+z);
				}
			}

		}
	}
}