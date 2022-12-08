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
double selection;	// strength of selection (num. of pairs)
#define mutation 0.001	// mutation rate
#define mutationrange 0.1	// mutation range for stepping stone
#define iter 1	// repeat num
#define dis 10

#define steps 1500	// num of steps
#define nsteps 15000	// steps/moment

// individual parameters
double sM = 1.0;
double sF = 1.0;

// result
std::string FileName;
int endtime;
double sumMmyu, sumFmyu;

// functions
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

	int simtype, mink, maxk, n, n1, n2, generations, pregene=0;
	char bs=8;
	double mindis, karidis;

	// settiing
	int OK;

	cout << "Select the Population size" << endl;
    cout << "100, 144, 196, 256, 324, 400, 484, 576, 676, 784, 900, 1024, 1156, 1296, 1444, 1600, 1764, 1936, 2116, 2304, 2500" << endl;

    cin >> n;
	n1 = n2 = n/2;
	int nroot = sqrt(n);

	cout << "Selection rate" << endl;
	cout << "0¨0.05, 5¨0.30, 10¨0.55, 15¨0.8" << endl;

    cout << "Enter the min selection rate" << endl;
	cin >> mink;

    cout << "Enter the max selection rate" << endl;
	cin >> maxk;

	cout << "Enter the num of generations" << endl;
	cin >> generations;
	
	cout << "create new: 1, update: 2" << endl;
	cin >> simtype;

	cout << "input : " << n << " " << mink << " " << maxk << " " << simtype << endl;

	cout << "Are you OK?" << endl;  
	cout << "Yes: 1, No: 2" << endl; 
	cin >> OK;

	if(OK == 2){ return(2); }

	vector<double> MlocX(n1);
	vector<double> FlocX(n2);
	vector<double> MlocY(n1);
	vector<double> FlocY(n2);
	vector<double> Mangle(n1);
	vector<double> Fangle(n2);
	vector<int> Mstatus(n1);	// encounter = 0, still = 1
	vector<int> Fstatus(n2);

	// Levy
	int initialmyu = 19;	// 1.1 + 19*0.1 = 3
	vector<int> Mmyu(n1);
	vector<int> Fmyu(n2);
	vector<double> Mnokori(n1);
	vector<double> Fnokori(n2);
	vector<int> id(n);

	// selection 
	vector<double> Mreproduce(n1);
	vector<double> Freproduce(n2);
	int pairs;

	vector<double> Mmindis(n1);

	int i,j,k,p,q,z;
	double R;
	int ofspring;
	double xdis, ydis;
	
	for(k = mink; k <= maxk; k++){
		selection = n1 * 0.05 * (k+1);

		cout << "n: " << n << ", selection rate: " <<  selection/n1 << endl;

		for(z = 0; z < iter; z++){

			FileName = "res_EV_Levy_n";
			FileName += std::to_string(n) ;
			FileName += "_s";
			FileName += std::to_string(selection/n1).substr(0, 4);
			//FileName += "_mr";
			//FileName += std::to_string(mutationrange).substr(0, 3) ;
			FileName = std::regex_replace(FileName, std::regex("\\."), "");
			FileName += ".csv";
			
			if(simtype == 1){
				std::ofstream ofs( FileName );
				ofs << "generation, male, female, time, pairs";
				for(p=0; p<20; p++){
					ofs << ",M" << 1.1+0.1*p;
				}
				for(p=0; p<20; p++){
					ofs << ",F" << 1.1+0.1*p;
				}
				ofs << std::endl;
				fill(Mmyu.begin(), Mmyu.end(), initialmyu);
				fill(Fmyu.begin(), Fmyu.end(), initialmyu);

			} else if (simtype == 2){
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
						p = 0;
						q = 0;
						vector<string> strvec = split(line, ',');
						pregene = stoi(strvec.at(0))+1;
						for (i=5;i<strvec.size();i++){
							if(stoi(strvec.at(i)) > 0){
								if(i < 25){
									for(j = 0; j < (stoi(strvec.at(i))); j++){
										Mmyu[j+p] = (i-5);
									}
									p += stoi(strvec.at(i));
								} else {
									for(j = 0; j < (stoi(strvec.at(i))); j++){
										Fmyu[j+q] = (i-25);
									}
									q += stoi(strvec.at(i));
								}
							}
				       }
					}
				}

			}

			iota(id.begin(),id.end(), 0);
			shuffle(id.begin(), id.end(), mt19937());

			// new world
			fill(Mstatus.begin(), Mstatus.end(), 1);
			fill(Fstatus.begin(), Fstatus.end(), 1);

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
			// new generation
			for(j = 0; j < generations; j++){
				/*time_t timer;
				struct tm *t_st;
				time(&timer);
				t_st = localtime(&timer);
				cout << put_time(t_st, "%c") << endl;*/
				pairs = 0;
				endtime = nsteps;
				for(int time = 0; time < nsteps; time ++){
					if(pairs < selection){
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
											//if( pow((pow((MlocX[p]-FlocX[q]),2) + pow((MlocY[p]-FlocY[q]),2)),0.5) < 1){
											karidis = sqrt(((MlocX[p]-FlocX[q])*(MlocX[p]-FlocX[q])) + ((MlocY[p]-FlocY[q])*(MlocY[p]-FlocY[q])));
											if(karidis < mindis){
												mindis = karidis;
												if( mindis < 1 ){
													Mstatus[p] = 0;
													Fstatus[q] = 0;
													Mreproduce[pairs] = Mmyu[p];
													Freproduce[pairs] = Fmyu[q];
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
					} else {
						endtime = time;
						break;
					}
				}

				// next generation
				sumFmyu = 0;
				sumMmyu = 0;
				shuffle(id.begin(), id.end(), mt19937());
				fill(Mmindis.begin(), Mmindis.end(), 0);
				for(i = 0; i < n1; i++){
					MlocX[i] =  dis * (id[i] % nroot);
					MlocY[i] =  dis * (id[i] / nroot);
					FlocX[i] =  dis * (id[i+n1] % nroot);
					FlocY[i] =  dis * (id[i+n1] / nroot);

					Mangle[i] = rnd () * nipi;
					ofspring = (int)(rnd() * 0.99999 * (pairs));
					Mmyu[i] = Mreproduce[ofspring];
					Mstatus[i] = 1;

					Fangle[i] = rnd () * nipi;
					ofspring = (int)(rnd() * 0.99999 * (pairs));
					Fmyu[i] = Freproduce[ofspring];
					Fstatus[i] = 1;

					// mutation
					R = rnd();
					if( R < mutation){
						if(Mmyu[i] > 18){ Mmyu[i] -= 1;
						}else if(Mmyu[i] < 1){ Mmyu[i] += 1;
						}else if( R < mutation/2){
							Mmyu[i] -= 1;
						}else{
							Mmyu[i] += 1;
						}
					}
					Mnokori[i] = rbeki(1.1 + 0.1*Mmyu.at(i));

					R = rnd();
					if( R < mutation){
						if(Fmyu[i] > 18){ Fmyu[i] -= 1;
						}else if(Fmyu[i] < 1){ Fmyu[i] += 1;
						}else if( R < mutation/2){
							Fmyu[i] -= 1;
						}else{
							Fmyu[i] += 1;
						}
					}
					Fnokori[i] = rbeki(1.1 + 0.1*Fmyu.at(i));

					sumMmyu += (1.1 + 0.1*Mmyu.at(i));
					sumFmyu += (1.1 + 0.1*Fmyu.at(i));
				}
				std::ofstream ofs( FileName, std::ios::out | std::ios::app );
				ofs << pregene+j << "," << sumMmyu/n1 << "," << sumFmyu/n1 << "," << endtime*moment << "," << pairs;
					for(p=0; p<20; p++){
						ofs << "," << count(Mmyu.cbegin(), Mmyu.cend(), p );
					}
					for(p=0; p<20; p++){
						ofs << "," << count(Fmyu.cbegin(), Fmyu.cend(), p );
					}
				ofs << std::endl;

				/*cout << "gene: " << setfill('0') << setw(5) << j << ", pairs:" << setfill('0') << setw(4) << pairs;
				cout << ", time: " << setfill('0') << setw(4) << endtime;
				cout << ", male: " << setprecision(2) << sumMmyu/n1 << ", female: " << fixed << setprecision(2) << sumFmyu/n1 << "\r"; */
			
				printf("generation:%05d, pairs:%04d, time:%05d, male: %.1f, female: %.1f   \r",pregene+j,pairs, endtime, sumMmyu/n1, sumFmyu/n1);
			}
			cout << endl;
			cout << endl;
		}
	}
	cout << "Enter to end the simulation" << endl;
	cin >> OK;
	return(OK);
}