// Code for simulations for "The optimal movement patterns for mating encounters with sexually asymmetric detection ranges"
// 161219 N. Mizumoto created
// 180113 N. Mizumoto commented and checked

// manually imput: L2, p, Er
// output: number of pairs created
// * to change the parameters (n, vreceiver, vsender, µreceiver amd µsender), please change the code

// The simulation program was implemented in Microsoft Visual Studio C++ 2012 (Windows 10).

//// libraries ////
#include "stdafx.h"
#include <iostream>	// cout; cin
#include <sstream>	// to_string
#include <fstream>	// writing
#include <regex>	// std::regex_replace
using namespace std;

// for random sampling
#include <random>
std::uniform_real_distribution<double> distribution(0.0, 1.0) ;
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
auto rnd = [&]{ return distribution(engine) ; } ;
#define pi 3.1415926535
#define nipi 3.1415926535*2

/*// opencv (for plot)
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
using namespace cv;
*/

//// abbreviation ////
// Sender: S
// Receiver: R
// µ: myu

//// setting ////
// exploring ranges myu values
#define R_myu_min 0						 // 0 = 1.1, 19 = 3.0
#define R_myu_max 19
#define S_myu_min 0
#define S_myu_max 19

// parameters for individuals
#define R 1								// R = 1
#define n1 25							// number of receivers
#define n2 25							// number of senders
#define vS 1
#define vR 1

// replication
#define step_scale 7					// number of time steps (10^step_scale)
#define steps pow(10,step_scale)+1
#define moment 10						// each time step is divided into 10 substeps to check encounters

// variables
int n = (n1 + n2);
int L2, Er;
double p;

//// containers ////
vector<double> R_loc_X(n1), S_loc_X(n2);	// location (x-axis)
vector<double> R_loc_Y(n1), S_loc_Y(n2);	// location (y-axis)
vector<double> R_angle(n1), S_angle(n2);	// direction of movement
vector<double> R_move(n1), S_move(n1);		// persistence of stratight motion (for Lévy walk)
vector<int> id(n1+n2);						// ID: name of each individual
vector<double> Dis_to_Nearest(n1);			// Distance from each receiver to its nearest sender
vector<int> Nearest_S(n1);					// ID of the nearest sender
int pairs;
int i, j, h;

//// functions ////
// create random variables from power-law distribution
double U, OBJ;
double rPL(double Myu){ 
	U = 1 - rnd();
	OBJ = 1 * pow(U, 1/(1-Myu));
	return(OBJ);
}

// modification for PBC (periodic boundary conditions)
double aPBC(double Loc, double Length){
	if(Loc > Length){ Loc -= Length; } 
	else if(Loc < 0){ Loc += Length; }
	return(Loc);
}

// determination of moving direction during navigation (considering PBC)
double direction, dx, dy;
double NearSearch(double Rx, double Ry, double Sx, double Sy, double Length){
	dx = Sx-Rx;	dy = Sy-Ry;
	if(dx < -Length/2){
		Sx += Length;
	} else if (dx > Length/2) {
		Sx -= Length; 
	}
	if(dy < -Length/2){
		Sy += Length;
	} else if (dy > Length/2) {
		Sy -= Length; 
	}
	dx = Sx-Rx;	dy = Sy-Ry;
	direction = atan( dy/dx );
	if(dx < 0){	direction += pi; }
	return(direction);
}

//// main simulations ////
double dis, xdis, ydis, prob;
int mom, Time;
int MS(double RMYU, double SMYU, double L, double ER, double p) {

	// initialization
	pairs = 0;
	fill(Dis_to_Nearest.begin(), Dis_to_Nearest.end(), 0);
	for(i = 0; i < n1; i++){	// all individuals are randomly distributed in initial condition
		R_loc_X[i] =  L * rnd();
		R_loc_Y[i] =  L * rnd();
		R_angle[i] = nipi * rnd();
		R_move[i] = rPL(RMYU);
		S_loc_X[i] =  L * rnd();
		S_loc_Y[i] =  L * rnd();				
		S_angle[i] = nipi * rnd();
		S_move[i] = rPL(SMYU);
	}
	// redistribute individuals which encountering at the initial condition
	// & determining the nearest sender for each receiver
	while(*min_element(Dis_to_Nearest.begin(),Dis_to_Nearest.end())<1){	
		for(j=0; j<n1; j++){	// search for encountering pairs
			if(Dis_to_Nearest[j] < ER){	// skip non encountering individuals (also see the part of checking encounters L206-234)
				Dis_to_Nearest[j] = 99999;	// reset
				for(h=0; h<n2; h++){        // search for the nearest sender
					xdis = min(abs(R_loc_X[j]-S_loc_X[h]), abs(L-abs(R_loc_X[j]-S_loc_X[h])));
					ydis = min(abs(R_loc_Y[j]-S_loc_Y[h]), abs(L-abs(R_loc_Y[j]-S_loc_Y[h])));
					dis = sqrt((xdis*xdis) + (ydis*ydis));
					if(dis < Dis_to_Nearest[j]){
						Dis_to_Nearest[j] = dis;
						Nearest_S[j] = h;
						if( Dis_to_Nearest[j] < 1 ){ // redistribution
							R_loc_X[j] =  L * rnd();
							R_loc_Y[j] =  L * rnd();
							S_loc_X[h] =  L * rnd();
							S_loc_Y[h] =  L * rnd();
							Dis_to_Nearest[h] = 0;
							break;
						}
					}
				}
			}
		}
	}

	// begining of simulations
	for(Time=0; Time<steps; Time++){
		// calicuration of angle during navigation
		for(i = 0; i < n1; i ++){
			if(Dis_to_Nearest[i] < ER){
				prob = rnd();
				if(prob <= p){
					R_angle[i] = NearSearch(R_loc_X[i], R_loc_Y[i], S_loc_X[Nearest_S[i]], S_loc_Y[Nearest_S[i]], L);
				} else {
					R_angle[i] = nipi * rnd();
				}
			}
		}
		
		// move and check whether encounter
		for(mom = 0; mom < moment; mom++){	// divide each time step into 10 substeps
			// movements
			for(i = 0; i < n1; i ++){	// for receivers
				R_loc_X[i] += cos(R_angle[i]) * 0.1 * vR;	// move
				R_loc_Y[i] += sin(R_angle[i]) * 0.1 * vR;
				Dis_to_Nearest[i] -= 0.1 * (vR + vS);	// theoretically possible shortest distance to the nearest sender
				R_move[i] -= 0.1 * vR;
				if(R_move[i] <= 0){	// the end of moving steps in Lévy walk
					R_loc_X[i] += cos(R_angle[i]) * R_move[i];	// back excessive movement
					R_loc_Y[i] += sin(R_angle[i]) * R_move[i];
					R_angle[i] = rnd() * nipi;	// redirection
					R_loc_X[i] -= cos(R_angle[i]) * R_move[i];
					R_loc_Y[i] -= sin(R_angle[i]) * R_move[i];	
					R_move[i] += rPL(RMYU);	// move
				}
				R_loc_X[i] = aPBC(R_loc_X[i], L);	// adjustment for PBC
				R_loc_Y[i] = aPBC(R_loc_Y[i], L);
			}

			for(i = 0; i < n2; i ++){	// for senders
				S_loc_X[i] += cos(S_angle[i]) * 0.1 * vS;
				S_loc_Y[i] += sin(S_angle[i]) * 0.1 * vS;
				S_move[i] -= 0.1 * vS;
				if(S_move[i] <= 0){
					S_loc_X[i] += cos(S_angle[i]) * S_move[i];	
					S_loc_Y[i] += sin(S_angle[i]) * S_move[i];
					S_angle[i] = rnd() * nipi;
					S_loc_X[i] -= cos(S_angle[i]) * S_move[i];	
					S_loc_Y[i] -= sin(S_angle[i]) * S_move[i];	
					S_move[i] += rPL(SMYU);
				}
				S_loc_X[i] = aPBC(S_loc_X[i], L);
				S_loc_Y[i] = aPBC(S_loc_Y[i], L);
			}

			// check for encounters
			for(j=0; j<n1; j++){
				if(Dis_to_Nearest[j] < ER){	// not check when the senders cannot encounter from the theoretically possible shortest distance to the nearest sender
					Dis_to_Nearest[j] = 99999;	// reset Dis_to_Nearest
					for(h=0; h<n2; h++){
						xdis = min(abs(R_loc_X[j]-S_loc_X[h]), abs(L-abs(R_loc_X[j]-S_loc_X[h])));
						ydis = min(abs(R_loc_Y[j]-S_loc_Y[h]), abs(L-abs(R_loc_Y[j]-S_loc_Y[h])));
						dis = sqrt((xdis*xdis) + (ydis*ydis));
						if(dis < Dis_to_Nearest[j]){
							Dis_to_Nearest[j] = dis;
							Nearest_S[j] = h;
							if( Dis_to_Nearest[j] < 1 ){	// check encounter
								// when encounter, redistribute new individuals
								R_loc_X[j] =  L * rnd();
								R_loc_Y[j] =  L * rnd();
								S_loc_X[h] =  L * rnd();
								S_loc_Y[h] =  L * rnd();
								Dis_to_Nearest[j] = ER+0.1;
								R_angle[j] = nipi * rnd();
								S_angle[h] = nipi * rnd();
								R_move[j] = rPL(RMYU);
								S_move[j] = rPL(SMYU);
								pairs ++;	// count up
								break;
							}
						}
					}
				}
			}
		}
		/* // comment out when showing the time development of individual locations and encounter events
		cv::Mat Img( cv::Size(L, L), CV_8UC3, cv::Scalar(255, 255, 255));
		cv::namedWindow("search image", cv::WINDOW_AUTOSIZE);
		for(j=0; j<n2; j++){
			cv::circle(Img, cv::Point(S_loc_X[j], S_loc_Y[j]), ER, cv::Scalar(0,200,200), -1, CV_AA);
		}
		for(j=0; j<n1; j++){
			cv::circle(Img, cv::Point(R_loc_X[j], R_loc_Y[j]), 1, cv::Scalar(200,0,0), -1, CV_AA);
		}
		for(j=0; j<n2; j++){
			cv::circle(Img, cv::Point(S_loc_X[j], S_loc_Y[j]), 1, cv::Scalar(0,0,200), -1, CV_AA);
		}
		cv::imshow("search image", Img);
		cv::imwrite("perce0.png", Img);
		cv::waitKey(1);
		*/
	}
	return(pairs);
}

//// main ////
int main(){
	int R_myu, S_myu, res1;
	double length;
	std::string FileName;

	// messages
	cout << "Simulations for Searching Efficiency (Levy ver)" << endl;
	cout << "probability of correct navigation (p)? ==> ";
	cin >> p;
	cout << "L2? ==> ";
	cin >> L2;
	cout << "Er? (1, 10) ==> ";
	cin >> Er;

	// create an file for output (csv file)
	FileName = "res_n-";
	FileName += std::to_string(n) ;
	if(p < 1){
		FileName += "_p-";
		FileName += std::to_string(p);
	}	
	FileName += "_L2-";
	FileName += std::to_string(L2);
	FileName += "_Er-";
	FileName += std::to_string(Er);
	FileName += "_vR-";
	FileName += std::to_string(vR);
	FileName += "_vS-";
	FileName += std::to_string(vS);
	FileName += "_step-10^";
	FileName += std::to_string(step_scale);
	FileName += ".csv";
	std::ofstream ofs( FileName );
	ofs << "receiver_myu, sender_myu, pair";
	ofs << std::endl;

	length = sqrt(L2);

	res1;
	for(R_myu = R_myu_min; R_myu < R_myu_max +1; R_myu++){
		for(S_myu = S_myu_min; S_myu < S_myu_max +1; S_myu++){
			printf("receiver_myu:%.2f, sender_myu:%.2f,        \r", 1.1+0.1*R_myu, 1.1+0.1*S_myu);
			res1 = MS(1.1+0.1*R_myu, 1.1+0.1*S_myu, length, Er, p);
			ofs << 1.1+0.1*R_myu << "," << 1.1+0.1*S_myu << "," << res1 << endl;
		}
	}
}

