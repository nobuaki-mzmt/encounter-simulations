// DensitySepSearch.cpp

// Encounter sim
// More generalized version
// 191105 N Mizumoto created
// 191220 Major update

// Summary
// 2D conditions in PBC including multiple panels, focused on the center panel
// a male and a female search for each other

// 1 step = 0.2 sec

// The simulation program was implemented in Microsoft Visual Studio C++ 2017 (Windows 10).

// library
#include "stdafx.h"
#include <iostream>	// cout; cin
#include <sstream>	// to_string
#include <fstream>	// writing
using namespace std;

///*// opencv
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
using namespace cv;
//*/

// Random sampling
#include <random>
#include <time.h>
std::uniform_real_distribution<double> distribution(0.0, 1.0);
std::random_device seed_gen;
std::mt19937 engine(seed_gen());
auto rnd = [&] { return distribution(engine); };

#define pi 3.1415926535
#define nipi 3.1415926535*2


//// settings ////
// movement parameters
// {female, male(low density), male(high density), others}
double speed1[4] = { 6.77, 17.05, 27.02, 12.88 };	// speed: speed parameter v (mm/sec)
double speed2[4] = { 6.77, 29.38, 29.38, 12.88 };	// speed: speed parameter v (mm/sec)
													// {female, male, male, others}
double sinuous1[4] = { 0.73, 0.64, 0.64, 0.78 };			// angle: sinuousity parameter 
double sinuous2[4] = { 0.69, 0.72, 0.72, 0.78 };			// angle: sinuousity parameter 

															// environmental parameters
int steps = 200;	// number of time steps (*0.2 = seconds)
int detection = 10;	// the distance between two becomes smaller than φ, they encounter
int iter;			// number of replicates
double L[4] = { 10000, 124, 124, 124 };	// the size of PBC for conditions before encounters 
int n[4] = { 1, 2, 4, 8 };	// number of pairs
double d = 15;	// female and male are initially separated by the distance d for conditions after separation
#define num_ind 8

				//// variables ////
double MLocX[num_ind];			// location (x-axis)
double MLocY[num_ind];			// location (y-axis)
double FLocX[num_ind];			// location (x-axis)
double FLocY[num_ind];			// location (y-axis)
double MAngle[num_ind];		// direction of movement
double FAngle[num_ind];		// direction of movement
vector<int> MSearch(num_ind);		// encounter or searching
vector<int> FSearch(num_ind);		// encounter or searching
std::string FileName;	// result name
vector<int> MovementPattern(4);		// put the sets of movement pattern parameters
double Dis_to_Nearest;			// Distance from each receiver to its nearest sender

								//// functions ////
								// random angle
double dir;
double allocate() {
	dir = rnd()*nipi;
	return(dir);
}

// Wrapped Cauchy distribution (Bartumeus and Levin 2008)
double U, s, Y, Z;
double WrappedCauchy(double rho) {
	U = 1 - rnd();
	Y = 2 * atan(((1 - rho) / (1 + rho)) * tan(pi*(U - 0.5)));
	return(Y);
}

// PBC
// modification for PBC
double aPBC(double Loc, double Length) {
	if (Loc >= Length) {
		Loc -= Length;
	}
	else if (Loc < 0) {
		Loc += Length;
	}
	return(Loc);
}

// Random search
double xdis, ydis, sepdis;
int i, j, k, simtime, endtime, endtimeF, endtimeM;
double Ran_search(vector<int>& MovementPattern, int& iter, int& steps, vector<int>& res_f_reunion, vector<int>& res_f_change,
	vector<int>& res_m_reunion, vector<int>& res_m_change, double& L, int& n, int emp) {
	for (i = 0; i < iter; i++) {
		cout << "\r" << "female:" << MovementPattern[0] << ", male:" << MovementPattern[1] << ", L:" << L << ", n:" << n << ", rep:" << i + 1;

		// initialization
		fill(MSearch.begin(), MSearch.end(), 1);
		fill(FSearch.begin(), FSearch.end(), 1);
		for (j = 1; j < (8 - n + 1); j++) {	// search for encountering pairs
			FSearch[j] = 0;
			MSearch[j] = 0;
		}

		MLocX[0] = L / 2 - (0.25 * 1.41421356*d);
		MLocY[0] = L / 2 + 0.25 * 1.41421356*d;
		FLocX[0] = L / 2 + 0.25 * 1.41421356*d;
		FLocY[0] = L / 2 - 0.25 * 1.41421356*d;


		// redistribute individuals which encountering at the initial condition
		for (j = 1; j < num_ind; j++) {	// search for encountering pairs
			Dis_to_Nearest = 0;
			while (Dis_to_Nearest < d) {
				MLocX[j] = rnd() * L;
				MLocY[j] = rnd() * L;
				xdis = min(abs(MLocX[j] - FLocX[0]), abs(L - abs(MLocX[j] - FLocX[0])));
				ydis = min(abs(MLocY[j] - FLocY[0]), abs(L - abs(MLocY[j] - FLocY[0])));
				Dis_to_Nearest = sqrt((xdis*xdis) + (ydis*ydis));
			}
			Dis_to_Nearest = 0;
			while (Dis_to_Nearest < d) {
				FLocX[j] = rnd() * L;
				FLocY[j] = rnd() * L;
				xdis = min(abs(MLocX[0] - FLocX[j]), abs(L - abs(MLocX[0] - FLocX[j])));
				ydis = min(abs(MLocY[0] - FLocY[j]), abs(L - abs(MLocY[0] - FLocY[j])));
				Dis_to_Nearest = sqrt((xdis*xdis) + (ydis*ydis));
			}
		}

		MAngle[0] = 0.75 * pi;
		FAngle[0] = -0.25 * pi;
		for (j = 1; j < num_ind; j++) {
			MAngle[j] = rnd()*nipi;
			FAngle[j] = rnd()*nipi;
		}

		endtime = steps;
		endtimeF = steps;
		endtimeM = steps;
		for (simtime = 0; simtime < steps + 1; simtime++) {

			if (simtime <= 10) {
				// change in position
				FLocX[0] += cos(FAngle[0])*0.2*speed1[MovementPattern[0]];
				FLocY[0] += sin(FAngle[0])*0.2*speed1[MovementPattern[0]];
				MLocX[0] += cos(MAngle[0])*0.2*speed1[MovementPattern[1]];
				MLocY[0] += sin(MAngle[0])*0.2*speed1[MovementPattern[1]];

				for (j = 1; j < num_ind; j++) {
					FLocX[j] += cos(FAngle[j])*0.2*speed1[3];
					FLocY[j] += sin(FAngle[j])*0.2*speed1[3];
					MLocX[j] += cos(MAngle[j])*0.2*speed1[3];
					MLocY[j] += sin(MAngle[j])*0.2*speed1[3];
				}

				for (j = 0; j < num_ind; j++) {
					FLocX[j] = aPBC(FLocX[j], L);
					MLocX[j] = aPBC(MLocX[j], L);
					FLocY[j] = aPBC(FLocY[j], L);
					MLocY[j] = aPBC(MLocY[j], L);
				}

				//if (simtime % 5 == 4) {
				MAngle[0] += WrappedCauchy(sinuous1[MovementPattern[1]]);
				FAngle[0] += WrappedCauchy(sinuous1[MovementPattern[0]]);
				for (j = 1; j < num_ind; j++) {
					MAngle[j] += WrappedCauchy(sinuous1[3]);
					FAngle[j] += WrappedCauchy(sinuous1[3]);
				}
				//}
			}
			else {
				// change in position
				FLocX[0] += cos(FAngle[0])*0.2*speed2[MovementPattern[0]];
				FLocY[0] += sin(FAngle[0])*0.2*speed2[MovementPattern[0]];
				MLocX[0] += cos(MAngle[0])*0.2*speed2[MovementPattern[1]];
				MLocY[0] += sin(MAngle[0])*0.2*speed2[MovementPattern[1]];

				for (j = 1; j < num_ind; j++) {
					FLocX[j] += cos(FAngle[j])*0.2*speed2[3];
					FLocY[j] += sin(FAngle[j])*0.2*speed2[3];
					MLocX[j] += cos(MAngle[j])*0.2*speed2[3];
					MLocY[j] += sin(MAngle[j])*0.2*speed2[3];
				}

				for (j = 0; j < num_ind; j++) {
					FLocX[j] = aPBC(FLocX[j], L);
					MLocX[j] = aPBC(MLocX[j], L);
					FLocY[j] = aPBC(FLocY[j], L);
					MLocY[j] = aPBC(MLocY[j], L);
				}

				//if (simtime % 5 == 4) {
				MAngle[0] += WrappedCauchy(sinuous2[MovementPattern[1]]);
				FAngle[0] += WrappedCauchy(sinuous2[MovementPattern[0]]);
				for (j = 1; j < num_ind; j++) {
					MAngle[j] += WrappedCauchy(sinuous2[3]);
					FAngle[j] += WrappedCauchy(sinuous2[3]);
				}
				//}
			}

			// encounter determination
			if (FSearch[0] == 1) {
				if (MSearch[0] == 1) {
					xdis = min(abs(FLocX[0] - MLocX[0]), abs(L - abs(FLocX[0] - MLocX[0])));
					ydis = min(abs(FLocY[0] - MLocY[0]), abs(L - abs(FLocY[0] - MLocY[0])));
					sepdis = sqrt(xdis*xdis + ydis*ydis);
					if (sepdis <= detection) {
						endtime = simtime + 1;
						FSearch[0] = 0;
						MSearch[0] = 0;
					}
				}
			}
			if (FSearch[0] == 1) {
				for (j = 0; j < num_ind; j++) {
					if (MSearch[j] == 1) {
						xdis = min(abs(FLocX[0] - MLocX[j]), abs(L - abs(FLocX[0] - MLocX[j])));
						ydis = min(abs(FLocY[0] - MLocY[j]), abs(L - abs(FLocY[0] - MLocY[j])));
						sepdis = sqrt(xdis*xdis + ydis*ydis);
						if (sepdis <= detection) {
							endtimeF = simtime + 1;
							FSearch[0] = 0;
							MSearch[j] = 0;
							break;
						}
					}
				}
			}
			if (MSearch[0] == 1) {
				for (j = 0; j < num_ind; j++) {
					if (FSearch[j] == 1) {
						xdis = min(abs(FLocX[j] - MLocX[0]), abs(L - abs(FLocX[j] - MLocX[0])));
						ydis = min(abs(FLocY[j] - MLocY[0]), abs(L - abs(FLocY[j] - MLocY[0])));
						sepdis = sqrt(xdis*xdis + ydis*ydis);
						if (sepdis <= detection) {
							endtimeM = simtime + 1;
							FSearch[j] = 0;
							MSearch[0] = 0;
							break;
						}
					}
				}
			}
			for (j = 0; j < num_ind; j++) {
				if (FSearch[j] == 1) {
					for (k = 0; k < num_ind; k++) {
						if (MSearch[k] == 1) {
							xdis = min(abs(FLocX[j] - MLocX[k]), abs(L - abs(FLocX[j] - MLocX[k])));
							ydis = min(abs(FLocY[j] - MLocY[k]), abs(L - abs(FLocY[j] - MLocY[k])));
							sepdis = sqrt(xdis*xdis + ydis*ydis);
							if (sepdis <= detection) {
								FSearch[j] = 0;
								MSearch[k] = 0;
								break;
							}
						}
					}
				}
			}

			/*///*
			cv::Mat Img(cv::Size(L, L), CV_8UC3, cv::Scalar(255, 255, 255));
			cv::namedWindow("search image", cv::WINDOW_AUTOSIZE);
			if (FSearch[0] == 1) {
			cv::circle(Img, cv::Point(FLocX[0], FLocY[0]), 7, cv::Scalar(0, 0, 200), 1, CV_AA);
			}
			if (MSearch[0] == 1) {
			cv::circle(Img, cv::Point(MLocX[0], MLocY[0]), 7, cv::Scalar(200, 0, 0), 1, CV_AA);
			}
			for (j = 1; j < num_ind; j++) {
			if (FSearch[j] == 1) {
			cv::circle(Img, cv::Point(FLocX[j], FLocY[j]), 7, cv::Scalar(0, 0, 200), -1, CV_AA);
			}
			if (MSearch[j] == 1) {
			cv::circle(Img, cv::Point(MLocX[j], MLocY[j]), 7, cv::Scalar(200, 0, 0), -1, CV_AA);
			}
			}

			cv::imshow("search image", Img);
			cv::waitKey(1);
			//*///*///*/

			if (MSearch[0] == 0 && FSearch[0] == 0) {
				break;
			}
		}
		if (endtime < endtimeM) {
			for (j = endtime; j < steps; j++) {
				res_f_reunion[j] ++;
				res_m_reunion[j] ++;
			}
		}
		else {
			for (j = endtimeM; j < steps; j++) {
				res_m_change[j] ++;
			}
			for (j = endtimeF; j < steps; j++) {
				res_f_change[j] ++;
			}
		}
	}
	cout << endl;
	return(0);
}

int OK, search_mode, emp;
int main() {
	// chech the setting
	cout << "Density Sep Search Simulations" << endl;
	cout << "Brownian walk version (different speed)" << endl;

	cout << "Analysisng steps" << endl;
	cout << "1 step = 0.2 sec" << endl;
	cout << "60 sec = 300 steps, OK?" << endl;
	cin >> OK;
	if (OK != 1) {
		cout << "Enter steps (*0.2=sec)" << endl;
		cin >> steps;
	}
	cout << endl;

	cout << "Enter the num of rep (enter number)" << endl;
	cin >> iter;
	cout << endl;

	cout << "Detection range = " << detection << endl;
	cout << "The separated distance = " << d << endl;
	cout << "The size of system (PBC)" << endl;
	cout << "L = " << L[0] << "and" << L[1] << endl;
	cout << "The speed (<2s) = " << speed1[0] << ", " << speed1[1] << ", " << speed1[2] << " and " << speed1[3] << endl;
	cout << "The speed (>2s) = " << speed2[0] << ", " << speed2[1] << ", " << speed2[2] << " and " << speed2[3] << endl;

	cout << "Are you OK?" << endl;
	cout << "Yes: 1, No: 2" << endl;
	cin >> OK;
	if (OK == 2) { return(2); }

	// output data
	vector<int> res_f_reunion(steps);
	vector<int> res_m_reunion(steps);
	vector<int> res_f_change(steps);
	vector<int> res_m_change(steps);
	FileName = "EnSim-Brownian";
	FileName += "_rep";
	FileName += std::to_string(iter);
	FileName += "_sec";
	FileName += std::to_string(int(steps*0.2));

	FileName += "_L";
	FileName += std::to_string(int(L[0]));
	FileName += "-";
	FileName += std::to_string(int(L[1]));

	FileName += ".csv";

	std::ofstream ofs(FileName);
	ofs << "Sex,Encounter,Fspeed,Mspeed,L,n";
	for (i = 0; i < steps; i++) {
		ofs << "," << i + 1;
	}
	ofs << endl;

	// simulations are performed in order of
	// "observed movement after separation", "observed movement before encounter",
	// "virtual movements with females after separation", "virtual movements with males after separation"
	//int MaleMovementPattern[2] = { 1,2 };	// movement patterns for males
	//int FemaleMovementPattern[2] = { 0,0 };	// movement patterns for females
	int MaleMovementPattern[2] = { 2,2 };	// movement patterns for males
	int FemaleMovementPattern[2] = { 0,3 };	// movement patterns for females

	int p, q;
	for (p = 0; p < 2; p++) {
		// run simulations
		// put movement patterns (0 for female, 1 for male)
		MovementPattern[1] = MaleMovementPattern[p];
		MovementPattern[0] = FemaleMovementPattern[p];

		for (q = 0; q < 4; q++) {
			//if (q == 1 || q == 2) { continue; }
			fill(res_f_reunion.begin(), res_f_reunion.end(), 0);
			fill(res_m_reunion.begin(), res_m_reunion.end(), 0);
			fill(res_f_change.begin(), res_f_change.end(), 0);
			fill(res_m_change.begin(), res_m_change.end(), 0);
			// simulation function for random search (movement patterns, replications, num of steps, container for results, size of PBC)
			//
			Ran_search(MovementPattern, iter, steps, res_f_reunion, res_f_change, res_m_reunion, res_m_change, L[q], n[q], emp);

			// for result outputs
			for (j = 0; j < 2; j++) {
				for (k = 0; k < 2; k++) {
					if (j == 0) {
						ofs << "female" << ",";
					}
					else {
						ofs << "male" << ",";
					}

					if (k == 0) {
						ofs << "reunion" << ",";
					}
					else {
						ofs << "change" << ",";
					}

					for (i = 0; i < 2; i++) {
						ofs << speed1[MovementPattern[i]] << ",";
					}
					ofs << L[q] << "," << n[q] << ",";

					if (k == 0) {
						if (j == 0) {
							for (i = 0; i < steps - 1; i++) {
								ofs << res_f_reunion[i] << ",";
							}
							ofs << res_f_reunion[steps - 1] << endl;
						}
						else {
							for (i = 0; i < steps - 1; i++) {
								ofs << res_m_reunion[i] << ",";
							}
							ofs << res_m_reunion[steps - 1] << endl;
						}
					}
					else {
						if (j == 0) {
							for (i = 0; i < steps - 1; i++) {
								ofs << res_f_change[i] << ",";
							}
							ofs << res_f_change[steps - 1] << endl;
						}
						else {
							for (i = 0; i < steps - 1; i++) {
								ofs << res_m_change[i] << ",";
							}
							ofs << res_m_change[steps - 1] << endl;
						}
					}
				}

				if (j == 0) {
					ofs << "female" << ",";
				}
				else {
					ofs << "male" << ",";
				}

				ofs << "total" << ",";

				for (i = 0; i < 2; i++) {
					ofs << speed1[MovementPattern[i]] << ",";
				}
				ofs << L[q] << "," << n[q] << ",";

				for (i = 0; i < steps - 1; i++) {
					if (j == 0) {
						ofs << res_f_reunion[i] + res_f_change[i] << ",";
					}
					else {
						ofs << res_m_reunion[i] + res_m_change[i] << ",";
					}
				}
				if (j == 0) {
					ofs << res_f_change[steps - 1] << endl;
				}
				else {
					ofs << res_m_change[steps - 1] << endl;
				}
			}
		}
	}
	return 0;
}

