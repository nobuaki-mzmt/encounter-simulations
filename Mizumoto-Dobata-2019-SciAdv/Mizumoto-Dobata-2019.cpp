// Simulations of termite searches
// 170924 N Mizumoto created
// 180605 N Mizumoto commented
// 180919-180927 N Mizumoto updated
// 190215 N Mizumoto updated

// summary
// a male and a female search for each other
// movement from experimental data of termites (before/after - male/female)
// two conditions (before: PBC, after: NBC with distance d separated)

// 1 step = 0.2 sec
// CRW (resample every 0.2 sec)
// pause/move cut in 0.2*n sec

// The simulation program was implemented in Microsoft Visual Studio C++ 2017 (Windows 10).

//// libraries ////
#include "stdafx.h"
#include <iostream>					// cout; cin
#include <sstream>					// to_string
#include <fstream>					// writing
using namespace std;

/*// opencv (for plot)
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
// movement parameters (order is "female after, female before, male after, male before")
//*// for R. speratus
double MoveRho[4] = { 0.84, 0.87, 0.78, 0.86 };	// rho: sinuousity parameter (during moving)
double PauseRho[4] = { 0.73, 0.68, 0.68, 0.69 };	// rho: sinuousity parameter (after pause)
double speed[4] = { 5.4, 14.8, 15.0, 14.7 };	// speed: speed parameter v
double MovePattern[4] = { 0, 0, 1, 1 }; // move patterns, 0:Truncated power-law, 1:Stretched exponential
double MoveParam1[4] = { 2.145, 1.120, 0.499, 0.479 }; // myu for TP, lambda for SE
double MoveParam2[4] = { 8.4*5, 27.8*5, 0.705, 0.851 }; // Xmax for TP, beta for SE
double PausePattern[4] = { 0, 0, 0, 1 }; // move patterns, 0:Truncated power-law, 1:Stretched exponential
double PauseParam1[4] = { 1.488, 1.894, 2.563, 2.795 }; // myu for TP, lambda for SE
double PauseParam2[4] = { 661*5, 22*5, 2.4*5, 0.431 }; // Xmax for TP, beta for SE
double PauseProb[4] = {0.90796296, 0.2638889, 0.08719298, 0.2189474};	// Probability begining from pause
//*/

/*// for C. formosanus
double MoveRho[4] = { 0.86, 0.86, 0.80, 0.84 };	// rho: sinuousity parameter (during moving)
double PauseRho[4] = { 0.74, 0.68, 0.67, 0.75 };	// rho: sinuousity parameter (after pause)
double speed[4] = { 10.1, 22.3, 19.0, 20.7 };
double MovePattern[4] = { 0, 1, 0, 1 }; // move patterns, 0:Truncated power-law, 1:Stretched exponential
double MoveParam1[4] = { 1.666, 0.527, 1.361, 0.948 }; // myu for TP, lambda for SE
double MoveParam2[4] = { 49.2*5, 0.501, 64.2*5, 0.387 }; // Xmax for TP, beta for SE
double PausePattern[4] = { 0, 0, 0, 0 }; // move patterns, 0:Truncated power-law, 1:Stretched exponential
double PauseParam1[4] = { 1.791, 2.215, 2.204, 1.965 }; // myu for TP, lambda for SE
double PauseParam2[4] = { 41.4*5, 6*5, 6.6*5, 7.6*5 }; // Xmax for TP, beta for SE
double PauseProb[4] = { 0.5035, 0.0805, 0.1260606, 0.1171212 };	// Probability begining from pause
																//*/

																//other parameters
int steps = 900;	// number of time steps (*0.2 = seconds)
int detection = 7;	// the distance between two becomes smaller than φ, they encounter (7 for R.speratus and 10 for C.formosanus)
int iter;			// number of replicates
double L = 223.606;	// the size of PBC for conditions before encounters
double d = 16.767;	// female and male are initially separated by the distance d for conditions after separation (or 22.97364 for C.formosanus)

					//// variables ////
double LocX[2];			// location (x-axis)
double LocY[2];			// location (y-axis)
double Angle[2];		// direction of movement
int Pnokori[2];			// persistence of Pause duration
int Wnokori[2];			// persistence of Move duration
vector<int> MovementPattern(4);		// put the sets of movement pattern parameters
std::string FileName;	// result name


						/*// Empirical data of Pausing/Moving time for Reticulitermes speratus
						double CumProb_P[4][100] = { { 0.69822485, 0.56804734, 0.46745562, 0.42011834, 0.39053254, 0.33727811, 0.33136095, 0.30769231, 0.28402367, 0.26627219, 0.26035503, 0.25443787, 0.23076923, 0.22485207, 0.21301775, 0.20710059, 0.18934911, 0.17751479, 0.16568047, 0.15384615, 0.14792899, 0.14201183, 0.13609467, 0.13017751, 0.12426036, 0.11834320, 0.11242604, 0.10650888, 0.10059172, 0.09467456, 0.08875740, 0.08284024, 0.07692308, 0.07100592, 0.06508876, 0.05917160, 0.05325444, 0.04733728, 0.04142012, 0.03550296, 0.02958580, 0.02366864, 0.01775148, 0.01183432, 0.00591716, 0 },
						{ 0.559670782, 0.349794239, 0.259259259, 0.205761317, 0.172839506, 0.144032922, 0.119341564, 0.106995885, 0.098765432, 0.094650206, 0.090534979, 0.082304527, 0.078189300, 0.074074074, 0.069958848, 0.061728395, 0.057613169, 0.041152263, 0.037037037, 0.032921811, 0.028806584, 0.024691358, 0.020576132, 0.016460905, 0.012345679, 0.008230453, 0.004115226, 0 },
						{ 0.343434343, 0.148148148, 0.084175084, 0.053872054, 0.037037037, 0.020202020, 0.013468013, 0.010101010, 0.003367003, 0 },
						{ 0.641104294, 0.444785276, 0.297546012, 0.217791411, 0.187116564, 0.147239264, 0.125766871, 0.098159509, 0.088957055, 0.076687117, 0.067484663, 0.058282209, 0.052147239, 0.042944785, 0.036809816, 0.030674847, 0.024539877, 0.018404908, 0.015337423, 0.012269939, 0.009202454, 0.006134969, 0.003067485, 0 } };
						double StepLabel_P[4][100] = { { 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 4.0, 4.2, 4.4, 5.2, 5.4, 5.6, 5.8, 6.0, 6.6, 8.0, 8.8, 11.2, 14.0, 21.2, 30.6, 33.6, 34.4, 41.2, 49.8, 61.6, 89.2, 95.6, 109.8, 112.2, 284.8, 288.2, 309.2, 313.6, 331.0, 373.2, 661.0 },
						{ 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.0, 3.6, 3.8, 4.2, 4.4, 4.8, 5.8, 6.0, 8.4, 11.4, 12.2, 14.8, 19.6, 20.6, 22.0 },
						{ 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.2, 2.4 },
						{ 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 4.2, 4.8, 7.8, 8.0, 13.8, 14.6 } };
						double CumProb_M[4][100] = { { 0.443820225, 0.280898876, 0.179775281, 0.134831461, 0.095505618, 0.078651685, 0.061797753, 0.056179775, 0.050561798, 0.044943820, 0.039325843, 0.028089888, 0.022471910, 0.016853933, 0.011235955, 0.005617978, 0 },
						{ 0.839506173, 0.761316872, 0.703703704, 0.633744856, 0.604938272, 0.567901235, 0.543209877, 0.514403292, 0.485596708, 0.448559671, 0.427983539, 0.403292181, 0.399176955, 0.382716049, 0.366255144, 0.358024691, 0.349794239, 0.341563786, 0.337448560, 0.325102881, 0.316872428, 0.304526749, 0.288065844, 0.279835391, 0.263374486, 0.242798354, 0.234567901, 0.222222222, 0.209876543, 0.197530864, 0.193415638, 0.189300412, 0.185185185, 0.172839506, 0.156378601, 0.152263374, 0.148148148, 0.144032922, 0.139917695, 0.135802469, 0.131687243, 0.127572016, 0.119341564, 0.106995885, 0.098765432, 0.094650206, 0.086419753, 0.082304527, 0.078189300, 0.069958848, 0.065843621, 0.057613169, 0.053497942, 0.049382716, 0.045267490, 0.041152263, 0.037037037, 0.032921811, 0.028806584, 0.024691358, 0.016460905, 0.012345679, 0.008230453, 0 },
						{ 0.875000000, 0.800675676, 0.760135135, 0.709459459, 0.668918919, 0.638513514, 0.601351351, 0.574324324, 0.540540541, 0.510135135, 0.472972973, 0.456081081, 0.432432432, 0.398648649, 0.378378378, 0.364864865, 0.358108108, 0.337837838, 0.334459459, 0.317567568, 0.307432432, 0.287162162, 0.266891892, 0.246621622, 0.243243243, 0.233108108, 0.222972973, 0.219594595, 0.212837838, 0.195945946, 0.189189189, 0.185810811, 0.168918919, 0.158783784, 0.141891892, 0.138513514, 0.125000000, 0.118243243, 0.104729730, 0.101351351, 0.097972973, 0.094594595, 0.087837838, 0.084459459, 0.077702703, 0.074324324, 0.070945946, 0.067567568, 0.064189189, 0.060810811, 0.057432432, 0.054054054, 0.047297297, 0.043918919, 0.037162162, 0.033783784, 0.027027027, 0.023648649, 0.020270270, 0.016891892, 0.013513514, 0.010135135, 0.006756757, 0.003378378, 0 },
						{ 0.90372671, 0.85403727, 0.80434783, 0.74534161, 0.68012422, 0.62422360, 0.57453416, 0.51242236, 0.46894410, 0.41614907, 0.38198758, 0.35093168, 0.32298137, 0.28260870, 0.25776398, 0.23602484, 0.21739130, 0.20807453, 0.19565217, 0.19254658, 0.18633540, 0.18322981, 0.17391304, 0.15838509, 0.14906832, 0.14596273, 0.13975155, 0.12422360, 0.11490683, 0.11180124, 0.10869565, 0.09937888, 0.09627329, 0.09316770, 0.09006211, 0.08385093, 0.08074534, 0.07453416, 0.07142857, 0.06521739, 0.06211180, 0.05590062, 0.05279503, 0.04347826, 0.04037267, 0.03416149, 0.03105590, 0.02795031, 0.02484472, 0.02173913, 0.01863354, 0.01552795, 0.01242236, 0.00931677, 0.00621118, 0.00310559, 0 } };

						double StepLabel_M[4][100] = { { 0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.6, 2.4, 2.6, 2.8, 3.0, 3.4, 3.6, 3.8, 4.0, 5.6, 8.4 },
						{ 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.6, 6.8, 7.0, 7.2, 7.6, 7.8, 8.0, 8.4, 8.8, 9.4, 9.6, 10.0, 10.2, 10.4, 10.6, 11.6, 11.8, 12.8, 13.2, 13.4, 13.6, 13.8, 14.8, 16.6, 16.8, 17.0, 19.6, 20.4, 22.0, 25.0 },
						{ 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.4, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.8, 9.4, 9.6, 10.0, 10.2, 10.6, 11.2, 11.4, 11.6, 11.8, 12.2, 13.0, 13.8, 14.0, 16.0, 16.6, 17.0, 17.6, 18.6, 21.4, 22.2, 22.8, 24.6, 27.0, 27.4 },
						{ 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 7.2, 7.4, 7.6, 7.8, 8.0, 8.4, 8.6, 8.8, 9.0, 9.2, 10.0, 10.2, 10.8, 11.2, 11.4, 11.6, 12.4, 12.6, 14.0, 14.2, 17.0, 18.2, 18.8, 27.2 } };
						//*/
						//*// Empirical data of Pausing/Moving time for Coptotermes formosanus
double CumProb_P[4][100] = { { 0.5250, 0.3675, 0.3200, 0.2925, 0.2700, 0.2350, 0.2000, 0.1775, 0.1600, 0.1550, 0.1425, 0.1250, 0.1200, 0.1125, 0.1025, 0.0925, 0.0900, 0.0875, 0.0825, 0.0775, 0.0750, 0.0725, 0.0700, 0.0675, 0.0650, 0.0625, 0.0575, 0.0550, 0.0525, 0.0500, 0.0450, 0.0425, 0.0400, 0.0375, 0.0350, 0.0325, 0.0275, 0.0250, 0.0225, 0.0200, 0.0175, 0.0150, 0.0125, 0.0100, 0.0075, 0.0050, 0.0025, 0 },
{ 0.418367347, 0.219387755, 0.132653061, 0.117346939, 0.091836735, 0.086734694, 0.071428571, 0.066326531, 0.056122449, 0.051020408, 0.035714286, 0.030612245, 0.025510204, 0.020408163, 0.015306122, 0.005102041, 0 },
{ 0.426934097, 0.266475645, 0.189111748, 0.131805158, 0.091690544, 0.063037249, 0.057306590, 0.054441261, 0.040114613, 0.037249284, 0.034383954, 0.031518625, 0.020057307, 0.017191977, 0.014326648, 0.008595989, 0.005730659, 0.002865330, 0 },
{ 0.501976285, 0.316205534, 0.233201581, 0.193675889, 0.154150198, 0.130434783, 0.118577075, 0.090909091, 0.083003953, 0.075098814, 0.063241107, 0.047430830, 0.039525692, 0.035573123, 0.023715415, 0.019762846, 0.015810277, 0.011857708, 0.007905138, 0.003952569, 0 } };
double StepLabel_P[4][100] = { { 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.8, 5.2, 5.4, 5.6, 6.2, 6.4, 6.6, 7.6, 8.0, 8.2, 8.4, 8.8, 9.2, 9.6, 9.8, 11.2, 11.6, 12.2, 14.6, 17.4, 19.6, 19.8, 24.2, 28.4, 28.6, 34.0, 36.2, 41.4 },
{ 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 1.8, 2.0, 2.2, 2.6, 3.0, 5.0, 5.2, 5.4, 5.8, 6.0 },
{ 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 3.2, 3.6, 4.0, 4.6, 5.6, 6.6 },
{ 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.4, 3.8, 5.4, 6.0, 6.2, 6.6, 7.6 } };
double CumProb_M[4][100] = { { 0.593908629, 0.467005076, 0.395939086, 0.345177665, 0.304568528, 0.276649746, 0.246192893, 0.218274112, 0.197969543, 0.180203046, 0.172588832, 0.159898477, 0.154822335, 0.149746193, 0.134517766, 0.124365482, 0.121827411, 0.109137056, 0.106598985, 0.101522843, 0.096446701, 0.093908629, 0.088832487, 0.083756345, 0.081218274, 0.078680203, 0.073604061, 0.068527919, 0.063451777, 0.060913706, 0.058375635, 0.055837563, 0.053299492, 0.050761421, 0.045685279, 0.043147208, 0.040609137, 0.038071066, 0.035532995, 0.032994924, 0.030456853, 0.027918782, 0.025380711, 0.022842640, 0.020304569, 0.017766497, 0.015228426, 0.012690355, 0.010152284, 0.007614213, 0.002538071, 0 },
{ 0.880, 0.825, 0.800, 0.775, 0.720, 0.695, 0.670, 0.620, 0.595, 0.570, 0.550, 0.540, 0.535, 0.510, 0.495, 0.490, 0.485, 0.475, 0.460, 0.445, 0.430, 0.425, 0.420, 0.415, 0.410, 0.405, 0.400, 0.385, 0.375, 0.370, 0.360, 0.350, 0.345, 0.335, 0.320, 0.310, 0.305, 0.290, 0.285, 0.280, 0.275, 0.270, 0.255, 0.245, 0.240, 0.235, 0.215, 0.210, 0.205, 0.200, 0.185, 0.180, 0.175, 0.165, 0.160, 0.155, 0.150, 0.145, 0.130, 0.125, 0.120, 0.115, 0.110, 0.105, 0.100, 0.095, 0.090, 0.085, 0.075, 0.070, 0.065, 0.055, 0.050, 0.045, 0.040, 0.030, 0.025, 0.020, 0.015, 0.010, 0.005, 0 },
{ 0.731988473, 0.616714697, 0.579250720, 0.536023055, 0.501440922, 0.472622478, 0.449567723, 0.429394813, 0.406340058, 0.386167147, 0.365994236, 0.348703170, 0.334293948, 0.314121037, 0.308357349, 0.299711816, 0.291066282, 0.270893372, 0.256484150, 0.250720461, 0.247838617, 0.236311239, 0.224783862, 0.210374640, 0.201729107, 0.195965418, 0.187319885, 0.184438040, 0.175792507, 0.172910663, 0.170028818, 0.158501441, 0.155619597, 0.152737752, 0.149855908, 0.141210375, 0.138328530, 0.132564841, 0.129682997, 0.126801153, 0.118155620, 0.115273775, 0.109510086, 0.103746398, 0.097982709, 0.089337176, 0.086455331, 0.083573487, 0.080691643, 0.077809798, 0.072046110, 0.066282421, 0.063400576, 0.060518732, 0.057636888, 0.054755043, 0.051873199, 0.048991354, 0.046109510, 0.040345821, 0.037463977, 0.034582133, 0.031700288, 0.028818444, 0.025936599, 0.023054755, 0.020172911, 0.017291066, 0.014409222, 0.011527378, 0.008645533, 0.002881844, 0 },
{ 0.803149606, 0.728346457, 0.681102362, 0.649606299, 0.618110236, 0.598425197, 0.578740157, 0.562992126, 0.515748031, 0.492125984, 0.488188976, 0.468503937, 0.440944882, 0.417322835, 0.397637795, 0.381889764, 0.370078740, 0.362204724, 0.334645669, 0.326771654, 0.322834646, 0.314960630, 0.311023622, 0.307086614, 0.299212598, 0.295275591, 0.291338583, 0.279527559, 0.275590551, 0.271653543, 0.259842520, 0.251968504, 0.244094488, 0.240157480, 0.232283465, 0.224409449, 0.220472441, 0.208661417, 0.204724409, 0.196850394, 0.192913386, 0.181102362, 0.173228346, 0.169291339, 0.165354331, 0.161417323, 0.153543307, 0.149606299, 0.145669291, 0.137795276, 0.129921260, 0.125984252, 0.122047244, 0.114173228, 0.110236220, 0.106299213, 0.102362205, 0.098425197, 0.094488189, 0.090551181, 0.086614173, 0.082677165, 0.078740157, 0.074803150, 0.070866142, 0.062992126, 0.059055118, 0.055118110, 0.051181102, 0.047244094, 0.043307087, 0.039370079, 0.035433071, 0.031496063, 0.027559055, 0.023622047, 0.019685039, 0.015748031, 0.011811024, 0.007874016, 0.003937008, 0 } };
double StepLabel_M[4][100] = { { 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 5.0, 5.4, 5.8, 6.0, 6.4, 6.6, 7.4, 8.2, 8.4, 9.4, 11.0, 11.2, 13.0, 14.4, 14.6, 15.6, 16.2, 16.4, 16.8, 18.8, 19.4, 20.4, 21.6, 22.0, 23.6, 26.6, 27.0, 30.2, 49.2 },
{ 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 6.2, 6.4, 6.6, 6.8, 7.2, 7.4, 7.6, 7.8, 8.4, 8.6, 8.8, 9.4, 9.6, 9.8, 10.0, 10.2, 10.4, 10.8, 11.0, 11.2, 11.4, 11.8, 12.0, 12.2, 12.6, 12.8, 13.0, 13.2, 13.4, 13.8, 14.0, 14.2, 14.4, 14.6, 14.8, 15.4, 15.6, 16.0, 16.2, 16.4, 17.0, 17.2, 17.6, 18.8, 19.0, 24.4, 24.6, 24.8, 28.2, 32.0, 32.6, 53.2, 62.4, 63.2, 75.2, 483.8 },
{ 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.6, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 7.0, 7.2, 7.4, 7.6, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0, 10.2, 10.6, 10.8, 11.0, 11.4, 11.6, 12.6, 13.0, 13.8, 14.2, 14.4, 15.0, 15.2, 15.4, 16.0, 16.2, 21.0, 21.2, 21.8, 25.0, 27.8, 28.0, 30.6, 45.0, 45.8, 53.2, 64.2 },
{ 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.4, 6.6, 6.8, 7.0, 7.2, 7.6, 8.0, 8.2, 8.4, 8.6, 8.8, 9.2, 9.4, 9.6, 9.8, 10.0, 10.2, 10.4, 10.6, 11.0, 11.4, 11.8, 12.0, 12.2, 12.6, 13.0, 13.8, 14.0, 14.2, 14.6, 15.4, 15.6, 16.2, 16.8, 17.2, 17.6, 18.8, 19.0, 20.0, 21.2, 24.2, 26.2, 27.6, 31.0, 33.0, 37.4, 38.8, 45.4, 51.6, 75.0, 133.0, 328.2 } };
//*/

//// functions ////
// random angle
double dir;
double allocate() {
	dir = rnd()*nipi;
	return(dir);
}

// Truncated power-law (Bartumeus et al 2010)
double U, OBJ;
double TP(double Myu, double step_max, double step_min = 1) {
	U = 1 - rnd();
	OBJ = pow(pow(step_max, 1 - Myu) - U * (pow(step_max, 1 - Myu) - pow(step_min, 1 - Myu)), 1 / (1 - Myu));
	return(OBJ);
}

// Stretched exponention (Bartumeus et al 2010)
double SE(double Rambda, double Beta, double step_min = 1) {
	U = 1 - rnd();
	OBJ = pow(pow(step_min, Beta) - 1 / Rambda*log(U), 1 / Beta);
	return(OBJ);
}

// Wrapped Cauchy distribution (Bartumeus and Levin 2008)
double s, Y, Z;
double WrappedCauchy(double rho) {
	U = 1 - rnd();
	Y = 2 * atan(((1 - rho) / (1 + rho)) * tan(pi*(U - 0.5)));
	return(Y);
}



// For empirical data simulations
int StepLength, check;
double EmpData(int pattern, int State) {
	U = rnd();
	if (State == 1) {
		for (check = 0; check < 100; check++) {
			if (U >= CumProb_P[pattern - 1][check]) {
				OBJ = StepLabel_P[pattern - 1][check];
				break;
			}
		}
	}
	else {
		for (check = 0; check < 100; check++) {
			if (U >= CumProb_M[pattern - 1][check]) {
				OBJ = StepLabel_M[pattern - 1][check];
				break;
			}
		}
	}
	StepLength = OBJ * 5;
	return(StepLength);
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

double Rn;
// initialization (1:f-after, 2:f-before, 3:m-after, 4:m-before)
// determine direction and move/pause status of individuals at the initial condition
int initial(int id, int pattern, int emp) {
	Rn = rnd();
	Angle[id] = allocate();
	if (emp != 1) {
		if (Rn <= PauseProb[pattern - 1]) {
			if (MovePattern[pattern - 1] == 0) {
				Wnokori[id] = 0; Pnokori[id] = TP(MoveParam1[pattern - 1], MoveParam2[pattern - 1]);
			}
			else {
				Wnokori[id] = 0; Pnokori[id] = SE(MoveParam1[pattern - 1], MoveParam2[pattern - 1]);
			}
		}
		else {
			if (MovePattern[pattern - 1] == 0) {
				Pnokori[id] = 0; Wnokori[id] = TP(MoveParam1[pattern - 1], MoveParam2[pattern - 1]);
			}
			else {
				Pnokori[id] = 0; Wnokori[id] = SE(MoveParam1[pattern - 1], MoveParam2[pattern - 1]);
			}
		}
	}
	else {
		if (Rn <= PauseProb[pattern - 1]) {
			Pnokori[id] = EmpData(pattern, 1);
			Wnokori[id] = 0;
		}
		else {
			Wnokori[id] = EmpData(pattern, 2);
			Pnokori[id] = 0;
		}
	}
	return(0);
}

// change position (1:f-after, 2:f-before, 3:m-after, 4:m-before)
int change_position(int &j, int pattern, int emp) {
	if (Pnokori[j] > 0) {
		Pnokori[j] -= 1;
		if (Pnokori[j] <= 0) {
			if (emp != 1) {
				if (MovePattern[pattern - 1] == 0) {
					Wnokori[j] = TP(MoveParam1[pattern - 1], MoveParam2[pattern - 1]);
				}
				else {
					Wnokori[j] = SE(MoveParam1[pattern - 1], MoveParam2[pattern - 1]);
				}
			}
			else {
				Wnokori[j] = EmpData(pattern, 2);
			}
			Pnokori[j] = 0;
			Angle[j] += WrappedCauchy(PauseRho[pattern - 1]);
		}
	}
	else {
		LocX[j] += cos(Angle[j]) * 0.2 * speed[pattern - 1];
		LocY[j] += sin(Angle[j]) * 0.2 * speed[pattern - 1];
		Wnokori[j] -= 1;
		if (Wnokori[j] <= 0) {
			if (emp != 1) {
				if (PausePattern[pattern - 1] == 0) {
					Pnokori[j] = TP(PauseParam1[pattern - 1], PauseParam2[pattern - 1]);
				}
				else {
					Pnokori[j] = SE(PauseParam1[pattern - 1], PauseParam2[pattern - 1]);
				}
			}
			else {
				Pnokori[j] = EmpData(pattern, 1);
			}
			Wnokori[j] = 0;
		}
		Angle[j] += WrappedCauchy(MoveRho[pattern - 1]);
	}
	return(0);
}


// Reunion search (after separation)
double xdis, ydis, sepdis;
int i, j, simtime, endtime;
double Sep_search(vector<int>& MovementPattern, int& iter, int& steps, vector<int>& res, double& d, int emp) {
	for (i = 0; i < iter; i++) {
		printf("female:%d, male:%d, rep:%d        \r", MovementPattern[0], MovementPattern[1], i + 1);
		// initialization
		LocX[0] = 0;
		LocY[0] = 0;
		LocX[1] = d;
		LocY[1] = 0;

		// begins from pause or walk randomly
		for (j = 0; j < 2; j++) {
			initial(j, MovementPattern[j], emp);
		}
		endtime = steps;	// record the encountered time for result output
		for (simtime = 0; simtime < steps + 1; simtime++) {
			// change in position
			for (j = 0; j < 2; j++) {
				change_position(j, MovementPattern[j], emp);
			}
			// encounter determination
			xdis = LocX[1] - LocX[0];
			ydis = LocY[1] - LocY[0];
			sepdis = sqrt(xdis*xdis + ydis*ydis);
			if (sepdis <= detection) {
				endtime = simtime + 1;
				break;
			}
			else if (sepdis > (steps - simtime) / 2) {
				break;
			}
			/*///*
			cv::Mat Img(cv::Size(500, 500), CV_8UC3, cv::Scalar(255, 255, 255));
			cv::namedWindow("search image", cv::WINDOW_AUTOSIZE);
			cv::circle(Img, cv::Point(LocX[0] + (500 / 2), LocY[0] + (500 / 2)), 7, cv::Scalar(0, 0, 200), -1, CV_AA);
			cv::circle(Img, cv::Point(LocX[1] + (500 / 2), LocY[1] + (500 / 2)), 7, cv::Scalar(200, 0, 0), -1, CV_AA);
			cv::imshow("search image", Img);
			cv::waitKey(1);
			//*///*/
		}

		for (j = endtime; j < steps; j++) {
			res[j] ++;
		}
	}
	return(0);
}

// Random search (before encounters)
double Ran_search(vector<int>& MovementPattern, int& iter, int& steps, vector<int>& res, double& L, int emp) {
	for (i = 0; i < iter; i++) {
		printf("female:%d, male:%d, rep:%d        \r", MovementPattern[0], MovementPattern[1], i + 1);
		// initialization
		LocX[0] = 0;
		LocY[0] = 0;
		LocX[1] = rnd()*L;
		LocY[1] = rnd()*L;

		// begins from pause or walk randomly
		for (j = 0; j < 2; j++) {
			initial(j, MovementPattern[j], emp);
		}
		endtime = steps;
		for (simtime = 0; simtime < steps + 1; simtime++) {
			// change in position
			for (j = 0; j < 2; j++) {
				change_position(j, MovementPattern[j], emp);
			}
			LocX[0] = aPBC(LocX[0], L);
			LocX[1] = aPBC(LocX[1], L);
			LocY[0] = aPBC(LocY[0], L);
			LocY[1] = aPBC(LocY[1], L);


			// encounter determination
			xdis = min(abs(LocX[1] - LocX[0]), abs(L - abs(LocX[1] - LocX[0])));
			ydis = min(abs(LocY[1] - LocY[0]), abs(L - abs(LocY[1] - LocY[0])));
			sepdis = sqrt(xdis*xdis + ydis*ydis);
			if (sepdis <= detection) {
				endtime = simtime + 1;
				break;
			}
			/*///*
			cv::Mat Img(cv::Size(L, L), CV_8UC3, cv::Scalar(255, 255, 255));
			cv::namedWindow("search image", cv::WINDOW_AUTOSIZE);
			cv::circle(Img, cv::Point(LocX[0], LocY[0]), 7, cv::Scalar(200, 0, 0), -1, CV_AA);
			cv::circle(Img, cv::Point(LocX[1], LocY[1]), 7, cv::Scalar(0, 0, 200), -1, CV_AA);
			cv::imshow("search image", Img);
			cv::waitKey(1);
			//*///*///*/
		}

		for (j = endtime; j < steps; j++) {
			res[j] ++;
		}
	}
	return(0);
}


//// main simulations ////
int OK, search_mode, emp;
int main() {

	// chech the setting
	cout << "Encounter simulation for Female (exp-surr) with Brownian Walers" << endl;
	cout << "Analysisng steps" << endl;
	cout << "1 step = 0.2 sec" << endl;
	cout << "180 sec = 900 steps, OK?" << endl;
	cin >> OK;
	if (OK != 1) {
		cout << "Enter steps (*0.2=sec)" << endl;
		cin >> steps;
	}
	cout << endl;
	cout << "Detection range" << endl;
	cout << "Press 1: 7 mm for R. speratus" << endl;
	cout << "Press 2: 10 mm for C. formosanus" << endl;
	cin >> OK;
	if (OK == 2) { detection = 10; }
	else if (OK != 1) {
		cout << "Enter the ditection radius (mm)" << endl;
		cin >> detection;
	}
	cout << endl;
	cout << "Enter the num of rep (enter number)" << endl;
	cin >> iter;
	cout << endl;
	cout << "Enter the search condition" << endl;
	cout << "0: random search, 1: reunion search" << endl;
	cin >> search_mode;
	cout << endl;
	if (search_mode == 0) {
		cout << "You selected random search" << endl;
		cout << "Enter the size of system (PBC)" << endl;
		cout << "Default: 100*5^0.5 = 223.606 mm, OK?" << endl;
		cin >> OK;
		if (OK != 1) {
			cout << "Enter (e.g., High: 50*5^0.5 = 111.803 mm)" << endl;
			cin >> L;
		}
	}
	else {
		cout << "You selected separation search" << endl;
		cout << "The separated distance (NBC)" << endl;
		cout << "Press 1: 16.767 mm for R.speratus" << endl;
		cout << "Press 2: 22.97364 mm for C.formosanus" << endl;
		cin >> OK;
		if (OK == 2) { d = 22.97364; }
		else if (OK != 1) {
			cout << "Enter the separation distance (mm)" << endl;
			cin >> d;
		}
	}

	cout << endl;
	cout << "Use empirical data?" << endl;
	cout << "Yes: 1, No: 2" << endl;
	cin >> OK;
	if (OK == 1) { emp = 1; }
	else { emp = 0; }
	cout << endl;
	cout << "Are you OK?" << endl;
	cout << "Yes: 1, No: 2" << endl;
	cin >> OK;
	if (OK == 2) { return(2); }

	// output data
	vector<int> res(steps);
	FileName = "EnSim";
	FileName += "_rep";
	FileName += std::to_string(iter);
	if (search_mode == 0) {
		FileName += "_RanSearch_L";
		FileName += std::to_string(int(L));
	}
	else {
		FileName += "_SepSearch_SepDis";
		FileName += std::to_string(int(d));
	}
	FileName += "_sec";
	FileName += std::to_string(int(steps*0.2));
	if (emp == 1) {
		FileName += "_empdata";
	}
	FileName += ".csv";
	std::ofstream ofs(FileName);
	ofs << "Female,Male";
	for (i = 0; i < steps; i++) {
		ofs << "," << i + 1;
	}
	ofs << endl;

	// simulations are performed in order of
	// "observed movement after separation", "observed movement before encounter",
	// "virtual movements with females after separation", "virtual movements with males after separation"
	// number indicates: 1:female after, 2:female before, 3:male after, 4:male before
	int MaleMovementPattern[4] = { 3,4,1,3 };	// movement patterns for males
	int FemaleMovementPattern[4] = { 1,2,1,3 };	// movement patterns for females
	int p;
	for (p = 0; p < 4; p++) {
		// run simulations
		// put movement patterns (0 for female, 1 for male)
		MovementPattern[1] = MaleMovementPattern[p];
		MovementPattern[0] = FemaleMovementPattern[p];
		fill(res.begin(), res.end(), 0);
		if (search_mode == 0) {
			// simulation function for random search (movement patterns, replications, num of steps, container for results, size of PBC)
			Ran_search(MovementPattern, iter, steps, res, L, emp);
		}
		else if (search_mode == 1) {
			// simulation function for random search (movement patterns, replications, num of steps, container for results, separated distance d)
			Sep_search(MovementPattern, iter, steps, res, d, emp);
		}

		// for result outputs
		cout << endl;
		for (j = 0; j < 2; j++) {
			switch (MovementPattern[j]) {
			case 1:
				ofs << "female_after" << ",";
				break;
			case 2:
				ofs << "female_before" << ",";
				break;
			case 3:
				ofs << "male_after" << ",";
				break;
			case 4:
				ofs << "male_before" << ",";
				break;
			}
		}
		for (i = 0; i < steps - 1; i++) {
			ofs << res[i] << ",";
		}
		ofs << res[steps - 1] << endl;

	}
}

