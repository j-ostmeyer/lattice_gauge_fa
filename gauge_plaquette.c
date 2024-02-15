#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "gauge_flags.h"
#include "gauge_aux.h"
#include "gauge_alg2group.h"

#include "gauge_plaquette.h"

double plaquette_av(double complex *u, unsigned *nnt, unsigned ns, unsigned nn, gauge_flags *mode){
	const unsigned links = nn/2;
	const unsigned nd = mode->space_dim;
	double plaquettes = 0;

	for(unsigned i = 0; i < ns; i++){
		for(unsigned nu = 1; nu < links; nu++){
			for(unsigned mu = 0; mu < nu; mu++){
				plaquettes += plaquette_fields(u, nnt, ns, nn, i, mu, nu, mode);
			}
		}
	}

	return plaquettes * 2/nd/(nd-1) / ns;
}

double plaquette_fields(double complex *u, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, unsigned nu, gauge_flags *mode){
	const unsigned n = mode->gauge_dim, mat_dim = n*n;
	const unsigned *nnl = nnt + ns*nn;
	double complex *z = mode->zdummy;

	copy_mat(u + nnl[pos*nn + mu], z + 3*mat_dim, n, 0);
	copy_mat(u + nnl[nnt[pos*nn + mu]*nn + nu], z + 2*mat_dim, n, 0);
	copy_mat(u + nnl[nnt[pos*nn + nu]*nn + mu], z + 1*mat_dim, n, 1);
	copy_mat(u + nnl[pos*nn + nu], z, n, 1);

	return trace_prod(z, 4, n);
}

void sum_of_plaquettes(double complex *u, double complex *pl, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, int direction, gauge_flags *mode){
	// calculates the sum of all plaquettes around U_mu(n)
	// direction:  1 only for U_{mu,nu}(n)
	// direction: -1 only for U_{mu,nu}(n)^+
	// direction:  0 all for i/2*(U_{mu,nu}(n) - U_{mu,nu}(n)^+)
	const unsigned n = mode->gauge_dim, mat_dim = n*n;
	double complex *z = pl + mat_dim;
	double complex *st = z + 6*mat_dim;
	const unsigned *nnl = nnt + ns*nn;

	sum_of_staples(u, st, z, nnt, ns, nn, pos, mu, mode);
	copy_mat(u + nnl[pos*nn + mu], z + 3*mat_dim, n, 0);

	mat_mul(z + 3*mat_dim, st, pl, n);

	if(direction == -1) dagger(pl, n);
	if(direction == 0) dagger_asym(pl, n);
}

void sum_of_staples(double complex *u, double complex *st, double complex *z, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, gauge_flags *mode){
	// calculates the sum of all staples around U_mu(n)
	const unsigned n = mode->gauge_dim, mat_dim = n*n, links = nn/2;
	double complex *prod = z + 3*mat_dim;

	for(unsigned i = 0; i < mat_dim; i++) st[i] = 0;
	
	for(unsigned nu = 0; nu < nn; nu++){
		if(nu % links == mu % links) continue;

		staple_fields(u, z, nnt, ns, nn, pos, mu, nu, mode);

		for(unsigned i = 0; i < mat_dim; i++) st[i] += prod[i];
	}
}

void staple_fields(double complex *u, double complex *z, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, unsigned nu, gauge_flags *mode){
	const unsigned n = mode->gauge_dim, mat_dim = n*n, links = nn/2;
	const unsigned *nnl = nnt + ns*nn;

	const int turn = nu >= links;

	copy_mat(u + nnl[nnt[pos*nn + mu]*nn + nu], z + 2*mat_dim, n, turn);
	copy_mat(u + nnl[nnt[pos*nn + nu]*nn + mu], z + 1*mat_dim, n, 1);
	copy_mat(u + nnl[pos*nn + nu], z, n, 1 - turn);

	mat_prod(z, 3, n);
}

double strong_coupling_plaquette(double beta, gauge_flags *mode){
	const unsigned d = mode->space_dim;
	switch(mode->gauge_group){
		case 1:
			//return 2*strong_coupling_pl_u1(beta, d);
			return 2*strong_coupling_pl_pade_u1(beta, d);
		case 2:
			//return strong_coupling_pl_su2(beta, d);
			return strong_coupling_pl_pade_su2(beta, d);
		case 3:
			//return strong_coupling_pl_su3(beta/2, d)/3;
			return strong_coupling_pl_pade_su3(beta/2, d)/3;
		default:
			return 0;
	}
}

double strong_coupling_pl_u1(double b, unsigned d){
	double p = b / 8.;
	p -= pow(b, 3) / 256.;
	p += pow(b, 5) * (d / 1024. - 11 / 6144.);
	p -= pow(b, 7) * (d / 4096. - 757 / 1572864.);
	p += (0.0001179218292236328-0.0001351038614908854*d+0.00003814697265625*d*d)*pow(b, 9);
	p += (-0.0000549834359575201+0.00005514174699783325*d-0.00001382827758789062*d*d)*pow(b, 11);
	p += (-0.00001882212729286895+0.00002578772190544341*d-0.00001236051321029663*d*d + 2.086162567138672e-6*d*d*d)*pow(b, 13);
	p += (8.28879524594832e-6-0.00001163224331793134*d+5.532056093215942e-6*d*d-8.94069671630859e-7*d*d*d)*pow(b, 15);

	return p;
}

double strong_coupling_pl_pade_u1(double b, double d){
	// 7th order Pade approximation using 16th order expansion
	// has extended convergence radius, but gives only 4-8 significant digits due to rounding errors
	// numerical stability gets worse with larger beta
	return (b/8. + ((-12974646990677 + 25354264107816*d - 19632612358080*pow(d,2) + 
					7550201272320*pow(d,3) - 1448268595200*pow(d,4) + 111476736000*pow(d,5))*
				pow(b,3))/
			(1792.*(16833264719 - 29034944040*d + 18798474240*pow(d,2) - 5439052800*pow(d,3) + 
					597196800*pow(d,4))) + ((-16817142476005 + 33556458843282*d - 
					26174544921360*pow(d,2) + 9950911810560*pow(d,3) - 1842849792000*pow(d,4) + 
					133772083200*pow(d,5))*pow(b,5))/
			(86016.*(16833264719 - 29034944040*d + 18798474240*pow(d,2) - 5439052800*pow(d,3) + 
					 597196800*pow(d,4))) + ((-85690276124445601 + 206530276810616400*d - 
					 205471846249588800*pow(d,2) + 108095412123763200*pow(d,3) - 
					 31762666465689600*pow(d,4) + 4953245810688000*pow(d,5) - 321052999680000*pow(d,6))
				 *pow(b,7))/
			(1.6515072e8*(16833264719 - 29034944040*d + 18798474240*pow(d,2) - 
						  5439052800*pow(d,3) + 597196800*pow(d,4))))/
		(1 + (3*(-1071401178137 + 2095918291628*d - 1625085253200*pow(d,2) + 
				 626010658560*pow(d,3) - 120340684800*pow(d,4) + 9289728000*pow(d,5))*pow(b,2))/
		 (56.*(16833264719 - 29034944040*d + 18798474240*pow(d,2) - 5439052800*pow(d,3) + 
			   597196800*pow(d,4))) + (5*(-6702008183149 + 13079522494806*d - 
			   10018435829328*pow(d,2) + 3760483539456*pow(d,3) - 692026675200*pow(d,4) + 
			   50164531200*pow(d,5))*pow(b,4))/
		 (10752.*(16833264719 - 29034944040*d + 18798474240*pow(d,2) - 5439052800*pow(d,3) + 
				  597196800*pow(d,4))) + ((-6625616580454381 + 15993559077377430*d - 
				  15937021694837520*pow(d,2) + 8398161262375200*pow(d,3) - 
				  2471998362624000*pow(d,4) + 386204184576000*pow(d,5) - 25082265600000*pow(d,6))*
			  pow(b,6))/
		 (1.29024e6*(16833264719 - 29034944040*d + 18798474240*pow(d,2) - 5439052800*pow(d,3) + 
					 597196800*pow(d,4))));
}

double strong_coupling_pl_su2(double b, unsigned d){
	double p = b / 4.;
	p -= pow(b, 3) / 96.;
	p += pow(b, 5) * (d / 512. - 5 / 1536.);
	p -= pow(b, 7) * (d / 1536. - 29 / 23040.);

	return p;
}

double strong_coupling_pl_pade_su2(double b, double d){
	// 7th order Pade approximation using 16th order expansion
	// has extended convergence radius and is numerically very stable
	return (b/4. - (5*(-4885734200 + 9494876445*d - 6023960973*pow(d,2) + 1765299690*pow(d,3) - 
					350320950*pow(d,4) + 45927000*pow(d,5))*pow(b,3))/
			(192.*(-624270496 + 1079133400*d + 163731525*pow(d,2) - 508863600*pow(d,3) + 
				   117369000*pow(d,4))) + ((-21907990656 + 213950313976*d + 607383501865*pow(d,2) - 
				   1034993172150*pow(d,3) + 263147081400*pow(d,4) - 65583756000*pow(d,5))*pow(b,5))
			/(64512.*(-624270496 + 1079133400*d + 163731525*pow(d,2) - 508863600*pow(d,3) + 
					117369000*pow(d,4))) + ((148537739520512 - 357822682005800*d + 
					278360263900775*pow(d,2) - 82393047985500*pow(d,3) + 26156806436250*pow(d,4) - 
					11231218215000*pow(d,5) + 578680200000*pow(d,6))*pow(b,7))/
			(1.548288e7*(-624270496 + 1079133400*d + 163731525*pow(d,2) - 508863600*pow(d,3) + 
						 117369000*pow(d,4))))/
		(1 + ((23180130008 - 45316115425*d + 30447267915*pow(d,2) - 9844225650*pow(d,3) + 
			   1986342750*pow(d,4) - 229635000*pow(d,5))*pow(b,2))/
		 (48.*(-624270496 + 1079133400*d + 163731525*pow(d,2) - 508863600*pow(d,3) + 
			   117369000*pow(d,4))) + ((171517025296 - 115199205478*d + 932058064525*pow(d,2) - 
			   1300303859400*pow(d,3) + 379720183500*pow(d,4) - 83587140000*pow(d,5))*pow(b,4))
		 /(16128.*(-624270496 + 1079133400*d + 163731525*pow(d,2) - 508863600*pow(d,3) + 
				 117369000*pow(d,4))) + ((12450521980528 - 29898858364050*d + 
				 23725766219155*pow(d,2) - 7556476979310*pow(d,3) + 2054993210550*pow(d,4) - 
				 825101518500*pow(d,5) + 48223350000*pow(d,6))*pow(b,6))/
		 (258048.*(-624270496 + 1079133400*d + 163731525*pow(d,2) - 508863600*pow(d,3) + 
				   117369000*pow(d,4))));
}

double strong_coupling_pl_su3(double b, unsigned d){
	double p = b / 2.;
	p += pow(b, 2) / 8.;
	p -= pow(b, 4) * 5 / 384.;
	p += pow(b, 5) * (d / 1296. - 113 / 20736.);
	p += pow(b, 6) * (d * 7 / 5184. - 931 / 414720.);
	p += pow(b, 7) * (d * 5 / 5184. - 1069 / 829440.);

	return p;
}

double strong_coupling_pl_pade_su3(double b, double d){
	// 7th order Pade approximation using 16th order expansion
	// has extended convergence radius, but gives only around 4 significant digits due to rounding errors
	return (-2.*b*(3.844516909690548e40*pow(d,12)*pow(b,6) - 
				3.754411044619676e36*pow(d,11)*pow(b,2)*
				(1.65888e6 - 2.592e6*b + 1.070784e6*pow(b,2) + 3528.*pow(b,3) + 550079.*pow(b,4)) - 
				1.466566814304561e33*pow(d,10)*
				(3.4064105472e11 - 8.8844967936e11*b + 9.05872896e11*pow(b,2) - 
				 5.734478592e11*pow(b,3) + 4.13795366784e11*pow(b,4) - 2.81405846448e11*pow(b,5) + 
				 6.6891211501e10*pow(b,6)) + 
				8.183966597681701e28*pow(d,9)*(3.9931714955575296e17 - 1.1364400773708595e18*b + 
					1.4732573319897754e18*pow(b,2) - 1.191377282758053e18*pow(b,3) + 
					6.759585093127357e17*pow(b,4) - 3.279743979409614e17*pow(b,5) + 
					7.343942379399352e16*pow(b,6)) - 
				3.22256764928e18*pow(d,8)*(9.329279697109007e28 - 4.549649271326095e29*b + 
					6.697218687472394e29*pow(b,2) - 5.761944574810583e29*pow(b,3) + 
					3.0222090603958893e29*pow(b,4) - 1.3790628454727891e29*pow(b,5) + 
					2.976555136562208e28*pow(b,6)) + 
				2.2478848e13*pow(d,7)*(-2.3325621377337923e35 - 5.785805818557892e34*b + 
					4.168104238898142e35*pow(b,2) - 5.03276328592795e35*pow(b,3) + 
					2.6537708136410914e35*pow(b,4) - 1.4080373884541189e35*pow(b,5) + 
					3.16648320125498e34*pow(b,6)) - 
				1.0976e9*pow(d,6)*(-1.0809390714410109e41 + 1.2196455633541351e41*b - 
					6.713103935384783e40*pow(b,2) + 7.377774617667407e39*pow(b,3) + 
					1.2563178005396512e39*pow(b,4) - 8.70365825295142e39*pow(b,5) + 
					2.4547301958106215e39*pow(b,6)) + 
				1.96e7*pow(d,5)*(-4.900523197982517e43 + 6.635083688486213e43*b - 
						5.044848402674635e43*pow(b,2) + 2.197666650724452e43*pow(b,3) - 
						8.355960688593772e42*pow(b,4) + 2.3585526661232162e41*pow(b,5) + 
						2.06371414287494e41*pow(b,6)) + 
				91875.*pow(d,4)*(4.5595670242247544e46 - 6.444729018552185e46*b + 
						5.239525079411576e46*pow(b,2) - 2.570466618243503e46*pow(b,3) + 
						1.0114593752780085e46*pow(b,4) - 1.4015622374385224e45*pow(b,5) + 
						7.141393601297785e43*pow(b,6)) - 
				84000.*pow(d,3)*(1.2818554432475188e47 - 1.8272193208015583e47*b + 
						1.518343485340916e47*pow(b,2) - 7.710811656393412e46*pow(b,3) + 
						3.0304259249801957e46*pow(b,4) - 5.127704938470678e45*pow(b,5) + 
						4.8359012949188694e44*pow(b,6)) - 
				6720.*d*(2.0232940687143836e48 - 2.854963556694176e48*b + 
						2.3912566997942844e48*pow(b,2) - 1.2307538759653781e48*pow(b,3) + 
						4.751197903248125e47*pow(b,4) - 8.905594208477374e46*pow(b,5) + 
						1.0152565920133657e46*pow(b,6)) + 
				1400.*pow(d,2)*(1.1678751391971075e49 - 1.6596517096165056e49*b + 
						1.3888198256146435e49*pow(b,2) - 7.132403494391353e48*pow(b,3) + 
						2.780357900275714e48*pow(b,4) - 5.049989004476572e47*pow(b,5) + 
						5.454709793005624e46*pow(b,6)) + 
				144.*(3.320749232195959e49 - 4.6472823288861413e49*b + 3.8851940246948103e49*pow(b,2) - 
						1.9961145531141983e49*pow(b,3) + 7.630479766349323e48*pow(b,4) - 
						1.4519733662236799e48*pow(b,5) + 1.6982514866819356e47*pow(b,6))))/
						(2.4028230685565924e39*pow(d,12)*pow(b,6)*(-80. + 21.*b) - 
						 1.5017644178478703e37*pow(d,11)*pow(b,2)*
						 (-1.65888e6 + 3.00672e6*b - 1.61712e6*pow(b,2) + 129996.*pow(b,3) - 
						  599965.*pow(b,4) + 146560.*pow(b,5)) - 
						 1.0475477245032577e32*pow(d,10)*(-1.907589906432e13 + 5.452215681024e13*b - 
							 6.435942137856e13*pow(b,2) + 4.770616725504e13*pow(b,3) - 
							 3.581286049152e13*pow(b,4) + 2.6136440450064e13*pow(b,5) - 
							 1.0112553889214e13*pow(b,6) + 1.761681801031e12*pow(b,7)) + 
						 5.114979123551063e27*pow(d,9)*(-2.555629757156819e19 + 7.912123934462706e19*b - 
							 1.1406877908350239e20*pow(b,2) + 1.0409981228479807e20*pow(b,3) - 
							 6.786720481049653e19*pow(b,4) + 3.6959365325259784e19*pow(b,5) - 
							 1.2450874049549892e19*pow(b,6) + 2.0270271716620598e18*pow(b,7)) - 
						 1.61128382464e18*pow(d,8)*(-7.463423757687205e29 + 3.826305011003056e30*b - 
							 6.314351202728679e30*pow(b,2) + 6.168707461161659e30*pow(b,3) - 
							 3.848347842282634e30*pow(b,4) + 1.9743660852622065e30*pow(b,5) - 
							 6.354808289918345e29*pow(b,6) + 1.0175007213206423e29*pow(b,7)) + 
						 4.4957696e13*pow(d,7)*(4.6651242754675845e35 - 9.119905155317672e32*b - 
							 8.333928501507454e35*pow(b,2) + 1.2270496308573067e36*pow(b,3) - 
							 8.161656448265873e35*pow(b,4) + 4.7208528016176926e35*pow(b,5) - 
							 1.604812349832729e35*pow(b,6) + 2.736941759804287e34*pow(b,7)) - 
						 1.0976e9*pow(d,6)*(4.3237562857640435e41 - 5.959521324857551e41*b + 
							 4.175121905368301e41*pow(b,2) - 1.2262936411069995e41*pow(b,3) + 
							 2.0226721606331807e40*pow(b,4) + 3.68656134150205e40*pow(b,5) - 
							 1.992276948607639e40*pow(b,6) + 4.452662051439754e39*pow(b,7)) + 
						 980000.*pow(d,5)*(3.920418558386014e45 - 6.288171590385475e45*b + 
								 5.6079216197360764e45*pow(b,2) - 3.0580194922222785e45*pow(b,3) + 
								 1.338341970151851e45*pow(b,4) - 2.560630744982035e44*pow(b,5) + 
								 8.328956758304963e42*pow(b,6) + 9.287470148082609e42*pow(b,7)) + 
						 6125.*pow(d,4)*(-2.7357402145348528e48 + 4.550772464765024e48*b - 
								 4.281408163838202e48*pow(b,2) + 2.5413887771521406e48*pow(b,3) - 
								 1.1643814192071336e48*pow(b,4) + 3.0020072120297123e47*pow(b,5) - 
								 4.5581693867296186e46*pow(b,6) + 3.4988219044989203e44*pow(b,7)) - 
						 7000.*pow(d,3)*(-6.152906127588091e48 + 1.0308879271744503e49*b - 
								 9.865268547572522e48*pow(b,2) + 6.007274801556028e48*pow(b,3) - 
								 2.7694409884057035e48*pow(b,4) + 7.651747351634882e47*pow(b,5) - 
								 1.3419033298468803e47*pow(b,6) + 6.761541261678656e45*pow(b,7)) + 
						 1400.*pow(d,2)*(-4.67150055678843e49 + 7.80648197766313e49*b - 
								 7.506899796874357e49*pow(b,2) + 4.608032686642098e49*pow(b,3) - 
								 2.1177667853528852e49*pow(b,4) + 5.995240359517597e48*pow(b,5) - 
								 1.1015025302740487e48*pow(b,6) + 7.048765463016341e46*pow(b,7)) - 
						 1120.*d*(-4.85590576491452e49 + 8.065888977294652e49*b - 7.755488323829946e49*pow(b,2) + 
								 4.766225503979745e49*pow(b,3) - 2.1773543231814385e49*pow(b,4) + 
								 6.219495928458045e48*pow(b,5) - 1.1641938332201265e48*pow(b,6) + 
								 8.161713823133586e46*pow(b,7)) + 
						 16.*(-1.1954697235905452e51 + 1.971889069296647e51*b - 1.8916421162142936e51*pow(b,2) + 
								 1.1603797441228477e51*pow(b,3) - 5.264702583576784e50*pow(b,4) + 
								 1.5075116150452165e50*pow(b,5) - 2.8428282032800036e49*pow(b,6) + 
								 2.0833536977474147e48*pow(b,7)));
}
