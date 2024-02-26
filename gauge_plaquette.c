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
				plaquettes += creal(plaquette_tr(u, nnt, ns, nn, i, mu, nu, mode));
			}
		}
	}

	return plaquettes * 2/nd/(nd-1) / ns;
}

double complex plaquette_tr(double complex *u, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, unsigned nu, gauge_flags *mode){
	const unsigned n = mode->gauge_dim, mat_dim = n*n, links = nn/2;
	const unsigned *nnl = nnt + ns*nn;
	double complex *z = mode->zdummy;

	const int turnM = mu >= links;

	copy_mat(u + nnl[pos*nn + mu], z + 3*mat_dim, n, turnM);
	copy_staple(u, z, nnt, ns, nn, pos, mu, nu, mode);

	return trace_prod(z, 4, n);
}

void plaquette_mat(double complex *u, double complex *z, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, unsigned nu, int direction, gauge_flags *mode){
	const unsigned n = mode->gauge_dim, mat_dim = n*n, links = nn/2;
	const unsigned *nnl = nnt + ns*nn;

	const int turnM = mu >= links;

	copy_mat(u + nnl[pos*nn + mu], z + 3*mat_dim, n, turnM);
	copy_staple(u, z, nnt, ns, nn, pos, mu, nu, mode);

	mat_prod(z, 4, n);

	if(direction == 1) dagger(z + 4*mat_dim, n);
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
	const unsigned n = mode->gauge_dim;

	copy_staple(u, z, nnt, ns, nn, pos, mu, nu, mode);

	mat_prod(z, 3, n);
}

void copy_staple(double complex *u, double complex *z, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, unsigned nu, gauge_flags *mode){
	const unsigned n = mode->gauge_dim, mat_dim = n*n, links = nn/2;
	const unsigned *nnl = nnt + ns*nn;

	const int turnM = mu >= links;
	const int turnN = nu >= links;

	copy_mat(u + nnl[nnt[pos*nn + mu]*nn + nu], z + 2*mat_dim, n, turnN);
	copy_mat(u + nnl[nnt[pos*nn + nu]*nn + mu], z + 1*mat_dim, n, 1 - turnM);
	copy_mat(u + nnl[pos*nn + nu], z, n, 1 - turnN);
}

double topo_charge(double complex *u, unsigned *nnt, unsigned ns, unsigned nn, gauge_flags *mode){
	const unsigned nd = mode->space_dim;

	switch(nd){
		case 2:
			return topo_charge_2d(u, nnt, ns, nn, mode);
		case 4:
			return topo_charge_4d(u, nnt, ns, nn, mode);
		default:
			return 0;
	}
}

double topo_charge_2d(double complex *u, unsigned *nnt, unsigned ns, unsigned nn, gauge_flags *mode){
	const unsigned links = nn/2;
	double charge = 0;

	for(unsigned i = 0; i < ns; i++){
		for(unsigned nu = 1; nu < links; nu++){
			for(unsigned mu = 0; mu < nu; mu++){
				charge += cimag(plaquette_tr(u, nnt, ns, nn, i, mu, nu, mode));
			}
		}
	}

	return charge / (2*M_PI);
}

double topo_charge_4d(double complex *u, unsigned *nnt, unsigned ns, unsigned nn, gauge_flags *mode){
	const unsigned links = nn/2;
	double charge = 0;

	for(unsigned i = 0; i < ns; i++){
		const unsigned mu = 0; // 8 equivalent permutations can be summarised with same mu
		for(unsigned nu = 1; nu < links; nu++){
			charge += creal(clover_field_tr(u, nnt, ns, nn, i, mu, nu, mode));
		}
	}

	return charge / (8*M_PI*M_PI);
}

double complex clover_field_tr(double complex *u, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, unsigned nu, gauge_flags *mode){
	// only works for mu = 0, nu > 0 in 4D
	unsigned rho = 0, sigma = 0;
	switch(nu){
		case 1:
			rho = 2; sigma = 3;
			break;
		case 2:
			rho = 3; sigma = 1;
			break;
		case 3:
			rho = 1; sigma = 2;
			break;
	}

	const unsigned n = mode->gauge_dim, mat_dim = n*n;
	double complex *cl = mode->zdummy;
	double complex *z = cl + 2*mat_dim;

	clover_mat(u, cl + mat_dim, z, nnt, ns, nn, pos, mu, nu, mode);
	clover_mat(u, cl, z, nnt, ns, nn, pos, rho, sigma, mode);

	return trace_prod(cl, 2, n);
}

void clover_mat(double complex *u, double complex *cl, double complex *z, unsigned *nnt, unsigned ns, unsigned nn, unsigned pos, unsigned mu, unsigned nu, gauge_flags *mode){
	// calculates the clover matrix, i.e. the field strength tensor F_{mu,nu}(n)
	const unsigned n = mode->gauge_dim, mat_dim = n*n, links = nn/2;
	double complex *prod = z + 4*mat_dim;

	for(unsigned i = 0; i < mat_dim; i++) cl[i] = 0;
	
	for(unsigned dir1 = 0; dir1 < 2; dir1++){
		for(unsigned dir2 = 0; dir2 < 2; dir2++){

			plaquette_mat(u, z, nnt, ns, nn, pos, mu + dir1*links, nu + dir2*links, (dir1 + dir2) % 2, mode);

			for(unsigned i = 0; i < mat_dim; i++) cl[i] += prod[i];
		}
	}
	
	dagger_asym(cl, n);
}

double strong_coupling_plaquette(double beta, gauge_flags *mode){
	const unsigned d = mode->space_dim;
	switch(mode->gauge_group){
		case 1:
			//return 2*strong_coupling_pl_u1(2*beta, d);
			return 2*strong_coupling_pl_pade_u1(beta, d);
		case 2:
			//return strong_coupling_pl_su2(beta, d);
			return strong_coupling_pl_pade_su2(beta, d);
		case 3:
			//return strong_coupling_pl_su3(beta/2, d)/3;
			return strong_coupling_pl_pade_su3(beta, d) / 3;
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

double strong_coupling_pl_pade_u1(double b, long d){
	// 7th order Pade approximation in Horner scheme using 16th order expansion
	// has extended convergence radius and is numerically very stable
	// might lead to integer overflow for d > 4, not tested
	return (b*(-5429737867760640 + d*(9365511549542400 + d*(-6063635850854400 + (1754420871168000 - 192631799808000*d)*d)) + 
				pow(b,2)*(74733966666299520 + d*(-146040561261020160 + 
						d*(113083847182540800 + d*(-43489159328563200 + (8342027108352000 - 642105999360000*d)*d))) + 
					pow(b,2)*(8072228388482400 + d*(-16107100244775360 + 
							d*(12563781562252800 + d*(-4776437669068800 + (884567900160000 - 64210599936000*d)*d))) + 
						pow(b,2)*(85690276124445601 + d*(-206530276810616400 + 
								d*(205471846249588800 + d*(-108095412123763200 + d*(31762666465689600 + d*(-4953245810688000 + 321052999680000*d))))))))))/
		(-21718951471042560 + d*(37462046198169600 + d*(-24254543403417600 + (7017683484672000 - 770527199232000*d)*d)) + 
		 pow(b,2)*(296220997731317760 + d*(-579479489269309440 + 
				 d*(449303570804736000 + d*(-173079426878668800 + (33271792533504000 - 2568423997440000*d)*d))) + 
			 pow(b,2)*(64339278558230400 + d*(-125563415950137600 + 
					 d*(96176983961548800 + d*(-36100641978777600 + (6643456081920000 - 481579499520000*d)*d))) + 
				 pow(b,2)*(424039461149080384 + d*(-1023587780952155520 + 
						 d*(1019969388469601280 + d*(-537482320792012800 + d*(158207895207936000 + d*(-24717067812864000 + 1605264998400000*d)))))))));
}

double strong_coupling_pl_su2(double b, unsigned d){
	double p = b / 4.;
	p -= pow(b, 3) / 96.;
	p += pow(b, 5) * (d / 512. - 5 / 1536.);
	p -= pow(b, 7) * (d / 1536. - 29 / 23040.);

	return p;
}

double strong_coupling_pl_pade_su2(double b, long d){
	// 7th order Pade approximation in Horner scheme using 16th order expansion
	// has extended convergence radius and is numerically very stable
	// might lead to integer overflow for d > 4, not tested
	return (b*(-2416376294277120 + d*(4177023234048000 + d*(633758888448000 + d*(-1969668513792000 + 454302535680000*d))) + 
				pow(b,2)*(1969928029440000 + d*(-3828334182624000 + 
						d*(2428861064313600 + d*(-711768835008000 + (141249407040000 - 18517766400000*d)*d))) + 
					pow(b,2)*(-5257917757440 + d*(51348075354240 + 
							d*(145772040447600 + d*(-248398361316000 + (63155299536000 - 15740101440000*d)*d))) + 
						pow(b,2)*(148537739520512 + d*(-357822682005800 + 
								d*(278360263900775 + d*(-82393047985500 + d*(26156806436250 + d*(-11231218215000 + 578680200000*d))))))))))/
		(-9665505177108480 + d*(16708092936192000 + d*(2535035553792000 + d*(-7878674055168000 + 1817210142720000*d))) + 
		 pow(b,2)*(7476982735380480 + d*(-14617166191488000 + 
				 d*(9821070738662400 + d*(-3175353425664000 + (640714717440000 - 74071065600000*d)*d))) + 
			 pow(b,2)*(164656344284160 + d*(-110591237258880 + 
					 d*(894775741944000 + d*(-1248291705024000 + (364531376160000 - 80243654400000*d)*d))) + 
				 pow(b,2)*(747031318831680 + d*(-1793931501843000 + 
						 d*(1423545973149300 + d*(-453388618758600 + d*(123299592633000 + d*(-49506091110000 + 2893401000000*d)))))))));
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
	// 7th order Pade approximation in Horner scheme using 16th order expansion
	// has extended convergence radius and is numerically acceptable, typically 5 - 10 correct digits
	return (b*(-6.971979427980059e54 + d*(1.9823749694687036e55 + d*(-2.3838667341291357e55 + 
						d*(1.5699139984541014e55 + d*(-6.107711012712467e54 + 
								d*(1.400412713241068e54 + d*(-1.7298276607783069e53 + 
										d*(7.64477656077332e51 + d*(4.38336545468768e50 + d*(-4.764741395787104e49 + 7.283772392725176e47*d))))))))) + 
				b*(3.252354065047677e54 + d*(-9.324082579078643e54 + 
						d*(1.1292270232230704e55 + d*(-7.459440155240282e54 + 
								d*(2.877652065896283e54 + d*(-6.320315318304427e53 + 
										d*(6.505999235840243e52 + d*(6.320830928269892e50 + 
												d*(-7.125514542870058e50 + (4.5200855898663824e49 - 6.3324381654144724e47*d)*d)))))))) + 
					b*(-9.063380620808054e53 + d*(2.60321769366405e54 + 
							d*(-3.149843364494011e54 + d*(2.0661618148519184e54 + 
									d*(-7.798378140069205e53 + d*(1.60184026481725e53 + 
											d*(-1.1936650664754908e52 + d*
												(-1.517845742476422e51 + d*
												 (3.496322925749124e50 + d*(-1.9532483847566252e49 + d*(2.152207466147538e47 + 1.0089550177791875e45*d))))))))
							  )) + b*(1.5521786765016005e53 + d*(-4.466159665103164e53 + 
								  d*(5.392097041759863e53 + d*(-3.4976241673400517e53 + 
										  d*(1.275272750976058e53 + d*(-2.32601038312676e52 + 
												  d*(4.372836526989942e50 + d*
													  (6.109058929915167e50 + 
													   d*(-1.0026858339213654e50 + d*(5.265103619157495e48 + (-4.5413978401982573e46 - 5.254974050933268e44*d)*d))
													  ))))))) + b*(-1.9778203554377446e52 + 
								  d*(5.747048983768932e52 + d*(-7.006501908694799e52 + 
										  d*(4.5820039985700556e52 + d*(-1.6727009418660066e52 + 
												  d*(2.9479829309358826e51 + 
													  d*(2.482081952170178e49 + 
														  d*(-1.0737667934401395e50 + 
															  d*(1.7530691665306857e49 + 
																  d*(-9.95763935094146e47 + d*(1.0923453951091167e46 + 7.236293896803663e43*d)))))))))) + 
								  b*(1.2545049884172595e51 + d*(-3.5907355848580776e51 + 
										  d*(4.241990763760321e51 + d*(-2.5843632889892215e51 + 
												  d*(7.726111833879854e50 + d*
													  (-2.7736579353609023e49 + 
													   d*(-5.731881179063687e49 + 
														   d*(1.8990635060026257e49 + 
															   d*(-2.6664739872867804e48 + 
																   d*(1.6104789105861567e47 + d*(-2.4762028545115307e45 + 7.94733729925093e40*d)))))))))) + 
									  b*(-4.890964281643975e49 + d*(1.3645048596659634e50 + 
											  d*(-1.5273187420415746e50 + d*
												  (8.1243141754637e49 + d*(-1.312231074238468e49 + 
																		   d*(-8.089759440069764e48 + 
																			   d*(5.3886237258434764e48 + 
																				   d*(-1.423577891511282e48 + 
																					   d*(1.918430057876717e47 + 
																						   d*(-1.2020515825660675e46 + 
																							   d*(1.9620086191198836e44 + (4.1304453460266934e42 - 7.689033819381096e40*d)*d)))))))))))))))))
																							   )/(-4.1831876567880356e55 + d*(1.1894249816812223e56 + 
																									   d*(-1.4303200404774814e56 + d*(9.419483990724608e55 + 
																											   d*(-3.6646266076274805e55 + d*(8.402476279446409e54 + 
																													   d*(-1.0378965964669841e54 + d*(4.586865936463992e52 + 
																															   d*(2.630019272812608e51 + d*(-2.8588448374722627e50 + 4.370263435635106e48*d))))))))) + 
																								   b*(2.3000114104276092e55 + d*(-6.585637032181537e55 + 
																										   d*(7.967295506402991e55 + d*(-5.26062109237122e55 + 
																												   d*(2.031976790173393e55 + d*(-4.4923955476031905e54 + 
																														   d*(4.768513371893299e53 + d*(-2.9889723424724984e49 + 
																																   d*(-4.494476998456419e51 + (2.950288423709185e50 - 4.163651518884942e48*d)*d)))))))) + 
																									   b*(-7.354704547841173e54 + d*(2.110733702213558e55 + 
																											   d*(-2.553847310896656e55 + d*(1.6780821799420862e55 + 
																													   d*(-6.372340875852684e54 + d*(1.3354704545239493e54 + 
																															   d*(-1.113575154209736e53 + d*(-9.104583644573137e51 + 
																																	   d*(2.4723335053208427e51 + d*(-1.417806399496407e50 + d*(1.6382954395956014e48 + 6.053730106675124e45*d)))))))))
																											   ) + b*(1.5038521483832105e54 + d*(-4.323919777210425e54 + 
																													   d*(5.2255090666521396e54 + d*(-3.406124812482268e54 + 
																															   d*(1.2608465070646059e54 + d*(-2.4274558729260446e53 + 
																																	   d*(1.0902437193880246e52 + d*
																																		   (4.468391266760596e51 + d*
																																			(-8.051026226395748e50 + d*(4.312993769478846e49 + (-4.0479334431759076e47 - 3.6574619394495545e45*d)*d)))))
																																 )))) + b*(-2.2743515161051707e53 + 
																														 d*(6.58431947330067e53 + d*(-8.005158448633906e53 + 
																																 d*(5.23424346808678e53 + d*(-1.925595772013797e53 + 
																																		 d*(3.5412528530217974e52 + d*
																																			 (-5.9942294014796436e50 + 
																																			  d*(-9.907090275354576e50 + 
																																				  d*(1.6742107700617287e50 + d*(-9.372762066045259e48 + d*(1.0129233739277988e47 + 6.557039843553399e44*d)))
																																				))))))) + b*(2.170816725665112e52 + 
																															d*(-6.269251895885709e52 + d*(7.554002852992172e52 + 
																																	d*(-4.820600831529975e52 + d*(1.654856475631379e52 + 
																																			d*(-2.2584763170741546e51 + 
																																				d*(-3.641732755589385e50 + 
																																					d*(1.9101479860428887e50 + 
																																						d*(-2.863137723390713e49 + 
																																							d*(1.7014174385256077e48 + (-2.4641251848071503e46 - 1.7570103053629658e43*d)*d))))))))) + 
																															b*(-1.3645575375744017e51 + d*(3.911691279619625e51 + 
																																	d*(-4.626310627151004e51 + d*(2.8179969926784485e51 + 
																																			d*(-8.375636248115674e50 + 
																																				d*(2.448713286941659e49 + 
																																					d*(6.560169536375233e49 + 
																																						d*(-2.1644599728247644e49 + 
																																							d*(3.071819941870083e48 + 
																																								d*(-1.9105788250023414e47 + 
																																									d*(3.1780148446688085e45 + (2.7030182668622923e43 - 5.766775364535822e41*d)*d)))))))))) + 
																																b*(3.3333659163958636e49 + d*(-9.141119481909617e49 + 
																																		d*(9.868271648222878e49 + d*(-4.733078883175059e49 + 
																																				d*(2.1430284165055886e48 + 
																																					d*(9.101720745120957e48 + 
																																						d*(-4.887241867660274e48 + 
																																							d*(1.2304659560698615e48 + 
																																								d*(-1.6394824538234835e47 + 
																																									d*(1.0368201665922193e46 + 
																																										d*(-1.845445761968825e44 + d*(-2.2009859307978387e42 + 5.045928443968844e40*d)))))))))))))))))
																																										));
}
