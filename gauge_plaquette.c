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
			//return 2*strong_coupling_pl_u1(2*beta, d);
			return strong_coupling_pl_pade_u1(beta, d);
		case 2:
			//return strong_coupling_pl_su2(beta, d);
			return strong_coupling_pl_pade_su2(beta, d);
		case 3:
			//return strong_coupling_pl_su3(beta/2, d)/3;
			return strong_coupling_pl_pade_su3(beta, d);
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
	return (b*(-1.0893717856218843e53 + d*(3.0974608897948494e53 + d*
					(-3.7247917720767746e53 + d*(2.4529906225845334e53 + 
												 d*(-9.54329845736323e52 + d*(2.188144864439169e52 + 
														 d*(-2.7028557199661045e51 + d*(1.1944963376208312e50 + 
																 d*(6.8490085229495e48 + d*(-7.44490843091735e47 + 1.1380894363633088e46*d))))))))) + 
				b*(1.0163606453273991e53 + d*(-2.913775805962076e53 + 
						d*(3.528834447572095e53 + d*(-2.331075048512588e53 + 
								d*(8.992662705925884e52 + d*(-1.9750985369701335e52 + 
										d*(2.033124761200076e51 + d*(1.9752596650843413e49 + 
												d*(-2.2267232946468931e49 + (1.4125267468332445e48 - 1.9788869266920226e46*d)*d)))))))) + 
					b*(-5.664612888005033e52 + d*(1.6270110585400312e53 + 
							d*(-1.968652102808757e53 + d*(1.291351134282449e53 + 
									d*(-4.873986337543253e52 + d*(1.0011501655107813e52 + 
											d*(-7.460406665471817e50 + d*(-9.486535890477637e49 + 
													d*(2.1852018285932025e49 + d*(-1.2207802404728908e48 + d*(1.3451296663422113e46 + 6.305968861119922e43*d))))))
									  )))) + b*(1.9402233456270006e52 + 
							  d*(-5.582699581378955e52 + d*(6.740121302199829e52 + 
									  d*(-4.3720302091750647e52 + d*(1.5940909387200725e52 + 
											  d*(-2.90751297890845e51 + d*(5.466045658737428e49 + 
													  d*(7.636323662393959e49 + 
														  d*(-1.2533572924017067e49 + d*(6.581379523946869e47 + (-5.6767473002478217e45 - 6.568717563666585e43*d)*d))
														))))))) + b*(-4.9445508885943614e51 + 
									d*(1.436762245942233e52 + d*(-1.7516254771736998e52 + 
											d*(1.1455009996425139e52 + d*(-4.1817523546650165e51 + 
													d*(7.3699573273397065e50 + 
														d*(6.205204880425445e48 + 
															d*(-2.6844169836003487e49 + 
																d*(4.3826729163267144e48 + 
																	d*(-2.489409837735365e47 + d*(2.7308634877727917e45 + 1.8090734742009157e43*d)))))))))) + 
									b*(6.2725249420862975e50 + d*(-1.7953677924290388e51 + 
											d*(2.1209953818801605e51 + d*(-1.2921816444946107e51 + 
													d*(3.863055916939927e50 + d*
														(-1.3868289676804512e49 + 
														 d*(-2.8659405895318434e49 + 
															 d*(9.495317530013128e48 + 
																 d*(-1.3332369936433902e48 + 
																	 d*(8.052394552930783e46 + d*(-1.2381014272557653e45 + 3.973668649625465e40*d)))))))))) + 
										b*(-4.890964281643975e49 + d*(1.3645048596659634e50 + 
												d*(-1.5273187420415746e50 + d*
													(8.1243141754637e49 + d*(-1.312231074238468e49 + 
																			 d*(-8.089759440069764e48 + 
																				 d*(5.3886237258434764e48 + 
																					 d*(-1.423577891511282e48 + 
																						 d*(1.918430057876717e47 + 
																							 d*(-1.2020515825660675e46 + 
																								 d*(1.9620086191198836e44 + (4.1304453460266934e42 - 7.689033819381096e40*d)*d)))))))))))))))))
																								 )/(-3.268115356865653e53 + d*(9.29238266938455e53 + d*
																										 (-1.1174375316230324e54 + d*(7.3589718677536e53 + 
																																	  d*(-2.862989537208969e53 + d*(6.564434593317507e52 + 
																																			  d*(-8.108567159898313e51 + d*(3.583489012862494e50 + 
																																					  d*(2.05470255688485e49 + d*(-2.2334725292752052e48 + 3.4142683090899266e46*d))))))))) + 
																									 b*(3.5937678287931393e53 + d*(-1.0290057862783652e54 + 
																											 d*(1.2448899228754673e54 + d*(-8.219720456830031e53 + 
																													 d*(3.1749637346459266e53 + d*(-7.019368043129985e52 + 
																															 d*(7.45080214358328e51 + d*(-4.670269285113279e47 + 
																																	 d*(-7.022620310088154e49 + (4.609825662045601e48 - 6.505705498257722e46*d)*d)))))))) + 
																										 b*(-2.2983451712003664e53 + d*(6.596042819417369e53 + 
																												 d*(-7.98077284655205e53 + d*(5.244006812319019e53 + 
																														 d*(-1.9913565237039637e53 + d*(4.1733451703873415e52 + 
																																 d*(-3.479922356905425e51 + d*(-2.845182388929105e50 + 
																																		 d*(7.726042204127633e49 + d*(-4.430644998426272e48 + d*(5.119673248736254e46 + 1.8917906583359764e44*d))))))))))
																											 + b*(9.399075927395066e52 + d*(-2.7024498607565155e53 + 
																													 d*(3.2659431666575873e53 + d*(-2.1288280078014176e53 + 
																															 d*(7.880290669153787e52 + d*(-1.5171599205787779e52 + 
																																	 d*(6.814023246175154e50 + d*(2.7927445417253725e50 + 
																																			 d*(-5.031891391497343e49 + d*(2.695621105924279e48 + (-2.529958401984942e46 - 2.2859137121559715e44*d)*d)))))
																															   )))) + b*(-2.8429393951314634e52 + 
																													   d*(8.230399341625837e52 + d*(-1.0006448060792383e53 + 
																															   d*(6.542804335108475e52 + d*(-2.4069947150172464e52 + 
																																	   d*(4.426566066277247e51 + d*(-7.492786751849554e49 + 
																																			   d*(-1.238386284419322e50 + 
																																				   d*(2.0927634625771608e49 + 
																																					   d*(-1.1715952582556574e48 + d*(1.2661542174097485e46 + 8.196299804441749e43*d)))))))))) + 
																													   b*(5.42704181416278e51 + d*(-1.5673129739714272e52 + 
																															   d*(1.888500713248043e52 + d*(-1.2051502078824938e52 + 
																																	   d*(4.1371411890784476e51 + d*
																																		   (-5.6461907926853865e50 + 
																																			d*(-9.104331888973462e49 + 
																																				d*(4.775369965107222e49 + 
																																					d*(-7.157844308476783e48 + 
																																						d*(4.253543596314019e47 + (-6.160312962017876e45 - 4.3925257634074145e42*d)*d))))))))) + 
																														   b*(-6.822787687872008e50 + d*(1.9558456398098125e51 + 
																																   d*(-2.313155313575502e51 + d*(1.4089984963392243e51 + 
																																		   d*(-4.187818124057837e50 + 
																																			   d*(1.2243566434708295e49 + 
																																				   d*(3.2800847681876167e49 + 
																																					   d*(-1.0822299864123822e49 + 
																																						   d*(1.5359099709350416e48 + 
																																							   d*(-9.552894125011707e46 + 
																																								   d*(1.5890074223344043e45 + (1.3515091334311461e43 - 2.883387682267911e41*d)*d)))))))))) + 
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
