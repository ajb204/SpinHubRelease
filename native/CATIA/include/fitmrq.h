struct Fitmrq {

	friend class Catia;

	/*
	 Based on the fitmrq.h found in Numerical Recepies 3.

	 Modified by Flemming on Oct 5 to adapt the catia structure.
	 - remove the pointer to a function, but rather use the
	 CalcR2drv
	 - Changed so that no loop over the points are performed,
	 until after the function call.
	 - Changed the convergence criteria
	 - added a few functions to set some parameters from the outside
	 of the structure.

	 Modified by Flemming on Oct 12 2007
	 - Change the vector arrays ai and a, from bool/double to
	 *int / *double, which points to the address of
	 LocalParamF[AtomNumber][ParamName],LocalParam[AtomNumber][ParamName],
	 .. we remember that int and bool is essential the same, but
	 we can keep some extra information in the int

	 */

	static const Int NDONE = 5;
	Int ndat, ma, mfit;
	VecDoub_I x, y, sig;
	Doub tol;
	std::vector<double*> a;
	std::vector<int*> ia;
	std::vector<int> AtomNumber;
	std::vector<std::string> ParamName;
	MatDoub alpha;
	MatDoub covar;
	Doub chisq;
	Doub RedChiSq;
	bool print;
	Catia* cpmg;
	Int ITMAX;
	bool conv; // is the fit converged ?
	//
	// Constructor
	Fitmrq(VecDoub_I xx, VecDoub_I yy, VecDoub_I ssig, const std::vector<double*> aa, const std::vector<int*> aai, const std::vector<int> aaa, const std::vector<std::string> aan,
			Catia &ccpmg, const Doub TOL = 1.e-4) :
		ndat(xx.size()), ma(aa.size()), x(xx), y(yy), sig(ssig), tol(TOL), a(aa), ia(aai), AtomNumber(aaa), ParamName(aan), alpha(ma, ma), covar(ma, ma), cpmg(&ccpmg),
				ITMAX(1000), conv(false) {
		// Do some cleaning:
		//std::cerr<<"# we are now fitting "<<a.size()<<" parameters"<<std::endl;
		//std::cerr<<"# DataPoints "<<ndat<<std::endl;
		//for (Int i=0;i<ma;i++) ia[i] = true;

	}

	//void catia::CalcR2drv(VecDoub,VecDoub&,MatDoub&,VecDoub)
	void hold(const Int i, const double val) {
		*(ia[i]) = false;
		*(a[i]) = val;
	}
	void free(const Int i) {
		*(ia[i]) = true;
	}
	void fix(const Int i) {
		*(ia[i]) = false;
	}

	void SetTol(double NewTol) {
		tol = NewTol;
	}
	void SetPrint(bool p) {
		print = p;
	}
	void SetMaxIter(int m) {
		ITMAX = m;
	}
	void FetchParam(double* aa) {
		Int j;
		for (j = 0; j < ma; j++) {
			aa[j] = *(a[j]);
		}
	}

	inline bool fetchia(int* i) {
		if (*i > 0) {
			return true;
		} else {
			return false;
		}
	}

	bool converged() {
		return conv;
	}

	void RaiseStopFittingFlag(/*int signal*/) {
		(*cpmg)._stopFitting = true;
		std::cerr << "\n I have caught a CTRL-C, fitting will be aborted after this iteration has finished" << std::endl;
		return;
	}
	//
	void fit() {
		//Start recording the ctrl-c flag
		//signal(SIGINT,(*this).RaiseStopFittingFlag);
		Int j, k, l, iter, mfitinit, done = 0;
		Doub alamda = .001, ochisq;
		double dummy;
		VecDoub atry(ma), beta(ma), da(ma);
		mfit = 0;
		for (j = 0; j < ma; j++)
			if (fetchia(ia[j]))
				mfit++;
		mfitinit = mfit;
		mrqcof(a, alpha, beta);
		//
		// allocate space (mfit might have changed!!)
		MatDoub oneda(mfit, 1), temp(mfit, mfit);
		ochisq = chisq;
		for (iter = 0; iter < ITMAX; iter++) {
			RedChiSq = chisq / (x.size() - mfit);
			if (print) {
				fprintf(stdout, "# lamda=%6g ChiSq=%13.6e ChiSq/(DF)=%13.6e Done=%d\n", alamda, chisq, RedChiSq, done);
			}
			if (done == NDONE) {
				// The Hansen criteria: change 0.-> something
				// ensures that we try all different steps around the minimum
				if (alamda > 0.) {
					alamda = 0.;
					conv = true;
				} else {
					done--;
				}
			}
			for (j = 0; j < mfit; j++) {
				for (k = 0; k < mfit; k++)
					covar[j][k] = alpha[j][k];
				covar[j][j] = alpha[j][j] * (1.0 + alamda);
				for (k = 0; k < mfit; k++)
					temp[j][k] = covar[j][k];
				oneda[j][0] = beta[j];
			}
			if (getenv("I_LOVE_CATIA") || true) {
				gaussj(temp, oneda);
			}
			for (j = 0; j < mfit; j++) {
				for (k = 0; k < mfit; k++)
					covar[j][k] = temp[j][k];
				da[j] = oneda[j][0];
			}
			if (done == NDONE || (*cpmg)._stopFitting) {
				covsrt(covar);
				covsrt(alpha);
				signal(SIGINT, SIG_DFL);
				return;
			}

			for (j = 0; j < ma; j++)
				atry[j] = *(a[j]); //store the parameters in atry

			for (j = 0, l = 0; l < ma; l++)
				if (fetchia(ia[l]))
					*(a[l]) = *(a[l]) + da[j++];

			//mrqcof should be called with 'a' (pointers to global catia)
			mrqcof(a, covar, da);
			// swap atry and a:
			for (l = 0; l < ma; l++) {
				dummy = *(a[l]);
				*(a[l]) = atry[l];
				atry[l] = dummy;
			}
			//fprintf(stderr,"chi-ochi = %15.8e \n",chisq-ochisq);
			//fprintf(stderr,"MAX(tol,tol*chisq)= %15.8e \n", MAX(tol,tol*chisq));

			if (fabs(chisq - ochisq) < MAX(tol, tol * chisq))
				done++;
			if (chisq < ochisq) {
				alamda = (0.1 * alamda > 1E-13 ? 0.1 * alamda : 1E-13);
				ochisq = chisq;
				for (j = 0; j < mfit; j++) {
					for (k = 0; k < mfit; k++)
						alpha[j][k] = covar[j][k];
					beta[j] = da[j];
				}
				for (l = 0; l < ma; l++) {
					*(a[l]) = atry[l];
				}
			} else {
				alamda *= 10.0;
				chisq = ochisq;
			}
		}
		//throw("Fitmrq too many iterations");
		covsrt(covar);
		covsrt(alpha);
		signal(SIGINT, SIG_DFL);
		return;
	}
	void mrqcof(const std::vector<double*> mra, MatDoub_O &alpha, VecDoub_O &beta) {
		Int i, j, k, l, m, mfitnew;
		Doub /*ymod,*/wt, sig2i, dy;
		//VecDoub dyda(ma);
		MatDoub dyda(ndat, ma);
		VecDoub ycalc(ndat);
		(*cpmg).CalcR2drv(x, ycalc, dyda, mra, ia, AtomNumber, ParamName);
		mfitnew = 0;
		for (j = 0; j < ma; j++)
			if (fetchia(ia[j]))
				mfitnew++;
		mfit = mfitnew;
		//
		// change the
		for (j = 0; j < mfit; j++) {
			for (k = 0; k <= j; k++)
				alpha[j][k] = 0.0;
			beta[j] = 0.;
		}
		chisq = 0.;
		for (i = 0; i < ndat; i++) {
			//funcs(x[i],a,ymod,dyda);
			sig2i = 1.0 / (sig[i] * sig[i]);
			//dy=y[i]-ymod;
			dy = y[i] - ycalc[i];
			for (j = 0, l = 0; l < ma; l++) {
				if (fetchia(ia[l])) {
					wt = dyda[i][l] * sig2i;
					for (k = 0, m = 0; m < l + 1; m++)
						if (fetchia(ia[m]))
							alpha[j][k++] += wt * dyda[i][m];
					beta[j++] += dy * wt;
				}
			}
			chisq += dy * dy * sig2i;
		}
		for (j = 1; j < mfit; j++)
			for (k = 0; k < j; k++)
				alpha[k][j] = alpha[j][k];
	}

	void covsrt(MatDoub_IO &covar) {
		Int i, j, k;
		for (i = mfit; i < ma; i++)
			for (j = 0; j < i + 1; j++)
				covar[i][j] = covar[j][i] = 0.0;
		k = mfit - 1;
		for (j = ma - 1; j >= 0; j--) {
			if (fetchia(ia[j])) {
				for (i = 0; i < ma; i++)
					SWAP(covar[i][k], covar[i][j]);
				for (i = 0; i < ma; i++)
					SWAP(covar[k][i], covar[j][i]);
				k--;
			}
		}
	}

	void gaussj(MatDoub_IO &a, MatDoub_IO &b) {
		Int i, icol, irow, j, k, l, ll, n = a.nrows(), m = b.ncols();
		Doub big, dum, pivinv;
		VecInt indxc(n), indxr(n), ipiv(n);
		for (j = 0; j < n; j++)
			ipiv[j] = 0;
		for (i = 0; i < n; i++) {
			big = 0.0;
			for (j = 0; j < n; j++)
				if (ipiv[j] != 1)
					for (k = 0; k < n; k++) {
						if (ipiv[k] == 0) {
							if (std::abs(a[j][k]) >= big) {
								big = std::abs(a[j][k]);
								irow = j;
								icol = k;
							}
						}
					}
			++(ipiv[icol]);
			if (irow != icol) {
				for (l = 0; l < n; l++)
					SWAP(a[irow][l], a[icol][l]);
				for (l = 0; l < m; l++)
					SWAP(b[irow][l], b[icol][l]);
			}
			indxr[i] = irow;
			indxc[i] = icol;
			if (a[icol][icol] == 0.0)
				throw("gaussj: Singular Matrix");
			pivinv = 1.0 / a[icol][icol];
			a[icol][icol] = 1.0;
			for (l = 0; l < n; l++)
				a[icol][l] *= pivinv;
			for (l = 0; l < m; l++)
				b[icol][l] *= pivinv;
			for (ll = 0; ll < n; ll++)
				if (ll != icol) {
					dum = a[ll][icol];
					a[ll][icol] = 0.0;
					for (l = 0; l < n; l++)
						a[ll][l] -= a[icol][l] * dum;
					for (l = 0; l < m; l++)
						b[ll][l] -= b[icol][l] * dum;
				}
		}
		for (l = n - 1; l >= 0; l--) {
			if (indxr[l] != indxc[l])
				for (k = 0; k < n; k++)
					SWAP(a[k][indxr[l]], a[k][indxc[l]]);
		}
	}

	void gaussj(MatDoub_IO &a) {
		MatDoub b(a.nrows(), 0);
		gaussj(a, b);
	}

};

