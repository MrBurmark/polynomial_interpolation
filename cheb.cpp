
#include <iostream>
#include <cmath>
#include <mgl2/mgl.h>

typedef enum {
	evenly_spaced,
	circle_even,
	circle_cheb,
} PointChooser;

const double pi = std::acos(-1);

// Class representing polynomials with all zeros == 0 on Domain -1 to 1
class ZeroPoly
{
public:
	const std::size_t order = 0;
	ZeroPoly(const std::size_t param_order)
	:	order(param_order)
	{
	}
	double operator()(const double x) const
	{
		double y = 1.0;
		for(std::size_t i = 0; i < order; ++i)
		{
			y *= x;
		}
		return y;
	}
};

// Class representing Chebyshev polynomials on Domain -1 to 1
class Cheb
{
public:
	const std::size_t order = 0;
	Cheb(const std::size_t param_order)
	:	order(param_order)
	,	mag(param_order > 1 ? std::pow(2.0, param_order - 1) : 1.0)
	,	zeros(new double[param_order])
	{
		for (std::size_t i = 1; i <= order; ++i)
		{
			zeros[i-1] = std::cos(pi*(2*i-1)/(double)(2*order));
		}
	}
	double operator()(const double x) const
	{
		double y = mag;
		for(std::size_t i = 0; i < order; ++i)
		{
			y *= (x - zeros[i]);
		}
		return y;
	}
	~Cheb()
	{
		delete zeros;
	}
private:
	const double mag = 1.0;
	double * const zeros = NULL;
};

inline std::size_t index(std::size_t i, std::size_t j, std::size_t n)
{
	return i * n + j;
}

void printmat(double * A, std::size_t n)
{
	std::cout << "mat (" << n << ", " << n << ")\n";
	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t iind = index(i, 0, n);
		for (std::size_t j = 0; j < n; ++j)
		{
			if (j != 0) std::cout << ",  ";
			std::cout << A[iind];
			++iind;
		}
		std::cout << "\n";
	}
}

void printvec(double * a, std::size_t n)
{
	std::cout << "vec (" << n << ")\n";
	for (std::size_t i = 0; i < n; ++i)
	{
		if (i != 0) std::cout << ",  ";
		std::cout << a[i];
	}
	std::cout << "\n";
}


void gaussianElim(double * A, double * y, std::size_t n)
{
	// printmat(A, n);
	// printvec(y, n);
	// forward elimination
	for (std::size_t i = 0; i < n; ++i)
	{
		std::size_t imax = i;
		// find pivot
		for (std::size_t ii = i + 1; ii < n; ++ii)
		{
			imax = A[index(ii, i, n)] > A[index(imax, i, n)] ? ii : imax;
		}

		// printmat(A, n);
		// printvec(y, n);

		// swap
		if (i != imax)
		{
			std::size_t iind = index(i, i, n);
			std::size_t imaxind = index(imax, i, n);
			for (std::size_t j = i; j < n; ++j)
			{
				double tmp = A[iind];
				A[iind] = A[imaxind];
				A[imaxind] = tmp;
				++iind; ++imaxind;
			}

			double tmp = y[i];
			y[i] = y[imax];
			y[imax] = tmp;
		}

		// printmat(A, n);
		// printvec(y, n);

		// eliminate
		std::size_t iind = index(i, i, n);
		const double aiinv = 1.0 / A[iind];
		A[iind++] = 1.0;
		for (std::size_t j = i + 1; j < n; ++j)
		{
			A[iind] *= aiinv;
			++iind;
		}
		const double yi = (y[i] *= aiinv);

		// printmat(A, n);
		// printvec(y, n);


		for (std::size_t ii = i + 1; ii < n; ++ii)
		{

			std::size_t iind = index(i, i + 1, n);
			std::size_t iiind = index(ii, i, n);
			const double aiiinv = 1.0 / A[iiind];
			A[iiind++] = 0.0;
			for (std::size_t j = i + 1; j < n; ++j)
			{
				A[iiind] = A[iind] - A[iiind] * aiiinv;
				++iiind; ++iind;
			}
			y[ii] = yi - y[ii] * aiiinv;
		}

		// printmat(A, n);
		// printvec(y, n);
	}

	// backwards elimination
	std::size_t i = n - 1;
	do
	{
		std::size_t iind = index(i, i + 1, n);
		for (std::size_t j = i + 1; j < n; ++j)
		{
			y[i] -= y[j] * A[iind];
			++iind;
		}
	}
	while (i-- > 0);

	// printmat(A, n);
	// printvec(y, n);
}

template <class P>
class Poly
{
public:
	const std::size_t order = 0;
	// takes order + 1 (x,y) pairs
	Poly(const std::size_t param_order, const double * const x, const double * const y)
	:	order(param_order)
	,	center(x[0])
	,	halfwidthinv(1.0)
	,	coefficients(new double[param_order + 1])
	,	polys(new P*[param_order + 1])
	{

		double minx = x[0], maxx = x[0];
		for (std::size_t i = 1; i < order + 1; ++i)
		{
			minx = std::min(minx, x[i]);
			maxx = std::max(maxx, x[i]);
		}

		center = order > 0 ? (maxx + minx) / 2.0 : x[0];
		halfwidthinv = order > 0 ? 2.0 / (maxx - minx) : 1.0;

		setPolyCoef(x, y);
	}
	double operator()(double x) const
	{
		double y = 0.0;
		x = (x - center) * halfwidthinv;
		for (std::size_t i = 0; i < order + 1; ++i)
		{
			y += coefficients[i] * polys[i][0](x);
		}
		return y;
	}
	~Poly()
	{
		delete coefficients;

		for (std::size_t i = 0; i < order + 1; ++i)
		{
			delete polys[i];
		}

		delete polys;
	}
private:
	double center;
	double halfwidthinv;
	double * const coefficients = NULL;
	P ** const polys = NULL;

	void setPolyCoef(const double * const x, const double * const y) 
	{
		// make Chebychev polynomials and copy 
		for (std::size_t i = 0; i < order + 1; ++i)
		{
			coefficients[i] = y[i];
			polys[i] = new P(i);
		}

		double * A = new double[(order + 1) * (order + 1)];

		for (std::size_t i = 0; i < order + 1; ++i)
		{
			std::size_t iind = index(i, 0, order + 1);
			double xx = (x[i] - center) * halfwidthinv;

			for (std::size_t j = 0; j < order + 1; ++j)
			{
				A[iind] = polys[j][0](xx);
				// std::cout << xx << ", " << A[iind] << "\n";
				++iind;
			}
		}

		gaussianElim(A, coefficients, order + 1);

		delete A; A = NULL;
	}
};

template <class P, PointChooser pct>
void mgls_prepare1d(mglData *x, mglData *y, double (*f) (double),
	const std::size_t n=101, const std::size_t maxorder = 3,
	const double viewlbound=-1.0, const double viewrbound=1.0,
	const double interplbound=-1.0, const double interprbound=1.0)
{
	if(x && y) {

		x->Create(n, maxorder + 2);
		y->Create(n, maxorder + 2); // creates plots

		for(std::size_t i = 0; i < n; i++)
		{
			double xx = (viewlbound * (n - 1 - i) + viewrbound * i) / (double)(n - 1);
			x->a[i+(maxorder + 1)*n] = xx;
			y->a[i+(maxorder + 1)*n] = f(xx);
		}

		for (std::size_t k = 0; k < maxorder + 1; k++)
		{
			double px[k+1];
			double py[k+1];

			for(std::size_t i = 0; i < k + 1; i++)
			{
				switch (pct)
				{
					case PointChooser::evenly_spaced:
					{
						// choose points evenly spaced
						px[i] = k > 0 ? (interplbound * (k - i) + interprbound * i) / (double)k
									  : (interplbound + interprbound) / 2.0;
						break;
					}
					case PointChooser::circle_even:
					{
						// choose points closer to edges by projecting points evenly distributed on a circle
						px[i] = k > 0 ?
								std::cos(pi*(i)/(double)(k))
								* (interprbound - interplbound) / 2.0
								+ (interplbound + interprbound) / 2.0
								: (interplbound + interprbound) / 2.0;
						break;
					}
					case PointChooser::circle_cheb:
					{
						// choose points as zeros of next higher order chebyshev polynomial
						// doesn't interpolate end points
						px[i] = k > 0 ?
								std::cos(pi*(2*i+1)/(double)(2*(k+1)))
								* (interprbound - interplbound) / 2.0
								+ (interplbound + interprbound) / 2.0
								: (interplbound + interprbound) / 2.0;
						break;
					}
				}

				py[i] = f(px[i]);
			}

			Poly<P> p(k, px, py);

			double maxerr = -1.0;
			double errx = 1.0/0.0;

			for(std::size_t i = 0; i < n; i++)
			{
				double xx = (viewlbound * (n - 1 - i) + viewrbound * i) / (double)(n - 1);
				x->a[i+k*n] = xx;
				y->a[i+k*n] = p(xx);

				if (interplbound <= xx && xx <= interprbound)
				{
					double curerr = std::abs(y->a[i+k*n] - y->a[i+(maxorder + 1)*n]);
					if (curerr > maxerr)
					{
						maxerr = curerr;
						errx = xx;
					}
				}
				// std::cout << x->a[i+k*n] << ", " << y->a[i+k*n] << "\n";
			}
			std::cout << "Max err in P" << k << " ( " << errx << ", " << maxerr << " )\n";
		}
	}
}

double f(double x)
{
	return std::tan(x*x);
}

int main (int argc, char ** argv) {

	mglGraph * gr = new mglGraph(0, 1000, 1000);

	const std::size_t n = 1001;

	const std::size_t maxorder = 45;

	const double viewlbound = -std::acos(-1.0);
	const double viewrbound = std::acos(-1.0);

	const double viewlowbound = -std::acos(-1.0);
	const double viewuprbound = std::acos(-1.0);

	const double interplbound = -1.0;
	const double interprbound = 1.0;

	mglData x;
	mglData y;
	mgls_prepare1d<Cheb, PointChooser::circle_even>
		(&x, &y, f, n, maxorder, viewlbound, viewrbound, interplbound, interprbound);
	gr->SetRanges(viewlbound, viewrbound, viewlowbound, viewuprbound);
	gr->SetOrigin((viewlbound + viewrbound) / 2.0 , 0, 0);
	gr->Title("Plot cheb");
	gr->Box();
	gr->Plot(x, y);

    gr->WriteFrame("sample.png");

    delete gr;

	return 0;
}
