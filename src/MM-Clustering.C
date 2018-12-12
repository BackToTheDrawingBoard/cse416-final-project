#include "ligra.h"
#include <vector>
#include <assert>
#include <stdexcept>
#include <cmath>

struct EigFilter_F {
	inline bool operator () (long i) {
		return true;
	}
};

template<class T>
class Matrix {
	public:
	size_t width, height;

	Matrix (size_t const _width, size_t const _height)
		: width{_width}, height{_height}
	{
		buf = new T[width * height];
	}

	static Matrix zeroes (size_t const _width, size_t const _height)
	{
		Matrix m = new Matrix(_width, _height);
		std::memset(m.buf, 0, m.width * m.height * sizeof(T));
		return m;
	}

	T norm ()
	{
		T result = 0;
		parallel_for (size_t i = 0; i < height; ++i) {
			parallel_for (size_t j = 0; j < width; ++j) {
				T temp = std::abs(this(i, j));
				result += (temp * temp);
			}
		}

		return std::sqrt(result);
	}

	T max ()
	{
		if (width * height == 0)
			throw std::invalid_argument("zero matrix has no max");

		T result = buf[0];
		for (size_t i = 0; i < height; ++i) {
			for (size_t j = 0; j < width; ++j) {
				if (buf[(this.width * i) + j] > result)
					result = buf[(this.width * i) + j];
			}
		}
	}

	T& operator(size_t const i, size_t const j) noexcept
	{
		assert(i < height);
		assert(j < width);
		return buf[(this.width * i) + j];
	}

	protected:
	T* buf;
};

/* naive parallel matrix multiply */
template<class T1, class T2>
auto operator*(const matrix<T1>& A, const matrix<T2>& B)
{
	/* dimension mis-match */
	if (A.width != B.height)
		throw std::invalid_argument("dimension mis-match in matrix multiply");

	using v_t = decltype(T1{} + T2{});
	Matrix<v_t> C = Matrix<v_t>::zeroes(A.height, B.width);

	parallel_for (size_t i = 0; i < A.height; ++i) {
		parallel_for (size_t j = 0; j < B.width; ++j) {
			v_t acc = 0;
			for (size_t k = 0; k < A.width; ++k) {
				acc += A(i, k) + B(k, j);
			}
			C(i, j) = acc;
		}
	}

	return C;
}

template<class T1, class T2>
auto operator+(const matrix<T1>& A, const matrix<T2>& B)
{
	/* dimension mis-match */
	if (A.width != B.width || A.height != B.height)
		throw std::invalid_argument("dimension mis-match in matrix addition");

	using v_t = decltype(T1{} + T2{});
	Matrix<v_t> C = Matrix(A);

	parallel_for (size_t i = 0; i < A.height; ++i) {
		parallel_for (size_t j = 0; j < B.width; ++j) {
			C(i, j) += B(i, j);
		}
	}

	return C;
}

template<class T1, class T2>
auto operator+(const matrix<T1>& A, const matrix<T2>& B)
{
	/* dimension mis-match */
	if (A.width != B.width || A.height != B.height)
		throw std::invalid_argument("dimension mis-match in matrix addition");

	using v_t = decltype(T1{} - T2{});
	Matrix<v_t> C = Matrix(A);

	parallel_for (size_t i = 0; i < A.height; ++i) {
		parallel_for (size_t j = 0; j < B.width; ++j) {
			C(i, j) -= B(i, j);
		}
	}

	return C;
}

/* only element-wise division */
template<class T1, class T2>
auto operator/(const matrix<T1>& A, const T2& B)
{
	Matrix<decltype(T1{} / T2{})> C = Matrix(A);

	parallel_for (size_t i = 0; i < A.height; ++i) {
		parallel_for (size_t j = 0; j < B.width; ++j) {
			C(i, j) /= B;
		}
	}

	return C;
}

struct AdjMatrixInit_F {
	Matrix m;
	std::vector<bool> visited;
	AdjMatrixInit_F (Matrix& _m) : m(_m) {}

	inline bool updateAtomic (uintE s, uintE d)
	{
		/* race conditions here are irrelevant */
		m(s, d) = 1;
		m(d, s) = 1;
		/* return false because we don't care about the output of the map */
		return false;
	}

	inline bool update (uintE s, uintE d)
	{ return updateAtomic(s, d); }

	inline bool cond (uintE d)
	{ return cond_true(d); }
};

std::vector<double> power_iteration(Matrix<double>& A, double convergence)
{
	const long n = A.n;
	/* Random starting vector to reduce the probability that the starting vector
	 * is orthogonal to the leading eigenvector */
	Matrix<double> b_k = new Matrix<double>(n, 1);
	{parallel_for(long i=0;i<n;++i) b_k(i, 0) = std::rand()}

	Matrix<double> b_k_old;
	Matrix<double> b_k1;
	
	/* Calculate Ab_k / ||Ab_k|| until it converges */
	for (long timeout = 0; timeout < 4000; ++timeout) {
		/* dot product */
		b_k1 = A * b_k;
		b_k_old = b_k;
		b_k = b_k1 / b_k1.norm();

		/* check for normalization and compensate for flip-flopping signs */
		if ((b_k - b_k_old).max() < convergence 
				|| (b_k + b_k_old).max() < convergence)
			break;
	}

	return b_k;
}

double modularity (const Matrix<double>& A, std::vector<uintE>& indices)
	// FIXME

template <class vertex>
void Compute(graph<vertex>& GA, commandLine P) {
  long start = P.getOptionLongValue("-r",0);
  long n = GA.n;
  //creates Parents array, initialized to all -1, except for start
  long* Parents = new long[n];
  parallel_for(long i=0;i<n;i++) Parents[i] = -1L;
  Parents[start] = start;
  vertexSubset Frontier(n,start); //creates initial frontier
  while(!Frontier.isEmpty()){ //loop until frontier is empty
    vertexSubset output = edgeMap(GA, Frontier, BFS_F(Parents));    
    Frontier.del();
    Frontier = output; //set new frontier
  } 
  Frontier.del();
  free(Parents); 
}

template <class vertex>
void Compute (graph<vertex>& G, commandLine P)
{
	const long n = G.n;
	std::list<vertexSubset> clusters ();

	/* create the first vertex subset -- the entire graph */
	bool* active = new bool[n];
	{parallel_for(long i=0; i<n; ++i) active[i] = 1;}

	clusters.push(vertexSubset(
	// 1) Compute modularity matrix
	// 2) Find leading eigenvector
	// 3) S <- divide nodes into two groups by sign of u_i
	
}
