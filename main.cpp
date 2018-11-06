#include <iostream>
#include <bitset>
#include <array>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <assert.h>

typedef uint64_t u64;
typedef int64_t i64;


//x^n mod q
u64 exp_mod(u64 x, u64 n, u64 q) {
	u64 result = 1;
	for (u64 i = 0; i < n; i ++) {
		result *= x;
		result %= q;
	}
	return static_cast<u64>(result);
}

// n mod q
// except it handles negative correctly
u64 pos_mod(i64 n, u64 q) {
	while (n < 0) {
		n += q;
	}
	return (n % q);
}

template<int Nm>
std::vector<u64> encrypt(u64 p, u64 q, u64 x_0, std::bitset<Nm> m) {
	u64 n = p * q;
	u64 k = (u64)floor(log2(n));
	u64 h = (u64)floor(log2(k));
	u64 t = (u64)ceil((float)Nm / (float)h);

	//Lowest h bits
	u64 h_mask = ~0ull >> (64 - h);

	std::vector<u64> c;
	std::vector<u64> x;
	x.push_back(x_0);
	for (i64 i = 1; i <= t; i ++) {
		u64 x_i = (x[i - 1] * x[i - 1]) % n;
		x.push_back(x_i);

		u64 p_i = x_i & h_mask;
		u64 m_i = (m >> ((t - i) * h)).to_ullong() & h_mask;
		u64 c_i = p_i ^ m_i;

		c.push_back(c_i);
	}
	u64 x_t_plus_1 = (x[t] * x[t]) % n;
	c.push_back(x_t_plus_1);

	return c;
}

template<int Nm>
std::bitset<Nm> decrypt(u64 p, u64 q, i64 a, i64 b, std::vector<u64> c) {
	u64 n = p * q;
	u64 k = (u64)floor(log2(n));
	u64 h = (u64)floor(log2(k));
	u64 t = (u64)ceil((float)Nm / (float)h);

	//Lowest h bits
	u64 h_mask = ~0ull >> (64 - h);

	u64 d_1 = exp_mod((p + 1)/4, t + 1, p - 1);
	u64 d_2 = exp_mod((q + 1)/4, t + 1, q - 1);

	u64 x_t_plus_1 = c[t];
	u64 u = exp_mod(x_t_plus_1, d_1, p);
	u64 v = exp_mod(x_t_plus_1, d_2, q);

	u64 x_0 = pos_mod(((i64)v * a * (i64)p + (i64)u * b * (i64)q), n);

	std::bitset<Nm> m;
	std::vector<u64> x;
	x.push_back(x_0);

	for (i64 i = 1; i <= t; i ++) {
		u64 x_i = (x[i - 1] * x[i - 1]) % n;
		x.push_back(x_i);

		u64 p_i = x_i & h_mask;
		u64 c_i = c[i - 1];
		u64 m_i = (p_i ^ c_i) & h_mask;

		m |= m_i << ((t - i) * h);
	}

	return m;
}

int main(int argc, const char **argv) {
	u64 p = 499;
	u64 q = 547;
	i64 a = -57;
	i64 b = 52;
	i64 x_0 = 159201;


	auto m = std::bitset<20>{"10011100000100001100"};
	std::cout << "m = " << m << std::endl;


	auto c = encrypt<20>(p, q, x_0, m);
	std::cout << "C(m) = {";
	for (i64 i = 0; i < c.size(); i ++) {
		std::cout << c[i];
		if (i < (c.size() - 1))
			std::cout << ", ";
	}
	std::cout << "};" << std::endl;


	auto dec = decrypt<20>(p, q, a, b, c);
	std::cout << "D(C(m)) = " << dec << std::endl;
	assert(dec == m);


	return 0;
}
