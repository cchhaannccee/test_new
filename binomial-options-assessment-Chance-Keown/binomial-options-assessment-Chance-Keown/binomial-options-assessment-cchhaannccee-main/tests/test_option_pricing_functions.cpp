#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <math.h>

extern "C"
{
#include "OptionPricingFunctions.h"
}

// error tolerance for floating point comparisons
#define EPSILON 1e-14

//error between european call price methods
void european_call_error(double S0, double volatility, double strike_price, double risk_free_rate, double dividend_yield, double expiration_time, unsigned int depth)
{
	double delta_t = expiration_time / (double) depth;

    REQUIRE(delta_t < (volatility * volatility) / ((risk_free_rate - dividend_yield) * (risk_free_rate - dividend_yield)));

    // get induction and recursion answers for the same option pricing problem
	double induction_answer = price_vanilla_option_european_induction_call
	  (depth, S0, risk_free_rate,
		volatility, strike_price, dividend_yield,
		expiration_time);

	double recursion_answer = price_vanilla_option_european_recursion_call(
		depth, S0, risk_free_rate,
		volatility, strike_price, dividend_yield,
		expiration_time);

    // check answers are equal within the tolerance epsilon
    REQUIRE(induction_answer == Approx(recursion_answer).epsilon(EPSILON));
}

//error between american call price methods
void american_call_error(double S0, double volatility, double strike_price, double risk_free_rate, double dividend_yield, double expiration_time, unsigned int depth)
{
	double delta_t = expiration_time / (double) depth;

    REQUIRE(delta_t < (volatility * volatility) / ((risk_free_rate - dividend_yield) * (risk_free_rate - dividend_yield)));

    // get induction and recursion answers for the same option pricing problem
	double induction_answer = price_vanilla_option_american_induction_call
	  (depth, S0, risk_free_rate,
		volatility, strike_price, dividend_yield,
		expiration_time);

	double recursion_answer = price_vanilla_option_american_recursion_call(
		depth, S0, risk_free_rate,
		volatility, strike_price, dividend_yield,
		expiration_time);

    // check answers are equal within the tolerance epsilon
    REQUIRE(induction_answer == Approx(recursion_answer).epsilon(EPSILON));
}

//error between european put price methods
void european_put_error(double S0, double volatility, double strike_price, double risk_free_rate, double dividend_yield, double expiration_time, unsigned int depth)
{
	double delta_t = expiration_time / (double) depth;

    REQUIRE(delta_t < (volatility * volatility) / ((risk_free_rate - dividend_yield) * (risk_free_rate - dividend_yield)));

    // get induction and recursion answers for the same option pricing problem
	double induction_answer = price_vanilla_option_european_induction_put
	  (depth, S0, risk_free_rate,
		volatility, strike_price, dividend_yield,
		expiration_time);

	double recursion_answer = price_vanilla_option_european_recursion_put(
		depth, S0, risk_free_rate,
		volatility, strike_price, dividend_yield,
		expiration_time);

    // check answers are equal within the tolerance epsilon
    REQUIRE(induction_answer == Approx(recursion_answer).epsilon(EPSILON));
}

//error between american put price methods
void american_put_error(double S0, double volatility, double strike_price, double risk_free_rate, double dividend_yield, double expiration_time, unsigned int depth)
{
	double delta_t = expiration_time / (double) depth;

    REQUIRE(delta_t < (volatility * volatility) / ((risk_free_rate - dividend_yield) * (risk_free_rate - dividend_yield)));

    // get induction and recursion answers for the same option pricing problem
	double induction_answer = price_vanilla_option_american_induction_put
	  (depth, S0, risk_free_rate,
		volatility, strike_price, dividend_yield,
		expiration_time);

	double recursion_answer = price_vanilla_option_american_recursion_put(
		depth, S0, risk_free_rate,
		volatility, strike_price, dividend_yield,
		expiration_time);

    // check answers are equal within the tolerance epsilon
    REQUIRE(induction_answer == Approx(recursion_answer).epsilon(EPSILON));
}

//error between american call price methods
void exact_call_error(double S0, double volatility, double strike_price, double risk_free_rate, double dividend_yield, double expiration_time, unsigned int depth)
{
	double delta_t = expiration_time / (double) depth;

    REQUIRE(delta_t < (volatility * volatility) / ((risk_free_rate - dividend_yield) * (risk_free_rate - dividend_yield)));

    // get induction and recursion answers for the same option pricing problem
	double induction_answer = price_vanilla_option_american_induction_call
	  (depth, S0, risk_free_rate,
		volatility, strike_price, dividend_yield,
		expiration_time);

	double recursion_answer = price_vanilla_option_american_recursion_call(
		depth, S0, risk_free_rate,
		volatility, strike_price, dividend_yield,
		expiration_time);

    // check answers are equal within the tolerance epsilon
    REQUIRE(induction_answer == Approx(recursion_answer).epsilon(EPSILON));
}

//test cases to compare error between european call price methods with varied inputs
//potential edge cases and unwanted inputs considered
//test case combinations of these parameters where delta t >= vol^2/(r-dividend)^2 not considered
TEST_CASE("european call error varied params", "[tests]")
{
	double c_S0[] = {100, 0, -10};
	double c_vol[] = {.001, .05, 1.3};
	double c_K[] = {95, 110, 0, -10};
	double c_r[] = {.02, -.02, 0, 1.1};
	double c_q[] = {.05, -.05, 0, 1.2};
	double c_T [] = {1, -1, 0, M_SQRT1_2};
	unsigned int c_depth[] = {1, 5, 10};


	for(int i = 0; i <= 2; i++){
		for(int j = 0; j <= 2; j++){
			for(int k = 0; k <= 3; k++){
				for(int l = 0; l <= 3; l++){
					for(int m = 0; m <= 3; m++){
						for(int n = 0; n <= 3; n++){
							for(int o = 0; o <= 2; o++){
								if (c_T[n]/(double) c_depth[o] >= c_vol[j]*c_vol[j] / ((c_r[l]-c_q[m]) * (c_r[l]-c_q[m]))) {
									break;
								} else {
									european_call_error(c_S0[i], c_vol[j], c_K[k], c_r[l], c_q[m], c_T[n], c_depth[o]);
								}
							}
						}
					}
				}
			}
		}
	}
}

//test cases to compare error between american call price methods with varied inputs
//potential edge cases and unwanted inputs considered
//test case combinations of these parameters where delta t >= vol^2/(r-dividend)^2 not considered
TEST_CASE("american call error varied params", "[tests]")
{
	double c_S0[] = {100, 0, -10};
	double c_vol[] = {.001, .05, 1.3};
	double c_K[] = {95, 110, 0, -10};
	double c_r[] = {.02, -.02, 0, 1.1};
	double c_q[] = {.05, -.05, 0, 1.2};
	double c_T [] = {1, -1, 0, M_SQRT1_2};
	unsigned int c_depth[] = {1, 5, 10};


	for(int i = 0; i <= 2; i++){
		for(int j = 0; j <= 2; j++){
			for(int k = 0; k <= 3; k++){
				for(int l = 0; l <= 3; l++){
					for(int m = 0; m <= 3; m++){
						for(int n = 0; n <= 3; n++){
							for(int o = 0; o <= 2; o++){
								if (c_T[n]/(double) c_depth[o] >= c_vol[j]*c_vol[j] / ((c_r[l]-c_q[m]) * (c_r[l]-c_q[m]))) {
									break;
								} else {
									american_call_error(c_S0[i], c_vol[j], c_K[k], c_r[l], c_q[m], c_T[n], c_depth[o]);
								}
							}
						}
					}
				}
			}
		}
	}
}

//test cases to compare error between european put price methods with varied inputs
//potential edge cases and unwanted inputs considered
//test case combinations of these parameters where delta t >= vol^2/(r-dividend)^2 not considered
TEST_CASE("european put error varied params", "[tests]")
{
	double c_S0[] = {100, 0, -10};
	double c_vol[] = {.001, .05, 1.3};
	double c_K[] = {95, 110, 0, -10};
	double c_r[] = {.02, -.02, 0, 1.1};
	double c_q[] = {.05, -.05, 0, 1.2};
	double c_T [] = {1, -1, 0, M_SQRT1_2};
	unsigned int c_depth[] = {1, 5, 10};


	for(int i = 0; i <= 2; i++){
		for(int j = 0; j <= 2; j++){
			for(int k = 0; k <= 3; k++){
				for(int l = 0; l <= 3; l++){
					for(int m = 0; m <= 3; m++){
						for(int n = 0; n <= 3; n++){
							for(int o = 0; o <= 2; o++){
								if (c_T[n]/(double) c_depth[o] >= c_vol[j]*c_vol[j] / ((c_r[l]-c_q[m]) * (c_r[l]-c_q[m]))) {
									break;
								} else {
									european_put_error(c_S0[i], c_vol[j], c_K[k], c_r[l], c_q[m], c_T[n], c_depth[o]);
								}
							}
						}
					}
				}
			}
		}
	}
}

//test cases to compare error between american put price methods with varied inputs
//potential edge cases and unwanted inputs considered
//test case combinations of these parameters where delta t >= vol^2/(r-dividend)^2 not considered
TEST_CASE("american put error varied params", "[tests]")
{
	double c_S0[] = {100, 0, -10};
	double c_vol[] = {.001, .05, 1.3};
	double c_K[] = {95, 110, 0, -10};
	double c_r[] = {.02, -.02, 0, 1.1};
	double c_q[] = {.05, -.05, 0, 1.2};
	double c_T [] = {1, -1, 0, M_SQRT1_2};
	unsigned int c_depth[] = {1, 5, 10};


	for(int i = 0; i <= 2; i++){
		for(int j = 0; j <= 2; j++){
			for(int k = 0; k <= 3; k++){
				for(int l = 0; l <= 3; l++){
					for(int m = 0; m <= 3; m++){
						for(int n = 0; n <= 3; n++){
							for(int o = 0; o <= 2; o++){
								if (c_T[n]/(double) c_depth[o] >= c_vol[j]*c_vol[j] / ((c_r[l]-c_q[m]) * (c_r[l]-c_q[m]))) {
									break;
								} else {
									american_put_error(c_S0[i], c_vol[j], c_K[k], c_r[l], c_q[m], c_T[n], c_depth[o]);
								}
							}
						}
					}
				}
			}
		}
	}
}


