#define ERROR -99999.0
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>

#include "OptionPricingFunctions.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ADD FUNCTIONS HERE IF REQUIRED
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double normal_CDF(double x)
{
	return 0.5 * erfc(-x * M_SQRT1_2);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ADD FUNCTIONS HERE IF REQUIRED (END)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// IMPLEMENT THE FUNCTIONS IN THIS SECTION
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void price_vanilla_option_european_bs(double S0, double r,
									  double volatility, double strike_price, double dividend_yield,
									  double expiration_time, double *call_price, double *put_price)
{
	// Analytical pricing of vanilla European option using Black-Scholes pricing forumula
	double K = strike_price;
	double T = expiration_time;
	double q = dividend_yield;
	double d1 = (1/(volatility * sqrt(T))) * (log(S0/K) + (r-q+ 0.5 * pow(volatility, 2))*T);
	double d2 = d1 - volatility * sqrt(T);
	
	//printf("%f\n", normal_CDF(d1));
	//printf("%f\n", normal_CDF(d2));
	*call_price = S0*exp(-q*T)*normal_CDF(d1) - K*exp(-r*T)*normal_CDF(d2);
	*put_price = -S0*exp(-q*T)*normal_CDF(-d1) + K*exp(-r*T)*normal_CDF(-d2);
	//printf("%f", *call_price);

}

double price_vanilla_option_european_recursion(unsigned int i, unsigned int j,
											   unsigned int depth, double S0, double r,
											   double volatility, double strike_price, double dividend_yield,
											   double expiration_time, option_fxn exercise_profit_calculator)
{
	// GENERAL NOTE about the option_fxn parameter - the code for this function will be the same whether pricing a call or a put option,
	// achieve this by calling exercise_profit_calculator variable like:
	// exercise_profit_calculator(strike_price, spot_price) whenever you need to calculate the exercise price.
	// For example if using for call option functions
	// -> exercise_profit_calculator(strike_price, spot_price) is equivalent to max(spot_price - strike_price, 0)

	//////////////// YOUR CODE HERE ////////////////
	//parameters
	double K = strike_price;
	double T = expiration_time;
	double q = dividend_yield;
	double dt = expiration_time / (double) depth;
	double u = exp(volatility * sqrt(dt));
	double d = 1/u;
	double p = (exp(dt*(r-q)) - d) / (u - d);

	//accounting for poor inputs of expiration time
	if (T == 0) {
		return exercise_profit_calculator(K, S0);
	}
	if (T < 0) {
		return 0;
	}
	
	//calculate option value if at the end of the tree
	//otherwise, use the recursion equation which will iterate until the end of the tree is reached
	if (i == depth) {
		//creates signed integer values so they can be subtracted
		//allows the power function to be run without error in subtraction
		int ii = (int) i;
		int jj = (int) j;
		double m = pow(u, 2*jj - ii);
		double option_value = exercise_profit_calculator(K, S0*m);
		//printf("%d, Prob: %f, Option value: %f .\n", j, p, option_value);
		return(exercise_profit_calculator(K, S0*m));
	} else {
		return exp(-r*dt) * (p * price_vanilla_option_european_recursion(i + 1, j +1, depth, S0, r, volatility, K, q,
												T, exercise_profit_calculator)
							+ (1-p) * price_vanilla_option_european_recursion(i + 1, j, depth, S0, r, volatility, K, q,
												T, exercise_profit_calculator) );
	}

	///////////////////////////////////////////////
	return 0.0;	
}

double price_vanilla_option_american_recursion(unsigned int i, unsigned int j,
											   unsigned int depth, double S0, double r,
											   double volatility, double strike_price, double dividend_yield,
											   double expiration_time, option_fxn exercise_profit_calculator)
{
	//////////////// YOUR CODE HERE ////////////////
	//parameters
	double K = strike_price;
	double T = expiration_time;
	double q = dividend_yield;
	double dt = expiration_time / (double) depth;
	double u = exp(volatility * sqrt(dt));
	double d = 1/u;
	double p = (exp(dt*(r-q)) - d) / (u - d);

	//accounting for poor inputs of expiration time
	if (T == 0) {
		return exercise_profit_calculator(K, S0);
	}
	if (T < 0) {
		return 0;
	}

	//creates signed integer values so they can be subtracted
	//allows the power function to be run without error in subtraction
	int ii = (int) i;
	int jj = (int) j;
	double m = pow(u, 2*jj - ii);
	double exercise_value = exercise_profit_calculator(K, S0*m);
	//printf("\nExercise: %f", exercise_value);

	//calculate option value if at the end of the tree
	//otherwise, use the recursion equation which will iterate until the end of the tree is reached
	if (i == depth) {
		//printf("%d, Prob: %f, Option value: %f .\n", j, p, exercise_value);
		return(exercise_profit_calculator(K, S0*m));
	} else {
		double hold_value = exp(-r*dt) * (p * price_vanilla_option_american_recursion(i + 1, j +1, depth, S0, 
						r, volatility, K, q, T, exercise_profit_calculator)
						+ (1-p) * price_vanilla_option_american_recursion(i + 1, j, depth, S0, r, volatility, 
						K, q, T, exercise_profit_calculator) );
		//printf("\nHOLD, %f", hold_value);
		return(fmax(hold_value, exercise_value));
	}
	////////////////////////////////////////////////
	return 0.0;
}

double price_vanilla_option_european_induction(
	unsigned int depth, double S0, double r,
	double volatility, double strike_price, double dividend_yield,
	double expiration_time, option_fxn exercise_profit_calculator)
{
	// array to store the option values at each timestep
	double *depth_values = (double *)malloc((depth + 1) * sizeof(double));

	// return error code if memory allocation failed
	if (depth_values == NULL)
		return ERROR;

	//////////////// YOUR CODE HERE ////////////////
	//parameters
	double K = strike_price;
	double T = expiration_time;
	double q = dividend_yield;
	double dt = expiration_time / (double) depth;
	double u = exp(volatility * sqrt(dt));
	double d = 1/u;
	double p = (exp(dt*(r-q)) - d) / (u - d);

	//accounting for poor inputs of expiration time
	if (T == 0) {
		return exercise_profit_calculator(K, S0);
	}
	if (T < 0) {
		return 0;
	}

	//runs through all tree depths and correspondig node heights
	//starts with last depth where option value calculated
	for ( int deep = (int) depth; deep >= 0; deep--) {
		for (int height = 0; height <= deep; height++) {
			if(deep == depth) {
				depth_values[height] = exercise_profit_calculator(K, S0*pow(u, 2*height - deep));
				//printf("%d Option Value: %f\n", deep, depth_values[height]);
			} else {
				depth_values[height] = exp(-r * dt) * (p*depth_values[height + 1] + (1-p)*depth_values[height]);
				//printf("%d Option Value: %f\n", deep, depth_values[height]);
			}
		}
	}

	///////////////////////////////////////////////

	//0, 0 node of tree
	double ret = depth_values[0];

	// free the memory allocated for the array
	free(depth_values);

	return ret;
}

double price_vanilla_option_american_induction(
	unsigned int depth, double S0, double r,
	double volatility, double strike_price, double dividend_yield,
	double expiration_time, option_fxn exercise_profit_calculator)
{
	// array to store the option values at each timestep
	double *depth_values = (double *)malloc((depth + 1) * sizeof(double));

	// return error code if memory allocation failed
	if (depth_values == NULL)
		return ERROR;

	//////////////// YOUR CODE HERE ////////////////
	//parameters
	double K = strike_price;
	double T = expiration_time;
	double q = dividend_yield;
	double dt = expiration_time / (double) depth;
	double u = exp(volatility * sqrt(dt));
	double d = 1/u;
	double p = (exp(dt*(r-q)) - d) / (u - d);

	//accounting for poor inputs of expiration time
	if (T == 0) {
		return exercise_profit_calculator(K, S0);
	}
	if (T < 0) {
		return 0;
	}

	//runs through all tree depths and correspondig node heights
	//starts with last depth where option value calculated
	for (int deep = (int) depth; deep >= 0; deep--) {
		for (int height = 0; height <= deep; height++) {
			if(deep == depth) {
				depth_values[height] = exercise_profit_calculator(K, S0*pow(u, 2*height - deep));
				//printf("%d Option Value: %f\n", deep, depth_values[height]);
			} else {
				depth_values[height] = fmax(exercise_profit_calculator(K, S0*pow(u, 2*height - deep)),
											exp(-r * dt) * (p*depth_values[height + 1] + (1-p)*depth_values[height]));
				//printf("%d Option Value: %f\n", deep, depth_values[height]);
			}
		}
	}
	///////////////////////////////////////////////

	//0, 0 node of tree
	double ret = depth_values[0];

	// free the memory allocated for the array
	free(depth_values);

	return ret;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// DO NOT MODIFY ANYTHING BELOW THIS LINE
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double call_exercise_profit(double strike_price, double spot_price)
{
	return fmax(spot_price - strike_price, 0);
}

double put_exercise_profit(double strike_price, double spot_price)
{
	return fmax(strike_price - spot_price, 0);
}

double price_vanilla_option_european_recursion_call(
	unsigned int depth, double S0, double r,
	double volatility, double strike_price, double dividend_yield,
	double expiration_time)
{
	return price_vanilla_option_european_recursion(0, 0,
												   depth, S0, r,
												   volatility, strike_price, dividend_yield,
												   expiration_time, &call_exercise_profit);
}

double price_vanilla_option_european_recursion_put(
	unsigned int depth, double S0, double r,
	double volatility, double strike_price, double dividend_yield,
	double expiration_time)
{
	return price_vanilla_option_european_recursion(0, 0,
												   depth, S0, r,
												   volatility, strike_price, dividend_yield,
												   expiration_time, &put_exercise_profit);
}

double price_vanilla_option_american_recursion_call(
	unsigned int depth, double S0, double r,
	double volatility, double strike_price, double dividend_yield,
	double expiration_time)
{
	return price_vanilla_option_american_recursion(0, 0,
												   depth, S0, r,
												   volatility, strike_price, dividend_yield,
												   expiration_time, &call_exercise_profit);
}

double price_vanilla_option_american_recursion_put(
	unsigned int depth, double S0, double r,
	double volatility, double strike_price, double dividend_yield,
	double expiration_time)
{
	return price_vanilla_option_american_recursion(0, 0,
												   depth, S0, r,
												   volatility, strike_price, dividend_yield,
												   expiration_time, &put_exercise_profit);
}

double price_vanilla_option_european_induction_call(
	unsigned int depth, double S0, double r,
	double volatility, double strike_price, double dividend_yield,
	double expiration_time)
{
	return price_vanilla_option_european_induction(
		depth, S0, r,
		volatility, strike_price, dividend_yield,
		expiration_time, &call_exercise_profit);
}

double price_vanilla_option_european_induction_put(
	unsigned int depth, double S0, double r,
	double volatility, double strike_price, double dividend_yield,
	double expiration_time)
{
	return price_vanilla_option_european_induction(
		depth, S0, r,
		volatility, strike_price, dividend_yield,
		expiration_time, &put_exercise_profit);
}

double price_vanilla_option_american_induction_call(
	unsigned int depth, double S0, double r,
	double volatility, double strike_price, double dividend_yield,
	double expiration_time)
{
	return price_vanilla_option_american_induction(
		depth, S0, r,
		volatility, strike_price, dividend_yield,
		expiration_time, &call_exercise_profit);
}

double price_vanilla_option_american_induction_put(
	unsigned int depth, double S0, double r,
	double volatility, double strike_price, double dividend_yield,
	double expiration_time)
{
	return price_vanilla_option_american_induction(
		depth, S0, r,
		volatility, strike_price, dividend_yield,
		expiration_time, &put_exercise_profit);
}