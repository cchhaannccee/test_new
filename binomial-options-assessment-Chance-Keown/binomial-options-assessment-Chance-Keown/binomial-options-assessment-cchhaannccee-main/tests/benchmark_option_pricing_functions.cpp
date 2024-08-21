#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"

extern "C"
{
#include "OptionPricingFunctions.h"
}

TEST_CASE("Benchmark recursion vs induction", "[benchmark]")
{
    double S0 = 100.0;
    double volatility = 0.05;
    double strike_price = 95;
    double dividend_yield = 0.05;
    double risk_free_rate = 0.02;
    double expiration_time = 1;

    //////////////////////////////////////////
    // Induction
    //////////////////////////////////////////
    //EUROPEAN
    BENCHMARK("induction european depth 5")
    {
        return price_vanilla_option_european_induction_call(5, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };

    BENCHMARK("induction european depth 10")
    {
        return price_vanilla_option_european_induction_call(10, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };
    BENCHMARK("induction european depth 15")
    {
        return price_vanilla_option_european_induction_call(15, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };
    BENCHMARK("induction european depth 20")
    {
        return price_vanilla_option_european_induction_call(20, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };
    BENCHMARK("induction european depth 50")
    {
        return price_vanilla_option_european_induction_call(50, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };
    BENCHMARK("induction european depth 200")
    {
        return price_vanilla_option_european_induction_call(200, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };
    BENCHMARK("induction european depth 1000")
    {
        return price_vanilla_option_european_induction_call(1000, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };

    //AMERICAN
    BENCHMARK("induction american depth 5")
    {
        return price_vanilla_option_american_induction_call(5, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };

    BENCHMARK("induction american depth 10")
    {
        return price_vanilla_option_american_induction_call(10, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };
    BENCHMARK("induction american depth 15")
    {
        return price_vanilla_option_american_induction_call(15, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };
    BENCHMARK("induction american depth 20")
    {
        return price_vanilla_option_american_induction_call(20, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };
    BENCHMARK("induction american depth 50")
    {
        return price_vanilla_option_american_induction_call(50, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };
    BENCHMARK("induction american depth 200")
    {
        return price_vanilla_option_american_induction_call(200, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };
    BENCHMARK("induction american depth 1000")
    {
        return price_vanilla_option_american_induction_call(1000, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };
    //////////////////////////////////////////
    // recursion
    //////////////////////////////////////////
    //EUROPEAN
    BENCHMARK("recursion european depth 5")
    {
        return price_vanilla_option_european_recursion_call(5, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };

    BENCHMARK("recursion european depth 10")
    {
        return price_vanilla_option_european_recursion_call(10, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };
    BENCHMARK("recursion european depth 15")
    {
        return price_vanilla_option_european_recursion_call(15, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };
    BENCHMARK("recursion european depth 18")
    {
        return price_vanilla_option_european_recursion_call(18, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };
    BENCHMARK("recursion european depth 20")
    {
        return price_vanilla_option_european_recursion_call(20, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };
    //AMERICAN
    BENCHMARK("recursion american depth 5")
    {
        return price_vanilla_option_american_recursion_call(5, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };

    BENCHMARK("recursion american depth 10")
    {
        return price_vanilla_option_american_recursion_call(10, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };
    BENCHMARK("recursion american depth 15")
    {
        return price_vanilla_option_american_recursion_call(15, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };
    BENCHMARK("recursion american depth 18")
    {
        return price_vanilla_option_american_recursion_call(18, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };
    BENCHMARK("recursion american depth 20")
    {
        return price_vanilla_option_american_recursion_call(20, S0, risk_free_rate,
                                                            volatility, strike_price, dividend_yield,
                                                            expiration_time);
    };
};