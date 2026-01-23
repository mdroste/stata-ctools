adopath ++ "./build"
webuse lifeexp, clear
describe country
capture drop country_code
cencode country, generate(country_code) verbose
describe country_code
sum country_code
list country country_code in 1/10
label list country_code
