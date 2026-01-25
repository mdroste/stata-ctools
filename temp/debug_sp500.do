sysuse sp500, clear
describe
tostring date, generate(date_str) format(%tdCCYY-NN-DD)
list date date_str in 1/5
