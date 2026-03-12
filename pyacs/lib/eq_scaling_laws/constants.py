
MU_CRUST = 30e9
MU_SUBDUCTION = 40e9
MU_SCR = 33e9
BETA = 2/3

# Leonard 2014 BSSA canonical coefficients
LEONARD14_COEFFICIENTS = {
    "crustal-strike":  {"Cw":1.17,"Cd":1.50e-5,"mu":MU_CRUST},
    "crustal-reverse": {"Cw":1.44,"Cd":1.80e-5,"mu":MU_CRUST},
    "crustal-normal":  {"Cw":1.02,"Cd":2.00e-5,"mu":MU_CRUST},
    "subduction":      {"Cw":1.90,"Cd":1.20e-5,"mu":MU_SUBDUCTION},
    "scr":             {"Cw":0.90,"Cd":2.50e-5,"mu":MU_SCR},
}
