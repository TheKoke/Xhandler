MASS_EXCESSES = {
    (1, 1)  :  7.28900, (1, 2)  :  13.1357, (1, 3)  :  14.9498,
    (2, 3)  :  14.9312, (2, 4)  :  2.42490,
    (3, 6)  :  14.0869, (3, 7)  :  14.9071, (3, 8)  :  20.9458,
    (4, 8)  :  4.94170, (4, 9)  :  11.3485, (4, 10) :  12.6075, (4, 11) :  20.1772,
    (5, 8)  :  22.9216, (5, 10) :  12.0506, (5, 11) :  8.66770, (5, 12) :  13.3694,
    (6, 11) :  10.6494, (6, 12) :  0.00000, (6, 13) :  3.12500, (6, 14) :  3.01990,
    (7, 13) :  5.34550, (7, 14) :  2.86340, (7, 15) :  0.10140, (7, 16) :  5.68390,
    (8, 14) :  8.00780, (8, 15) :  2.85560, (8, 16) :  -4.7370, (8, 17) :  -0.8088, (8, 18) :  -0.7828, (8, 19) : 3.3329,
    (9, 17) :  1.95170, (9, 18) :  0.87310, (9, 19) :  -1.4874, (9, 20) :  -0.0175, (9, 21) :  -0.0476, (9, 22) : 2.7934,
    (10, 19) : 1.75210, (10, 20) : -0.7419, (10, 21) : -5.7318, (10, 22) : -8.0247, (10, 23) : -5.1540
}


PROTON_STATES = [0]
DEUTERON_STATES = [0]
TRITON_STATES = [0]
HE_3_STATES = [0]
HE_4_STATES = [0, 23.64, 24.25, 25.95, 27.42]
LI_6_STATES = [0, 2.186, 3.563, 4.312, 5.366, 5.65]
LI_7_STATES = [0, 0.478, 4.63, 6.68, 7.46, 9.67]
LI_8_STATES = [0, 0.981, 2.255, 3.21, 5.4, 6.53]
BE_8_STATES = [0, 3.03, 11.35, 16.626, 16.922, 17.64]
BE_9_STATES = [0, 1.684, 2.429, 3.049, 4.704, 6.38, 11.282]
BE_10_STATES = [0, 3.368, 5.958, 6.179, 6.263, 7.371, 7.542]
BE_11_STATES = [0, 0.320, 1.783, 2.654, 3.400, 3.889, 3.955]
B_8_STATES = [0, 0.769, 2.320, 3.500, 10.619]
B_10_STATES = [0, 0.718, 1.74, 2.154, 3.587, 4.774, 5.11]
B_11_STATES = [0, 2.125, 4.445, 5.02, 6.742, 6.792, 7.286]
B_12_STATES = [0, 0.953, 1.674, 2.621, 2.723, 3.389, 3.76, 4.302]
C_11_STATES = [0, 2.0, 4.319, 4.804, 6.34, 6.478, 6.905]
C_12_STATES = [0, 4.444, 7.654, 9.641, 10.847, 11.836]
C_13_STATES = [0, 3.089, 3.684, 3.854, 6.864, 7.492, 7.547]
C_14_STATES = [0, 6.094, 6.589, 6.728, 6.903, 7.012, 7.341]
N_13_STATES = []
N_14_STATES = []
N_15_STATES = []
N_16_STATES = []
O_14_STATES = []
O_15_STATES = []
O_16_STATES = []
O_17_STATES = []
O_18_STATES = []
O_19_STATES = []
F_17_STATES = []
F_18_STATES = []
F_19_STATES = []
F_20_STATES = []
F_21_STATES = []
F_22_STATES = []
NE_19_STATES = []
NE_20_STATES = []
NE_21_STATES = []
NE_22_STATES = []
NE_23_STATES = []


STATES = {
    (1, 1) :  PROTON_STATES, (1, 2) : DEUTERON_STATES, (1, 3):    TRITON_STATES,
    (2, 3) :    HE_3_STATES, (2, 4) :     HE_4_STATES,
    (3, 6) :    LI_6_STATES, (3, 7) :     LI_7_STATES, (3, 8) :     LI_8_STATES,
    (4, 8) :    BE_8_STATES, (4, 9) :     BE_9_STATES, (4, 10) :   BE_10_STATES, (4, 11) :  BE_11_STATES,
    (5, 8) :     B_8_STATES, (5, 10) :    B_10_STATES, (5, 11) :    B_11_STATES, (5, 12) :   B_12_STATES,
    (6, 11) :   C_11_STATES, (6, 12) :    C_12_STATES, (6, 13) :    C_13_STATES, (6, 14) :   C_14_STATES, 
    (7, 13) :   N_13_STATES, (7, 14) :    N_14_STATES, (7, 15) :    N_15_STATES, (7, 16) :   N_16_STATES,
    (8, 14) :   O_14_STATES, (8, 15) :    O_15_STATES, (8, 16) :    O_16_STATES, (8, 17) :   O_17_STATES, (8, 18) :   O_18_STATES, (8, 19) : O_19_STATES,
    (9, 17) :   F_17_STATES, (9, 18) :    F_18_STATES, (9, 19) :    F_19_STATES, (9, 20) :   F_20_STATES, (9, 21) :   F_21_STATES, (9, 22) : F_22_STATES,
    (10, 19) : NE_19_STATES, (10, 20) :  NE_20_STATES, (10, 21) :  NE_21_STATES, (10, 22) : NE_22_STATES, (10, 23) : NE_23_STATES,
}