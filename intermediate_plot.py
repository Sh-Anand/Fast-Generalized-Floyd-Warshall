import matplotlib.pyplot as plt

x = [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
x_O0 = [4, 8, 16, 32, 64, 128, 256]

y_min_plus_baseline = [0.606635,
                       0.993403,
                       1.378475,
                       1.568875,
                       1.690234,
                       1.644125,
                       1.693905,
                       1.230458,
                       0.753746,
                       0.693780,
                       0.709773]

y_min_plus_baseline_O0 = [0.17302,
                          0.227283,
                          0.216871,
                          0.202217,
                          0.22129,
                          0.179687,
                          0.143159
                          ]

y_min_plus_opt = [0.807062,
                 0.968963,
                 1.158469,
                 1.209638,
                 1.237534,
                 0.990288,
                 0.720366,
                 0.609435,
                 0.598687,
                 0.535071,
                 0.529993]

y_min_plus_opt_O0 = [0.156174,
                     0.202252,
                     0.163699,
                     0.118519,
                     0.115895,
                     0.114807,
                     0.116]

# plt.plot(x, y_min_plus_baseline, label="baseline", linewidth=2.0)
# plt.plot(x, y_min_plus_opt, label="basic optimization", linewidth=2.0)
plt.plot(x_O0, y_min_plus_baseline_O0, label="baseline, flag -O0", linewidth=2.0)
plt.plot(x_O0, y_min_plus_opt_O0, label="basic optimization, flag -O0", linewidth=2.0)
plt.legend()
plt.title("min plus base optimization")
plt.xscale('log')
plt.xlabel("N")
plt.ylabel("flops/cycle")
plt.savefig("./plots/min_plus_base_opt")

#==============================================================
y_max_min_baseline = [0.644512,
                      1.192917,
                      1.340972,
                      1.616624,
                      1.735605,
                      1.653051,
                      1.697671,
                      1.754933,
                      1.338690,
                      0.792861,
                      0.672087]

y_max_min_baseline_O0 = [0.090395,
                         0.091512,
                         0.086091,
                         0.083615,
                         0.079356,
                         0.081099,
                         0.080396]

y_max_min_opt = [0.523303,
                 0.573927,
                 0.619274,
                 0.690119,
                 0.721193,
                 0.671105,
                 0.693432,
                 0.690820,
                 0.631541,
                 0.540784,
                 0.534826]

y_max_min_opt_O0 = [0.1,
                    0.1067,
                    0.108453,
                    0.104666,
                    0.097379,
                    0.102803,
                    0.101958]


# # plt.plot(x, y_max_min_baseline, label="baseline", linewidth=2.0)
# # plt.plot(x, y_max_min_opt, label="basic optimization", linewidth=2.0)
# plt.plot(x_O0, y_max_min_baseline_O0, label="baseline, flag -O0", linewidth=2.0)
# plt.plot(x_O0, y_max_min_opt_O0, label="basic opt, flag -O0", linewidth=2.0)
# plt.legend()
# plt.title("max min base optimization")
# plt.xscale('log')
# plt.xlabel("N")
# plt.ylabel("flops/cycle")
# plt.savefig("./plots/max_min_base_opt")

#==============================================================

y_or_and_baseline = [0.416667,
                     0.578858,
                     0.738271,
                     0.913595,
                     1.012184,
                     0.954492,
                     0.988117,
                     0.990364,
                     0.888463,
                     0.728210,
                     0.694072]

y_or_and_baseline_O0 = [0.155794,
                        0.15987,
                        0.166228,
                        0.171097,
                        0.17807,
                        0.169556,
                        0.165158
                        ]

y_or_and_opt = [0.50513,
                0.548356,
                0.606734,
                0.708268,
                0.768382,
                0.735261,
                0.720346,
                0.724475,
                0.659623,
                0.552817,
                0.549917]

y_or_and_opt_O0 = [0.107149,
                   0.126507,
                   0.134511,
                   0.143284,
                   0.143353,
                   0.1455,
                   0.143868]

# # plt.plot(x, y_or_and_baseline, label="baseline", linewidth=2.0)
# # plt.plot(x, y_or_and_opt, label="basic optimization", linewidth=2.0)
# plt.plot(x_O0, y_or_and_baseline_O0, label="baseline, flag -O0", linewidth=2.0)
# plt.plot(x_O0, y_or_and_opt_O0, label="basic optimization, flag -O0", linewidth=2.0)
# plt.legend()
# plt.title("or and base optimization")
# plt.xscale('log')
# plt.xlabel("N")
# plt.ylabel("flops/cycle")
# plt.savefig("./plots/or_and_base_opt")
