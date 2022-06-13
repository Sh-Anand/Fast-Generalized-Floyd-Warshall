import matplotlib.pyplot as plt

x = [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]

y_min_plus_baseline = [0.633663,
                       0.841413,
                       1.002325,
                       1.09367,
                       1.175905,
                       1.131389,
                       1.14603,
                       1.147604,
                       1.045012,
                       0.85719,
                       0.816123
                       ]

y_min_plus_opt = [0.684492,
                  0.773414,
                  0.777156,
                  0.827611,
                  0.837963,
                  0.803499,
                  0.819047,
                  0.835346,
                  0.75816,
                  0.660645,
                  0.658655
                  ]

y_min_plus_tiled = [0.043956,
                    0.185675,
                    0.083806,
                    2.230557,
                    0.278412,
                    0.292648,
                    1.316907,
                    0.073953,
                    0.721084,
                    1.187293,
                    0.268666
                    ]

y_min_plus_tiled_vect= [0.035526,
                        0.128369,
                        0.064637,
                        2.76325,
                        0.210248,
                        0.197561,
                        1.4197,
                        0.04444,
                        2.132966,
                        1.306201,
                        0.157309
                        ]

plt.plot(x, y_min_plus_baseline, label="baseline", marker='o')
plt.plot(x, y_min_plus_opt, label="basic optimization", marker='o')
plt.plot(x, y_min_plus_tiled, label="tiled", marker='o')
plt.plot(x, y_min_plus_tiled_vect, label="tiled vectorized", marker='o')
plt.legend(bbox_to_anchor=(1.38, 1.02))
plt.suptitle("(min, +) performances", fontweight="bold", fontsize=12, x=0.285)
plt.title("flags: -march=native -O3 -ffast-math", loc='left')
plt.xscale('log', base=2)
plt.xlabel("N")
plt.ylabel("flops/cycle")
plt.savefig("./plots/min_plus_performances.pdf", bbox_inches="tight")

#==============================================================
y_max_min_baseline = [0.61244,
                      0.855472,
                      1.006759,
                      1.143495,
                      1.237267,
                      1.138247,
                      1.146389,
                      1.21808,
                      1.162767,
                      0.879274,
                      0.83987
                      ]

y_max_min_opt = [0.649746,
                 0.720619,
                 0.76148,
                 0.832764,
                 0.844169,
                 0.806621,
                 0.833565,
                 0.837667,
                 0.784712,
                 0.645265,
                 0.678994
                 ]

y_max_min_tiled = [0.056313,
                   0.266251,
                   0.124268,
                   1.845201,
                   0.287868,
                   0.333735,
                   1.125903,
                   0.119227,
                   0.721661,
                   1.220366,
                   0.444199
                   ]

y_max_min_tiled_vect = [0.036607,
                        0.136719,
                        0.063257,
                        2.659935,
                        0.212769,
                        0.190021,
                        1.40029,
                        0.040739,
                        2.022769,
                        1.163224,
                        0.128313
                        ]

# plt.plot(x, y_max_min_baseline, label="baseline", marker='o')
# plt.plot(x, y_max_min_opt, label="basic optimization", marker='o')
# plt.plot(x, y_max_min_tiled, label="tiled", marker='o')
# plt.plot(x, y_max_min_tiled_vect, label="tiled vectorized", marker='o')
# plt.legend(bbox_to_anchor=(1.38, 1.02))
# plt.suptitle("(max, min) performances", fontweight="bold", fontsize=12, x=0.31)
# plt.title("flags: -march=native -O3 -ffast-math", loc='left')
# plt.xscale('log', base=2)
# plt.xlabel("N")
# plt.ylabel("flops/cycle")
# plt.savefig("./plots/max_min_performances.svg", bbox_inches="tight")

#==============================================================

y_or_and_baseline = [0.719101,
                     0.905393,
                     1.073235,
                     1.252552,
                     1.295312,
                     1.219713,
                     1.244866,
                     1.255767,
                     1.163521,
                     0.876437,
                     0.841047
                     ]

y_or_and_opt = [0.688172,
                0.747991,
                0.752043,
                0.693011,
                0.851258,
                0.800325,
                0.81615,
                0.82648,
                0.776472,
                0.641221,
                0.645115
                ]

y_or_and_tiled = [0.05699,
                  0.257222,
                  0.114975,
                  1.840485,
                  0.253245,
                  0.312642,
                  1.113622,
                  0.115247,
                  0.719165,
                  1.22158,
                  0.425342
                  ]

y_or_and_tiled_vect = [0.042314,
                       0.12776,
                       0.073844,
                       1.913907,
                       0.242062,
                       0.187179,
                       1.192032,
                       0.044487,
                       1.983928,
                       1.313406,
                       0.146329
                       ]

# plt.plot(x, y_or_and_baseline, label="baseline", marker='o')
# plt.plot(x, y_or_and_opt, label="basic optimization", marker='o')
# plt.plot(x, y_or_and_tiled, label="tiled", marker='o')
# plt.plot(x, y_or_and_tiled_vect, label="tiled vectorized", marker='o')
# plt.legend(bbox_to_anchor=(1.38, 1.02))
# plt.suptitle("(or, and) performances", fontweight="bold", fontsize=12, x=0.29)
# plt.title("flags: -march=native -O3 -ffast-math", loc='left')
# plt.xscale('log', base=2)
# plt.xlabel("N")
# plt.ylabel("flops/cycle")
# plt.savefig("./plots/or_and_performances.svg", bbox_inches="tight")
