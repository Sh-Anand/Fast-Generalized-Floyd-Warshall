import matplotlib.pyplot as plt

x = [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
y_min_plus = [0.807062,
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

#
# plt.plot(x, y_min_plus, linewidth=2.0)
# plt.title("min plus base optimization")
# plt.xscale('log')
# plt.xlabel("N")
# plt.ylabel("flops/cycle")
# plt.savefig("./plots/min_plus_base_opt")

y_max_min = [0.523303,
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


# plt.plot(x, y_max_min, linewidth=2.0)
# plt.title("max min base optimization")
# plt.xscale('log')
# plt.xlabel("N")
# plt.ylabel("flops/cycle")
# plt.savefig("./plots/max_min_base_opt")

y_or_and = [0.416667,
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


plt.plot(x, y_or_and, linewidth=2.0)
plt.title("or and base optimization")
plt.xscale('log')
plt.xlabel("N")
plt.ylabel("flops/cycle")
plt.savefig("./plots/or_and_base_opt")