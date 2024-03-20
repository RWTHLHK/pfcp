import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":
  data = pd.read_csv("/home/rwth1393/projects/pfcp/problems/update_method_test_out.csv")
  plt.plot(data['e_zz'].values, data['uncracked_cauchy_stress'].values, label='uncracked_cauchy_stress')
  plt.plot(data['e_zz'].values, data['cracked_cauchy_stress'].values, label='cracked_cauchy_stress')
  plt.xlabel('Strain')
  plt.ylabel('Stress(MPa)')
  plt.legend()
  plt.savefig("cracked_stress.png")

