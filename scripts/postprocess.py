import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":
  data = pd.read_csv("/home/rwth1393/projects/pfcp/problems/update_method_test_out.csv")
  plt.plot(data['e_zz'].values, data['uncracked_pk2'].values, label='uncracked_pk2')
  plt.plot(data['e_zz'].values, data['cracked_pk2'].values, label='cracked_pk2')
  plt.xlabel('Strain')
  plt.ylabel('Stress(MPa)')
  plt.legend()
  plt.savefig("cracked_stress.png")

