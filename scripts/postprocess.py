import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":
  data = pd.read_csv("/home/rwth1393/projects/pfcp/problems/update_method_test_out.csv")
  plt.plot(data['e_zz'].values, data['d'].values, label='damage evolution')
  plt.xlabel('Strain')
  plt.ylabel('damage')
  plt.legend()
  plt.savefig("damage.png")

