import numpy as np
import matplotlib.pyplot as plt
plt.style.use('D:/E/E/BIT_conference/experiment/CsBiBr/CsBiBr_kc') # this can be commented out

sqefilepath = "SQE_300K.dat"
sqedata = np.log10(np.loadtxt(sqefilepath))
print("sqedata shape:", sqedata.shape)

nqh = 100
deltaH = 0.005

ne = 100
deltaE = 0.5
E_array = np.arange(0, (ne + 1) * deltaE, deltaE)
print("E_array shape:", E_array.shape)

H04 = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5])
weneed = [int(i / deltaH) - 1 for i in H04]  # Corrected calculation method
print("weneed:", weneed)

TA2energyrange = np.array([[0, 15],
                          [3, 9],
                          [8, 14],
                          [12, 22],
                          [18, 24],
                          [20, 26]])

TA2energyrange_index = (TA2energyrange / deltaE).astype(int)  # Using deltaE=0.5
print("TA2energyrange_index:")
print(TA2energyrange_index)

qlabel = ['q = 0.05', 'q = 0.1', 'q = 0.2', 'q = 0.3', 'q = 0.4', 'q = 0.5']

fig = plt.figure(figsize=(4.5,2.8))
ax = fig.gca()

for i in range(len(H04)):
    print(f"\nPlotting curve {i}: {qlabel[i]}")
    
    # Ensure indices are within valid range
    start_idx = max(0, TA2energyrange_index[i, 0])
    end_idx = min(len(E_array), TA2energyrange_index[i, 1])
    
    print(f"Indices: {start_idx} to {end_idx}")
    print(f"Energy range: {E_array[start_idx]:.2f} to {E_array[end_idx-1]:.2f} meV")
    
    # Check if sqedata index is valid
    if weneed[i] >= sqedata.shape[0]:
        print(f"Warning: weneed[{i}] = {weneed[i]} exceeds sqedata first dimension")
        continue
        
    # Extract data
    energy_slice = E_array[start_idx:end_idx]
    data_slice = sqedata[weneed[i], start_idx:end_idx]
    
    print(f"Data points: {len(data_slice)}")
    
    # Plot the curve
    ax.plot(energy_slice, data_slice, ls='-', marker='o', markersize=3, 
            linewidth=1.5, label=qlabel[i])

# Add legend and labels
ax.legend(loc='best', ncol=2, fontsize=10)
ax.set_xlabel('Energy (meV)', fontsize=10)
ax.set_ylabel('Intensity (a.u.)', fontsize=10)
ax.tick_params(axis='x', labelsize=10)
ax.tick_params(axis='y', labelsize=10)
# ax.set_title('S(Q,E) of TA2 in GaN along H04', fontsize=16)
# ax.grid(True, alpha=0.3)

# Set x-axis range
ax.set_xlim([0, 30])
ax.set_ylim([-3, 3])

plt.tight_layout()
# plt.show()
plt.savefig("TA2_in_GaN_along_H04.jpg",dpi=650)