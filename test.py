# import digital_lab, ntpms
import main
import numpy as np
import matplotlib.pyplot as plt
import digital_lab, ntpms


# main.meshConvergence(0, 4, 0.131075844, 0.3, 0.2, 0.05, 100)
CSV_file = "tpms_list_"+main.now.strftime("%Y_%m_%d_%H%M")+".csv"
Header = ['Count', 'Lattice Type', 'Cell Size', 'Thickness', 'Porosity', 'Mesh Size', 'Stress', 'Strain', 'Stiffness']
output_file = open(CSV_file, mode='w')
output_file.write(",".join(Header)+'\n')
output_file.close()

thickness = [0.04392770758, 0.06046035723, 0.07747550122, 0.09383708602, 0.1066880127, 0.137524254, 0.1330446372]
cell = [1.5, 2, 2.5, 3, 3.5, 4, 4.5]

count = 0

for i in range(0,len(thickness)):
    main.generate(CSV_file,0,cell[i],thickness[i],1,100,count)
    count += 1

# print(ntpms.mesh(0,4,0.12355126121064888,0.75, 0))
# _, strain = digital_lab.compressionV2("tpms_mesh.cdb",100,1)
# stress = 6.25 / (4 * 4 * 0.000001)     # convert back to m^2
# stiffness = stress / strain
# print("stiffness")
# print(stiffness)
# _, strain = digital_lab.compressionV2("tpms_mesh.cdb",100,1)
# stress = 6.25 / (4 * 4 * 0.000001)     # convert back to m^2
# stiffness = stress / strain
# print("stiffness")
# print(stiffness)
# print(ntpms.mesh(0,5,0.125,1))
# digital_lab.compression("tpms_mesh.cdb",100,1)
# print(ntpms.mesh(0,5,0.3,12))
# digital_lab.compression("tpms_mesh.cdb",100,12)
# print(ntpms.mesh(0,5,0.2,13))
# digital_lab.compression("tpms_mesh.cdb",100,13)
# print(ntpms.mesh(0,5,0.1,14))
# digital_lab.compression("tpms_mesh.cdb",100,14)
