import numpy as np
import os
import shutil

def update_beta_and_twist(parameters_file, beta_change, twist, volume):
    with open(parameters_file, 'r') as file:
        lines = file.readlines()
    
    with open(parameters_file, 'w') as file:
        for line in lines:
            if line.startswith("beta"):
                file.write(f"beta                           {beta_change}\n")
            elif line.startswith("twist_coeff"):
                file.write(f"twist_coeff                    {twist}\n") 
            elif line.startswith("lattice size"):
                file.write(f"lattice size                        {volume[0]},{volume[1]},{volume[2]},{volume[3]}\n")
            else:
                file.write(line)

def run_script(parameters_file):
    betas = np.linspace(10.85,30, 10)
    twists = [1,2]
    volume = [44,44,64,6]
    for beta in betas:
        for twist in twists:
            update_beta_and_twist(parameters_file, beta, twist, volume)
            
            folder_name = f"beta-{beta:.3f}-twist-{twist}-{volume[0]}-{volume[1]}-{volume[2]}-{volume[3]}"
            if not os.path.exists(folder_name):
                os.makedirs(folder_name, exist_ok=True)
                
                os.system(f"./build/suN_gauge_twist > {folder_name}/out.txt")
                
                for file_name in os.listdir('.'):
                    if file_name.startswith("fourier_profile_") or file_name.startswith("surface_smooth_"):
                        shutil.move(file_name, os.path.join(folder_name, file_name))
            else:
                print(f"Folder {folder_name} already exists. Skipping...")
            
    
    
if __name__ == '__main__':
    run_script("/home/haaaaron/HILA/applications/suN_gauge_twist/parameters")