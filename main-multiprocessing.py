import sys
import os
import time
import math
import csv
import click
import numpy as np
import logging
import matplotlib.pyplot as plt
from ising import run_ising #import run_ising function from ising.py
import multiprocessing as mp

def run_simulation(index,temp,n,num_steps,num_burnin,num_analysis,flip_prop,j,b,data_filename,corr_filename,data_listener,corr_listener,plots):
    print("Working on Temp {0}".format(round(temp,3)))
    if plots:
        #initialize vars for plotting values
        temp_arr, M_mean_arr, E_mean_arr, M_std_arr, E_std_arr, CV_arr, X_arr, CV_err_arr, X_err_arr = [],[],[],[],[],[],[],[],[]
    try:
        #run the Ising model
        Msamp, Esamp, spin = run_ising(n,temp,num_steps,num_burnin,flip_prop,j,b,disable_tqdm=True)

        try:
            K_b = 1
            B = 1/(K_b*temp)
            M_mean = np.average(Msamp[-num_analysis:])
            E_mean = np.average(Esamp[-num_analysis:])
            M2 = [x**2 for x in Msamp]
            E2 = [x**2 for x in Esamp]
            M2_mean = np.average(M2[-num_analysis:])
            E2_mean = np.average(E2[-num_analysis:])
            M_std = np.std(Msamp[-num_analysis:])
            E_std = np.std(Esamp[-num_analysis:])
            M2_std = np.std(M2[-num_analysis:])
            E2_std = np.std(E2[-num_analysis:])
            CV = (B/temp)*(E2_mean-(E_mean**2))
            X = B*(M2_mean-(M_mean**2))
            CV_err = math.sqrt((E2_std**2)+(2*E_mean*E_std)**2)
            X_err = math.sqrt(M2_std**2+(2*M_mean*M_std)**2)

            data_array = [np.abs(M_mean),M_std,E_mean,E_std,CV,X,CV_err,X_err]
            data_listener.put([temp]+data_array)

            corr = compute_autocorrelation(spin)
            [corr_listener.put([temp]+corr_value) for corr_value in corr]

            print("Done with Temp {0}".format(round(temp,3)))
            return True
			
            if plots:
                #for plotting
                M_mean, E_mean, M_std, E_std, CV, X, CV_err, X_err = get_plot_values(temp,Msamp,Esamp,num_analysis)
                temp_arr.append(temp)
                M_mean_arr.append(M_mean)
                E_mean_arr.append(E_mean)
                M_std_arr.append(M_std)
                E_std_arr.append(E_std)
                CV_arr.append(CV)
                X_arr.append(X)
                CV_err_arr.append(CV_err)
                X_err_arr.append(X_err)
        except:
            logging.error("Temp="+str(round(temp,3))+": Statistical Calculation Failed. No Data Written.")
            return False

    except KeyboardInterrupt:
        print("\n\nProgram Terminated. Good Bye!")
        data_listener.put('kill')
        corr_listener.put('kill')
        sys.exit()

    except:
        logging.error("Temp="+str(round(temp,3))+": Simulation Failed. No Data Written")
    
    if plots:
        plot_graphs(temp_arr, M_mean_arr, E_mean_arr, M_std_arr, E_std_arr, CV_arr, X_arr, CV_err_arr, X_err_arr)
        print("plotting")
		
#simulation options (enter python main.py --help for details)
@click.command()
@click.option('--t_min', default=2.0, prompt='Minimum Temp', help='Minimum Temperature (inclusive)', type=float)
@click.option('--t_max', default=2.6, prompt='Maximum Temp', help='Maximum Temperature (inclusive)', type=float)
@click.option('--t_step', default=0.1, prompt='Temp Step Size', help='Temperature Step Size', type=float)

@click.option('--n', prompt='Lattice Size', help='Lattice Size (NxN)',type=int)
@click.option('--num_steps', default=100000, help='Total Number of Steps',type=int)
@click.option('--num_analysis', default=50000, help='Number of Steps used in Analysis',type=int)
@click.option('--num_burnin', default=0, help='Total Number of Burnin Steps',type=int)

@click.option('--j', default=1.0, help='Interaction Strength',type=float)
@click.option('--b', default=0.0, help='Applied Magnetic Field',type=float)
@click.option('--flip_prop', default=0.1, help='Proportion of Spins to Consider Flipping per Step',type=float)

@click.option('--output', default='data', help='Directory Name for Data Output',type=str)
@click.option('--plots', default=True, help='Turn Automatic Plot Creation Off or On',type=bool)

def main(t_min,t_max,t_step,n,num_steps,num_analysis,num_burnin,j,b,flip_prop,output,plots):
    data_filename, corr_filename = initialize_simulation(n,num_steps,num_analysis,num_burnin,output,j,b,flip_prop)
    run_processes(t_min,t_max,t_step,n,num_steps,num_burnin,num_analysis,flip_prop,j,b,data_filename,corr_filename,plots)
    print('\n\nSimulation Finished! Data written to '+ data_filename)

def initialize_simulation(n,num_steps,num_analysis,num_burnin,output,j,b,flip_prop):
    check_step_values(num_steps, num_analysis, num_burnin)
    data_filename, corr_filename = get_filenames(output)
    write_sim_parameters(data_filename,corr_filename,n,num_steps,num_analysis,num_burnin,j,b,flip_prop)
    print('\nSimulation Started! Data will be written to ' + data_filename + '\n')
    return data_filename, corr_filename
	
def get_plot_values(temp,Msamp,Esamp,num_analysis): #only for plotting at end
    try:
        K_b = 1
        B = 1/(K_b*temp)
        M_mean = np.average(Msamp[-num_analysis:])
        E_mean = np.average(Esamp[-num_analysis:])
        M2_mean = np.average((Msamp[-num_analysis:])**2)
        E2_mean = np.average((Esamp[-num_analysis:])**2)
        M_std = np.std(Msamp[-num_analysis:])
        E_std = np.std(Esamp[-num_analysis:])
        M2_std = np.std((Msamp[-num_analysis:])**2)
        E2_std = np.std((Esamp[-num_analysis:])**2)
        CV = (B/temp)*(E2_mean-(E_mean**2))
        X = B*(M2_mean-(M_mean**2))
        CV_err = math.sqrt(E2_std+2*E_mean*E_std)
        X_err = math.sqrt(M2_std+2*M_mean*M_std)
        return M_mean, E_mean, M_std, E_std, CV, X, CV_err, X_err
    except:
        logging.error("Temp={0}: Error getting plot values".format(temp))
        return False

def plot_graphs(temp_arr,M_mean_arr,M_std_arr,E_mean_arr,E_std_arr, CV_arr, X_arr, CV_err_arr, X_err_arr): #plot graphs at end
    plt.figure(1)
    plt.ylim(0,1)
    plt.errorbar(temp_arr, np.absolute(M_mean_arr), yerr=M_std_arr, uplims=True, lolims=True,fmt='o')
    plt.xlabel('Temperature')
    plt.ylabel('Magnetization')
    plt.figure(2)
    plt.errorbar(temp_arr, E_mean_arr, yerr=E_std_arr, fmt='o')
    plt.xlabel('Temperature')
    plt.ylabel('Energy')
    plt.figure(3)
    plt.errorbar(temp_arr, CV_arr, yerr=CV_err_arr, fmt='o')
    plt.xlabel('Temperature')
    plt.ylabel('CV')
    plt.figure(4)
    plt.errorbar(temp_arr, X_arr, yerr=X_err_arr, fmt='o')
    plt.xlabel('Temperature')
    plt.ylabel('X')
    plt.show()
	
def check_step_values(num_steps,num_analysis,num_burnin): #simulation size checks and exceptions
    if (num_burnin > num_steps):
        raise ValueError('num_burning cannot be greater than available num_steps. Exiting simulation.')
        sys.exit()

    if (num_analysis > num_steps - num_burnin):
        raise ValueError('num_analysis cannot be greater than available num_steps after burnin. Exiting simulation.')
        sys.exit()

def get_filenames(dirname): #make data folder if doesn't exist, then specify filename
    try:
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        data_filename = os.path.join(dirname,'data_'+str(time.strftime("%Y%m%d-%H%M%S"))+".csv")
        corr_filename = os.path.join(dirname,'corr_'+str(time.strftime("%Y%m%d-%H%M%S"))+".csv")
        #Write simulation parameters to file
        return data_filename, corr_filename
    except:
        raise ValueError('Directory name not valid. Exiting simulation.')
        sys.exit()

def get_temp_array(t_min,t_max,t_step):
    if (t_min > t_max):
        raise ValueError('T_min cannot be greater than T_max. Exiting Simulation')
        sys.exit()
    try:
        T = np.arange(t_min,t_max,t_step).tolist()
        return T
    except:
        raise ValueError('Error creating temperature array. Exiting simulation.')
        sys.exit()

def write_sim_parameters(data_filename,corr_filename,n,num_steps,num_analysis,num_burnin,j,b,flip_prop):
    try:
        with open(data_filename,'w') as csv_file:
            writer = csv.writer(csv_file, delimiter=',', lineterminator='\n')
            #Write simulations parameters to CSV file
            writer.writerow(['Lattice Size (NxN)','Total Steps','Steps Used in Analysis','Burnin Steps','Interaction Strength','Applied Mag Field','Spin Prop'])
            writer.writerow([n,num_steps,num_analysis,num_burnin,j,b,flip_prop])
            writer.writerow([])
            writer.writerow(['Temp','Mean Mag','SD Mag','Mean Enegy','SD Energy','CV','X','CV_err','X_err'])
        with open(corr_filename,'w') as csv_file:
            writer = csv.writer(csv_file, delimiter=',', lineterminator='\n')
            #Write simulations parameters to CSV file
            writer.writerow(['Lattice Size (NxN)','Total Steps','Steps Used in Analysis','Burnin Steps','Interaction Strength','Applied Mag Field','Spin Prop'])
            writer.writerow([n,num_steps,num_analysis,num_burnin,j,b,flip_prop])
            writer.writerow([])
            writer.writerow(['Temp','Mean Mag','SD Mag','Mean Enegy','SD Energy','CV','X','CV_err','X_err'])
    except:
        logging.error('Could not save simulation parameters. Exiting simulation')
        sys.exit()

def compute_autocorrelation(spin):
    n = len(spin)
    corr_array = []
    for k in range(1,int(n/2)):
        col_mean, row_mean = spin.mean(axis=0),spin.mean(axis=1)
        #compute r values for rows and cols
        r_col = [np.multiply(spin[j,:]-col_mean,spin[(j+k)%n,:]-col_mean) for j in range(1,n)]
        r_row = [np.multiply(spin[:,j]-row_mean,spin[:,(j+k)%n]-row_mean) for j in range(1,n)]
        #normalize r values
        r_col = np.divide(r_col,float(n))
        r_row = np.divide(r_row,float(n))
        #calculate corr for k and add it to array
        corr = (r_col.mean() + r_row.mean())/2.0
        corr_array.append([k,corr])
    return corr_array

def listener(q, fn):
    '''listens for messages on the q, writes to file. '''
    f = open(fn, 'a') 
    writer = csv.writer(f, delimiter=',', lineterminator='\n')
    while 1:
        m = q.get()
        if m == 'kill':
            break
        writer.writerow(m)
        f.flush()
    f.close()

def run_processes(t_min,t_max,t_step,n,num_steps,num_burnin,num_analysis,flip_prop,j,b,data_filename,corr_filename,plots):
    
    T = get_temp_array(t_min, t_max, t_step)
    
    #must use Manager queue here, or will not work
    manager = mp.Manager()
    data_listener = manager.Queue()
    corr_listener = manager.Queue()    
    pool = mp.Pool(mp.cpu_count())

    #put listener to work first
    data_watcher = pool.apply_async(listener, args=(data_listener, data_filename,))
    corr_watcher = pool.apply_async(listener, args=(corr_listener, corr_filename,))

    #fire off workers 
    jobs = [pool.apply_async(run_simulation, args=(index,temp,n,num_steps,num_burnin,num_analysis,flip_prop,j,b,data_filename,corr_filename,data_listener,corr_listener,plots)) for index,temp in enumerate(T)]

    # collect results from the workers through the pool result queue   
    [job.get() for job in jobs]

    #now we are done, kill the listener
    data_listener.put('kill')
    corr_listener.put('kill')
    pool.close()

if __name__ == "__main__":
   main()