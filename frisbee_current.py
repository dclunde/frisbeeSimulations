# -*- coding: utf-8 -*-
"""
Created on Thu Sep 29 14:44:42 2016

@author: dlunde
"""
    
def Simulation(overwrite = [],**kwargs):
#     global x,y,z,vx,vy,vz,deltavx,deltavy,deltavz,delta_t,mass,\
#        delta_omega_x,delta_omega_y,delta_omega_z,Ixx,Iyy,Izz,omega_x,omega_y,omega_z,\
#        alpha,phi,mat_x,mat_y,mat_z,mat_vz,mat_vy,mat_vz,mat_a,mat_p,mat_t,mat_wx,mat_wy,mat_wz
    

    from frisbee_current_header import *
    from frisbee_current_def import *

    PRINT = False         #If true print outputs

    #@profile for timing things
#    from frisbee_current_header import grab_from_file,initialize_variables
    if 'choose_throw_type' in kwargs:    
        choose_throw_type()
    if 'edit_initial_conditions' in kwargs:
        edit_initial_conditions()
    if 'change_temperature' in kwargs:
        change_temperature()
    if 'which_planet' in kwargs:
        which_planet()
    if 'grab_from_file' in kwargs:
        grab_from_file()
        
#    choose_throw_type()
#    edit_initial_conditions()
    #which_planet()
    initial_time = time.clock()
    if PRINT:
        print "Beginning Frisbee Simulation"    


    if not overwrite==[]:
#        x,y,z,vx,vy,vz,omega_x,omega_y,omega_z,alpha,phi=overwrite
        initialize_variables(overwrite=overwrite)
    else:
        initialize_variables()

#    else:
    x,y,z,vx,vy,vz,omega_x,omega_y,omega_z,alpha,phi=\
        mat_x[0],mat_y[0],mat_z[0],mat_vx[0],mat_vy[0],mat_vz[0],\
        mat_wx[0],mat_wy[0],mat_wz[0],mat_a[0],mat_p[0],
    
    #    frisbee_current_header.initialize_variables()
    
#    phi = phi + 90
#    omega_x = omega_x * 2 * math.pi
#    omega_y = omega_y * 2 * math.pi
#    omega_z = omega_z * 2 * math.pi
    
#    print "in current",x,y,z,vx,vy,vz,omega_x,omega_y,omega_z,alpha,phi
     
    
    plt.close("Frisbee Simulation") #Close Old Plots
    plt.close('Frisbee Analysis')   #Close Old Plots
    plt.close('Frisbee Energy')     #Close Old Plots
    CONTINUE = True                 #Boolean for simulation continuation
    
    i=0     #Time Step variable
    update_vectors(0)               #Initialize vectors
    #    cl,cd = lift_drag(i)             #Initialize lift and drag
        
    i = 1                           #Start time step
    while CONTINUE: 
        
        #velocity = vx*vx + vy*vy + vz*vz
        forces_sheet(forces(i),i)    #Send forces into the output file
        
        forces_out = forces(i)       #Create a list because it doesn't work otherwise
        
        deltavx = forces_out[0] + forces_out[3] + forces_out[6] + forces_out[9]#  + forces_out[12] + forces_out[15]
        deltavy = forces_out[1] + forces_out[4] + forces_out[7] + forces_out[10]# + forces_out[13] + forces_out[16]
        deltavz = forces_out[2] + forces_out[5] + forces_out[8] + forces_out[11]# + forces_out[14] + forces_out[17]
    
        # Delta V * mass = Force * time
        deltavx = deltavx * delta_t / mass
        deltavy = deltavy * delta_t / mass    
        deltavz = deltavz * delta_t / mass
           
        # Update Velocity and Position
        vx = vx + deltavx
        vy = vy + deltavy
        vz = vz + deltavz
        x  = x  + vx * delta_t
        y  = y  + vy * delta_t
        z  = z  + vz * delta_t
        
        
        # MOMENTS 
        
        moments_sheet(moments(i,delta_omega_x,delta_omega_y,delta_omega_z),i)    #Send forces into the output file
    
        moments_out= moments(i,delta_omega_x,delta_omega_y,delta_omega_z)
    
        delta_omega_x = moments_out[0] * delta_t / Ixx
        detla_omega_y = moments_out[1] * delta_t / Iyy
        delta_omega_z = moments_out[2] * delta_t / Izz
    #    delta_omega_y = W_mag + -rho * (np.linalg.norm(velocity_vec[i]))**2 *\
    #        area/5 * cd *  (1 / (2 * mass)) * -W_mag     
        
        # Update angular velocity 
        omega_x = omega_x + delta_omega_x
        omega_y = omega_y + delta_omega_y
        omega_z = omega_z + delta_omega_z
    
        #update alpha and phi
        alpha = alpha + omega_z * delta_t
        phi   = phi   + omega_x * delta_t
        
        #alpha = math.degrees(math.atan(vy/vx))
        #alpha = alpha + 1
         
        # Record Data
        mat_x[i],mat_y[i],mat_z[i]    = x,y,z
        mat_vx[i],mat_vy[i],mat_vz[i] = vx,vy,vz
        mat_a[i],mat_p[i],mat_t[i]    = alpha,phi,delta_t*i
        mat_wx[i],mat_wy[i],mat_wz[i] = omega_x,omega_y,omega_z
    
        # Update Vectors and Cooeficients    
        update_vectors(i)
    #    cl,cd = lift_drag(i)
        
        if PRINT:
            if i*delta_t*100%10 == 0:
                print i*delta_t,"sec\r",
        
        #Bounce
    #    if y < diameter and 0<angle_between(frisbee_vec[i],[0,1,0],degrees=True)<10:
    #        vy = - 0.5 * vy
    #        vx = 0.8 * vx
    #        alpha = - 0.3 * alpha 
    #        print "BOUNCE!"
        
        i = i + 1    
        
        #Unrealisitic stop
        if absmax([vx,vy,vz]) > 100:
            print "Frisbee unrealisitic %f,%f,%f, ending simulation" % (vx,vy,vz)
            CONTINUE = False
        
        #Hit Ground
        if y <= 0.0:
            print "Frisbee hit the ground, ending simulation"
            CONTINUE = False
            
        #Run out of time
        if i >= total_t/delta_t-1:
            print "Ran out of time, ending simulation"
            CONTINUE = False
            
            
        if i > 3/delta_t:
            CONTINUE = False

    
    if PRINT:
        print "Finished simulation in",time.clock()-initial_time,"seconds"                

    #plot_frisbee(i,step=True)

    if 'plot_frisbee' in kwargs:
        plot_frisbee(i,**kwargs)
    if 'analyze_plots' in kwargs:
        analyze_plots(i)
    if 'excel_sheet' in kwargs:
        excel_sheet(i)
    if 'calculate_energy' in kwargs:
        calculate_energy(i,plot=True,save_data=True)
    if 'analyze' in kwargs:
        analyze(i)
    if 'clear_old' in kwargs:
        if kwargs['clear_old']:
            clear_mats(i)
    
    if PRINT:
        print "Finished everything  in",time.clock()-initial_time,"seconds"
    if 'returnLast' in kwargs:
        return i
    else:
        return
    
    

#Simulation(analyze_plots=1,calculate_energy=1)
