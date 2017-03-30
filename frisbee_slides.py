# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 12:12:24 2017

@author: dlunde
"""

def Sliders():
    #from frisbee_current import Simulation
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider,Button
    
    input_vec = [0,1,0,14,3,0,0,10,0,8,0]
    plt.close("Frisbee Sliders")
    
    fig_sliders = plt.figure("Frisbee Sliders")
    fig_sliders.canvas.set_window_title('Frisbee Sliders') 
    
    heights = [0.92]*11
    for hi in range (1,11,1):
        heights[hi]=heights[hi-1]-0.08
    
    ax1= plt.axes([0.2, heights[0], 0.7, 0.05])
    ax2= plt.axes([0.2, heights[1], 0.7, 0.05])
    ax3= plt.axes([0.2, heights[2], 0.7, 0.05])
    ax4= plt.axes([0.2, heights[3], 0.7, 0.05])
    ax5= plt.axes([0.2, heights[4], 0.7, 0.05])
    ax6= plt.axes([0.2, heights[5], 0.7, 0.05])
    ax7= plt.axes([0.2, heights[6], 0.7, 0.05])
    ax8= plt.axes([0.2, heights[7], 0.7, 0.05])
    ax9= plt.axes([0.2, heights[8], 0.7, 0.05])
    ax10=plt.axes([0.2, heights[9], 0.7, 0.05])
    ax11=plt.axes([0.2,heights[10], 0.7, 0.05])
    
    sli1 = Slider(ax1, 'X Pos (m)', -10, 10, valinit=input_vec[0])
    sli2 = Slider(ax2, 'Y Pos (m)', 0, 5, valinit=input_vec[1])
    sli3 = Slider(ax3, 'Z Pos (m)', -3, 3, valinit=input_vec[2])
    sli4 = Slider(ax4, 'X Vel (m/s)', 0, 30, valinit=input_vec[3])
    sli5 = Slider(ax5, 'Y Vel (m/s)', -5, 25, valinit=input_vec[4])
    sli6 = Slider(ax6, 'Z Vel (m/s)', -15, 15, valinit=input_vec[5])
    sli7 = Slider(ax7, 'X A Vel (rot/sec)', -20, 20, valinit=input_vec[6])
    sli8 = Slider(ax8, 'Y A Vel (rot/sec)', -20, 20, valinit=input_vec[7])
    sli9 = Slider(ax9, 'X A Vel (rot/sec)', -20, 20, valinit=input_vec[8])
    sli10 = Slider(ax10, 'Alpha (degrees)', -70, 180, valinit=input_vec[9])
    sli11 = Slider(ax11, 'Phi (degrees)', -150, 150, valinit=input_vec[10])
    
    startax = plt.axes([0.7,0.01,0.13,0.1])
    start_button = Button(startax, 'Start', color = 'g', hovercolor='0.75')
    
    resetax = plt.axes([0.55,0.01,0.13,0.1])
    reset_button = Button(resetax, 'Reset', color = 'b', hovercolor='0.75')
    
    printax = plt.axes([0.4,0.01,0.13,0.1])
    print_button = Button(printax, 'Print', color = 'r', hovercolor='0.75')
    
    gifax = plt.axes([0.85,0.01,0.07,0.05])
    gif_button = Button(gifax, 'Gif it', color = 'y', hovercolor='0.75')
    
    #def change_val(val):
        
    #print xsli.__getattribute__('val')
    #xsli.on_changed(change_val)
    
    def write_input():
        ff = open('frisbee_input.txt','w')
        print >> ff,"X Pos,Y Pos,Z Pos,X Vel,Y Vel,Z Vel,X A Vel,Y A Vel,Z A Vel,alpha,phi"
        for each in range(0,10,1):
            print >> ff,input_vec[each],",",
        print >> ff,input_vec[10]
        ff.close()
        return
    
    def get_vals():
        vals_vec = [0.0]*11
        vals_vec[0] = sli1.__getattribute__('val')
        vals_vec[1] = sli2.__getattribute__('val')
        vals_vec[2] = sli3.__getattribute__('val')
        vals_vec[3] = sli4.__getattribute__('val')
        vals_vec[4] = sli5.__getattribute__('val')
        vals_vec[5] = sli6.__getattribute__('val')
        vals_vec[6] = sli7.__getattribute__('val')
        vals_vec[7] = sli8.__getattribute__('val')
        vals_vec[8] = sli9.__getattribute__('val')
        vals_vec[9] =sli10.__getattribute__('val')
        vals_vec[10]=sli11.__getattribute__('val')
        return vals_vec    
        
    def set_inputs(event):
    #    print "from sliders",input_vec
    #    write_input()
    #    execfile("frisbee_current.py")
        from frisbee_current import Simulation
        Simulation(CHARLIE=False,overwrite=get_vals(),clear_old=True)
        return
    
    
    start_button.on_clicked(set_inputs)
    
    def reset_inputs(event):
        input_vec = [0,1,0,14,3,0,0,10,0,8,0]
        sli1.set_val(input_vec[0])
        sli2.set_val(input_vec[1])
        sli3.set_val(input_vec[2])
        sli4.set_val(input_vec[3])
        sli5.set_val(input_vec[4])
        sli6.set_val(input_vec[5])
        sli7.set_val(input_vec[6])
        sli8.set_val(input_vec[7])
        sli9.set_val(input_vec[8])
        sli10.set_val(input_vec[9])
        sli11.set_val(input_vec[10])
        return
    reset_button.on_clicked(reset_inputs)
    
    def print_inputs(event):
        for each in get_vals():
            print int(each),
        print ""
        return
        
    print_button.on_clicked(print_inputs)
    
    def gif_run(event):
        from frisbee_current import Simulation
        Simulation(CHARLIE=False,gif=True,overwrite=get_vals())
        return
    
    gif_button.on_clicked(gif_run)
    
    plt.show(block=True) #this doesn't allow it to be closed with close
    #fig_sliders.canvas.manager.window.attributes('-topmost', 1)
    return
    
Sliders()